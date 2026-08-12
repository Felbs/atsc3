#!/usr/bin/env python3
"""m11_watch.py -- live ATSC 3.0: tune, decode and WATCH, continuously.

    python -m atsc3 watch --rf 33

M9 measured the pieces and left the plumbing as a design.  This is the
plumbing, wired as M9 specified it:

    reader  ->  front end  ->  decoder  ->  transport  ->  player
    (SDR)       (Stage 1)      (Stage 2)    (Stage 3/4)   (ffplay on stdin)

each on its own thread with a bounded queue in front of it, because NumPy and
SciPy release the GIL inside their kernels -- which is the same reason M9's
front end is threaded rather than forked, and it matters twice over here: a
live pipeline is ONE stream, and process parallelism buys throughput by
spending latency.

WHAT IS MEASURED AND REPORTED, ALWAYS
-------------------------------------
Sustained x-real-time, instantaneous x-real-time, every queue's occupancy, the
player's buffer in seconds, and the underrun count.  A player that stutters is
not "working", so the stutter has to have a number.  The occupancy meter is
the honest kind: seconds of media pushed minus seconds of wall clock since
playback began.  It does not ask the player how it is doing, because a pipe
cannot answer.

The radio side is fleet-standard and stated in `m11_stream.SdrSource`:
priority 60, polite wait, heartbeat on a TIMER inside the read loop, yield to
anything that outranks us, release in `finally`, and never inside the meteor
window.
"""
from __future__ import annotations

import argparse
import collections
import datetime as dt
import hashlib
import json
import os
import queue
import signal
import sys
import threading
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
sys.path.insert(0, HERE)

import m11_stream as ST                                           # noqa: E402
import m6_cells as C                                              # noqa: E402
import m6_payload as P6                                           # noqa: E402

_sib = os.environ.get("ATSC3_SIBLING_TOOLS")            # optional sibling toolbox
if _sib:
    sys.path.insert(0, _sib)
try:
    import radio_lock
except Exception:                                                 # noqa: BLE001
    radio_lock = None

FRAME_SEC = ST.FRAME_SEC


def _decode_proc_main(shm_name, nslots, slot_len, in_q, out_q, cfg):
    """One decode worker PROCESS (E52).

    Python threads cannot scale this decoder past ~1.2x on any core count:
    the chain is gather/scatter-heavy and NumPy holds the GIL for advanced
    indexing, so 1->6 threads bought 36% (measured, TR 2990WX).  Processes
    sidestep the GIL entirely.  Frame windows travel by shared memory (a
    5.7 MB memcpy, ~1 ms); decoded packets travel back by queue (bytes
    pickle as a memcpy).  Each worker owns a full FrameDecoder; every cache
    it builds is a pure function of the code tables, so any worker produces
    byte-identical output for a given (window, t0) -- the main process
    re-imposes FIFO order, making the whole arrangement byte-identical to
    the serial decoder.  The parent gates that, per run, end to end.
    """
    import traceback as tb
    try:
        os.environ.setdefault("M9_NO_TORCH", "1")
        from multiprocessing import shared_memory
        import numpy as np
        import m9_fast
        shm = shared_memory.SharedMemory(name=shm_name)
        buf = np.ndarray((nslots, slot_len), dtype=np.complex128,
                         buffer=shm.buf)
        fd = m9_fast.FrameDecoder(threads=cfg["threads"], backend="cpu",
                                  iters=cfg["iters"],
                                  cpu_fast=cfg["cpu_fast"])
        fd.prewarm()
        out_q.put(("ready", os.getpid()))
        while True:
            msg = in_q.get()
            if msg is None:
                break
            slot, idx, n, t0, coh = msg
            t = time.perf_counter()
            r16, pkts, diag = fd.decode_frame(buf[slot, :n], t0)
            out_q.put(("ok", (idx, slot, pkts,
                              dict(converged=diag["converged"],
                                   bch_ok=diag["bch_ok"]),
                              time.perf_counter() - t, coh)))
        out_q.put(("done", dict(fd.tm)))
        fd.close()
        shm.close()
    except Exception as e:                                     # noqa: BLE001
        try:
            out_q.put(("err", f"{type(e).__name__}: {e}\n{tb.format_exc()}"))
        except Exception:                                      # noqa: BLE001
            pass


def in_meteor_window(now=None):
    now = now or dt.datetime.now()
    mins = now.hour * 60 + now.minute
    return 2 * 60 + 10 <= mins <= 5 * 60 + 40


def player_cmd(kind, title, extra=None):
    if kind == "none":
        return None
    if kind == "ffplay":
        c = ["ffplay", "-hide_banner", "-loglevel", "error", "-nostats",
             "-autoexit", "-window_title", title, "-i", "pipe:0"]
    elif kind == "mpv":
        c = ["mpv", "--title=" + title, "--profile=low-latency", "-"]
    else:
        c = kind.split()
    return c + (extra.split() if extra else [])


class Watch:
    def __init__(self, a):
        self.a = a
        self.stop = threading.Event()
        # LOSSLESS REPLAY (the port to slower machines made this load-bearing):
        # every shed point below exists to protect the SDR reader, whose only
        # answer to backpressure is dropping SAMPLES.  A capture replay that is
        # not throttled to real time has no radio to protect -- backpressure is
        # free, and a gate that sheds frames measures the MACHINE, not the
        # code.  Measured before this: the m11 e2e gate lost 30-50% of its
        # datagrams to frame shedding on every box slower than the capture
        # feed, first divergence always at the startup warmup stall.
        self.lossless = bool(a.capture) and not a.realtime and not a.shed
        self.q_raw = queue.Queue(maxsize=a.raw_queue)
        self.q_frame = queue.Queue(maxsize=a.frame_queue)
        self.q_bb = queue.Queue(maxsize=a.bb_queue)
        self.src = None
        self.fe = None
        self.fd = None
        self.tr = None
        self.player = None
        self.live = None
        self.ex = None
        self.dg_fh = None
        self.t_start = None
        self.n_frames = 0
        self.n_conv = 0
        self.n_bch = 0
        self.bb_abs = 0
        self.bad_run = 0
        self.stats = collections.Counter()
        self.hist = collections.deque(maxlen=64)
        self.err = None
        self.lock = threading.Lock()
        # E82 -- LDM / CTI mode.  `plan` is None until the multiplex has been
        # identified FROM L1; the Frame length, the PLP layout and the
        # interleaver are all signalled, so none of them may be a constant
        # here.  `frame_sec` replaces the module-level FRAME_SEC in every
        # report, because a Frame is 247.1 ms on RF33 and 242.5 ms on RF25.
        self.ldm = False
        self.plan = None
        self.pipe = None
        self.frame_sec = FRAME_SEC
        self.blocks_per_frame = 74

    def _put_blocking(self, q, item):
        """Backpressure put for lossless replay: wait for the consumer, but
        keep watching the stop flag so shutdown can never deadlock a full
        queue against an exited consumer.  Returns False iff stopping."""
        while not self.stop.is_set():
            try:
                q.put(item, timeout=0.5)
                return True
            except queue.Full:
                pass
        return False

    # -- threads -----------------------------------------------------------
    def t_reader(self):
        try:
            # E82: the identification sniff already consumed the first blocks.
            # They are pushed here, IN ORDER, before anything new is read, so
            # the front end sees the same sample stream it would have seen if
            # no sniff had happened -- which is what makes the RF33
            # byte-identity claim survive the addition.
            pre = getattr(self, "_sniffed", None) or []
            self._sniffed = []
            if not pre:
                self.src.start()
            for x in pre:
                if self.stop.is_set():
                    break
                if isinstance(x, bytes):
                    self.q_raw.put("BREAK")
                    continue
                if self.lossless:
                    if not self._put_blocking(self.q_raw, x):
                        break
                else:
                    self.q_raw.put(x)
            while not self.stop.is_set():
                x = self.src.read()
                if x is None:
                    break
                if isinstance(x, bytes):            # a continuity break
                    self.q_raw.put("BREAK")
                    continue
                if self.lossless:
                    if not self._put_blocking(self.q_raw, x):
                        break                     # stopping; drain and leave
                    continue
                try:
                    self.q_raw.put(x, timeout=5.0)
                except queue.Full:
                    # A dropped block is a hole in the sample stream, and a
                    # hole is exactly what t0 tracking cannot survive.  Saying
                    # so costs a re-acquisition; NOT saying so would let the
                    # front end keep counting +FRAME_SAMPLES across a gap and
                    # report healthy Frames that are not where it thinks.
                    self.stats["raw_dropped"] += 1
                    self.q_raw.put("BREAK")
        except Exception as e:                                 # noqa: BLE001
            import traceback
            ST.log(f"  *** reader thread died ***\n{traceback.format_exc()}")
            self.err = self.err or f"reader: {type(e).__name__}: {e}"
        finally:
            self.q_raw.put(None)

    def t_front(self):
        try:
            while True:
                x = self.q_raw.get()
                if x is None:
                    break
                if isinstance(x, str):
                    ST.log("  *** sample continuity BREAK (SDR overflow) -- "
                           "re-acquiring the bootstrap ***")
                    self.fe.reacquire()
                    continue
                self.fe.push(x)
                for fr in self.fe.frames():
                    if self.a.frames and fr[0] >= self.a.frames:
                        raise StopIteration
                    # STEP 2 -- DROP A FRAME, NEVER A SAMPLE.
                    #
                    # The old `put(timeout=5.0)` stalled the front end for up
                    # to five seconds when the decoder fell behind, which let
                    # q_raw fill, which made the reader drop samples. That is
                    # the expensive trade backwards: a dropped Frame costs
                    # 247 ms of picture and nothing else, while a dropped
                    # SAMPLE breaks sample continuity, costs the bootstrap
                    # lock and ~15 s of re-acquisition. Shed here, instantly,
                    # so the pressure never reaches the radio.  On a lossless
                    # replay there is no radio: block instead.
                    if self.lossless:
                        self._put_blocking(self.q_frame, fr)
                    else:
                        drop_oldest(self.q_frame, fr, self.stats,
                                    "frame_dropped")
        except StopIteration:
            ST.log(f"  Frame limit {self.a.frames} reached")
            self.stop.set()
        except Exception as e:                                 # noqa: BLE001
            # E55: live1b died here with a bare one-liner ("TypeError:
            # 'NoneType' object is not subscriptable") and the site had to
            # be reasoned out statically.  Log the full traceback AT the
            # thread, where the stack still exists.
            import traceback
            ST.log(f"  *** front-end thread died ***\n"
                   f"{traceback.format_exc()}")
            self.err = self.err or f"front end: {type(e).__name__}: {e}"
        finally:
            self.q_frame.put(None)

    # -- process decode pool (E52) ----------------------------------------
    def _procs_start(self):
        import multiprocessing as mp
        from multiprocessing import shared_memory
        a = self.a
        n = a.decode_procs
        self.slot_len = ST.FRAME_WINDOW + 4096
        self.nslots = n + 4
        self.shm = shared_memory.SharedMemory(
            create=True, size=self.nslots * self.slot_len * 16)
        self.shm_buf = np.ndarray((self.nslots, self.slot_len),
                                  dtype=np.complex128, buffer=self.shm.buf)
        ctx = mp.get_context("spawn")
        self.d_in = ctx.Queue()
        self.d_out = ctx.Queue()
        cfg = dict(threads=a.threads, iters=a.iters,
                   cpu_fast=(a.accel == "cpu" and not a.exact_cpu))
        self.procs = [ctx.Process(
            target=_decode_proc_main,
            args=(self.shm.name, self.nslots, self.slot_len, self.d_in,
                  self.d_out, cfg), daemon=True) for _ in range(n)]
        t = time.time()
        for p in self.procs:
            p.start()
        ready = 0
        while ready < n:
            kind, payload = self.d_out.get(timeout=180)
            if kind == "err":
                raise RuntimeError(f"decode worker failed to start: {payload}")
            if kind == "ready":
                ready += 1
        ST.log(f"  decode pool: {n} processes x {a.threads} threads, "
               f"{self.nslots} shared {self.slot_len * 16 >> 20} MB frame "
               f"slots, warmed in {time.time() - t:.1f} s")

    def _procs_stop(self):
        if not getattr(self, "procs", None):
            return
        for _ in self.procs:
            try:
                self.d_in.put(None)
            except Exception:                                  # noqa: BLE001
                break
        for p in self.procs:
            p.join(timeout=8)
            if p.is_alive():
                p.terminate()
        try:
            self.shm.close()
            self.shm.unlink()
        except Exception:                                      # noqa: BLE001
            pass
        self.procs = []

    def _t_decode_procs(self):
        """q_frame -> worker processes -> STRICT frame order -> q_bb.

        Dispatch order is recorded in a deque; results park in a dict until
        the head of the deque is available, so hand-over order equals arrival
        order equals the serial decoder's order -- the byte stream cannot
        tell how many processes decoded it.
        """
        free = list(range(self.nslots))
        order = collections.deque()
        parked = {}
        done_tm = collections.Counter()

        def pump(block):
            try:
                kind, payload = (self.d_out.get(timeout=0.5) if block
                                 else self.d_out.get_nowait())
            except queue.Empty:
                return False
            if kind == "err":
                raise RuntimeError(f"decode worker: {payload}")
            if kind == "done":
                done_tm.update(payload)
                return False
            if kind == "ok":
                idx, slot, pkts, diag, dt_, coh = payload
                free.append(slot)
                parked[idx] = (pkts, diag, dt_, coh)
            while order and order[0] in parked:
                self._decode_done(parked.pop(order.popleft()))
            return True

        try:
            while True:
                item = self.q_frame.get()
                if item is None:
                    break
                idx, w, t0, coh = item
                # A set stop flag is NOT permission to shed: on a replay the
                # frame limit sets it while the tail Frames are still queued,
                # and bailing here dropped them (measured: 17 of 19 reached
                # the datagram file).  Wait for a slot; only a worker error
                # (or a worker that stopped answering) abandons the stream.
                wait_until = time.time() + 30
                while not free and not self.err:
                    pump(True)
                    if time.time() > wait_until:
                        self.err = self.err or "decoder: no free slot in 30 s"
                        break
                if not free:
                    break
                slot = free.pop()
                n = len(w)
                self.shm_buf[slot, :n] = w
                self.d_in.put((slot, idx, n, t0, coh))
                order.append(idx)
                while pump(False):
                    pass
            # Drain EVERYTHING in flight, stop flag or not: the frame-limit
            # path (and plain EOF) sets stop while the last Frames are still
            # inside the workers, and a drain that quits on the flag SHEDS
            # them -- the e2e gate caught exactly that as a shorter datagram
            # file.  The deadline covers only a genuinely dead worker.
            deadline = time.time() + 30
            while order and time.time() < deadline and not self.err:
                pump(True)
            # stop the workers and collect their stage timers
            for _ in self.procs:
                self.d_in.put(None)
            got_done = 0
            deadline = time.time() + 8
            while got_done < len(self.procs) and time.time() < deadline:
                try:
                    kind, payload = self.d_out.get(timeout=0.5)
                except queue.Empty:
                    continue
                if kind == "done":
                    done_tm.update(payload)
                    got_done += 1
            self._procs_tm = done_tm
        except Exception as e:                                 # noqa: BLE001
            import traceback
            self.err = self.err or (f"decoder: {type(e).__name__}: {e}\n"
                                    f"{traceback.format_exc()}")
        finally:
            self.q_bb.put(None)

    def t_decode(self):
        if self.ldm:
            return self._t_decode_ldm()
        if getattr(self, "use_procs", False):
            return self._t_decode_procs()
        # E52 -- OVERLAP FRAME DECODES, EMIT IN ORDER.  One decode thread
        # left cores idle whenever a stage's internal parallelism ran out;
        # Frames are independent after the front end (fine timing tracked
        # THERE, not here), so up to --decode-workers Frames decode
        # concurrently.  Results are handed downstream strictly FIFO -- the
        # deque is popped only at its head -- so the Baseband stream is
        # byte-identical to the serial order and everything after this
        # thread is unchanged.  The bad_run/reacquire bookkeeping also runs
        # at hand-over, in order; on a weak live signal a re-acquire can
        # land up to (workers-1) Frames later than the serial chain's, which
        # is inside the noise of an already-heuristic trigger.
        try:
            import concurrent.futures as cf
            nw = max(1, int(getattr(self.a, "decode_workers", 1)))
            dex = cf.ThreadPoolExecutor(max_workers=nw,
                                        thread_name_prefix="decode") \
                if nw > 1 else None
            pending = collections.deque()
            while True:
                item = self.q_frame.get()
                if item is None:
                    break
                if dex is None:
                    self._decode_done(self._decode_job(item))
                    continue
                pending.append(dex.submit(self._decode_job, item))
                while pending and (pending[0].done()
                                   or len(pending) >= nw):
                    self._decode_done(pending.popleft().result())
            while pending:
                self._decode_done(pending.popleft().result())
            if dex is not None:
                dex.shutdown(wait=False)
        except Exception as e:                                 # noqa: BLE001
            import traceback
            self.err = self.err or (f"decoder: {type(e).__name__}: {e}\n"
                                    f"{traceback.format_exc()}")
        finally:
            self.q_bb.put(None)

    # -- E82: the LDM / CTI decode thread ---------------------------------
    def _t_decode_ldm(self):
        """q_frame -> LdmPipeline -> q_bb.

        The RF33 thread hands Frames to workers and re-imposes FIFO order
        because its Frames are INDEPENDENT.  A CTI multiplex has no such
        freedom: a FEC Block is a diagonal across ~two Frames of one
        continuous index space, so the interleaver state is inherently
        serial and the parallelism lives INSIDE the pipeline (the Frame
        demodulator's per-pilot-class batches, the batched demap, and the
        LDPC batch) rather than across Frames.

        Everything this thread emits is `(bytes, bounds)` -- the same pair
        `_decode_done` builds -- so the transport, the lanes, the viewer and
        the datagram dump are the RF33 chain, unmodified.
        """
        try:
            while True:
                item = self.q_frame.get()
                if item is None:
                    break
                idx, w, t0, coh = item
                t = time.perf_counter()
                before = self.pipe.n_blocks
                out = self.pipe.push_frame(idx, w, t0, ps=self.fe.ps)
                dt_ = time.perf_counter() - t
                nb = self.pipe.n_blocks - before
                with self.lock:
                    self.n_frames += 1
                    self.n_conv = self.pipe.n_conv
                    self.n_bch = self.pipe.n_bch
                    self.hist.append((time.time(), dt_, nb, float(coh)))
                for stream, bounds in out:
                    # the pipeline owns the absolute Baseband offset (the ALP
                    # anchors are indexes into the WHOLE stream), so this
                    # thread must NOT advance it a second time
                    self.bb_abs = self.pipe.bb_abs
                    if self.lossless:
                        if not self._put_blocking(self.q_bb,
                                                  (stream, bounds)):
                            end = time.time() + 15
                            while time.time() < end and not self.err:
                                try:
                                    self.q_bb.put((stream, bounds),
                                                  timeout=0.25)
                                    break
                                except queue.Full:
                                    pass
                    else:
                        drop_oldest(self.q_bb, (stream, bounds), self.stats,
                                    "bb_dropped")
                # THE SAME GUARD THE RF33 PATH HAS, AGAINST A DIFFERENT
                # FAILURE.  There a weak Frame means fine timing has drifted;
                # here it can ALSO mean the commutator phase is wrong, and
                # the phase tracker is the thing that decides which -- but a
                # dead front end must still cost a re-acquisition, or a lost
                # bootstrap would look like a lost phase forever.
                if self.pipe.ph is not None and self.pipe.ph.state == "cold":
                    self.bad_run += 1
                    if self.bad_run >= max(self.a.bad_run, 8):
                        ST.log(f"  *** {self.bad_run} Frames with no CTI "
                               f"phase -- re-acquiring the front end ***")
                        self.bad_run = 0
                        self.fe.reacquire()
                        self.pipe.reset()
                else:
                    self.bad_run = 0
            # E85: drain the pool before signalling end of stream, or the tail
            # Frames are decoded inside a worker and then thrown away -- the
            # shed the RF33 e2e gate caught as a short datagram file.
            for stream, bounds in self.pipe.flush():
                self.bb_abs = self.pipe.bb_abs
                if self.lossless:
                    self._put_blocking(self.q_bb, (stream, bounds))
                else:
                    drop_oldest(self.q_bb, (stream, bounds), self.stats,
                                "bb_dropped")
        except Exception as e:                                 # noqa: BLE001
            import traceback
            self.err = self.err or (f"ldm decoder: {type(e).__name__}: {e}\n"
                                    f"{traceback.format_exc()}")
        finally:
            self.q_bb.put(None)

    def _decode_job(self, item):
        idx, w, t0, coh = item
        t = time.perf_counter()
        r16, pkts, diag = self.fd.decode_frame(w, t0)
        return pkts, diag, time.perf_counter() - t, coh

    def _decode_done(self, res):
        pkts, diag, dt_, coh = res
        stream, bounds = bytearray(), []
        for p in pkts:
            if p is None:
                continue
            ptr, pay = P6.bb_split(p)
            if ptr is not None:
                bounds.append(self.bb_abs + len(stream) + ptr)
            stream += pay
        self.bb_abs += len(stream)
        with self.lock:
            self.n_frames += 1
            self.n_conv += diag["converged"]
            self.n_bch += diag["bch_ok"]
            self.hist.append((time.time(), dt_, diag["converged"],
                              float(coh)))
        if diag["converged"] < self.a.min_blocks:
            self.bad_run += 1
            self.stats["weak_frame"] += 1
            # E53: a weak Frame is the timing tracker's REAL guard (the
            # coherence guard measured blind while FEC fell to 0/74) --
            # demand a full fine-timing rescan before waiting out bad_run.
            try:
                self.fe.ft.request_full()
            except AttributeError:
                pass
            if self.bad_run >= self.a.bad_run:
                ST.log(f"  *** {self.bad_run} consecutive Frames under "
                       f"{self.a.min_blocks}/74 FEC Blocks -- "
                       f"re-acquiring ***")
                self.bad_run = 0
                self.fe.reacquire()
        else:
            self.bad_run = 0
        # Same rule one stage down. This put had NO timeout at all, so a
        # stalled transport blocked the decoder outright and the
        # back-pressure reached the reader every time.
        if self.lossless:
            if not self._put_blocking(self.q_bb, (bytes(stream), bounds)):
                # stop is already set (frame limit / EOF) but the transport
                # is still draining towards its None: a put refused HERE
                # sheds an already-decoded Frame's Baseband on a LOSSLESS
                # replay -- measured as datagram files ~2 Frames short.
                # Bounded, error-aware, and only on the lossless path.
                end = time.time() + 15
                while time.time() < end and not self.err:
                    try:
                        self.q_bb.put((bytes(stream), bounds), timeout=0.25)
                        break
                    except queue.Full:
                        pass
        else:
            drop_oldest(self.q_bb, (bytes(stream), bounds),
                        self.stats, "bb_dropped")

    def t_transport(self):
        try:
            while True:
                item = self.q_bb.get()
                if item is None:
                    break
                stream, bounds = item
                for seg in self.tr.feed(stream, bounds):
                    self._emit(seg)
            for seg in self.tr.flush():
                self._emit(seg)
        except Exception as e:                                 # noqa: BLE001
            import traceback
            self.err = self.err or (f"transport: {type(e).__name__}: {e}\n"
                                    f"{traceback.format_exc()}")

    def _emit(self, seg):
        # STEP 1 -- THE FILE COMES FIRST, AND IT IS THE ONLY THING THAT MUST
        # NOT BE SKIPPED.
        #
        # Piping segments straight into a player's stdin couples the two: a
        # player that stalls back-pressures push() -> transport -> q_bb ->
        # decode -> q_frame -> front -> q_raw, and the reader answers a full
        # q_raw by dropping SAMPLES, which costs the bootstrap lock and ~15 s
        # of re-acquisition. A hiccup in the display was taking down the
        # radio. STVT solved this in June by writing live.ts and having the
        # player tail it ~25 s behind; the file is the cushion.
        if self.live is not None:
            # Hand the receive-path counters over on every emit so live.json
            # always carries them. `seq_wild` counting up is the direct,
            # on-air evidence that corrupted mpu_sequence_numbers arrive --
            # the thing E20 could only be argued about after the fact.
            st = dict(self.tr.stats)
            for fl in getattr(self.tr, "mmtp", {}).values():
                for k, v in fl.stats.items():
                    st[k] = st.get(k, 0) + v
            for rs in getattr(self.tr, "routes", {}).values():
                for k, v in rs.stats.items():
                    st["route_" + k] = st.get("route_" + k, 0) + v
            self.live.stats = st
            self.live.write(self.tr.init, seg)
            dump_route_info(self.tr, self.live.dir)
        if isinstance(seg, dict) and seg.get("transport") == "route":
            # A ROUTE segment is another SERVICE.  It goes to its lanes and
            # nowhere else -- pushing it into the MMTP programme's player
            # stdin would interleave two services into one bitstream.
            self.stats["segment"] += 1
            return
        if self.a.player != "none":
            if self.player is None:
                self.player = ST.Player(
                    cmd=player_cmd(self.a.player,
                                   f"ATSC 3.0 -- RF{self.a.rf}",
                                   self.a.player_args),
                    prebuffer=self.a.prebuffer, out_path=self.a.record,
                    quiet=self.a.quiet_player)
                self.player.open(self.tr.init)
                ST.log(f"  player up ({self.a.player}), prebuffering "
                       f"{self.a.prebuffer:.1f} s")
            self.player.push(seg)
        elif self.a.record and self.player is None:
            self.player = ST.Player(cmd=None, prebuffer=0.0,
                                    out_path=self.a.record, quiet=True)
            self.player.open(self.tr.init)
            self.player.push(seg)
        elif self.a.record:
            self.player.push(seg)
        self.stats["segment"] += 1

    # -- E82: identify the multiplex before committing to a decoder -------
    SNIFF_SEC = 0.75

    def _sniff(self):
        """Read a little air, ask L1 what this is, and hand the samples back.

        Returns the source blocks consumed, so `run()` can push them into the
        front end in order.  Nothing here is allowed to have a side effect on
        the legacy path: if L1 does not verify, or verifies and says Hybrid
        time interleaving, this returns with `self.ldm` False and the chain
        that follows is the one that has always run.
        """
        import m44_ldm as M44
        blocks, total = [], 0
        need = int(self.SNIFF_SEC * self.rate)
        t0 = time.time()
        self.src.start()
        while total < need and time.time() - t0 < 30:
            x = self.src.read()
            if x is None:
                break
            # A continuity BREAK is KEPT, not swallowed.  It is replayed into
            # the front end with the samples so the bootstrap tracker learns
            # about it; dropping it here would splice two discontinuous
            # stretches into one apparently-continuous stream, which is
            # exactly what t0 tracking cannot survive.
            blocks.append(x)
            if isinstance(x, bytes):
                continue
            total += len(x)
        samples = [b for b in blocks if not isinstance(b, bytes)]
        if not samples:
            return blocks
        # identify on the LAST unbroken stretch, so a break inside the sniff
        # window cannot fabricate a multiplex out of two halves
        last = []
        for b in blocks:
            if isinstance(b, bytes):
                last = []
            else:
                last.append(b)
        plan, info = M44.sniff(np.concatenate(last or samples), self.rate,
                               log=ST.log)
        if plan is None:
            ST.log(f"  multiplex sniff: {info.get('why', 'unidentified')} "
                   f"-- taking the m9_fast (RF33-shaped) path")
            return blocks
        ST.log("  multiplex sniff (from L1, not from the channel number) --\n"
               "    " + plan.describe())
        if not plan.uses_cti:
            ST.log("  L1 signals the HYBRID time interleaver -- m9_fast owns "
                   "this multiplex; the LDM path stands down")
            return blocks
        if getattr(self.a, "ldm", "auto") == "off":
            return blocks
        self.ldm = True
        self.plan = plan
        self.frame_sec = plan.frame_sec
        self.blocks_per_frame = plan.plp_size / plan.ncell
        return blocks

    # -- run ---------------------------------------------------------------
    def build(self):
        a = self.a
        import m9_fast
        from concurrent.futures import ThreadPoolExecutor
        self.ex = ThreadPoolExecutor(max_workers=a.fe_threads)
        if a.capture:
            self.src = ST.FileSource(a.capture, a.rate, fmt=a.fmt,
                                     block=a.block, realtime=a.realtime,
                                     start=a.start)
            rate = a.rate
        else:
            self.src = ST.SdrSource(a.rf, rate=a.rate, ant=a.ant,
                                    ifgr=a.ifgr, rfgain=a.rfgain,
                                    wait=a.wait, block=a.block,
                                    radio_lock=radio_lock)
            rb = self.src.open()
            ST.log(f"  radio readback: {json.dumps(rb, default=str)}")
            if not rb["rate_exact"]:
                ST.log(f"  *** requested {a.rate:.0f} Hz, radio gave "
                       f"{rb['rate']:.0f} Hz -- the resampler stays in the "
                       f"pipeline ***")
            rate = rb["rate"]
        self.rate = rate
        # ---- E82: WHICH MULTIPLEX IS THIS?  Ask L1, not the channel number.
        #
        # The decision below reads exactly two things off the decoded
        # L1-Detail -- is there a Layer-1 PLP (LDM), and does the Core PLP
        # signal TI mode 1 (the Convolutional Time Interleaver) -- because
        # those are the two facts m9_fast cannot handle.  RF25 and RF30 both
        # answer yes and neither is named anywhere in this file.
        #
        # The sniff consumes source blocks and REPLAYS them into the front
        # end afterwards, in order, so a multiplex that takes the legacy path
        # sees byte-for-byte the same sample stream it always did.  A radio
        # cannot be rewound; this is the version that does not need to be.
        sniffed = []
        if getattr(a, "ldm", "auto") != "off":
            # THIS RUNS ON THE USER'S TELEVISION.  Identification is a
            # convenience; the legacy path is the product.  Any failure in
            # here -- a malformed L1, an unfamiliar Table H.1.1 row, a bug of
            # mine -- must cost a log line and nothing else, so the whole
            # thing is inside a net.  (The first build let a Reserved
            # CTI_depth on RF33 raise straight out of build(): the chain
            # refused to start at all, on the channel it was supposed to
            # leave untouched.)
            try:
                sniffed = self._sniff()
            except Exception as e:                             # noqa: BLE001
                import traceback
                self.ldm = False
                self.plan = None
                ST.log(f"  multiplex sniff FAILED ({type(e).__name__}: {e}) "
                       f"-- taking the m9_fast path\n"
                       f"{traceback.format_exc()}")
        self.fe = ST.FrontEnd(rate, ex=self.ex,
                              fast=(a.accel == "cpu" and not a.exact_cpu),
                              plan=self.plan if self.ldm else None)
        if not self.fe.rs.bypass:
            ST.log(f"  resampler ACTIVE: {rate/1e6:g} -> "
                   f"{ST.FS_POST/1e6:g} Msps (up={self.fe.rs.up} "
                   f"down={self.fe.rs.down}), phase carried across blocks")
        else:
            ST.log(f"  resampler BYPASSED: capturing at "
                   f"{ST.FS_POST/1e6:g} Msps directly")
        self.use_procs = (a.decode_procs > 0 and a.accel == "cpu"
                          and not self.ldm)
        if a.accel == "cpu":
            ST.log(f"  cpu mode: {'EXACT float64' if a.exact_cpu else 'fast float32 (gated against the exact path; --exact-cpu restores it)'}, "
                   + (f"{a.decode_procs} decode process(es)" if self.use_procs
                      else f"{a.decode_workers} decode worker thread(s)"))
            if not a.exact_cpu:
                import m16_margin as M16
                mg = M16.MarginCfg()
                ST.log(f"  E60 margin levers: {json.dumps(mg.asdict())}"
                       if mg.active() else
                       "  E60 margin levers: OFF (ATSC3_MARGIN=0)")
        if self.ldm:
            import m44_ldm as M44
            self.fd = None
            self._procs_tm = collections.Counter()
            ncpu = os.cpu_count() or 1
            nproc = a.ldm_procs
            if nproc is None:
                # MEASURED, not guessed (E85, rf25_fox_g4, 4 threads each):
                #   0 procs  0.4845x    4 procs  1.4902x
                #   6 procs  2.0142x    8 procs  2.1931x
                # Past 6 the parent's own serial work (front end + the CTI
                # gather) is the limit, so the extra processes buy little and
                # cost cores the audio workers and the viewer also want.
                nproc = (6 if ncpu >= 32 else 4 if ncpu >= 16
                         else 2 if ncpu >= 8 else 0)
            self.pipe = M44.LdmPipeline(
                plan=self.plan, iters=a.iters, threads=a.threads, ex=self.ex,
                accel=("gpu" if a.accel in ("gpu", "gpu-full") else "cpu"),
                log=ST.log, procs=nproc,
                proc_threads=a.ldm_proc_threads,
                fec_threads=a.ldm_fec_threads or a.threads)
            self.pipe.adopt(self.plan)
            t_warm = time.time()
            self.pipe.prewarm()
            ST.log(f"  LDM decoder prewarmed in {time.time() - t_warm:.2f} s "
                   f"(189 frequency-interleaver sequences and the 64K LDPC "
                   f"tables, built before the clock)")
        elif self.use_procs:
            self.fd = None
            self._procs_tm = collections.Counter()
            self._procs_start()
        else:
            self.fd = m9_fast.FrameDecoder(threads=a.threads, backend=a.accel,
                                           iters=a.iters,
                                           cpu_fast=(a.accel == "cpu"
                                                     and not a.exact_cpu))
            t_warm = time.time()
            self.fd.prewarm()
            ST.log(f"  decoder prewarmed in {time.time() - t_warm:.2f} s "
                   f"(tables built before the clock, not under the first "
                   f"Frame)")
        if a.dump_dg:
            self.dg_fh = open(a.dump_dg, "wb")
        want = {"video": (b"vide",),
                "av": (b"vide", b"soun"),
                "all": (b"vide", b"soun", b"subt", b"text")}.get(
                    getattr(a, "assets", "video"), (b"vide",))
        self.tr = ST.Transport(probe=a.flow_probe, pid=a.pid, want=want,
                               dg_fh=self.dg_fh, repair=not a.no_repair,
                               route=([a.route] if getattr(a, "route", None)
                                      else None))
        self._sniffed = sniffed

    def run(self):
        if getattr(self.a, "cpu_isolate", None):
            cores = apply_cpu_isolation(self.a.cpu_isolate)
            if cores:
                ST.log(f"  cpu isolation: pinned to {len(cores)} cores {cores}")
        if getattr(self.a, "live_dir", None):
            self.live = LiveWriter(self.a.live_dir)
            ST.log(f"  live file: {self.live.path}  "
                   f"(playback should tail this, not the decoder)")
        a = self.a
        self.build()
        self.t_start = time.time()
        ts = [threading.Thread(target=f, daemon=True, name=n) for f, n in (
            (self.t_reader, "reader"), (self.t_front, "front"),
            (self.t_decode, "decode"), (self.t_transport, "transport"))]
        for t in ts:
            t.start()
        deadline = self.t_start + a.secs if a.secs else None
        last = 0.0
        try:
            while any(t.is_alive() for t in ts[:3]):
                time.sleep(0.25)
                if self.err:
                    break
                now = time.time()
                if deadline and now >= deadline:
                    ST.log(f"  run length {a.secs:.0f} s reached")
                    break
                if self.stop.is_set():
                    break
                if getattr(self.src, "yield_reason", None):
                    ST.log(f"  YIELDING the radio: {self.src.yield_reason}")
                    break
                if now - last >= a.report:
                    last = now
                    self.report(now)
        except KeyboardInterrupt:
            ST.log("  interrupted")
        self.stop.set()
        try:
            self.src.close()
        except Exception:                                      # noqa: BLE001
            pass
        for t in ts:
            t.join(timeout=8)
        if self.player is not None:
            self.player.close()
        if self.live is not None:
            self.live.close()
        if self.dg_fh is not None:
            self.dg_fh.close()
        if self.fd is not None:
            self.fd.close()
        if self.pipe is not None:
            self.pipe.close()
        if getattr(self, "use_procs", False):
            self._procs_stop()
        if self.ex is not None:
            self.ex.shutdown(wait=False)
        return self.summary()

    # -- telemetry ---------------------------------------------------------
    def report(self, now=None):
        now = now or time.time()
        wall = now - self.t_start
        with self.lock:
            nf = self.n_frames
            h = list(self.hist)
        air = nf * self.frame_sec
        inst = ""
        if len(h) >= 2:
            span = h[-1][0] - h[0][0]
            if span > 0.5:
                inst = f"  inst {len(h) * self.frame_sec / span:5.2f}x"
        p = self.player
        occ = p.occupancy(now) if p else 0.0
        under = p.stats["underrun"] if p else 0
        # a repaired fragment is a fragment the picture SHOWS but does not show
        # whole, so it belongs on the running line, not only in the summary
        rep = self.tr.stats["segment_truncated"]
        rep = f" rep {rep}" if rep else ""
        # INTERVAL FEC, beside the cumulative one. The cumulative figure is a
        # lifetime average: after a fade it stays depressed for hours while
        # the signal is already perfect again, and reading it as "now" has
        # produced a false alarm about signal degradation twice (E16, and
        # again on 8/07 before catching it). The last ~64 Frames are what the
        # antenna is doing at this moment.
        recent = ""
        bpf = self.blocks_per_frame
        if self.ldm:
            tot_blocks = max(self.pipe.n_blocks, 1)
            recent = f" [phase {self.pipe.ph.state}]"
            ST.log(f"  {wall:6.1f}s  {nf:5d} Frames  "
                   f"{air/max(wall,1e-9):5.2f}x rt{inst}  "
                   f"FEC {self.n_bch}/{self.pipe.n_blocks} BCH0{recent}  "
                   f"q raw/frm/bb {self.q_raw.qsize()}/{self.q_frame.qsize()}/"
                   f"{self.q_bb.qsize()}  "
                   f"MPU {self.tr.stats['segment']} seg "
                   f"buf {occ:5.2f}s  underruns {under}")
            return
        if len(h) >= 8:
            got = sum(x[2] for x in h)
            recent = f" [{100.0 * got / (bpf * len(h)):5.1f}% now]"
        ST.log(f"  {wall:6.1f}s  {nf:5d} Frames  {air/max(wall,1e-9):5.2f}x rt"
               f"{inst}  FEC {self.n_conv}/{int(bpf)*max(nf,1)}{recent}  "
               f"q raw/frm/bb {self.q_raw.qsize()}/{self.q_frame.qsize()}/"
               f"{self.q_bb.qsize()}  "
               f"MPU {self.tr.stats['segment']} seg{rep} "
               f"buf {occ:5.2f}s  underruns {under}")

    def summary(self):
        wall = time.time() - self.t_start
        air = self.n_frames * self.frame_sec
        p = self.player
        mm = {}
        for k, v in self.tr.mmtp.items():
            for kk, vv in v.stats.items():
                mm[kk] = mm.get(kk, 0) + vv
        s = dict(
            wall_s=round(wall, 3), frames=self.n_frames,
            air_s=round(air, 3),
            x_real_time=round(air / wall, 4) if wall else 0.0,
            # ONE SOURCE OF TRUTH.  The LDM counters are the pipeline's;
            # this thread only mirrors them for the running line, and the
            # mirror stops when the decode loop exits -- so mixing a mirrored
            # numerator with the pipeline's denominator reported 7172/7662
            # (93.6%) for a run the pipeline itself scored 7662/7662.  A
            # number assembled from two clocks is a wrong number.
            fec_converged=(self.pipe.n_conv if self.ldm else self.n_conv),
            fec_total=(self.pipe.n_blocks if self.ldm
                       else 74 * self.n_frames),
            bch_zero=(self.pipe.n_bch if self.ldm else self.n_bch),
            front_end=dict(self.fe.stats), source=dict(self.src.stats),
            alp=dict(self.tr.walker.stats), ip=dict(self.tr.ip.stats),
            transport=dict(self.tr.stats), mmtp=mm,
            pipeline=dict(self.stats),
            rate=self.rate, resampler=not self.fe.rs.bypass,
            readback=getattr(self.src, "readback", None))
        nf = max(self.n_frames, 1)
        s["front_end_ms_per_frame"] = {
            k: round(1000 * v / nf, 2) for k, v in self.fe.tm.items()}
        if self.ldm:
            s["ldm"] = self.pipe.summary()
            dtm = collections.Counter(
                {k: v for k, v in self.pipe.tm.items()})
            dtm.update(self.pipe.fd.tm)
        else:
            dtm = (self._procs_tm if getattr(self, "use_procs", False)
                   else self.fd.tm)
        # E60: "n_*" entries are COUNTERS (weighted frames, SP rescues, CE
        # fallback symbols), not stage seconds -- report them as counts.
        s["decode_ms_per_frame"] = {
            k: round(1000 * v / nf, 2) for k, v in dtm.items()
            if not k.startswith("n_")}
        s["margin_counters"] = {k: int(v) for k, v in dtm.items()
                                if k.startswith("n_")}
        s["frame_budget_ms"] = round(1000 * self.frame_sec, 2)
        if p is not None:
            s["player"] = dict(
                segments=p.stats["segment"] + 1, media_s=round(p.media_s, 3),
                bytes=p.bytes, underruns=p.stats["underrun"],
                underrun_s=round(p.stats["underrun_s"], 3),
                min_buffer_s=(None if p.min_occ is None
                              else round(p.min_occ, 3)),
                prebuffer_s=p.prebuffer, player_died=p.dead)
        if self.err:
            s["error"] = self.err
        return s


def dump_route_info(tr, live_dir):
    """live_route.json: the ROUTE service description, written when the SLS
    (re)parses and never on the per-segment hot path (a version check is one
    dict lookup).  The viewer reads it to bind lanes to Representations."""
    routes = getattr(tr, "routes", {})
    if not routes:
        return
    seen = getattr(tr, "_route_info_seen", None)
    if seen is None:
        seen = tr._route_info_seen = {}
    for key, rs in routes.items():
        info = rs.route_info()
        if info is None or seen.get(key) == info["ver"]:
            continue
        seen[key] = info["ver"]
        p = os.path.join(live_dir, "live_route.json")
        try:
            tmp = p + ".tmp"
            with open(tmp, "w") as fh:
                json.dump(info, fh, indent=1)
            os.replace(tmp, p)
            ST.log(f"  live_route.json: serviceId {info['service_id']}, "
                   f"{len(info['reps'])} reps")
        except OSError:
            pass


def drop_oldest(q, item, stats, key):
    """Put without ever blocking the producer.

    A bounded queue offers two ways to handle a full queue and only one of
    them is safe here. Blocking propagates the stall upstream until it
    reaches the SDR reader, which can only answer by dropping samples --
    fatal, because sample continuity is what the bootstrap lock is made of.
    Shedding the OLDEST item costs one Frame of picture and stops the
    pressure dead. Newest-first also matters: the freshest Frame is the one
    the viewer is closest to needing.
    """
    while True:
        try:
            q.put_nowait(item)
            return
        except queue.Full:
            try:
                q.get_nowait()
                stats[key] += 1
            except queue.Empty:          # drained by the consumer; retry
                pass


class LiveWriter:
    """Append the fMP4 byte stream to a growing file, and nothing else.

    This is STVT's live.ts cushion, ported. The chain's only obligation is
    to keep writing; whoever watches tails the file however far behind it
    likes. A file append is microseconds and cannot back-pressure, which is
    exactly the property a player's stdin does not have.

    Also writes a small JSON heartbeat next to it, because a supervisor
    needs a cheap and HONEST liveness signal -- "the process exists" is not
    one (a livelocked chain kept a core at 98 % while emitting nothing for
    two minutes). Bytes written and their timestamp cannot be faked by a
    wedged thread.
    """

    NAME = {b"vide": "video", b"soun": "audio", b"subt": "subs",
            b"text": "subs"}

    def __init__(self, d, keep_mb=4096):
        self.dir = d
        os.makedirs(d, exist_ok=True)
        self.meta = os.path.join(d, "live.json")
        self.lanes = {}                     # pid -> lane dict
        self.stats = {}                     # receive-path counters, live
        self.t0 = time.time()
        self._last_meta = 0.0

    # `path` still names the video lane so single-asset callers read the same.
    # kind "video" (the MMTP programme) beats "route_video" (E47): when both
    # transports write lanes, the top-level path/segments fields must keep
    # describing the lane every existing consumer means by "the video".
    @property
    def path(self):
        alt = None
        for ln in self.lanes.values():
            if ln["handler"] == b"vide":
                if ln["kind"] == "video":
                    return ln["path"]
                alt = alt or ln["path"]
        return alt or os.path.join(self.dir, "live_video_pid0.m4s")

    @property
    def bytes(self):
        return sum(ln["bytes"] for ln in self.lanes.values())

    def _lane(self, pid, handler, init, kind=None):
        ln = self.lanes.get(pid)
        if ln is None:
            # The pid goes in the FILENAME, always. This multiplex carries
            # TWO `soun` assets (a second programme audio), and naming lanes
            # by handler alone made both open "live_audio.m4s" -- the second
            # open truncated the first and the two interleaved into garbage.
            # Lock order is not deterministic either, so which one won would
            # change between runs. live.json carries handler and pid so a
            # consumer can choose deliberately instead of by filename luck.
            # E47: a segment may carry its own lane kind ("route_video",
            # "route_audio", ...) so a second SERVICE's lanes can never be
            # mistaken for the MMTP programme's -- every existing consumer
            # matches kind by exact equality.
            kind = kind or self.NAME.get(handler, "data")
            name = f"{kind}_pid{pid}"
            path = os.path.join(self.dir, f"live_{name}.m4s")
            # ROLL, DO NOT TRUNCATE.
            #
            # The supervisor restarted a genuinely wedged chain on 8/07 and
            # the cure ate the patient: this open("wb") discarded 433.8 MB of
            # already-recorded television, and the .idx with it. Appending
            # instead would be WORSE -- a restarted chain re-acquires and
            # begins a fresh timeline (fragment 1, base_time 0, a new init
            # box), so "ab" splices two incompatible timelines into one file
            # that no player can read correctly.
            #
            # So neither: move the old lane aside and start clean. The
            # recording survives as live_video_pid12.0001.m4s, each file holds
            # exactly one continuous timeline, and the roll is VISIBLE to
            # downstream workers through `generation` in live.json rather than
            # being something they have to infer from a file that shrank.
            gen = 0
            if os.path.exists(path) and os.path.getsize(path) > 0:
                while True:
                    gen += 1
                    old = f"{os.path.splitext(path)[0]}.{gen:04d}.m4s"
                    if not os.path.exists(old):
                        break
                for src, dst in ((path, old),
                                 (os.path.splitext(path)[0] + ".idx",
                                  os.path.splitext(old)[0] + ".idx")):
                    try:
                        os.replace(src, dst)
                    except OSError:
                        pass
                ST.log(f"  rolled previous {kind} lane -> "
                       f"{os.path.basename(old)} (generation {gen})")
            fh = open(path, "wb", buffering=1 << 20)
            if init:
                fh.write(init)
            # A per-fragment index of MPU sequence numbers. Without it the
            # lanes cannot be ALIGNED, only concatenated: fragment N of the
            # audio lane is not necessarily the same 2.002 s slot as
            # fragment N of the video lane, because a lane that loses an MPU
            # has one fewer fragment from that point on. Concatenating then
            # SHIFTS everything after the hole instead of leaving a gap --
            # the same mistake m39 exists to avoid on the batch side.
            idx = open(os.path.splitext(path)[0] + ".idx", "w", buffering=1)
            ln = self.lanes[pid] = dict(fh=fh, path=path, handler=handler,
                                        bytes=len(init or b""), segments=0,
                                        media=0.0, name=name, kind=kind,
                                        pid=pid, idx=idx, first_seq=None,
                                        last_seq=None, init_bytes=len(init or b""),
                                        generation=gen)
            ST.log(f"  live lane: {kind} pid {pid} -> {os.path.basename(path)}")
        return ln

    def write(self, init, seg):
        # A segment from Transport is a dict -- it carries its own pid,
        # handler and init since Transport went multi-asset, so each asset
        # lands in its OWN file. One interleaved file would need a real
        # muxer and would couple the lanes; separate files let the audio
        # decoder (1.28x real time) run behind the video without either
        # waiting for the other -- the same decoupling, one level up.
        if not isinstance(seg, dict):
            return
        pid = seg.get("pid", 0)
        ln = self._lane(pid, seg.get("handler"), seg.get("init") or init,
                        seg.get("lane_kind"))
        body = seg["bytes"]
        off = ln["bytes"]
        ln["fh"].write(body)
        ln["fh"].flush()
        seq = seg.get("seq")
        ln["idx"].write(json.dumps(dict(seq=seq, off=off, len=len(body),
                                        dur=seg.get("dur", 0))) + "\n")
        if ln["first_seq"] is None:
            ln["first_seq"] = seq
        ln["last_seq"] = seq
        ln["bytes"] += len(body)
        ln["segments"] += 1
        ln["media"] += seg.get("dur", 0) / 90000.0
        now = time.time()
        if now - self._last_meta >= 1.0:
            self._last_meta = now
            self._dump_meta(now)

    def _dump_meta(self, now):
        lanes = {ln["name"]: dict(
                    bytes=ln["bytes"], segments=ln["segments"],
                    media_s=ln["media"], path=ln["path"],
                    kind=ln["kind"], pid=ln["pid"],
                    first_seq=ln["first_seq"], last_seq=ln["last_seq"],
                    init_bytes=ln["init_bytes"],
                    idx=os.path.splitext(ln["path"])[0] + ".idx",
                    generation=ln.get("generation", 0),
                    handler=(ln["handler"] or b"?").decode("ascii", "replace"))
                 for ln in self.lanes.values()}
        vid = next((l for l in self.lanes.values()
                    if l["kind"] == "video"), None) \
            or next((l for l in self.lanes.values()
                     if l["handler"] == b"vide"), None)
        try:
            tmp = self.meta + ".tmp"
            with open(tmp, "w") as fh:
                json.dump(dict(
                    bytes=self.bytes, updated=now, started=self.t0,
                    segments=vid["segments"] if vid else 0,
                    media_s=(vid["media"] if vid else 0.0),
                    path=self.path, lanes=lanes,
                    # Receive-path counters, CONTINUOUSLY. E20's lanes died
                    # of a corrupted MPU sequence number and the counters
                    # that would have named it only existed in the exit
                    # summary -- which a force-kill then threw away. A
                    # diagnostic you can read only after stopping the run
                    # is a diagnostic you do not have.
                    stats=dict(self.stats) if self.stats else {}), fh)
            os.replace(tmp, self.meta)
        except OSError:
            pass

    def close(self):
        # one final heartbeat so live.json describes the CLOSED state -- the
        # 1 s throttle otherwise leaves the last second of segments uncounted
        # (harmless live, wrong for anything that reads the dir afterwards)
        if self.lanes:
            self._dump_meta(time.time())
        for ln in self.lanes.values():
            for k in ("fh", "idx"):
                try:
                    ln[k].close()
                except (OSError, AttributeError):
                    pass


def apply_cpu_isolation(spec):
    """Pin this process to a core set so it is never preempted.

    stvt_run.sh:84 states the mechanism plainly: "preemption -> pipeline
    stall -> SDR overflow -> drought". Same chain here. taskset is
    Linux-only, so go through psutil, which works on Windows too.

    `spec` is a core list like "0-7" or "0,2,4"; empty disables.
    """
    if not spec:
        return None
    cores = []
    for part in str(spec).split(","):
        part = part.strip()
        if "-" in part:
            a, b = part.split("-", 1)
            cores.extend(range(int(a), int(b) + 1))
        elif part:
            cores.append(int(part))
    if not cores:
        return None
    try:
        import psutil
        p = psutil.Process()
        n = psutil.cpu_count()
        cores = sorted({c for c in cores if 0 <= c < n})
        p.cpu_affinity(cores)
        try:                      # best effort; not fatal if refused
            p.nice(psutil.HIGH_PRIORITY_CLASS if os.name == "nt" else -5)
        except (psutil.AccessDenied, PermissionError, OSError):
            pass
        return cores
    except Exception as e:                                     # noqa: BLE001
        ST.log(f"  cpu isolation unavailable ({type(e).__name__}: {e})")
        return None


def main(argv=None):
    ap = argparse.ArgumentParser(
        prog="atsc3 watch",
        description="Tune an ATSC 3.0 multiplex and play it live.")
    ap.add_argument("--rf", type=int, required=True, help="RF channel number (required)")
    ap.add_argument("--rate", type=float, default=ST.FS_POST,
                    help="capture rate; 6912000 makes the resampler vanish")
    ap.add_argument("--ant", default="Antenna B")
    ap.add_argument("--ifgr", type=int, default=None,
                    help="IFGR override. Default None = take it from the carrier gain table (E79).")
    ap.add_argument("--rfgain", type=int, default=None,
                    help="rfgain_sel override. Default None = look the "
                         "carrier up in lab/carrier_gain.json (E79): RF33 "
                         "needs 4, RF25 wants 2, and one global value is "
                         "wrong for one of them whichever it is.")
    ap.add_argument("--wait", type=float, default=120.0,
                    help="polite wait for the radio lock, seconds")
    ap.add_argument("--secs", type=float, default=0.0,
                    help="stop after N seconds (0 = until interrupted)")
    ap.add_argument("--frames", type=int, default=0,
                    help="stop after N Frames (0 = no limit).  Used by the "
                         "gate so both paths cover exactly the same air.")
    ap.add_argument("--capture", default=None,
                    help="replay a banked capture instead of opening the radio")
    ap.add_argument("--fmt", default=None)
    ap.add_argument("--start", type=float, default=0.0)
    ap.add_argument("--realtime", action="store_true",
                    help="throttle a replay to real time")
    ap.add_argument("--shed", action="store_true",
                    help="shed frames on a capture replay the way live air "
                         "must (default: an unthrottled replay applies "
                         "backpressure instead -- lossless, so slow machines "
                         "decode everything late rather than little on time)")
    ap.add_argument("--assets", default="video",
                    help="which assets to carry: video | av | all "
                         "(av = video+audio, all = +captions). Live was "
                         "video-only until 8/07.")
    ap.add_argument("--live-dir", default=None,
                    help="write the fMP4 stream to DIR/live.m4s continuously "
                         "so playback can tail it instead of being piped "
                         "straight from the decoder (see LiveWriter)")
    ap.add_argument("--route", default=None, metavar="IP:PORT",
                    help="ALSO carry this ROUTE/LCT flow as DASH lanes "
                         "(e.g. 239.255.32.1:8321 = serviceId 1, 132.1). "
                         "E47: STAGED -- the flow must ride a PLP this "
                         "decoder actually decodes; 132.1 rides PLP1 in "
                         "SUBFRAME 1, which the live decoder does NOT yet "
                         "decode (m9_fast is subframe-0 only), so on air "
                         "this currently yields nothing. Explicit flow "
                         "only: the Widevine-locked services stay untouched.")
    ap.add_argument("--cpu-isolate", default=None,
                    help="pin the chain to these cores, e.g. 0-7 or 0,2,4 -- "
                         "preemption stalls the pipeline and overflows the SDR")
    ap.add_argument("--player", default="ffplay",
                    choices=("ffplay", "mpv", "none"))
    ap.add_argument("--player-args", default=None)
    ap.add_argument("--quiet-player", action="store_true")
    ap.add_argument("--prebuffer", type=float, default=5.0,
                    help="seconds of media held before playback starts.  "
                         "Media arrives in indivisible 2.002 s MPUs, so a "
                         "buffer smaller than one delivery unit cannot absorb "
                         "even one late delivery -- measured: at 2.0 s the "
                         "occupancy sawtooths through zero and every lost MPU "
                         "is an underrun.  5.0 s holds three MPUs.")
    ap.add_argument("--record", default=None,
                    help="also write the played stream to this .mp4")
    ap.add_argument("--no-repair", action="store_true",
                    help="drop a holed MPU whole instead of playing the part "
                         "of it that arrived (the M11 behaviour; the control "
                         "for m12_repair_gate.py)")
    ap.add_argument("--dump-dg", default=None,
                    help="write decoded IP datagrams (m7 .dg format) -- the "
                         "gate artefact")
    ap.add_argument("--pid", type=int, default=None,
                    help="force an MMTP packet_id instead of picking the "
                         "asset whose moov handler is 'vide'")
    ap.add_argument("--ldm-procs", type=int, default=None,
                    help="E85: demodulate LDM Frames in N worker PROCESSES.  "
                         "The Frame demod is per-Frame independent and was "
                         "230 ms of a 242 ms budget on RF25 (189 OFDM symbols "
                         "against RF33's 35); the CTI and the commutator "
                         "phase are one continuous index space and stay in "
                         "the parent, which re-imposes Frame order before "
                         "anything reaches the interleaver.  Default: 4 on a "
                         "24+ thread box, 2 on 8+, else 0.  0 = serial.")
    ap.add_argument("--ldm-proc-threads", type=int, default=4)
    ap.add_argument("--ldm-fec-threads", type=int, default=None,
                    help="threads the PARENT gives the FEC once the demod has "
                         "left it.  Default: --threads.")
    ap.add_argument("--ldm", default="auto", choices=("auto", "off"),
                    help="E82: identify the multiplex from L1 at tune time "
                         "and, if it signals LDM + the Convolutional Time "
                         "Interleaver (RF25, RF30), decode its CORE layer "
                         "through m44_ldm instead of m9_fast.  `off` forces "
                         "the m9_fast path and skips the sniff entirely.  "
                         "There is no channel number in this decision.")
    ap.add_argument("--accel", default="gpu-full",
                    choices=("cpu", "gpu", "gpu-full"))
    ap.add_argument("--exact-cpu", action="store_true",
                    help="run --accel cpu in the HEAD-exact float64 path "
                         "instead of the gated float32 fast path (E52)")
    ap.add_argument("--no-margin", action="store_true",
                    help="disable the E60 margin levers (smoothed CE, "
                         "auto-weighted LLRs, exact-BP rescue) -- the chain "
                         "runs the pre-E60 code paths untouched.  Individual "
                         "levers: ATSC3_CE_W / ATSC3_WLLR / ATSC3_SP env")
    ap.add_argument("--decode-procs", type=int, default=None,
                    help="decode Frames in N worker PROCESSES (shared-memory "
                         "windows, strict-order hand-over).  The GIL caps "
                         "thread workers at ~1.2x on ANY core count; "
                         "processes actually scale.  Default: 2 when --accel "
                         "cpu on a machine with 8+ hardware threads, else 0. "
                         "0 forces the in-process thread path")
    ap.add_argument("--decode-workers", type=int, default=2,
                    help="Frames decoded concurrently; results are emitted "
                         "strictly in order, so downstream bytes are "
                         "identical to --decode-workers 1")
    ap.add_argument("--threads", type=int, default=16)
    ap.add_argument("--fe-threads", type=int, default=16)
    ap.add_argument("--iters", type=int, default=50)
    ap.add_argument("--block", type=int, default=854_000,
                    help="source block size in samples.  854000 is half a "
                         "Frame (123.6 ms) and measured fastest: smaller "
                         "blocks pay thread dispatch on every stage, larger "
                         "ones buy nothing and cost latency.")
    ap.add_argument("--raw-queue", type=int, default=24,
                    help="source blocks buffered ahead of the "
                         "front end; 24 x 123.6 ms = 3.0 s")
    ap.add_argument("--frame-queue", type=int, default=6)
    ap.add_argument("--bb-queue", type=int, default=32)
    ap.add_argument("--flow-probe", type=int, default=64)
    ap.add_argument("--min-blocks", type=int, default=60,
                    help="FEC Blocks per Frame below which a Frame is 'weak'")
    ap.add_argument("--bad-run", type=int, default=8,
                    help="consecutive weak Frames before re-acquiring")
    ap.add_argument("--report", type=float, default=5.0)
    ap.add_argument("--json", default=None, help="write the summary here")
    ap.add_argument("--force-meteor", action="store_true",
                    help=argparse.SUPPRESS)
    a = ap.parse_args(argv)

    if a.no_margin:
        # before any worker spawns: the config rides the environment
        os.environ["ATSC3_MARGIN"] = "0"
    if a.decode_procs is None:
        ncpu = os.cpu_count() or 1
        a.decode_procs = (0 if a.accel != "cpu" or ncpu < 8
                          else 4 if ncpu >= 24 else 2)
    if a.capture is None and in_meteor_window() and not a.force_meteor:
        ST.log("REFUSING: inside the 02:10-05:40 meteor window.")
        return 3
    if a.capture and not os.path.isabs(a.capture):
        a.capture = os.path.join(HERE, a.capture)

    ST.log(f"ATSC 3.0 live -- RF{a.rf}, accel={a.accel}, "
           f"{'replay ' + os.path.basename(a.capture) if a.capture else 'RADIO'}")
    w = Watch(a)
    if w.lossless:
        ST.log("  replay is LOSSLESS (backpressure, no frame shedding); "
               "--shed or --realtime restores live behaviour")

    def onsig(*_):
        w.stop.set()
    try:
        signal.signal(signal.SIGINT, onsig)
        signal.signal(signal.SIGTERM, onsig)
    except Exception:                                          # noqa: BLE001
        pass

    try:
        s = w.run()
    finally:
        try:
            w.src.close()
        except Exception:                                      # noqa: BLE001
            pass

    print("\n  === RUN SUMMARY ===")
    print(f"    air decoded        {s['air_s']:9.3f} s over "
          f"{s['wall_s']:.3f} s of wall")
    print(f"    SUSTAINED          {s['x_real_time']:9.4f}x REAL TIME  "
          f"({s['frames']} Frames)")
    print(f"    FEC Blocks         {s['fec_converged']}/{s['fec_total']} "
          f"converged, {s['bch_zero']} BCH zero")
    if "ldm" in s:
        L = s["ldm"]
        print(f"    LDM/CTI            phase {L['state']}, start_row "
              f"{L['start_row']}, C {L['C']}, origin Frame "
              f"{L['origin_frame']}  {L['phase']}  {L['stats']}")
    print(f"    ALP resyncs        {s['alp'].get('resync', 0)}")
    print(f"    IP datagrams       {s['ip'].get('udp', 0)}")
    print(f"    MMTP               {s['mmtp'].get('complete', 0)} MPUs "
          f"complete, {s['mmtp'].get('lost', 0)} lost, "
          f"{s['mmtp'].get('psn_break', 0)} sequence breaks")
    print(f"    front end          {s['front_end'].get('acquire', 0)} "
          f"acquisitions, {s['front_end'].get('reacquire', 0)} re-acquisitions")
    print(f"    SDR                {s['source'].get('overflow_reads', 0)} "
          f"overflow reads")
    if "player" in s:
        p = s["player"]
        print(f"    player             {p['segments']} segments, "
              f"{p['media_s']:.3f} s of media, {p['bytes']/1e6:.2f} MB")
        print(f"    *** UNDERRUNS      {p['underruns']}  "
              f"(min buffer {p['min_buffer_s']} s, prebuffer "
              f"{p['prebuffer_s']} s) ***")
    else:
        print("    player             NEVER STARTED -- no complete video MPU")
    b = s["frame_budget_ms"]
    print(f"\n    where the time goes, ms per Frame (budget {b} ms):")
    tots = []
    for tag, d in (("front end", s["front_end_ms_per_frame"]),
                   ("decode", s["decode_ms_per_frame"])):
        tot = sum(d.values())
        tots.append(tot)
        print(f"      {tag:10s} {tot:7.1f}  "
              + "  ".join(f"{k} {v:.1f}" for k, v in
                          sorted(d.items(), key=lambda kv: -kv[1])))
    if tots and max(tots) > 0:
        print(f"      HEADROOM   slowest stage {b/max(tots):5.2f}x, "
              f"stages summed {b/max(sum(tots), 1e-9):5.2f}x")
    if s.get("margin_counters"):
        print("      E60        " + "  ".join(
            f"{k[2:]} {v}" for k, v in sorted(s["margin_counters"].items())))
    if not a.capture:
        # On a radio, x-real-time cannot exceed 1.0: the air arrives at 1x and
        # that is the whole supply.  "Held real time" therefore means the
        # queues never backed up and the deficit is startup, not backlog --
        # which is what the headroom line above measures independently.
        print(f"      (live source: x-real-time is CAPPED at 1.000x by the "
              f"air itself; the {(1 - s['x_real_time']) * s['wall_s']:.1f} s "
              f"deficit is acquisition, not backlog)")
    if s.get("error"):
        print(f"\n    ERROR: {s['error']}")
    if a.dump_dg and os.path.exists(a.dump_dg):
        h = hashlib.sha256(open(a.dump_dg, "rb").read()).hexdigest()
        s["dg"] = dict(path=a.dump_dg, bytes=os.path.getsize(a.dump_dg),
                       sha256=h)
        print(f"    datagram file      {os.path.getsize(a.dump_dg)} bytes, "
              f"sha256 {h[:32]}...")
    if a.json:
        with open(a.json, "w") as f:
            json.dump(s, f, indent=1, default=float)
        print(f"\n  wrote {a.json}")
    return 0 if not s.get("error") else 1


if __name__ == "__main__":
    raise SystemExit(main())

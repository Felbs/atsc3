"""M46 -- automatic level control (ALC) for the decoded programme audio.

WHY. Services differ a lot in level. Broadcaster A's stereo decodes 14-17 dB
under its own loudness metadata (the stream, not the decoder: ac4-decoder's
README, "Known issues"), and the worker's fixed scale (0.9 / REF_PEAK) puts
every element at whatever level the decoder hands it. Switching services
therefore jumps by up to ~20 dB. This is a leveller AFTER the decoder: it
reads only PCM, so any decoder (libac4dec or m42) gets the same treatment,
and with it switched off the worker writes the same bytes as before.

WHAT IT DOES, per worker (= per service):

  * MEASURES the programme's loudness the way ITU-R BS.1770-4 / EBU R128
    do: K-weighting (the two 48 kHz biquads of BS.1770), channel weights
    1.0 for L R C, 1.41 for Ls Rs, LFE excluded; 400 ms blocks every
    100 ms; an absolute gate at -70 LUFS and a relative gate 10 LU under
    the (abs-gated) running level. The measurement is FEED-FORWARD: it
    reads the signal before the gain, so the gain never chases itself.
  * ESTIMATES the level as a running mean of the gated blocks, in the
    power domain. For the first `tau` seconds of programme it is the plain
    mean of every gated block so far (the integrated loudness of what has
    been heard: fast and unbiased), afterwards an exponential mean with
    time constant `tau` (default 15 s of programme). Gated blocks only:
    speech pauses and silence (and the worker's own pad silence) neither
    move the estimate nor the gain -- no pumping.
  * STEERS the gain toward target - estimate, clamped to [-max_cut,
    +max_boost], slew-limited: fast (`fast_db_s`, default 15 dB/s) for the
    first `fast_s` (3) seconds of gated programme, `mid_db_s` (4 dB/s)
    while the estimate is still the plain mean (the first `tau` seconds),
    then slow: at most
    `down_db_s` (default 2 dB/s) when the programme got louder and
    `up_db_s` (default 1 dB/s) when it got quieter -- a loud commercial is
    taken down within seconds, a quiet passage is not chased. The estimate
    itself (time constant `tau`) is what keeps the ride slow; the slew
    limits only bound how fast a step in it is followed. The gain is decided once per 100 ms sub-block from
    the sub-blocks BEFORE it and ramped linearly (in dB) across the
    sub-block, so it is smooth (no zipper) and the output does not depend
    on how the stream is split into passes.
  * LIMITS: a peak limiter after the gain, on the largest |sample| of any
    channel, with a ceiling (default -1 dBFS) -- the 16-bit wav never
    clips, whatever the boost. NO LOOKAHEAD and NO ADDED LATENCY: the
    attack is instant (the gain drops on the sample that would exceed the
    ceiling) and it releases linearly in dB (`release_db_s`, default
    20 dB/s). It is a safety net, not a loudness tool: at -24 LUFS, ATSC
    A/85 programme peaks sit at or under -2 dBTP, so it fires rarely (the
    gate counts how rarely). It is causal per sample, so it too is
    independent of the pass split.
  * ONE GAIN FOR ALL CHANNELS: the balance and the 5.1 image are kept.

COST: two second-order filters per channel (scipy sosfilt) and a few numpy
passes per call -- ~1 % of a core for 5.1 at 48 kHz.

CONFIGURATION (environment, read by `Alc.from_env`; the worker also has
--alc / --alc-target):
    ATSC3_ALC              1 = on, 0 = off (the worker's default: off;
                           the live launchers set 1)
    ATSC3_ALC_TARGET       target loudness, LUFS        (default -24)
    ATSC3_ALC_MAX_BOOST    largest gain, dB             (default 20)
    ATSC3_ALC_MAX_CUT      largest attenuation, dB      (default 20)
    ATSC3_ALC_CEILING      limiter ceiling, dBFS        (default -1)
    ATSC3_ALC_TAU          slow time constant, s        (default 15)
    ATSC3_ALC_TRACE        file: one line per 100 ms sub-block (gain trace)
"""
from __future__ import annotations

import math
import os

import numpy as np
from scipy.signal import sosfilt

FS = 48000
SUB = FS // 10                  # 100 ms sub-block
BLOCK_SUBS = 4                  # 400 ms gating block = 4 sub-blocks, hop 100 ms
ABS_GATE = -70.0                # LUFS
REL_GATE = -10.0                # LU under the abs-gated level

# BS.1770-4 K-weighting at 48 kHz: the high shelf (stage 1) and the RLB
# high-pass (stage 2), as printed in the Recommendation.
K_SOS = np.array([
    [1.53512485958697, -2.69169618940638, 1.19839281085285,
     1.0, -1.69065929318241, 0.73248077421585],
    [1.0, -2.0, 1.0,
     1.0, -1.99004745483398, 0.99007225036621],
])


def channel_weights(n_ch):
    """BS.1770 weights in the worker's channel order (L R C LFE Ls Rs)."""
    if n_ch == 6:
        return np.array([1.0, 1.0, 1.0, 0.0, 1.41, 1.41])
    return np.ones(n_ch)


def _lufs(p):
    return -0.691 + 10.0 * math.log10(p) if p > 0 else -math.inf


def loudness_integrated(x, fs=FS):
    """BS.1770-4 integrated loudness of a whole (n, ch) signal, LUFS.

    The offline reference the gate measures the ALC with -- independent of
    the class below except for the shared filter constants."""
    assert fs == FS
    x = np.asarray(x, dtype=float)
    if x.ndim == 1:
        x = x[:, None]
    w = channel_weights(x.shape[1])
    y = sosfilt(K_SOS, x, axis=0)
    n_sub = len(y) // SUB
    if n_sub < BLOCK_SUBS:
        return -math.inf
    ms = (y[:n_sub * SUB] ** 2).reshape(n_sub, SUB, -1).mean(axis=1)
    sub_p = ms @ w
    c = np.convolve(sub_p, np.ones(BLOCK_SUBS) / BLOCK_SUBS, mode="valid")
    lb = -0.691 + 10 * np.log10(np.maximum(c, 1e-30))
    c = c[lb > ABS_GATE]
    if not len(c):
        return -math.inf
    rel = _lufs(c.mean()) + REL_GATE
    lb = -0.691 + 10 * np.log10(c)
    c = c[lb > rel]
    return _lufs(c.mean()) if len(c) else -math.inf


def true_peak_dbtp(x, over=4):
    """Peak of the 4x-oversampled signal (BS.1770 Annex 2 style), dBTP."""
    from scipy.signal import resample_poly
    x = np.asarray(x, dtype=float)
    if x.ndim == 1:
        x = x[:, None]
    pk = 0.0
    for c in range(x.shape[1]):
        pk = max(pk, float(np.abs(resample_poly(x[:, c], over, 1)).max()))
    return 20 * math.log10(pk) if pk > 0 else -math.inf


class Alc:
    """A streaming loudness leveller; see the module docstring.

    process(y) takes (n, ch) float PCM at 48 kHz, full scale 1.0, and
    returns the same shape. Any split of a stream into calls gives the same
    output samples."""

    def __init__(self, n_ch, target=-24.0, max_boost=20.0, max_cut=20.0,
                 ceiling=-1.0, tau=15.0, fast_s=3.0, fast_db_s=15.0,
                 mid_db_s=4.0, down_db_s=2.0, up_db_s=1.0, release_db_s=20.0, trace=None):
        self.n_ch = int(n_ch)
        self.w = channel_weights(self.n_ch)
        self.target = float(target)
        self.max_boost = float(max_boost)
        self.max_cut = float(max_cut)
        self.ceiling = 10 ** (float(ceiling) / 20)
        self.ceiling_db = float(ceiling)
        self.n_slow = max(1, int(round(tau * 10)))     # blocks per tau
        self.n_fast = max(1, int(round(fast_s * 10)))
        self.fast_step = fast_db_s / 10.0              # dB per sub-block
        self.mid_step = mid_db_s / 10.0
        self.down_step = down_db_s / 10.0
        self.up_step = up_db_s / 10.0
        self.release = release_db_s / FS               # dB per sample
        self.zi = np.zeros((K_SOS.shape[0], 2, self.n_ch))
        self.acc = np.zeros(self.n_ch)                 # current sub-block
        self.acc_n = 0
        self.subs = []                                 # last 3 completed
        self.p_abs = 0.0                               # abs-gated mean power
        self.k_abs = 0
        self.p_g = 0.0                                 # gated mean power
        self.k_g = 0
        self.g_prev = 0.0                              # dB at sub-block start
        self.g_next = 0.0                              # dB at sub-block end
        self.lim_db = 0.0                              # limiter gain, dB <= 0
        self.n_sub = 0                                 # sub-blocks completed
        # counters for the log
        self.n_lim = 0                                 # samples limited
        self.lim_min = 0.0
        self.last_block = -math.inf
        self.trace = trace
        self._trace_fh = None

    # ---- configuration -------------------------------------------------
    @staticmethod
    def enabled_from_env(default="0"):
        return os.environ.get("ATSC3_ALC", default).strip().lower() in (
            "1", "on", "yes", "true")

    @classmethod
    def from_env(cls, n_ch, target=None):
        e = os.environ.get

        def f(name, dflt):
            try:
                return float(e(name, dflt))
            except ValueError:
                return float(dflt)
        return cls(n_ch,
                   target=(f("ATSC3_ALC_TARGET", -24.0) if target is None
                           else float(target)),
                   max_boost=f("ATSC3_ALC_MAX_BOOST", 20.0),
                   max_cut=f("ATSC3_ALC_MAX_CUT", 20.0),
                   ceiling=f("ATSC3_ALC_CEILING", -1.0),
                   tau=f("ATSC3_ALC_TAU", 15.0),
                   trace=e("ATSC3_ALC_TRACE") or None)

    def describe(self):
        return (f"target {self.target:.1f} LUFS, boost <= "
                f"+{self.max_boost:.0f} dB, cut <= -{self.max_cut:.0f} dB, "
                f"limiter {self.ceiling_db:.1f} dBFS (no lookahead), "
                f"tau {self.n_slow / 10:.0f} s")

    # ---- state for a resumed worker ------------------------------------
    def state(self):
        return dict(g=round(self.g_next, 4), p_abs=self.p_abs,
                    k_abs=self.k_abs, p_g=self.p_g, k_g=self.k_g)

    def restore(self, st):
        try:
            self.g_prev = self.g_next = float(st["g"])
            self.p_abs, self.k_abs = float(st["p_abs"]), int(st["k_abs"])
            self.p_g, self.k_g = float(st["p_g"]), int(st["k_g"])
        except (KeyError, TypeError, ValueError):
            pass

    # ---- read-outs -------------------------------------------------------
    @property
    def gain_db(self):
        return self.g_next

    @property
    def level(self):
        """The running gated input loudness, LUFS (-inf before any)."""
        return _lufs(self.p_g) if self.k_g else -math.inf

    def status(self):
        lv = self.level
        s = (f"alc {self.g_next:+.1f} dB (in "
             + (f"{lv:.1f}" if lv > -math.inf else "--") + " LUFS)")
        if self.n_lim:
            s += f", limited {self.n_lim} smp (min {self.lim_min:.1f} dB)"
        return s

    # ---- the controller --------------------------------------------------
    def _block_done(self):
        """A sub-block completed: gate, update the estimate, step the gain."""
        sub = self.acc / SUB
        self.subs.append(sub)
        self.acc = np.zeros(self.n_ch)
        self.acc_n = 0
        self.n_sub += 1
        lb = -math.inf
        if len(self.subs) >= BLOCK_SUBS:
            self.subs = self.subs[-BLOCK_SUBS:]
            p = float(np.mean(self.subs, axis=0) @ self.w)
            lb = _lufs(p)
            if lb > ABS_GATE:
                self.k_abs += 1
                self.p_abs += (p - self.p_abs) / min(self.k_abs, self.n_slow)
                if lb > _lufs(self.p_abs) + REL_GATE:
                    self.k_g += 1
                    self.p_g += (p - self.p_g) / min(self.k_g, self.n_slow)
        self.last_block = lb
        g = self.g_next
        if self.k_g:
            want = min(self.max_boost,
                       max(-self.max_cut, self.target - _lufs(self.p_g)))
            if self.k_g <= self.n_fast:
                g = g + min(self.fast_step, max(-self.fast_step, want - g))
            elif self.k_g < self.n_slow:
                g = g + min(self.mid_step, max(-self.mid_step, want - g))
            else:
                g = g + min(self.up_step, max(-self.down_step, want - g))
        self.g_prev, self.g_next = self.g_next, g
        if self.trace:
            try:
                if self._trace_fh is None:
                    self._trace_fh = open(self.trace, "a")
                self._trace_fh.write(
                    f"{self.n_sub * 0.1:.1f} {lb:.2f} {self.level:.2f} "
                    f"{g:.3f} {self.lim_db:.2f}\n")
                self._trace_fh.flush()
            except OSError:
                self.trace = None

    def process(self, y):
        y = np.asarray(y, dtype=float)
        n = len(y)
        if n == 0:
            return y
        if y.ndim == 1:
            y = y[:, None]
        k, self.zi = sosfilt(K_SOS, y, axis=0, zi=self.zi)
        k2 = k * k
        g_db = np.empty(n)
        i = 0
        while i < n:
            take = min(n - i, SUB - self.acc_n)
            # the gain across this sub-block: a linear ramp (in dB) from the
            # value decided before the previous sub-block to the value
            # decided before this one -- known before its first sample
            t = (self.acc_n + np.arange(take) + 1) / SUB
            g_db[i:i + take] = self.g_prev + (self.g_next - self.g_prev) * t
            self.acc += k2[i:i + take].sum(axis=0)
            self.acc_n += take
            i += take
            if self.acc_n == SUB:
                self._block_done()
        out = y * (10.0 ** (g_db / 20.0))[:, None]
        # the limiter: instant attack, linear-in-dB release, no lookahead.
        # lim[j] = min(need[j], lim[j-1] + release, 0) is a min-plus
        # recursion; its closed form is a running minimum.
        pk = np.abs(out).max(axis=1)
        need = np.zeros(n)
        hot = pk > self.ceiling
        if hot.any() or self.lim_db < 0.0:
            need[hot] = 20.0 * np.log10(self.ceiling / pk[hot])
            # (need <= 0 everywhere, so the 0 dB cap is implied); the
            # carried value of the previous call starts the recursion.
            r = self.release
            j = np.arange(n)
            lim = np.minimum(np.minimum.accumulate(need - r * j) + r * j,
                             self.lim_db + r * (j + 1))
            self.lim_db = float(lim[-1])
            act = lim < 0.0
            if act.any():
                out = out * (10.0 ** (lim / 20.0))[:, None]
                self.n_lim += int(hot.sum())
                self.lim_min = min(self.lim_min, float(lim.min()))
            # float rounding must never let a sample past the ceiling
            np.clip(out, -self.ceiling, self.ceiling, out=out)
        return out

    def close(self):
        if self._trace_fh is not None:
            self._trace_fh.close()
            self._trace_fh = None

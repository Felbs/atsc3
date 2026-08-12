#!/usr/bin/env python3
"""atsc3_doctor -- diagnose an ATSC 3.0 receiver install before you blame the sky.

    python tools/atsc3_doctor.py                 # the full check
    python tools/atsc3_doctor.py --radio         # ALSO enumerate the SDR
    python tools/atsc3_doctor.py --capture F.cs16 --rate 6.912e6
    python tools/atsc3_doctor.py --selftest      # prove the checks can FAIL

Every check here exists because it cost this project real hours.  The list is
not "things that could go wrong"; it is things that DID, with the entry number
that paid for them:

    E73  wrong SDR gain          0.0% FEC vs 100.0% FEC, same air, same minute
    E64  clipped captures        a clipped file made the L1 analyser lie
    E66  disconnected input      rms 354 -- a capture of nothing, banked as data
    E56  BLAS thread storm       480% CPU, 0.09x, on a worker that runs 1.59x
    E56  VLC flags per BUILD     an unknown option is FATAL; distro != Windows
    PORT missing spec tables     a fresh clone dies with KeyError: '2026'
    E54  LDPC kernel per box     never copy the binary, rebuild it
    E61  orphaned decode procs   42 orphans held the live chain at 0.41x
    E56  single-tenant SDR       two chains fight; one wedges in SoapySDR unmake

THE RADIO IS NOT TOUCHED unless you pass --radio.  This box may be running an
unattended soak, and a doctor that breaks the patient to take its temperature
is not a doctor.  Even with --radio the check only ENUMERATES; it never tunes.

Exit codes: 0 all pass (warnings allowed), 1 at least one FAIL, 2 selftest
failed (a check could not be made to fail, so a clean bill means nothing).
"""
from __future__ import annotations

import argparse
import importlib
import os
import shutil
import subprocess
import sys
import tempfile
import textwrap
from pathlib import Path

HERE = Path(__file__).resolve()
REPO = HERE.parent.parent

PASS, WARN, FAIL, SKIP = "PASS", "WARN", "FAIL", "SKIP"
RESULTS: list[tuple[str, str, str, str]] = []   # (status, check, detail, fix)


def record(status: str, check: str, detail: str, fix: str = "") -> str:
    RESULTS.append((status, check, detail, fix))
    mark = {PASS: "  ok  ", WARN: " warn ", FAIL: " FAIL ", SKIP: " skip "}[status]
    print("[%s] %-26s %s" % (mark, check, detail))
    if fix and status in (WARN, FAIL):
        for line in textwrap.wrap(fix, 92):
            print("                                     -> " + line)
    return status


# ---------------------------------------------------------------- python env
def check_python() -> None:
    v = sys.version_info
    detail = "%d.%d.%d  %s" % (v.major, v.minor, v.micro, sys.executable)
    if v >= (3, 10):
        record(PASS, "python", detail)
    else:
        record(FAIL, "python", detail, "Python 3.10+ is required (3.12 and 3.14 "
                                       "are both exercised by this project).")


def check_modules() -> None:
    required = {
        "numpy": "the entire DSP chain",
        "scipy": "filters and resampling",
    }
    optional = {
        "psutil": "process supervision (atsc3_run/atsc3_warden reap by pid+"
                  "create_time identity; without it orphan reaping is weaker)",
        "fitz": "pymupdf -- a RUNTIME import for the AC-4 scalefactor-band "
                "tables (m27_sfb parses the ETSI PDF on import)",
        # Reported on EVERY run, not just --radio. Live tuning is the headline
        # feature, so a default run that says "0 problems" and is then followed
        # by ModuleNotFoundError on `atsc3 watch` is a false all-clear.
        "SoapySDR": "LIVE TV. Without it `atsc3 watch --rf N` cannot open a "
                    "radio -- only `--capture FILE` replay works",
        "torch": "the OPTIONAL gpu-full accelerator; never required",
    }
    for mod, why in required.items():
        try:
            m = importlib.import_module(mod)
            record(PASS, "module " + mod, getattr(m, "__version__", "?"))
        except Exception as e:
            record(FAIL, "module " + mod, "MISSING (%s)" % type(e).__name__,
                   "pip install %s -- needed for %s" % (mod, why))
    for mod, why in optional.items():
        try:
            m = importlib.import_module(mod)
            record(PASS, "module " + mod, getattr(m, "__version__", "present"))
        except Exception:
            record(WARN, "module " + mod, "not installed", "optional: " + why)


# ------------------------------------------------------- the spec runtime law
SPEC_TABLES = ["A322_2026_tbl.txt", "A322_2024_tbl.txt", "A322_2026.txt",
               "A322_2024.txt", "A331_tbl.txt"]
SPEC_AC4 = ["ts_10319001v010301p.pdf"]


def check_spec(spec_lab: Path, spec_root: Path) -> None:
    """The failure this catches is `KeyError: '2026'` on a fresh clone.

    The NUMERIC tables now ship (lab/spec_bank/), so a fresh clone decodes
    without them.  The extracted spec text is optional: when it is present the
    decoder prefers it, because two independent editions cross-checking each
    other is stronger evidence than one banked file.  The standards DOCUMENTS
    are copyrighted by their publishers and are not redistributed here.
    """
    bank = spec_lab.parent / "spec_bank"
    banked = [f for f in ("nuc_a322.npz", "ac4_sfb_offsets.json")
              if (bank / f).exists()]
    missing = [f for f in SPEC_TABLES if not (spec_lab / f).exists()]
    if not missing:
        record(PASS, "spec tables (lab/spec)", "%d/%d present (spec text wins)"
               % (len(SPEC_TABLES), len(SPEC_TABLES)))
    elif len(banked) == 2:
        record(PASS, "spec tables", "banked numeric tables present "
               "(%d/2); extracted spec text optional" % len(banked))
    else:
        record(FAIL, "spec tables",
               "no banked tables in %s AND no extracted spec text in %s"
               % (bank, spec_lab),
               "The decoder needs the A/322 constellations. Normally they ship "
               "in lab/spec_bank/ and there is nothing to do -- if that "
               "directory is missing, re-clone. To re-derive them yourself, "
               "download A/322 from the ATSC site, run the extraction step, "
               "and then `python lab/bank_tables.py`.")

    missing4 = [f for f in SPEC_AC4 if not (spec_root / f).exists()]
    if not missing4:
        record(PASS, "spec AC-4 (spec/)", "present")
    else:
        record(PASS if (spec_lab.parent / "spec_bank" / "ac4_sfb_offsets.json").exists()
               else WARN,
               "spec AC-4 (spec/)",
               "using banked band offsets" if (spec_lab.parent / "spec_bank" / "ac4_sfb_offsets.json").exists()
               else "MISSING %s" % ", ".join(missing4),
               "Only the AUDIO path needs the document: m27_sfb reads the "
               "scalefactor-band tables straight out of ETSI TS 103 190-1 with "
               "pymupdf. Video decodes without it. ETSI publishes the document "
               "free but returns 403 to automated fetches -- download it by "
               "hand, do not script around the access control.")


# ------------------------------------------------------ the optional C kernel
def check_ldpc_kernel() -> None:
    src = REPO / "lab" / "ldpc_kernel.c"
    out = REPO / "lab" / ("ldpc_kernel.dll" if os.name == "nt"
                          else "ldpc_kernel.so")
    if os.environ.get("ATSC3_LDPC_KERNEL") == "0":
        record(WARN, "LDPC C kernel", "DISABLED by ATSC3_LDPC_KERNEL=0",
               "The numpy fallback is ~7x slower in the LDPC stage. Unset the "
               "variable unless you are deliberately gating against the "
               "reference path.")
        return
    if not src.exists():
        record(FAIL, "LDPC C kernel", "lab/ldpc_kernel.c missing", "Bad clone.")
        return
    if not out.exists():
        record(WARN, "LDPC C kernel", "not built",
               "OPTIONAL but worth ~7x on the LDPC stage: run "
               "`python lab/build_ldpc_kernel.py`, then `python lab/gate_e54.py`. "
               "Needs any C compiler. On Windows with MSYS2, "
               r"C:\msys64\mingw64\bin must LEAD your PATH or cc1 dies with no "
               "message at all.")
        return
    try:
        import ctypes
        ctypes.CDLL(str(out))
    except OSError as e:
        record(FAIL, "LDPC C kernel", "present but will not load: %s" % e,
               "This is what a COPIED binary looks like. The kernel is built "
               "per box and must never be copied between machines -- delete it "
               "and run `python lab/build_ldpc_kernel.py` here.")
        return
    stale = out.stat().st_mtime < src.stat().st_mtime
    if stale:
        record(WARN, "LDPC C kernel", "loads, but OLDER than ldpc_kernel.c",
               "Rebuild: python lab/build_ldpc_kernel.py")
    else:
        record(PASS, "LDPC C kernel", "built and loads (%s)" % out.name)


# ---------------------------------------------------------- the BLAS throttle
TINY_MATMUL = (
    "import os,time,numpy as np\n"
    "a=np.random.randn(64,2).astype(np.float32)\n"
    "b=np.random.randn(2,4).astype(np.float32)\n"
    "t=time.perf_counter()\n"
    "for _ in range(20000): a@b\n"
    "print('%.4f'%(time.perf_counter()-t))\n"
)


def check_blas() -> None:
    """A BLAS pool under a Python loop of tiny ops is a THROTTLE dressed as
    parallelism (E56: 480% CPU, ~0.09x, on a worker the big box runs at 1.59x;
    E82: the same shape DEADLOCKED eight workers, twice).

    So this does not just read the environment -- it MEASURES both ways and
    shows you the ratio on your own box.
    """
    pins = {k: os.environ.get(k) for k in
            ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS",
             "NUMEXPR_NUM_THREADS")}
    set_pins = {k: v for k, v in pins.items() if v}

    def timed(env_extra: dict) -> float | None:
        env = os.environ.copy()
        env.update(env_extra)
        try:
            r = subprocess.run([sys.executable, "-c", TINY_MATMUL],
                               capture_output=True, text=True, env=env,
                               timeout=180)
            return float(r.stdout.strip()) if r.returncode == 0 else None
        except Exception:
            return None

    unpinned = timed({k: "" for k in pins})       # let BLAS choose
    pinned = timed({k: "1" for k in pins})
    if unpinned is None or pinned is None or pinned <= 0:
        record(WARN, "BLAS tiny-op behaviour", "could not measure", "")
    else:
        ratio = unpinned / pinned
        detail = "20k tiny sgemms: unpinned %.2fs vs pinned-to-1 %.2fs (%.2fx)" \
                 % (unpinned, pinned, ratio)
        if ratio > 1.5:
            record(FAIL, "BLAS tiny-op behaviour", detail,
                   "Your BLAS fans every tiny matmul across a thread pool and "
                   "pays more in handoff than the arithmetic costs. Export "
                   "OPENBLAS_NUM_THREADS=1 (and OMP/MKL) for the AUDIO WORKERS, "
                   "and pin the chain to a small number (this project ships "
                   "workers=1, chain=4). Note the E82 corollary: pinning "
                   "through the environment only fixes the process that "
                   "remembered to -- code on a hot path should not call a "
                   "multi-threaded BLAS from inside a worker pool at all.")
        else:
            record(PASS, "BLAS tiny-op behaviour", detail)
    record(PASS if set_pins else WARN, "BLAS pins in env",
           str(set_pins) if set_pins else "none set",
           "" if set_pins else "Not fatal -- the launchers export their own "
           "pins. But if you run a worker by hand, export "
           "OPENBLAS_NUM_THREADS=1 first.")


# --------------------------------------------------------------- ffmpeg / VLC
def check_ffmpeg() -> None:
    exe = shutil.which("ffmpeg")
    if not exe:
        record(FAIL, "ffmpeg", "not on PATH",
               "Required for the mux and for every null-sink verification. "
               "Install ffmpeg and re-run.")
        return
    try:
        r = subprocess.run([exe, "-version"], capture_output=True, text=True,
                           timeout=30)
        ver = r.stdout.splitlines()[0] if r.stdout else "?"
    except Exception as e:
        record(WARN, "ffmpeg", "found but would not run: %s" % e, "")
        return
    record(PASS, "ffmpeg", ver[:70])
    # Not a fault -- a fact people trip over. ffmpeg DEMUXES AC-4 and has no
    # decoder for it, and MP4 has no tag to write it back out with.
    record(SKIP, "ffmpeg AC-4", "no AC-4 decoder in any ffmpeg build",
           "Expected. That is why this project has its own AC-4 decoder, and "
           "why proof clips are transcoded to AAC ('Could not find tag for "
           "codec ac4').")


VLC_FLAGS = ["--no-qt-privacy-ask", "--no-qt-updates-notif", "--dummy-quiet",
             "--play-and-exit", "--file-caching"]


def _vlc_help(exe: str) -> str | None:
    for args in (["--longhelp", "--advanced", "--help-verbose"],
                 ["--longhelp", "--advanced"], ["--longhelp"], ["-H"]):
        try:
            r = subprocess.run([exe] + args, capture_output=True, text=True,
                               timeout=60)
            blob = (r.stdout or "") + (r.stderr or "")
            if len(blob) > 2000:
                return blob
        except Exception:
            continue
    return None


def check_vlc(extra_flags: list[str] | None = None) -> None:
    """E56: option sets differ per BUILD, not per version -- and an unknown
    option is FATAL to VLC, so a flag that is merely cosmetic on one box kills
    the player one second after the window opens on another.
    """
    exe = shutil.which("vlc") or shutil.which("vlc.exe")
    if not exe:
        alt = shutil.which("mpv") or shutil.which("ffplay")
        record(WARN, "VLC", "not on PATH (fallback found: %s)" % (alt or "none"),
               "VLC is the reference viewer (audio-track picker + caption "
               "toggle). mpv/ffplay work but lose the track UI. Check the "
               "TOOL, not the intention: this project once logged 'mpv' "
               "windows for weeks that were really ffplay via a fallback "
               "branch.")
        return
    blob = _vlc_help(exe)
    if blob is None:
        record(WARN, "VLC", "found at %s, could not read its option list" % exe, "")
        return
    ver = ""
    try:
        r = subprocess.run([exe, "--version"], capture_output=True, text=True,
                           timeout=30)
        ver = (r.stdout or r.stderr or "").splitlines()[0][:50]
    except Exception:
        pass
    unsupported = [f for f in (VLC_FLAGS + (extra_flags or []))
                   if f.lstrip("-").split("=")[0] not in blob and f not in blob]
    if unsupported:
        record(WARN, "VLC flags on THIS build", "%s :: unsupported here: %s"
               % (ver, " ".join(unsupported)),
               "An unknown option is FATAL to VLC -- it exits about a second "
               "after the window opens, which reads exactly like a decode "
               "failure. Remove these from the player command line on this "
               "box. --no-qt-updates-notif is absent from distro builds (no "
               "update checker compiled in) and --dummy-quiet is Windows-only.")
    else:
        record(PASS, "VLC flags on THIS build", "%s :: all probed flags exist"
               % ver)


# ------------------------------------------------------------ capture hygiene
def check_capture(path: Path, rate: float | None) -> str:
    """E64/E66/E67: a clipped capture makes every downstream analyser lie, and
    a disconnected input produces a file that looks like data.

    The capture law: samples == wall x fs, or the file is VOID.
    """
    import numpy as np
    if not path.exists():
        return record(FAIL, "capture " + path.name, "no such file", "")
    n_bytes = path.stat().st_size
    if n_bytes % 4:
        return record(FAIL, "capture " + path.name,
                      "%d bytes is not a whole number of cs16 samples" % n_bytes,
                      "This reader assumes interleaved int16 I/Q (.cs16).")
    nsamp = n_bytes // 4
    x = np.fromfile(path, dtype=np.int16)
    if x.size == 0:
        return record(FAIL, "capture " + path.name, "empty file", "")
    peak = int(np.max(np.abs(x.astype(np.int32))))
    full = int(np.count_nonzero(np.abs(x.astype(np.int32)) >= 32767))
    frac = full / x.size
    rms = float(np.sqrt(np.mean((x.astype(np.float64)) ** 2)))
    dbfs = 20.0 * np.log10(max(peak, 1) / 32768.0)

    bits = ["%d samples" % nsamp, "max|s|=%d" % peak, "peak %.1f dBFS" % dbfs,
            "rms %.0f" % rms, "full-scale %.4f%%" % (100 * frac)]
    if rate:
        secs = nsamp / rate
        bits.append("%.3f s @ %.3f Msps" % (secs, rate / 1e6))

    if frac > 0.0:
        return record(FAIL, "capture " + path.name, "  ".join(bits),
                      "VOID -- CLIPPED. The tuner front end was saturated, so "
                      "every L1 and SNR number derived from this file is a "
                      "measurement of the clipping, not of the air. Lower "
                      "rfgain_sel (READ IT BACK -- an unverified writeSetting "
                      "is a hope) and re-capture. E66 retracted a whole L1 "
                      "verdict to this cause.")
    if rms < 1000:
        return record(FAIL, "capture " + path.name, "  ".join(bits),
                      "rms this low is a DISCONNECTED INPUT, not a weak "
                      "signal -- E66 caught an antenna mid-swap at rms 354. "
                      "Check the coax before you check the code, and rename "
                      "the file VOID rather than deleting it.")
    if dbfs > -3.0:
        return record(WARN, "capture " + path.name, "  ".join(bits),
                      "No headroom. Not clipped, but one fade away from it; "
                      "aim for about -20 dBFS peak.")
    return record(PASS, "capture " + path.name, "  ".join(bits))


# ------------------------------------------------------------- the radio, SDR
def check_radio(enumerate_devices: bool) -> None:
    if not enumerate_devices:
        record(SKIP, "SDR enumeration", "not attempted (pass --radio)",
               "The SDR is SINGLE-TENANT. If an unattended run is in progress, "
               "enumerating from a second process is how you get a chain "
               "wedged inside SoapySDR unmake.")
        return
    try:
        import SoapySDR
    except Exception as e:
        record(FAIL, "SDR enumeration", "SoapySDR not importable (%s)" % type(e).__name__,
               "The live path needs SoapySDR with a driver for your device. "
               "Replay from a banked capture works without it. Note that the "
               "interpreter that has SoapySDR may not be the one you are "
               "running the doctor in.")
        return
    try:
        devs = SoapySDR.Device.enumerate()
    except Exception as e:
        record(FAIL, "SDR enumeration", "enumerate() raised %s" % e,
               "A driver that raises on enumerate is usually a wedged USB "
               "device. On Linux try `systemctl restart sdrplay` after a "
               "kill -9 of any stale chain; on Windows a ghost-enumerated "
               "device is reboot-only -- no amount of replugging clears it.")
        return
    if not devs:
        record(FAIL, "SDR enumeration", "0 devices",
               "Zero devices with the vendor's own app working means the API "
               "directory is off PATH. Check the device is enumerated by the "
               "OS first, then the SoapySDR module path. Use a SHORT DIRECT "
               "USB 3.0 port -- not a hub, not a long cable -- and route the "
               "USB cable away from the antenna coax.")
        return
    record(PASS, "SDR enumeration", "; ".join(
        str(dict(d)).replace("'", "")[:60] for d in devs[:3]))


# -------------------------------------------------- who else is using the box
def check_single_tenant() -> None:
    """The SDR is single-tenant and orphaned decode workers are invisible
    weight (E61: 42 orphans held the live chain at 0.41x).

    SELF-MATCH LAW: a census that matches its own censor kills the wrong
    process.  This one only REPORTS, and it excludes itself by pid.
    """
    try:
        import psutil
    except Exception:
        record(SKIP, "other atsc3 processes", "psutil not installed", "")
        return
    me = os.getpid()
    mine = []
    for p in psutil.process_iter(["pid", "name", "cmdline"]):
        try:
            if p.info["pid"] == me:
                continue
            cl = " ".join(p.info["cmdline"] or "")
            if "atsc3" in cl and "atsc3_doctor" not in cl:
                mine.append((p.info["pid"], cl[:70]))
        except Exception:
            continue
    if not mine:
        record(PASS, "other atsc3 processes", "none")
        return
    record(WARN, "other atsc3 processes", "%d running" % len(mine),
           "Something is already using this repo and possibly the radio. Do "
           "NOT start a second chain, and do NOT kill from an unfiltered "
           "census. Listed: "
           + " | ".join("%d %s" % m for m in mine[:6]))


def check_disk(root: Path) -> None:
    try:
        usage = shutil.disk_usage(root)
    except Exception:
        return
    free_gb = usage.free / 1e9
    if free_gb < 20:
        record(FAIL, "free disk", "%.1f GB" % free_gb,
               "Live capture and the growing .ts sink eat tens of GB per "
               "session. A malformed header field once allocated 42 GB in one "
               "go, so keep real margin.")
    elif free_gb < 100:
        record(WARN, "free disk", "%.1f GB" % free_gb,
               "Enough to watch, tight for a long soak.")
    else:
        record(PASS, "free disk", "%.1f GB" % free_gb)


# ------------------------------------------------------- THE NEGATIVE CONTROL
def selftest() -> int:
    """Prove the checks can FAIL.

    A doctor that has never been seen to report illness is not evidence of
    health. Each leg below builds a KNOWN-BAD input and requires the matching
    check to fail on it, and a KNOWN-GOOD input and requires it to pass.
    """
    import numpy as np
    ok = True
    scratch = Path(tempfile.mkdtemp(prefix="atsc3_doctor_selftest_"))
    print("selftest scratch: %s\n" % scratch)

    # --- spec tables: an empty dir must FAIL, a populated one must PASS ------
    empty = scratch / "empty_spec"
    empty.mkdir()
    before = len(RESULTS)
    check_spec(empty, empty)
    got = [r[0] for r in RESULTS[before:]]
    if FAIL in got:
        print("  pass  spec check FAILS on an empty spec dir (the fresh-clone "
              "KeyError: '2026')")
    else:
        ok = False
        print("  FAIL  spec check did NOT fail on an empty spec dir")

    full = scratch / "full_spec"
    full.mkdir()
    for f in SPEC_TABLES:
        (full / f).write_text("stub", encoding="utf-8")
    for f in SPEC_AC4:
        (full / f).write_text("stub", encoding="utf-8")
    before = len(RESULTS)
    check_spec(full, full)
    got = [r[0] for r in RESULTS[before:]]
    if FAIL not in got:
        print("  pass  spec check PASSES when the tables are present")
    else:
        ok = False
        print("  FAIL  spec check failed on a populated spec dir")

    # --- capture integrity: clipped / dead / clean --------------------------
    rng = np.random.default_rng(7)
    n = 200000

    clean = (rng.normal(0, 3000, n * 2)).astype(np.int16)
    clean_f = scratch / "clean.cs16"
    clean.tofile(clean_f)
    before = len(RESULTS)
    check_capture(clean_f, 6.912e6)
    if RESULTS[before][0] == PASS:
        print("  pass  capture check PASSES a clean synthetic capture")
    else:
        ok = False
        print("  FAIL  capture check rejected a clean capture (%s)"
              % RESULTS[before][1])

    clipped = clean.copy()
    clipped[:400] = 32767
    clipped_f = scratch / "clipped.cs16"
    clipped.tofile(clipped_f)
    before = len(RESULTS)
    check_capture(clipped_f, 6.912e6)
    if RESULTS[before][0] == FAIL:
        print("  pass  capture check FAILS a clipped capture (E64/E66)")
    else:
        ok = False
        print("  FAIL  capture check did not catch clipping")

    dead = (rng.normal(0, 354, n * 2)).astype(np.int16)
    dead_f = scratch / "disconnected.cs16"
    dead.tofile(dead_f)
    before = len(RESULTS)
    check_capture(dead_f, 6.912e6)
    if RESULTS[before][0] == FAIL:
        print("  pass  capture check FAILS a disconnected-input capture (E66)")
    else:
        ok = False
        print("  FAIL  capture check did not catch a dead input")

    # --- VLC flag probe: an impossible flag must be reported ----------------
    exe = shutil.which("vlc") or shutil.which("vlc.exe")
    if exe:
        before = len(RESULTS)
        check_vlc(extra_flags=["--this-flag-cannot-exist-e83"])
        st, _, detail, _ = RESULTS[before]
        if st == WARN and "this-flag-cannot-exist" in detail:
            print("  pass  VLC probe reports a flag this build does not have")
        else:
            ok = False
            print("  FAIL  VLC probe missed an impossible flag (%s)" % detail)
    else:
        print("  skip  VLC not installed, cannot exercise the flag probe")

    shutil.rmtree(scratch, ignore_errors=True)
    print()
    print("DOCTOR SELFTEST %s" % ("PASSED -- the checks are proven to fail on "
                                  "known-bad input" if ok else
                                  "FAILED -- a clean bill of health from this "
                                  "build means nothing"))
    return 0 if ok else 2


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--radio", action="store_true",
                    help="also ENUMERATE the SDR (never tunes). Do not use "
                         "while an unattended run holds the radio.")
    ap.add_argument("--capture", action="append", default=[],
                    help="check a .cs16 capture for clipping / dead input")
    ap.add_argument("--rate", type=float, default=6.912e6)
    ap.add_argument("--skip-blas", action="store_true",
                    help="skip the BLAS measurement (it spawns two short "
                         "CPU-bound subprocesses)")
    ap.add_argument("--selftest", action="store_true")
    args = ap.parse_args()

    if args.selftest:
        return selftest()

    print("=" * 96)
    print("ATSC 3.0 RECEIVER -- INSTALL DOCTOR      repo: %s" % REPO)
    print("=" * 96)

    print("\n-- environment " + "-" * 81)
    check_python()
    check_modules()
    check_disk(REPO)

    print("\n-- inputs the decoder cannot run without " + "-" * 55)
    check_spec(REPO / "lab" / "spec", REPO / "spec")
    check_ldpc_kernel()

    print("\n-- the performance traps " + "-" * 71)
    if args.skip_blas:
        record(SKIP, "BLAS tiny-op behaviour", "--skip-blas", "")
    else:
        check_blas()

    print("\n-- the viewing stack " + "-" * 75)
    check_ffmpeg()
    check_vlc()

    print("\n-- the radio " + "-" * 83)
    check_radio(args.radio)
    check_single_tenant()

    if args.capture:
        print("\n-- captures " + "-" * 84)
        for c in args.capture:
            check_capture(Path(c), args.rate)

    print("\n" + "=" * 96)
    n = {s: sum(1 for r in RESULTS if r[0] == s) for s in (PASS, WARN, FAIL, SKIP)}
    print("SUMMARY   pass %d   warn %d   FAIL %d   skip %d"
          % (n[PASS], n[WARN], n[FAIL], n[SKIP]))
    if n[FAIL]:
        print("\nBlocking problems:")
        for st, check, detail, fix in RESULTS:
            if st == FAIL:
                print("  * %-26s %s" % (check, detail))
    print("=" * 96)
    return 1 if n[FAIL] else 0


if __name__ == "__main__":
    sys.exit(main())

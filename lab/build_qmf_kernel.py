#!/usr/bin/env python3
"""Build the E109 compiled AC-4 QMF bank -> lab/qmf_kernel.{dll,so}.

Same fleet laws as lab/build_ldpc_kernel.py, for the same reasons:
  * OPTIONAL: if this is never run (or fails), m9_fast falls back to the
    numpy cell_pool_fast path with a one-line log.  The numpy float64 EXACT path
    remains the gate reference either way.
  * Never copy the binary between machines -- run this script per box.
    Pure C ABI (no Python.h), so one build serves every CPython on that
    box, but not another box.
  * IEEE-pinned: no -ffast-math ever; -ffp-contract=off (gcc/clang) and
    /fp:precise (MSVC), so no FMA contraction can change a float32
    rounding relative to the numpy path this kernel replaces.

Run:  python lab/build_qmf_kernel.py [--check-only]
Then: python lab/gate_e109_qmf.py      (the identity gate -- run it)
"""
from __future__ import annotations

import argparse
import ctypes
import os
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, "qmf_kernel.c")
OUT = os.path.join(HERE, "qmf_kernel.dll" if os.name == "nt"
                   else "qmf_kernel.so")

GCC_FLAGS = ["-O3", "-ffp-contract=off", "-shared"]


def _run(cmd, **kw):
    print("  $", " ".join(cmd) if isinstance(cmd, list) else cmd)
    return subprocess.run(cmd, capture_output=True, text=True, **kw)


def build_windows():
    cands = []
    if shutil.which("gcc"):
        cands.append(("gcc (PATH)", shutil.which("gcc"), None))
    msys = r"C:\msys64\mingw64\bin"
    if os.path.exists(os.path.join(msys, "gcc.exe")):
        cands.append(("gcc (MSYS2 mingw64)", os.path.join(msys, "gcc.exe"),
                      msys))
    for name, gcc, prepend in cands:
        env = os.environ.copy()
        if prepend:
            env["PATH"] = prepend + os.pathsep + env["PATH"]
        r = _run([gcc] + GCC_FLAGS + ["-static-libgcc", "-o", OUT, SRC],
                 env=env)
        if r.returncode == 0 and os.path.exists(OUT):
            return name
        print(f"  {name} failed (rc {r.returncode}): "
              f"{(r.stderr or r.stdout).strip()[:300]}")
    vsw = (r"C:\Program Files (x86)\Microsoft Visual Studio\Installer"
           r"\vswhere.exe")
    if os.path.exists(vsw):
        r = _run([vsw, "-products", "*", "-latest", "-requires",
                  "Microsoft.VisualStudio.Component.VC.Tools.x86.x64",
                  "-property", "installationPath"])
        for root in r.stdout.splitlines():
            vcvars = os.path.join(root.strip(), "VC", "Auxiliary", "Build",
                                  "vcvars64.bat")
            if not os.path.exists(vcvars):
                continue
            cmd = (f'"{vcvars}" >nul 2>&1 && '
                   f'cl /nologo /O2 /fp:precise /LD /Fe:"{OUT}" "{SRC}"')
            r2 = _run(cmd, shell=True, cwd=HERE)
            if r2.returncode == 0 and os.path.exists(OUT):
                for ext in (".lib", ".exp", ".obj"):
                    p = os.path.join(HERE, "qmf_kernel" + ext)
                    if os.path.exists(p):
                        os.remove(p)
                return "MSVC (vcvars64)"
            print(f"  MSVC failed: {(r2.stderr or r2.stdout).strip()[:300]}")
    return None


def build_posix():
    for cc in ("cc", "gcc", "clang"):
        if not shutil.which(cc):
            continue
        r = _run([cc] + GCC_FLAGS + ["-fPIC", "-o", OUT, SRC])
        if r.returncode == 0 and os.path.exists(OUT):
            return cc
        print(f"  {cc} failed (rc {r.returncode}): "
              f"{(r.stderr or r.stdout).strip()[:300]}")
    return None


def smoke(path):
    """Round-trip a known signal through the kernel bank and check the SNR.

    The bank's own gate (m33_qmf.gate) analyses then synthesises noise and
    measures how well the output reproduces the input after the group delay.
    Reusing that here means the smoke test checks the DSP, not just that the
    symbols load.
    """
    import numpy as np
    dll = ctypes.CDLL(path)
    dll.qmf_kernel_abi.restype = ctypes.c_int32
    if int(dll.qmf_kernel_abi()) != 1:
        return "ABI mismatch"
    sys.path.insert(0, HERE)
    import m33_qmf as Q
    f64 = np.ctypeslib.ndpointer(np.float64, flags="C")
    dll.qmf_analyse_f64.restype = ctypes.c_int32
    dll.qmf_analyse_f64.argtypes = [f64, ctypes.c_int32, f64, f64, f64]
    dll.qmf_synthesise_f64.restype = ctypes.c_int32
    dll.qmf_synthesise_f64.argtypes = [f64, ctypes.c_int32, f64, f64, f64]
    w = Q.qwin()
    M = np.ascontiguousarray(Q.analysis_matrix(), np.complex128)
    # 255 is the constant the A-SPX renderer actually uses
    N = np.ascontiguousarray(Q.synthesis_matrix(255), np.complex128)
    rng = np.random.default_rng(7)
    x = rng.standard_normal(64 * 200)
    nts = len(x) // 64
    Qc = np.empty((64, nts), np.complex128)
    rc = dll.qmf_analyse_f64(np.ascontiguousarray(x), nts, w,
                             M.view(np.float64), Qc.view(np.float64))
    if rc != 0:
        return f"analyse returned {rc}"
    y = np.empty(nts * 64, np.float64)
    rc = dll.qmf_synthesise_f64(Qc.view(np.float64), nts, w,
                                N.view(np.float64), y)
    if rc != 0:
        return f"synthesise returned {rc}"
    best, lag = -1.0, 0
    for d in range(0, 1024):
        if d + 4096 > len(y):
            break
        c = float(np.corrcoef(x[:4096], y[d:d + 4096])[0, 1])
        if c > best:
            best, lag = c, d
    err = x[:4096] - y[lag:lag + 4096]
    snr = 10 * np.log10(np.sum(x[:4096] ** 2) / max(np.sum(err ** 2), 1e-30))
    if snr < 40.0:
        return f"round-trip SNR only {snr:.1f} dB at lag {lag}"
    # and agree with the reference implementation it replaces
    An = Q.analyse(x, w, M)
    Sn = Q.synthesise(An, w, N)
    da = float(np.abs(An - Qc).max())
    ds = float(np.abs(Sn - y).max())
    if da > 1e-10 or ds > 1e-10:
        return f"disagrees with m33_qmf: analyse {da:.3e}, synthesise {ds:.3e}"
    print(f"  round-trip SNR {snr:.1f} dB at group delay {lag}; "
          f"vs m33_qmf: analyse {da:.2e}, synthesise {ds:.2e}")
    return None


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--check-only", action="store_true",
                    help="only smoke-test an existing binary")
    a = ap.parse_args(argv)
    if not a.check_only:
        if not os.path.exists(SRC):
            print(f"FAIL: {SRC} missing")
            return 1
        tool = build_windows() if os.name == "nt" else build_posix()
        if tool is None:
            print("FAIL: no toolchain built the kernel; the numpy fast path "
                  "remains in use (a working state, just slower)")
            return 1
        print(f"built {OUT} with {tool}")
    err = smoke(OUT)
    if err:
        print(f"FAIL smoke test: {err}")
        return 1
    print(f"smoke test PASS: {OUT}")
    print("identity gate: python lab/gate_e109_qmf.py  (run it)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""Build the E103 compiled cell-pool equaliser -> lab/cellpool_kernel.{dll,so}.

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

Run:  python lab/build_cellpool_kernel.py [--check-only]
Then: python lab/gate_e103_cellpool.py      (the identity gate -- run it)
"""
from __future__ import annotations

import argparse
import ctypes
import os
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, "cellpool_kernel.c")
OUT = os.path.join(HERE, "cellpool_kernel.dll" if os.name == "nt"
                   else "cellpool_kernel.so")

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
                    p = os.path.join(HERE, "cellpool_kernel" + ext)
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
    """Equalise one hand-checkable symbol: a flat channel and unit pilots.

    With Yrow == refpk on the pilot bins the channel estimate is exactly 1,
    so every data cell must come back as its own bin value, moved to the slot
    hmat names.  That exercises the divide, the interpolation and the scatter
    together, and its expected answer needs no reference implementation.
    """
    import numpy as np
    dll = ctypes.CDLL(path)
    dll.cellpool_kernel_abi.restype = ctypes.c_int32
    if int(dll.cellpool_kernel_abi()) != 3:
        return "ABI mismatch (rebuild: the source gained derot_f64)"
    i32 = np.ctypeslib.ndpointer(np.int32, flags="C")
    f64 = np.ctypeslib.ndpointer(np.float64, flags="C")
    dll.cellpool_row_f64.restype = ctypes.c_int32
    dll.cellpool_row_f64.argtypes = [f64, ctypes.c_int32, i32, f64,
                                     ctypes.c_int32, i32, f64,
                                     ctypes.c_int32, i32, i32, i32,
                                     ctypes.c_int32, f64, f64]
    nfft, npk, nd = 64, 4, 6
    # +1 so no bin is 0+0j -- a zero reference pilot is a 0/0 divide,
    # which is a property of the fixture, not of the kernel.
    Y = ((np.arange(nfft) + 1.0)
         + 1j * (np.arange(nfft) + 2.0)).astype(np.complex128)
    bins_pk = np.array([0, 8, 16, 24], np.int32)
    refpk = Y[bins_pk].copy()                       # -> hp == 1 exactly
    jj = np.zeros(nd, np.int32)
    wt = np.zeros(nd, np.float64)                   # H = hp[0] == 1
    bins_d = np.array([1, 2, 3, 4, 5, 6], np.int32)
    dmap = np.arange(nd, dtype=np.int32)
    hmat = np.array([5, 4, 3, 2, 1, 0], np.int32)   # reversal
    hp = np.zeros(2 * npk, np.float64)
    xg = np.zeros(2 * nd, np.float64)
    rc = dll.cellpool_row_f64(Y.view(np.float64), nfft, bins_pk,
                              refpk.view(np.float64), npk, jj, wt, nd,
                              bins_d, dmap, hmat, nd, hp, xg)
    if rc != 0:
        return f"kernel returned {rc}"
    got = xg.view(np.complex128)
    want = Y[bins_d][::-1]
    if not np.allclose(got, want, rtol=0, atol=1e-12):
        return f"flat-channel identity failed: {got!r} != {want!r}"
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
    print("identity gate: python lab/gate_e103_cellpool.py  (run it)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

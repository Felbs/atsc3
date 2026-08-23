#!/usr/bin/env python3
"""Build the E102 compiled BICM demapper -> lab/demap_kernel.{dll,so}.

Same fleet laws as lab/build_ldpc_kernel.py, for the same reasons:
  * OPTIONAL: if this is never run (or fails), m9_fast falls back to the
    numpy fast path with a one-line log.  The numpy float64 EXACT path
    remains the gate reference either way.
  * Never copy the binary between machines -- run this script per box.
    Pure C ABI (no Python.h), so one build serves every CPython on that
    box, but not another box.
  * IEEE-pinned: no -ffast-math ever; -ffp-contract=off (gcc/clang) and
    /fp:precise (MSVC), so no FMA contraction can change a float32
    rounding relative to the numpy path this kernel replaces.

Run:  python lab/build_demap_kernel.py [--check-only]
Then: python lab/gate_e102_demap.py      (the identity gate -- run it)
"""
from __future__ import annotations

import argparse
import ctypes
import os
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, "demap_kernel.c")
OUT = os.path.join(HERE, "demap_kernel.dll" if os.name == "nt"
                   else "demap_kernel.so")

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
                    p = os.path.join(HERE, "demap_kernel" + ext)
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
    """Load the fresh binary and demap one hand-checkable QPSK cell.

    QPSK with points (+-1 +- 1j)/sqrt(2) and a received cell sitting exactly
    on point 0: every bit's '0' subset contains the nearest point, so both
    LLR differences must be positive and the d2 must be ~0.
    """
    import numpy as np
    dll = ctypes.CDLL(path)
    dll.demap_kernel_abi.restype = ctypes.c_int32
    abi = int(dll.demap_kernel_abi())
    if abi != 1:
        return f"ABI {abi} != 1"
    f32 = np.ctypeslib.ndpointer(np.float32, flags="C")
    dll.demap_llr_f32.restype = ctypes.c_int32
    dll.demap_llr_f32.argtypes = [f32, f32, ctypes.c_int32, ctypes.c_int32,
                                  ctypes.c_int32, ctypes.c_int32,
                                  f32, f32, f32, f32, f32]
    r = float(np.sqrt(0.5))
    pts = np.array([r + 1j * r, r - 1j * r, -r + 1j * r, -r - 1j * r],
                   np.complex128)
    pfx = (pts.real * -2.0).astype(np.float32)
    pfy = (pts.imag * -2.0).astype(np.float32)
    p2f = (pts.real ** 2 + pts.imag ** 2).astype(np.float32)
    zr = np.array([r], np.float32)
    zi = np.array([r], np.float32)
    out = np.zeros(2, np.float32)
    d2 = np.zeros(1, np.float32)
    rc = dll.demap_llr_f32(zr, zi, 1, 1, 4, 2, pfx, pfy, p2f, out, d2)
    if rc != 0:
        return f"kernel returned {rc}"
    if abs(float(d2[0])) > 1e-5:
        return f"d2 {float(d2[0])!r} should be ~0 on top of a point"
    if not (out[0] > 0 and out[1] > 0):
        return f"LLR signs wrong: {out!r}"
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
    print("identity gate: python lab/gate_e102_demap.py  (run it)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

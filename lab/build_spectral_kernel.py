#!/usr/bin/env python3
"""Build the E112 compiled AC-4 spectral Huffman parse -> lab/spectral_kernel.{dll,so}.

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

Run:  python lab/build_spectral_kernel.py [--check-only]
Then: python lab/gate_e112_spectral.py      (the identity gate -- run it)
"""
from __future__ import annotations

import argparse
import ctypes
import os
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, "spectral_kernel.c")
OUT = os.path.join(HERE, "spectral_kernel.dll" if os.name == "nt"
                   else "spectral_kernel.so")

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
                    p = os.path.join(HERE, "spectral_kernel" + ext)
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
    """Decode a tiny hand-built stream through a hand-built codebook.

    Book: two symbols, '0' -> 0 and '10' -> 1 (dim 2, mod 2, off 1, so every
    decoded pair is signed values in {-1, 0}); one section over one band.
    The expected bits and lines are computable by hand, so the test needs no
    reference implementation.
    """
    import numpy as np
    dll = ctypes.CDLL(path)
    dll.spectral_kernel_abi.restype = ctypes.c_int32
    if int(dll.spectral_kernel_abi()) != 1:
        return "ABI mismatch"
    i32 = np.ctypeslib.ndpointer(np.int32, flags="C")
    i64 = np.ctypeslib.ndpointer(np.int64, flags="C")
    u8 = np.ctypeslib.ndpointer(np.uint8, flags="C")
    dll.asf_spectral_i32.restype = ctypes.c_int32
    dll.asf_spectral_i32.argtypes = [u8, ctypes.c_int64, i64,
                                     i32, ctypes.c_int32, i32, ctypes.c_int32,
                                     i32, i32, i32, i32, i32,
                                     i32, ctypes.c_int32]
    # tree: root node 0: bit0 -> leaf sym0, bit1 -> node 1; node 1: bit0 ->
    # leaf sym1, bit1 -> invalid
    tree = np.array([-(0 + 1), 1, -(1 + 1), 0], np.int32)
    troot = np.full(12, -1, np.int32)
    troot[1] = 0                       # put the book at cb=1 (offset 1 = signed-in-word)
    cb_dim = np.zeros(12, np.int32); cb_dim[1] = 2
    cb_mod = np.zeros(12, np.int32); cb_mod[1] = 2
    cb_off = np.zeros(12, np.int32); cb_off[1] = 1
    # stream: codeword '10' (sym 1 -> vals (1//2-1, 1%2-1) = (-1, 0)) then '0'
    # (sym 0 -> (0-1, 0-1) = (-1, -1)); cb_off != 0 so NO sign bits.
    # bits: 1 0 0 -> byte 0b100_00000
    buf = np.array([0b10000000], np.uint8)
    sects = np.array([1, 0, 2], np.int32)          # cb=1, bands 0..2
    offsets = np.array([0, 2, 4], np.int32)        # band 0 -> lines 0..2, band 1 -> 2..4
    lines = np.zeros(4, np.int32)
    pos = np.zeros(1, np.int64)
    rc = dll.asf_spectral_i32(buf, 1, pos, sects, 1, offsets, 3,
                              cb_dim, cb_mod, cb_off, tree, troot, lines, 4)
    if rc != 0:
        return f"kernel returned {rc}"
    if pos[0] != 3:
        return f"bit position {pos[0]}, expected 3"
    if list(lines) != [-1, 0, -1, -1]:
        return f"lines {list(lines)}, expected [-1, 0, -1, -1]"
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
    print("identity gate: python lab/gate_e112_spectral.py  (run it)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

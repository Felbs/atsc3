#!/usr/bin/env python3
"""Build the E108 BCH syndrome kernel -> lab/demap_kernel.{dll,so}.

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

Run:  python lab/build_bch_kernel.py [--check-only]
Then: python lab/gate_e107_ldpc.py (its BCH legs) and the M17 byte gate      (the identity gate -- run it)
"""
from __future__ import annotations

import argparse
import ctypes
import os
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, "bch_kernel.c")
OUT = os.path.join(HERE, "bch_kernel.dll" if os.name == "nt"
                   else "bch_kernel.so")

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
                    p = os.path.join(HERE, "bch_kernel" + ext)
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
    """The kernel's parities equal numpy's on random words and masks."""
    import numpy as np
    dll = ctypes.CDLL(path)
    dll.bch_kernel_abi.restype = ctypes.c_int32
    if int(dll.bch_kernel_abi()) != 1:
        return "ABI mismatch"
    u64 = np.ctypeslib.ndpointer(np.uint64, flags="C")
    u8 = np.ctypeslib.ndpointer(np.uint8, flags="C")
    dll.bch_parity_u64.restype = ctypes.c_int32
    dll.bch_parity_u64.argtypes = [u64, ctypes.c_int32, ctypes.c_int32, u64, ctypes.c_int32, u8]
    rng = np.random.default_rng(108)
    words = rng.integers(0, 2 ** 63, (7, 101), dtype=np.int64).astype(np.uint64)
    mask = rng.integers(0, 2 ** 63, (13, 101), dtype=np.int64).astype(np.uint64)
    par = np.zeros((7, 13), np.uint8)
    if dll.bch_parity_u64(words, 7, 101, mask, 13, par) != 0:
        return "kernel returned non-zero"
    ref = (np.bitwise_count(words[:, None, :] & mask[None, :, :]).sum(2) & 1).astype(np.uint8)
    return None if np.array_equal(par, ref) else "parities differ from numpy"

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
    print("identity gate: python lab/gate_e107_ldpc.py (its BCH legs) and the M17 byte gate  (run it)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

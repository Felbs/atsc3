#!/usr/bin/env python3
"""Build the E107 layered LDPC kernel -> lab/demap_kernel.{dll,so}.

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

Run:  python lab/build_ldpc_layered_kernel.py [--check-only]
Then: python lab/gate_e107_ldpc.py      (the identity gate -- run it)
"""
from __future__ import annotations

import argparse
import ctypes
import os
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, "ldpc_layered_kernel.c")
OUT = os.path.join(HERE, "ldpc_layered_kernel.dll" if os.name == "nt"
                   else "ldpc_layered_kernel.so")

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
                    p = os.path.join(HERE, "ldpc_layered_kernel" + ext)
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
    """Load the fresh binary and decode two hand-checkable words on a 4-bit
    chain code (checks {0,1}, {1,2}, {2,3}): all-positive LLRs must retire at
    iteration 0; one weakly wrong bit between two strong neighbours must be
    corrected in one layered sweep."""
    import numpy as np
    dll = ctypes.CDLL(path)
    dll.ldpc_layered_abi.restype = ctypes.c_int32
    abi = int(dll.ldpc_layered_abi())
    if abi != 4:
        return f"ABI {abi} != 4"
    i32 = np.ctypeslib.ndpointer(np.int32, flags="C")
    f32 = np.ctypeslib.ndpointer(np.float32, flags="C")
    u8 = np.ctypeslib.ndpointer(np.uint8, flags="C")
    dll.ldpc_layered_f32.restype = ctypes.c_int32
    dll.ldpc_layered_f32.argtypes = [f32, ctypes.c_int32, ctypes.c_int32,
                                     ctypes.c_int32, ctypes.c_int32, i32, i32,
                                     ctypes.c_int32, ctypes.c_float, ctypes.c_int32, ctypes.c_int32, ctypes.c_int32, ctypes.c_int32,
                                     u8, u8, i32, i32, f32]
    n, cptr = 4, np.array([0, 2, 4, 6], np.int32)
    cidx = np.array([0, 1, 1, 2, 2, 3], np.int32)
    llr = np.array([[5, 5, 5, 5], [5, -0.5, 5, 5]], np.float32)
    bits = np.zeros((2, n), np.uint8); conv = np.zeros(2, np.uint8)
    its = np.zeros(2, np.int32); bad = np.zeros(2, np.int32)
    work = np.zeros((n + 6) * 8, np.float32)
    rc = dll.ldpc_layered_f32(llr, 2, n, 3, 6, cptr, cidx, 10, 0.75, 1, 0, 0, 0,
                              bits, conv, its, bad, work)
    if rc == -3:
        return "no AVX2 on this machine: the layered kernel cannot run here"
    if rc != 0:
        return f"kernel returned {rc}"
    if not (conv[0] == 1 and its[0] == 0 and bits[0].sum() == 0):
        return f"clean word: conv {conv[0]} its {its[0]} bits {bits[0]}"
    if not (conv[1] == 1 and its[1] == 1 and bits[1].sum() == 0):
        return f"one wrong bit: conv {conv[1]} its {its[1]} bits {bits[1]}"
    # E118: the lane-compacting entry point on 19 blocks (more than the eight
    # lanes): the two words above repeated, so lanes must be refilled, and
    # every block must come back exactly as the eight-at-a-time call gave it.
    dll.ldpc_layered_pack_f32.restype = ctypes.c_int32
    dll.ldpc_layered_pack_f32.argtypes = dll.ldpc_layered_f32.argtypes
    m = 19
    llr19 = np.ascontiguousarray(np.tile(llr, (10, 1))[:m])
    bits19 = np.zeros((m, n), np.uint8); conv19 = np.zeros(m, np.uint8)
    its19 = np.zeros(m, np.int32); bad19 = np.zeros(m, np.int32)
    rc = dll.ldpc_layered_pack_f32(llr19, m, n, 3, 6, cptr, cidx, 10, 0.75, 1, 0, 0, 0,
                                   bits19, conv19, its19, bad19, work)
    if rc != 0:
        return f"pack kernel returned {rc}"
    want_its = np.tile(its, 10)[:m]
    if not (conv19.all() and (its19 == want_its).all() and bits19.sum() == 0):
        return f"pack: conv {conv19} its {its19} (want {want_its})"
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
    print("identity gate: python lab/gate_e107_ldpc.py  (run it)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

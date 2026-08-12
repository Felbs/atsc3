#!/usr/bin/env python3
"""Build the E54 compiled min-sum LDPC kernel -> lab/ldpc_kernel.{dll,so}.

Fleet laws honoured here:
  * OPTIONAL: if this is never run (or fails), m3_ldpc falls back to the
    numpy batch decoder with a one-line log.  The numpy float64 path
    REMAINS the gate reference either way.
  * Never copy the binary between machines -- run this script per box
    (PORT.md).  The kernel is pure C ABI (no Python.h), so ONE build
    serves every CPython on that box, but not another box.
  * IEEE-pinned flags: no -ffast-math ever; -ffp-contract=off (gcc/clang)
    and /fp:precise (MSVC) so no FMA contraction or reassociation can
    change a float32 rounding vs the numpy fast path.  gcc and MSVC
    builds were measured output-IDENTICAL on the 222-block real-air set.

Windows toolchains tried in order: gcc on PATH, MSYS2 mingw64 gcc
(C:\\msys64\\mingw64\\bin must lead PATH or cc1 dies silently -- handled),
MSVC Build Tools via vswhere+vcvars64.  POSIX: cc, gcc, clang.

Run:  python lab/build_ldpc_kernel.py [--check-only]
"""
from __future__ import annotations

import argparse
import ctypes
import os
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, "ldpc_kernel.c")
OUT = os.path.join(HERE, "ldpc_kernel.dll" if os.name == "nt"
                   else "ldpc_kernel.so")

GCC_FLAGS = ["-O3", "-ffp-contract=off", "-shared"]


def _run(cmd, **kw):
    print("  $", " ".join(cmd) if isinstance(cmd, list) else cmd)
    return subprocess.run(cmd, capture_output=True, text=True, **kw)


def build_windows():
    # 1/2: gcc (PATH, then MSYS2's mingw64 with its bin dir prepended --
    # invoking that gcc without it on PATH makes cc1 fail with no message)
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
    # 3: MSVC via vswhere + vcvars64
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
                    p = os.path.join(HERE, "ldpc_kernel" + ext)
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
    """Load the fresh binary and decode a trivial all-erasure block."""
    import numpy as np
    dll = ctypes.CDLL(path)
    dll.ldpc_kernel_abi.restype = ctypes.c_int32
    abi = int(dll.ldpc_kernel_abi())
    if abi != 1:
        return f"ABI {abi} != 1"
    sys.path.insert(0, HERE)
    import m3_ldpc as LD
    checks = LD.parity_check("3/15", 16200)
    LD._KERNEL_STATE["lib"] = None            # force a re-probe of the file
    bits, conv, its, bad = LD._kernel_decode_checked(
        np.full((2, 16200), 5.0, np.float32), checks, 16200, 50, 1.0)
    if bits is None:
        return "kernel path declined the decode"
    if not (conv.all() and (its == 0).all() and (bad == 0).all()
            and not bits.any()):
        return "all-zero codeword did not retire at iteration 0"
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
            print("FAIL: no toolchain built the kernel; the numpy fallback "
                  "remains in use (that is a working state, just slower)")
            return 1
        print(f"built {OUT} with {tool}")
    err = smoke(OUT)
    if err:
        print(f"FAIL smoke test: {err}")
        return 1
    print(f"smoke test PASS: {OUT}")
    print("identity gate: python lab/gate_e54.py  (real air, run it)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

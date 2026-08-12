"""`python -m atsc3 ...` -- the front door.

    python -m atsc3 watch --rf 33          tune, decode and PLAY, live
    python -m atsc3 gate  --rf 33          prove the streaming path == the
                                           batch path, on air decoded now

The milestone tools themselves live in `lab/`, which is where every rung of
this ladder has lived since M1 (the repository ships the code and gitignores
the captures).  This module only dispatches, so `python -m atsc3 watch` and
`python lab/m11_watch.py` are the same program.
"""
from __future__ import annotations

import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
LAB = os.path.join(os.path.dirname(HERE), "lab")

USAGE = """usage: python -m atsc3 <command> [options]

  watch    tune an ATSC 3.0 multiplex and play it live
  gate     re-derive the batch reference and diff the streaming path against it

  python -m atsc3 watch --help   for the options
"""


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    if not argv or argv[0] in ("-h", "--help", "help"):
        print(USAGE)
        return 0
    cmd, rest = argv[0], argv[1:]
    if LAB not in sys.path:
        sys.path.insert(0, LAB)
    if cmd == "watch":
        import m11_watch
        return m11_watch.main(rest)
    if cmd == "gate":
        import m11_gate
        return m11_gate.main(rest)
    print(f"unknown command {cmd!r}\n")
    print(USAGE)
    return 2


if __name__ == "__main__":
    raise SystemExit(main())

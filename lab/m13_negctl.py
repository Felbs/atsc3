#!/usr/bin/env python3
"""M13 PLP1 negative controls -- prove the gate can FAIL.

Reuses m13_sf1's own pipeline. Runs the identical subframe-1 / PLP-1
extraction with the CORRECT params and then with three deliberately WRONG
ones. The correct run must converge; every control must collapse. A gate that
passes on wrong params is not a gate (lab burned 7 harnesses on 8/07 for
exactly this reason).
"""
import copy, os, sys, time
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import m13_sf1 as M13
import m6_cells as C
from m2_pilots import load, guess_format
from concurrent.futures import ThreadPoolExecutor

CAP = os.path.join(HERE, "long_rf33.cs16")
L1 = os.path.join(HERE, "m4_l1detail_hit_rf33.json")
import json
js = json.load(open(L1))
l1b, l1d = js["l1_basic"], js["l1_detail"]
g0 = C.Geometry.rf33_legacy()
g = M13.geometry_for_subframe(l1b, l1d, 1, g0)
plp = {p["id"]: p for p in g.plps}[1]
fmt = guess_format(CAP)

total = M13.frame_samples(g0) + (g.nfft + g.gi) * g.nsym
per = total / M13.FS_POST
span = per + 0.06
off = M13.subframe_offset(g0, 1)
y, _fs, _cfo, _h = load(CAP, 8e6, fmt=fmt, span_sec=span, start_sec=0.0)
t0, coh = M13.fine_t0(y, g0)
EX = ThreadPoolExecutor(max_workers=12)

NB = 12  # blocks per run; a collapse is obvious in the first dozen

def run(tag, p, fi):
    t = time.time()
    try:
        r = M13.decode_subframe(y, t0 + off, g, p, fi, iters=50,
                                max_blocks=NB, ex=EX, verbose=False)
    except Exception as exc:  # noqa: BLE001
        print(f"  {tag:34s} STRUCTURAL REJECT ({type(exc).__name__}) -- "
              f"wrong params cannot even form FEC blocks  ({time.time()-t:.0f}s)")
        return None
    print(f"  {tag:34s} conv {r.get('converged'):3d}/{r.get('n_blocks')}"
          f"  bch {r.get('bch_ok'):3d}  median-unsat {r.get('unsat_median')}"
          f"  ({time.time()-t:.0f}s)")
    return r

print(f"PLP1 negative controls  (t0 {t0}, coh {coh:.4f}, {NB} blocks each)")
print("=" * 72)

# CORRECT
run("CORRECT 256QAM/nti=3/fi=0", plp, 0)

# CONTROL A: wrong modulation (QPSK instead of 256QAM)
pA = copy.deepcopy(plp); pA["mod"] = "QPSK"
pA["cells_per_fec"] = pA["ninner"] // M13._BITS["QPSK"]
run("CTRL wrong mod = QPSK", pA, 0)

# CONTROL B: wrong TI-block count (nti=2 instead of 3)
pB = copy.deepcopy(plp); pB["nti"] = 2
run("CTRL wrong nti = 2", pB, 0)

# CONTROL C: wrong frequency-interleaver origin (no reset: fi=36)
run("CTRL wrong fi_offset = 36", plp, 36)
print("=" * 72)
print("gate is valid iff CORRECT converged and ALL controls collapsed to ~0")

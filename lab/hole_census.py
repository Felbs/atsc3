#!/usr/bin/env python3
"""hole_census.py -- count lane holes from a live *.idx (E117 referee).
A hole = consecutive-fragment seq delta above the modal delta; missing
slots = (delta - modal)/modal rounded. Prints per-lane summary."""
import json, sys, collections
for path in sys.argv[1:]:
    seqs=[json.loads(l)["seq"] for l in open(path)]
    if len(seqs)<3: print(f"{path}: only {len(seqs)} frags"); continue
    ds=[b-a for a,b in zip(seqs,seqs[1:])]
    modal=collections.Counter(ds).most_common(1)[0][0]
    holes=[(seqs[i],round(d/modal)-1) for i,d in enumerate(ds) if d>modal]
    miss=sum(h[1] for h in holes)
    total=len(ds)*1+miss
    print(f"{path}: {len(seqs)} frags, modal delta {modal}, "
          f"{len(holes)} holes, {miss} missing slots ({100*miss/max(total,1):.1f}%)")
    if holes: print("   last holes:", holes[-5:])

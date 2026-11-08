#!/usr/bin/env python3
"""atsc3_stress_summary.py OUTDIR -- one table from the files atsc3_stress.sh left."""
import glob, json, os, re, subprocess, sys


def lane_stats(path):
    r = subprocess.run(["ffprobe", "-v", "error", "-select_streams", "v:0", "-show_entries", "packet=pts_time",
                        "-of", "csv=p=0", path], capture_output=True, text=True)
    pts = sorted(float(x.rstrip(",")) for x in r.stdout.split() if x.rstrip(","))
    if len(pts) < 2:
        return f"{len(pts)} frames"
    gaps = sum(1 for a, b in zip(pts, pts[1:]) if b - a > 0.05)
    return f"{len(pts)} frames, {pts[-1] - pts[0]:.1f} s, {gaps} gaps"


def main(out):
    rows = []
    for js in sorted(glob.glob(os.path.join(out, "rf*_*.json")), key=lambda p: (int(re.search(r"rf(\d+)", p).group(1)), p)):
        name = os.path.basename(js)[:-5]
        log = open(js[:-5] + ".log").read() if os.path.exists(js[:-5] + ".log") else ""
        try:
            j = json.load(open(js))
        except Exception as e:                                 # noqa: BLE001
            rows.append(f"{name}: no result ({e})"); continue
        m = re.search(r"sniffed: (.*)", log); lay = re.search(r"layers: (.*)", log)
        cfo = re.findall(r"CFO ([+-][\d.]+) Hz", log)
        p = j.get("plp1") or {}
        ldm = "LDM" in (lay.group(1) if lay else "") or "LDM plan adopted" in log
        kern = "layered kernel unavailable" in log or "ABI mismatch" in log
        lanes = []
        d = js[:-5]
        for f in sorted(glob.glob(os.path.join(d, "live*video*.m4s"))):
            if os.path.getsize(f) > 0:
                lanes.append(f"{os.path.basename(f)}: {lane_stats(f)}")
        rows.append(
            f"{name}\n"
            f"   {m.group(1)[:80] if m else 'no L1 sniffed'}\n"
            + (f"   {lay.group(1)[:95]}\n" if lay else "")
            + f"   frames {j.get('frames')}  dropped {j.get('pipeline', {}).get('frame_dropped', 0)}  x rt {j.get('x_real_time', 0):.2f}  "
              f"reacq {j.get('front_end', {}).get('reacquire', 0)}  CFO {cfo[0] if cfo else '?'} Hz{'  [LDM path]' if ldm else ''}{'  [KERNEL FELL BACK]' if kern else ''}\n"
            f"   PLP 0 (core) {j.get('fec_converged')}/{j.get('fec_total')}   PLP 1 (M17) {p.get('fec_converged', '-')}/{p.get('fec_total', '-')}\n"
            + ("".join(f"   {x}\n" for x in lanes) if lanes else "   no video lane written\n"))
    print("\n".join(rows))


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else ".")

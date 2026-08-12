#!/usr/bin/env python3
"""gate_e50.py -- the clock-trim servo holds the cushion; open loop drains.

Simulates the three-clock reality: media appended at the broadcaster rate
(1.0 by definition), playhead advancing at (1 + ppm_err) * rate. The
controller samples every 5 s exactly as the viewer's poll loop does.
Negative control = the identical simulation with the servo OFF: cushion
must visibly drain, proving the gate can fail.
"""
import importlib.util, os, sys

spec = importlib.util.spec_from_file_location(
    "tv", os.path.join(os.path.dirname(__file__), "..", "tools",
                       "atsc3_tv.py"))
tv = importlib.util.module_from_spec(spec)
spec.loader.exec_module(tv)

def simulate(ppm_err, hours, servo, noise_seed=12345):
    """-> (max |cushion - target|, final rate). Deterministic noise."""
    trim = tv.ClockTrim(target=24.0)
    rng = noise_seed
    cushion, rate = 24.0, 1.0
    worst = 0.0
    t = 0.0
    while t < hours * 3600:
        t += 5.0
        # playhead consumes at (1+err)*rate; media arrives at 1.0
        cushion += 5.0 * (1.0 - (1.0 + ppm_err) * rate)
        # REALISTIC noise: VLC reports integer seconds -> the measured
        # cushion carries the playhead's fractional second, plus a little
        # jitter. Deterministic LCG for the jitter phase.
        rng = (1103515245 * rng + 12345) % (1 << 31)
        frac = (t * (1.0 + ppm_err) * rate) % 1.0
        meas = cushion + frac + ((rng / (1 << 31)) - 0.5) * 0.1
        if servo:
            trim.sample(t, meas)
            nr = trim.evaluate(t)
            if nr is not None:
                rate = nr
        worst = max(worst, abs(cushion - 24.0))
    return worst, rate

fails = 0
# CASE 1: -300 ppm (player fast), 2 h, servo ON -> cushion held within 2 s
w, r = simulate(-300e-6, 3, servo=True)
ok = w < 2.5
print(f"{'PASS' if ok else 'FAIL'} servo ON  -300ppm 3h: worst dev {w:.2f}s "
      f"(want <2.5), final rate {(r-1)*1e6:+.0f} ppm")
fails += not ok
# CASE 2: +500 ppm (player slow / mux fast), 2 h, servo ON
w, r = simulate(+500e-6, 3, servo=True)
ok = w < 2.5
print(f"{'PASS' if ok else 'FAIL'} servo ON  +500ppm 3h: worst dev {w:.2f}s "
      f"(want <2.5), final rate {(r-1)*1e6:+.0f} ppm")
fails += not ok
# CASE 3 (NEGATIVE CONTROL): -300 ppm, servo OFF -> must drain past 2 s
w, _ = simulate(-300e-6, 3, servo=False)
ok = w > 2.5
print(f"{'PASS' if ok else 'FAIL'} NEG-CTL servo OFF -300ppm 3h: worst dev "
      f"{w:.2f}s (must exceed 2.5 -- proves the gate can fail)")
fails += not ok
# CASE 4: zero error, servo ON -> deadband keeps rate at 1.0 (no hunting)
w, r = simulate(0.0, 3, servo=True)
ok = abs(r - 1.0) < 200e-6 and w < 1.5
print(f"{'PASS' if ok else 'FAIL'} servo ON  0ppm 3h: no hunting "
      f"(rate {(r-1)*1e6:+.0f} ppm, worst dev {w:.2f}s)")
fails += not ok
print(f"GATE E50: {4 - fails}/4")
sys.exit(1 if fails else 0)

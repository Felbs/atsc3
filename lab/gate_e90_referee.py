"""E90 gate: the referee self-heal, WITH a negative control.

Leg A  idle tuner (ch=none)          -> must retune and return a lock
Leg B  tuner parked on someone
       else's channel, not locked    -> must return None and NOT retune
"""
import os, sys, time, types, subprocess
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, r"tools")
import atsc3_judge as J

HB, ID, T = J.HDHR_BIN, J.HDHR_ID, J.HDHR_TUNER
def hd(*a):
    return subprocess.run([HB, ID, *a], capture_output=True, text=True,
                          timeout=15).stdout.strip()

a = types.SimpleNamespace(live_dir=r"data\e31")
j = J.Judge.__new__(J.Judge)          # no __init__: we only need the method
j.log = lambda m: print("   judge:", m.strip())

fails = []

# ---- LEG A: idle tuner must be reclaimed -------------------------------
hd("set", f"/tuner{T}/channel", "none"); time.sleep(2)
print("A pre :", hd("get", f"/tuner{T}/status"))
r = j.rf_referee()
print("A ret :", r)
if not (r and r[2] not in ("none", "")):
    fails.append("LEG A: idle tuner was not healed (got %r)" % (r,))

# ---- LEG B: someone else's tuner must be left alone --------------------
hd("set", f"/tuner{T}/channel", "auto:107000000"); time.sleep(3)
before = hd("get", f"/tuner{T}/channel")
print("B pre :", before, "|", hd("get", f"/tuner{T}/status"))
r = j.rf_referee()
after = hd("get", f"/tuner{T}/channel")
print("B ret :", r, "| channel after:", after)
if r is not None:
    fails.append("LEG B: reported a reading from an unlocked tuner: %r" % (r,))
if after != before:
    fails.append("LEG B: RETUNED a tuner it does not own (%s -> %s)"
                 % (before, after))

# ---- restore ----------------------------------------------------------
hd("set", f"/tuner{T}/channel", J.HDHR_CHANNEL); time.sleep(5)
print("restored:", hd("get", f"/tuner{T}/status"))

print()
if fails:
    print("GATE FAIL"); [print("  -", f) for f in fails]; sys.exit(1)
print("GATE PASS: heals an idle tuner, never touches a busy one")

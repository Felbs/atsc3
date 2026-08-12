# LAWS

Rules this receiver was taught by being wrong first. They are written to be
useful outside this project — none of them are really about television.

The lab notebook they were extracted from is not published: it is 364 KB of
unedited narrative, and the timestamps in it are worth nothing to a reader.
These are the transferable part.

---

## On measurement

**1. Instrument the output you actually consume, not a proxy for it.**
A player's clock ran forward 25 seconds through a completely frozen picture,
and the watchdog written that evening polled the clock. "Frames actually
displayed" is the honest primitive. Every proxy is a claim that two things move
together; that claim is exactly what fails when something breaks.

**2. Plausibility is the trap, not silence.**
Five separate instruments returned entirely reasonable numbers from a code path
that was not connected to the thing under test. A number that looks right is
not evidence that anything is right. A zero, from an instrument you have proven
fires, is a much stronger result than a plausible non-zero from one you haven't.

**3. A gate without a negative control is a decoration.**
Every detector here is required to demonstrate, on each run, that it still
fires on a planted specimen — and that it stays silent on a clean one. This
caught a refactor that had silently disabled a whole rule, before any scan
result from it was believed.

**4. Check the instrument before you believe the negative.**
Re-running an old control moved a published-looking result by 1.94 dB. The
*instrument* had drifted. Without the control, this project would have
announced an improvement that was partly its own measurement moving underneath
it.

**5. Profiling a starved pipeline measures the failure mode, not the code.**
Throughput numbers gathered while the input was signal-starved described the
starvation. When the antenna improved, the same code got 6× faster with no
change at all — and retracted a conclusion we had already drawn.

---

## On correctness

**6. A shared constant is invisible to cross-checking.**
Three independent subsystems agreed on a sample rate. All three were wrong by
4.27%, because all three read the same constant. Agreement between things that
share an input is not corroboration.

**7. Gate the symmetry, not the part.**
A stereo balance bug survived three separate quality gates: the left channel
was exemplary and the right was 20× too loud, and all three gates measured the
left. When two things are supposed to match, test the *relationship*, not each
side.

**8. "The parse closed" is never sufficient.**
Accepting short frames raised a bitstream's closure rate from 64.6% to 76.1%
while every correctness gate collapsed. A parser that consumes its input
without complaining has demonstrated only that it consumed its input.

**9. A decoder only knows the broadcasters it has met.**
Taking a second station's audio from 0% of frames to 100% took two fixes: one
core mode and one header field that the first broadcaster simply never set.
Every `if this flag is set` branch you have not implemented is a station you
cannot hear yet.

**10. Two roads, one destination.**
The strongest correctness claim available is an independent path arriving at a
byte-identical result. This receiver's live streaming output matches an offline
batch decode of the same signal, SHA-256 for SHA-256, on two operating systems
and two processors.

---

## On debugging

**11. A bug can wear the costume of a hardware limit.**
A transmitter written off as "needs a better antenna" decoded perfectly once a
one-frame interleaver offset was fixed — 242 ms of bookkeeping between "this
station is dead" and "this station is perfect." Meanwhile every *quality*
metric read excellent and every *correctness* metric read dead. When those two
disagree, believe the correctness metric.

**12. A crash truncates the evidence.**
A confident diagnosis of "monotone decay" turned out to be the acquisition
transient. The steady state did not exist in any log, because the crash always
arrived first. Fix the crash before theorising about the shape of the data.

**13. A missing spec clause is not a weak signal.**
Verification sat at 0 for 30 straight attempts because a repetition term — 48%
of the codeword — was simply absent from the implementation. The arithmetic had
never once matched the value the standard prints, and nothing was comparing
them.

**14. Near a threshold, a percentage is a cliff, not a scale.**
Three settings "ranked" 5.21%, 1.28% and 0.50% were within 0.38 dB of each
other. Ranking anything by a metric that is falling off a cliff produces a
confident, stable, meaningless ordering.

**15. The same invariant in two places must carry the same consequence.**
One code path treated a condition as recoverable and swept on; the other raised
on it. The receiver therefore died at exactly the moment it most needed to
re-acquire.

---

## On building things that read from the air

**16. Bound every field you read from a broadcast.**
ATSC 3.0's link-layer headers carry no integrity protection of their own. One
flipped bit killed a media lane permanently; another told the decoder to
allocate 42 GB. This is a *class* to audit, not two bugs to fix — treat the
wire as adversarial, because entropy alone is adversary enough.

**17. An open port still detects the loudest carrier.**
A receiver connected to the wrong antenna port still produced signal, still
locked, and still reported numbers. "It sees something" is not "it is connected
to what you think."

**18. A recorded failure is a reference, not an authority.**
Five entries in this project's notebook are formal retractions of earlier
entries in the same notebook. A negative result needs a stated cause and a
re-test trigger, or it silently becomes folklore that closes doors that are
actually open.

---

## On the ethics of a receiver

**19. Read what is in the clear; leave locked things locked.**
Encrypted services here are detected, enumerated and labeled — never attacked,
never circumvented. Where a locked service leaves its guide and captions in the
clear, those are read, because they are in the clear. This is enforced by a
build gate, not by a promise in a document.

**20. An audit report is a map of the exposure.**
The tool written to find identifiers in this repository was itself, in its
first draft, a tidy sorted index of every identifier it existed to find — worse
than the report, because the tool definitely ships. Rules became *shapes*; the
specific values live in a file that is never published.

# LLS is public by construction — a datacasting field note

*A receiver-project observation about ATSC 3.0 signaling. Contains no
credentials, no endpoints, and does not identify the operator involved. The
specific instance behind the security note was reported privately and resolved
before this was written.*

> **The lesson in one line:** ATSC 3.0 Low Level Signaling is transmitted
> unencrypted by design, so **every field an operator places in it is
> world-readable to any receiver in coverage.** Treat LLS as a public
> bulletin board, never as a place to put anything that should stay private.

## Background: what LLS is, and why it can't be secret

To find anything on an ATSC 3.0 multiplex, a receiver first has to read the
**Low Level Signaling (LLS)** — the small, fixed set of tables (A/331) that
say what services exist, where their components live, and how to tune them.
LLS is the bootstrap. It has to be readable by a receiver that knows nothing
yet, which means it is sent **in the clear**: it cannot be encrypted without
breaking the very step that lets a receiver join the multiplex.

Most LLS tables are standardized (the service list, the system time, region
ratings, and so on). One table id — **255** — is reserved for
*broadcaster-defined* use. Anything a station wants to signal that the standard
doesn't cover can ride there.

## What we observed

While bringing up an open ATSC 3.0 receiver, we decoded the LLS on several
real multiplexes. On some of them, table 255 carried a **datacasting overlay**
— a user-defined signaling structure describing a commercial *data* service
riding the broadcast (file delivery / edge distribution to receivers), rather
than a television service. This is one of NextGen TV's headline non-TV business
models, and it was plainly live on the air: real commercial traffic, not just
a spec promise.

Two details are worth stating for other receiver builders:

1. **It rides the same path you already decode.** A general ATSC 3.0 receiver
   parses LLS to find its TV channels; the datacast signaling sits right there
   in the same table set. You see it for free. "The signal" is much more than
   the picture.

2. **It behaved like a per-station deployment.** Multiplexes in one market
   carried the overlay; a multiplex examined identically in a neighboring
   market did not carry it at all — a clean negative control that the overlay
   is a deliberate, station-by-station choice rather than something baked into
   every ATSC 3.0 transmitter.

## The security angle, handled responsibly

Inside one such user-defined table, the datacasting structure included
**operational backend detail that should not be broadcast in the clear** —
session-scoped credential-like fields, in a standard UUID format, sitting in
plain XML. Because that XML rides LLS, it is recoverable by **any** ATSC 3.0
receiver in coverage: no decryption, no special access, a commodity SDR and an
open decoder reproduce it in seconds.

We did the boring, correct thing:

- We **never used, tested, or transmitted to** any of it. We only observed the
  signaling a receiver reads anyway.
- We reported it **privately**, through the operator's published security
  contact (the RFC 9116 `security.txt` channel), with the values withheld from
  the email and offered only over a secure channel.
- We set **no deadline** and published nothing until it was closed.

The operator responded promptly and cooperatively. Their assessment was that
these particular fields are **low-scope session tokens, unique per receiver
session and by design**, tied to a limited device-reporting API rather than to
any backend system, data, or privilege — and that the platform carrying them is
already scheduled for replacement. They confirmed the captured data could be
deleted and that no further action was needed.

We take that assessment at face value, and we want to be epistemically honest
about our own side of it: **we deliberately never probed the tokens**, so we
can't independently rank their impact — we defer to the owner on that. What we
*can* say with certainty is structural, and it does not depend on how much any
one token was worth:

> Anything placed in LLS is public the instant it is on the air. The fix for
> that class of exposure is never "encrypt the table" (you can't) — it is to
> assume the table is world-readable and put nothing in it that isn't safe to
> be. Scope any broadcast token to be useless on its own: short-lived,
> read-only, rate-limited, and worthless without a second factor the air
> doesn't carry.

That is the durable, vendor-neutral takeaway, and it is why this note names no
operator, no callsign, no endpoint, and no key: the *mechanism* is the science,
and the mechanism is what generalizes to every LLS field on every transmitter.

## Reproducing the class (not the instance)

You don't need our capture to learn the general lesson — you need any ATSC 3.0
signal and a decoder that surfaces LLS. Point a receiver at a multiplex, dump
the LLS table set, and read what's there. The standardized tables tell you the
services; table id 255, where present, tells you whatever that station chose to
add. Everything you can read, so can anyone else with an antenna.

The point of this note is not the one deployment we happened to see. It is the
invariant behind it, which every operator putting data on NextGen TV inherits
whether they think about it or not:

> **LLS is public by construction. Treat every field placed in it as
> world-readable.**

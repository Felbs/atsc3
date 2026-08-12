# Field observation: commercial datacasting on ATSC 3.0 (BroadSpan / Cast.era)

*Publishable. Contains no credentials. The specific keys observed are handled
privately under responsible disclosure and are not reproduced here.*

## What we saw

ATSC 3.0 reserves a Low Level Signaling table (`LLS_table_id 255`) for
broadcaster-defined use. On two Baltimore-market transmitters — WBFF's RF25
multiplex and one adjacent RF30 multiplex — that table carries a
**datacasting** overlay: a `UDST` (user-defined signaling table) describing
`broadSpanServices`, the commercial datacasting product Sinclair operates
(marketed as *BroadSpan*, technology by *Cast.era*).

This is one of NextGen TV's headline business models seen in the wild:
the station rents spare broadcast bandwidth to carry **data, not television** —
file distribution, software updates, and edge delivery to any ATSC 3.0
receiver (cars, phones, set-tops), offloading traffic from cellular and the
internet.

It is **not** present on the Washington DC multiplex (RF33) examined the same
way — consistent with datacasting being a per-host-station deployment.

## Why it is interesting to a receiver project

- It is direct evidence that the non-television side of ATSC 3.0 is live and
  carrying real commercial traffic, not just a spec promise.
- The signaling rides the same LLS path a receiver already decodes to find its
  TV services, so a general ATSC 3.0 receiver sees it for free — a reminder
  that "the signal" is much more than the picture.

## A security note, handled responsibly

The user-defined table was observed to contain operational backend detail that
should not be broadcast in the clear. Because LLS is unencrypted by design
(a receiver must read it to bootstrap), anything placed there is public to
every receiver in coverage. That specific issue was reported privately to the
operator through their published security contact; details are withheld here
until it is resolved. The general lesson is worth stating plainly:

> **LLS is public by construction. Treat every field placed in it as world-readable.**

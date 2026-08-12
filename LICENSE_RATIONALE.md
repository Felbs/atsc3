# Licence proposal and reasoning

reject or modify. It is engineering reasoning about licence choice, not legal
advice, and none of it has been reviewed by a lawyer.**

---

## The proposal

| | |
|---|---|
| **Licence** | **Apache License 2.0** (`SPDX-License-Identifier: Apache-2.0`) |
| **Applies to** | this repository's own code: the physical-layer receiver, the transport and media layers, the AC-4 decoder, the tools, the gates |
| **Does not apply to** | the ATSC and ETSI standards documents (not ours), captured broadcast content (not ours), any third-party code that may later be vendored |
| **Companion files** | `LICENSE` (verbatim), `NOTICE` (attribution + the patent and trademark statements) |

`release/LICENSE` is the canonical Apache-2.0 text, copied **verbatim** from a
known-good source and diffed against a second copy — not retyped. A licence with
a transcription error in it is not the licence you think you granted. Its
SHA-256 is `8173d5c2…707c90`; re-verify it against apache.org before publishing.

---

## Why Apache-2.0

### 1. The patent grant is the whole reason

This is the argument that decides it. **ATSC 3.0 and AC-4 are heavily
patent-encumbered technologies.** ATSC 3.0 has multiple patent pools; AC-4 is a
Dolby technology with its own licensing programme. Anyone who adopts this code
is walking into that landscape.

Apache-2.0 §3 grants an explicit patent licence from each contributor covering
their own contributions, and §3 terminates that grant for anyone who sues over
patents in the work. MIT and BSD grant patents only by implication, if at all —
which in a field this encumbered is precisely the wrong kind of silence.

**What it does not do, and this must be said out loud in `NOTICE`:** a
contributor can only license patents they hold. Nobody here holds any ATSC or
Dolby patents, so **§3 grants the reader nothing about third-party patents.**
Implementing a broadcast standard may require licences this project does not
have and cannot convey. Apache's grant makes *our* position clear and is honest
about its own limits; that is the best any licence can do here.

### 2. Trademarks are addressed explicitly (§6)

Apache-2.0 §6 says plainly that the licence grants no trademark rights. That
matters a lot for this project, which has already renamed a component once for
exactly this reason:

* **"AC-4" and "Dolby"** are Dolby Laboratories marks. The standalone decoder is
  named `ac4-decoder` — descriptive, nominative use, saying what format it
  decodes. **"Dolby" must not appear in any package, repo, or module name**, and
  the earlier `dolbyTuna` naming was correctly abandoned. See
  [AC4_SPLIT_PLAN.md](AC4_SPLIT_PLAN.md) for the naming discussion, including
  the open question of whether even `ac4-decoder` should carry a disclaimer in
  its package description.
* **"ATSC"** and **"NextGen TV"** are marks of ATSC and its licensors.
  `atsc3` as a repository name is descriptive use of a standard's designation;
  "NextGen TV" appears only in the sentence explaining what ATSC 3.0 is called
  in the marketplace, which is the safest possible use of it.
* **The commercial tuner in the scoreboard** is named as a measurement
  instrument. Nominative fair use — but the comparison must stay factual,
  measured, and must not imply endorsement or disparagement. The scoreboard says
  we lose three dimensions to one, which is not a comparative-advertising risk
  so much as the opposite.

The `NOTICE` draft below carries all of this.

### 3. Permissive matches the actual goal

The standing bar for this project is *a stranger installs it and watches
television in under 15 minutes*. That goal is adoption, not reciprocity.
Apache-2.0 lets a hobbyist vendor it, a distribution package it, and a
manufacturer evaluate it without a legal review that ends in "no".

### 4. It preserves downstream options

Apache-2.0 is one-way compatible with GPL-3.0: this code can be incorporated
into a GPLv3 project later. Choosing GPL now would foreclose the permissive
direction forever, and nothing about the project requires it.

### 5. It has a contribution clause (§5) and a NOTICE mechanism (§4d)

An inbound=outbound default for contributions, without needing a CLA, and a
standard place to put the patent and trademark statements this project actually
needs to make.

---

## Alternatives considered

| licence | why not |
|---|---|
| **MIT** | Simplest, and it is what one prior-art project uses — but **no patent grant and no defensive termination**, in a field defined by patents. The simplicity is bought with exactly the silence that matters here. |
| **BSD-2/3** | Same objection as MIT. |
| **GPL-3.0** | Would make combining with the GPL-3.0 prior-art SDR blocks straightforward, and would force reciprocity. Rejected for two reasons: it works against adoption, and **we did not read that code** — this implementation is from-spec, and its clean-room character is an asset worth keeping visible. Taking a GPL licence "just in case" would blur a distinction that was deliberately maintained. |
| **LGPL-3.0** | The library/application split does not map onto this codebase, which is a receiver, not a library. |
| **Dual MIT/Apache** | The Rust convention, and defensible. Adds explanatory burden for no benefit here — the reason for Apache is the patent clause, and offering MIT alongside lets an adopter opt out of the one clause that matters. |
| **Non-commercial / research-only** | Not open source, would put the project outside every distribution channel, and would not actually reduce anyone's patent exposure. |

---

## The parts that are NOT ours, and how to keep that straight

This is the section that matters more than the licence choice, because the
licence only governs what we wrote.

### The standards documents

`lab/spec/` and `spec/` hold extracted text and PDFs of ATSC A/322, A/331 and
ETSI TS 103 190-1. **They are gitignored and must stay that way.** Both bodies
publish them free of charge — which is what makes this project possible at all —
and neither grants redistribution rights. ETSI returns 403 to automated
fetching, and working around an access control is not on the table.

**The open question** (recorded in `PORT.md` as a pending decision, and it blocks
the 15-minute install goal):

> May the *parsed numeric tables* be committed instead of the documents?

The argument for: an LDPC parity-check matrix and a constellation table are
**facts**, and facts are weak candidates for copyright in themselves. Committing
the parsed arrays would make a clone self-contained and delete the worst step in
`INSTALL.md`.

The argument for caution: the *selection, arrangement and presentation* in a
standards document can attract protection even where individual numbers do not,
and both bodies assert copyright in the documents as a whole. The distance
between "these integers" and "this table as the standard lays it out" is a real
one and I am not qualified to place the line.

**Recommendation: ask, don't infer.** Both ATSC and ETSI have contacts for
exactly this question, and a written answer converts a licensing risk into a
citable permission. Until then the current arrangement — ship the *parser*, not
the *parsed* — is the conservative and defensible position, and `INSTALL.md`
says so plainly rather than hiding the cost.

### Captured broadcast content

Everything the receiver pulls off the air is somebody else's copyrighted work.
Eight such files are currently tracked in git (video, frames, programme guides).
**They must be untracked before publication** — see the project scrub audit. Private
research capture and publication are different acts, and no licence we choose
changes that.

### Prior art

This implementation is from-spec. `PRIOR_ART.md` records what else exists,
including a GPL-3.0 SDR implementation and an MIT-licensed library. **Nothing
has been copied from either.** If anything is ever vendored:

* from the MIT project — fine, retain its notice in `NOTICE`;
* from the GPL-3.0 project — **stop.** It would relicense the combined work and
  destroy the clean-room position. That is a decision to take deliberately in
  daylight, if ever, and never as a convenience.

### The third-party API key

the project scrub audit found a broadcaster's own session key inside captured
signalling that is currently tracked. That is not a licensing matter, it is a
"delete it" matter.

---

## Draft `NOTICE`

To be reviewed and adjusted; the copyright holder line is the operator's to fill.

```
atsc3 — an open-source ATSC 3.0 (NextGen TV) receiver
Copyright 2026 <copyright holder>

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

--------------------------------------------------------------------------
PATENTS

This software implements published broadcast standards. Those standards are
covered by patents held by third parties, and this project holds none of
them. The patent grant in Section 3 of the License extends only to patent
claims owned by the contributors to THIS work; it grants nothing with
respect to any third-party patent.

Implementing, distributing, or using an implementation of these standards
may require patent licences that this project does not have and cannot
convey. If you intend to use this software commercially, or to ship it in a
product, get your own advice. Nothing here is a representation that you are
free to do so.

--------------------------------------------------------------------------
TRADEMARKS

Section 6 of the License grants no trademark rights.

  "ATSC" and "NextGen TV" are trademarks of the Advanced Television Systems
  Committee, Inc. and/or its licensors.
  "AC-4" and "Dolby" are trademarks of Dolby Laboratories.

  All other product and company names mentioned are the trademarks of their
  respective owners.

They are used here only descriptively — to state which published standards
this software implements and which equipment it was measured against. This
project is not affiliated with, endorsed by, sponsored by, or certified by
any of them, and no such relationship should be inferred.

--------------------------------------------------------------------------
STANDARDS DOCUMENTS

This software reads tables from published standards documents at runtime.
Those documents are NOT included in this distribution and are not ours to
redistribute. They are published free of charge by their respective bodies;
see INSTALL.md for how to obtain them.

--------------------------------------------------------------------------
ENCRYPTED SERVICES

This software does not decrypt, circumvent, or attempt to defeat any content
protection system. Encrypted services are detected, enumerated and labeled as
locked, and that behaviour is enforced by a test in the build.
```

---

## Before publishing — the licence checklist

1. Decide the copyright holder line. A legal name is the usual choice; a
   pseudonym or handle is legally weaker but consistent with the rest of the
   scrub decisions. **This interacts with the privacy work — it is the one place
   where publishing a real name may be the technically correct choice and the
   operator may not want it.** Worth thinking about before, not after.
2. Verify `LICENSE` byte-for-byte against apache.org.
3. Add `NOTICE` at the repository root, with §4(d) of the licence in mind — it
   travels with derivative works.
4. Add `SPDX-License-Identifier: Apache-2.0` headers to source files, or state
   once in the README that the whole repository is under it. (The short SPDX
   header is preferred; it survives file-copying, which a README does not.)
5. Set the licence field in any packaging metadata, and check that GitHub's
   licence detector recognises the file.
6. Confirm the untracking in the project scrub audit is done. **A licence file does
   not launder other people's video.**
7. Decide the standards-table question above, or ship with the manual step and
   say why.

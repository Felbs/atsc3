#!/usr/bin/env python3
"""scrub_audit -- find every identifier in this repo that would locate or
identify the operator, BEFORE anything of it becomes public.

E83.  This is a GATE, not a grep.  A one-time grep proves nothing about the
next commit; a gate you can re-run, that proves on every run that it still
catches a planted identifier, is evidence.  So:

    python tools/scrub_audit.py                # scan the tracked tree
    python tools/scrub_audit.py --selftest     # NEGATIVE CONTROL: plant one
                                               # specimen per rule in a scratch
                                               # copy and prove each rule fires
    python tools/scrub_audit.py --selftest --keep-scratch   # leave it to look at
    python tools/scrub_audit.py --json OUT.json
    python tools/scrub_audit.py --severity HIGH   # only HIGH and FATAL

Exit codes:
    0   no findings at or above --fail-on (default FATAL)
    1   findings at or above --fail-on
    2   the selftest failed -- a rule did NOT fire on its own planted
        specimen, i.e. the instrument is broken and a clean run means nothing

DISPOSITIONS.  Not everything that names a place is a leak.  A receiver
project's decode gates need REAL captures off REAL transmitters, and a gate
that cannot say which transmitter is not a gate.  So each rule carries a
proposed disposition, and the operator decides:

    SCRUB       remove or redact before publish; no technical value that
                survives the removal
    PARAMETERIZE keep the capability, move the value out of the source into
                a flag / env var / config file with a neutral default
    CONSENT     legitimate technical evidence; publishing it is a deliberate
                choice the operator makes knowing it narrows their location

SELF-MATCH.  This file excludes itself from its own scan by default (the
fleet's pgrep/self-match law: an unfiltered census matches the censor).
--include-self overrides, and the report always says which mode it ran in.

THE RULE TABLE CARRIES NO REAL IDENTIFIERS, AND THIS IS LOAD-BEARING.  The
first draft of this file hardcoded the operator's hostname, login, device
serial, LAN address, antenna nickname and local station call signs as literal
patterns -- which made the AUDIT TOOL a tidy, sorted, published index of
exactly the things it exists to find.  An audit report is a dossier; an audit
tool with a literal table in it is a worse one, because the tool definitely
ships.  So every rule here is a SHAPE:

    /home/<any-login>        not /home/<the-login>
    HDHR_ID = "<8 hex>"      not HDHR_ID = "<the-serial>"
    <callsign> near <n.n>    not a list of the local stations

Site-specific literals live in `tools/scrub_local.json`, which is GITIGNORED
and never published.  With it, the audit is exact for this installation;
without it, the shapes still catch the classes.  The report says which mode it
ran in and how many local literals loaded -- the COUNT, never the values.

tools/scrub_local.json (gitignored) -- optional, per-installation literals:

    {
      "PRIVATE_REPO": ["some-private-repo"],
      "HOST_USER":    ["my-desktop", "mylogin"],
      "SITE_ANTENNA": ["The Big Yagi"],
      "CALLSIGN":     ["WXYZ", "WABC"],
      "LOCAL_PATH":   ["sibling-repo-name"],
      "BSID":         ["540"],
      "DEVICE_ID":    ["0123ABCD"]
    }

Each list is OR-ed into that rule as literal alternatives.  Values are escaped,
so they are matched literally and cannot inject regex.

ALLOWLIST.  A line ending in `# scrub-allow: RULE_ID reason` is suppressed for
that rule only, and the suppression is COUNTED and REPORTED -- an allowlist you
cannot see is an allowlist nobody reviews.
"""
from __future__ import annotations

import argparse
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
from dataclasses import dataclass, field
from pathlib import Path

HERE = Path(__file__).resolve()
REPO = HERE.parent.parent

SEVERITY_ORDER = {"INFO": 0, "MEDIUM": 1, "HIGH": 2, "FATAL": 3}


@dataclass
class Rule:
    rid: str
    severity: str
    disposition: str
    why: str
    fix: str
    pattern: str = ""
    path_pattern: str = ""      # matches the PATH, not the content
    specimen: str = ""          # what --selftest plants to prove the rule works
    specimen_path: str = ""     # filename the specimen must live under, if any
    flags: int = 0
    rx: re.Pattern = field(init=False, default=None)
    prx: re.Pattern = field(init=False, default=None)

    def __post_init__(self):
        self.rx = re.compile(self.pattern, self.flags) if self.pattern else None
        self.prx = re.compile(self.path_pattern, re.I) if self.path_pattern else None


# --------------------------------------------------------------------------
# THE RULES.  Ordered roughly by how badly a hit locates the operator.
# --------------------------------------------------------------------------
RULES: list[Rule] = [
    Rule(
        rid="PRIVATE_REPO",
        severity="FATAL",
        disposition="SCRUB",
        why="Names a private sibling repository. Some of these are SECRET by "
            "standing rule -- they must never appear in anything publishable, not "
            "even as a path in a comment. There is no generic SHAPE for a project "
            "name, so this rule is driven ENTIRELY by scrub_local.json and fires "
            "on nothing without it.",
        fix="Delete the reference. If the sentence needs a source, describe the "
            "tool ('the capture tool we use') without naming the repository.",
        pattern="",                       # local literals only -- by design
        specimen="",                      # supplied by the selftest's local file
    ),
    Rule(
        rid="HOST_USER",
        severity="FATAL",
        disposition="SCRUB",
        why="A hostname or login name is a direct identity handle -- it ties this "
            "repo to a named machine and a named person, and it correlates across "
            "every other public repo that carries the same string.",
        fix="Replace with a neutral placeholder: 'the Linux box', USER, <user>.",
        # SHAPES, not names. A GitHub handle used as a handle is fine -- the repo
        # is published under that account, so it is public by construction and
        # the noreply commit address is the CORRECT one. What must go is a handle
        # used as a LOGIN (user@host, /home/user) and every user-profile path.
        pattern=r"(?:ssh://)?\b[a-z_][a-z0-9_.-]{1,31}@(?:\d{1,3}\.){3}\d{1,3}\b"
                r"|/home/[a-z_][a-z0-9_.-]*"
                r"|[A-Za-z]:\\+Users\\+[A-Za-z0-9_.-]+"
                r"|/c/Users/[A-Za-z0-9_.-]+"
                r"|(?<![A-Za-z0-9])/Users/[A-Za-z0-9_.-]+",
        flags=re.I,
        specimen=r"ssh someone@10.0.0.5 and C:\Users\someone\Desktop",
    ),
    Rule(
        rid="EMAIL",
        severity="FATAL",
        disposition="SCRUB",
        why="A personal email address is the single strongest cross-service identity "
            "link. The fleet already has one pending history rewrite for exactly this.",
        fix="Remove, or use the GitHub noreply address.",
        # noreply addresses are the intended-public ones (GitHub's, and the
        # assistant trailer's) -- flagging them buries the one that matters.
        pattern=r"(?!noreply@)[A-Za-z0-9._%+-]+@(?!users\.noreply\.github\.com)"
                r"(?!anthropic\.com)[A-Za-z0-9.-]+\.(com|net|org|edu|io|dev)\b",
        specimen="astro.example.person@examplemail.com",
    ),
    Rule(
        rid="LAN_ADDR",
        severity="HIGH",
        disposition="PARAMETERIZE",
        why="A private LAN address is not routable, so it does not locate the "
            "operator by itself -- but it exposes the internal network layout and, "
            "paired with a device ID, makes a fingerprint. It is also simply WRONG "
            "for a stranger: their HDHomeRun is not at this address.",
        fix="Read from env/flag with no default, or discover the device. "
            "ATSC3_HDHR_HOST / --hdhr-host.",
        pattern=r"\b(192\.168\.\d{1,3}\.\d{1,3}|10\.\d{1,3}\.\d{1,3}\.\d{1,3}|"
                r"172\.(1[6-9]|2\d|3[01])\.\d{1,3}\.\d{1,3})\b",
        specimen="HDHR = \"192.168.0.99\"",
    ),
    Rule(
        rid="DEVICE_ID",
        severity="HIGH",
        disposition="PARAMETERIZE",
        why="The HDHomeRun device ID is a globally unique serial. It is queryable "
            "through the vendor's discovery service and is a durable identifier for "
            "one household's hardware.",
        fix="Env var with NO default (ATSC3_HDHR_ID), and auto-discovery when unset.",
        pattern=r"\b(?:hdhr|hdhomerun|device|tuner)[_-]?id\b\s*[:=(,]\s*"
                r"[\"']?[0-9A-F]{8}\b"
                r"|[\"']ATSC3_HDHR_ID[\"']\s*,\s*[\"'][0-9A-F]{8}[\"']"
                r"|\bdevice\s+[0-9A-F]{8}\b",
        flags=re.I,
        specimen='HDHR_ID = os.environ.get("ATSC3_HDHR_ID", "0123ABCD")',
    ),
    Rule(
        rid="GEO",
        severity="FATAL",
        disposition="SCRUB",
        why="Coordinates, a 6-character Maidenhead grid, or a street/ZIP is a direct "
            "location. Nothing in a receiver needs one.",
        fix="Delete. If a propagation claim needs a location, give the distance to "
            "the transmitter, not the position of the receiver.",
        pattern=r"\b[Ff][MNL]\d{2}[a-x]{2}\b|"
                r"\b(lat|latitude|lon|lng|longitude)\s*[:=]\s*-?\d{1,3}\.\d{3,}|"
                r"-?\b3[89]\.\d{4,}\s*,\s*-7[5-9]\.\d{4,}\b|"
                r"\bQTH\b",
        specimen="QTH FM19ab, latitude: 38.90000",
    ),
    Rule(
        rid="LOCAL_PATH",
        severity="MEDIUM",
        disposition="SCRUB",
        why="A hardcoded absolute path exposes the operator's drive layout and the "
            "names of sibling projects, and it is dead on arrival for anyone else.",
        fix="Make it relative to the repo root, or a flag. Keep sibling-repo names "
            "out entirely.",
        # Shape only: an absolute drive path. Sibling PROJECT NAMES have no shape,
        # so they come from scrub_local.json.
        pattern=r"\b[A-Za-z]:[\\/]{1,2}[A-Za-z0-9_.-]+[\\/]",
        specimen=r"python Z:\some-other-project\iq_capture.py --rf 27",
    ),
    Rule(
        rid="SITE_ANTENNA",
        severity="MEDIUM",
        disposition="CONSENT",
        why="Antenna model, mounting location (attic/shed/roof) and named ports "
            "describe the physical installation. Any one is harmless; the set of "
            "them plus a receivable-station list describes a specific house.",
        fix="Keep the ELECTRICAL facts (yagi, gain, splitter loss -- they carry the "
            "engineering) and drop the SITE facts (which room, which named port, "
            "the nickname). 'a UHF yagi' not '<nickname> on port B in the attic'.",
        # Shape: mounting locations and named switch ports. Antenna NICKNAMES and
        # model numbers have no shape -- scrub_local.json.
        pattern=r"\b[Pp]ort [ABCD]\b|\b(?:attic|shed|rooftop|roofline|crawlspace|"
                r"garage|balcony|back yard|backyard)\b",
        specimen="Antenna port B, mounted in the attic",
    ),
    Rule(
        rid="CALLSIGN",
        severity="MEDIUM",
        disposition="CONSENT",
        why="A broadcast call sign names one transmitter in one market. A handful of "
            "them narrows the operator to a metro area; the full receivable list "
            "narrows it much further, because reception is a function of position.",
        fix="Decide per file. In DECODE GATES the call sign is legitimate evidence "
            "(a gate that cannot name its specimen is not a gate) -- CONSENT. In "
            "prose that does not need it, generalise to 'the MMTP service' -- SCRUB. "
            "Never publish a COMPLETE receivable lineup; that is a position fix.",
        # Shape: a US call-sign-shaped token that appears near a virtual channel
        # number or a broadcast suffix -- which is how a call sign is actually
        # written in this repo. A bare \b[KW][A-Z]{3}\b would fire on every
        # acronym in a DSP codebase, so the context requirement IS the rule.
        # Known local call signs go in scrub_local.json for an exact match.
        pattern=r"\b[KW][A-Z]{2,3}\b(?=[^\n]{0,40}\b\d{1,3}\.\d\b)"
                r"|\b[KW][A-Z]{2,3}-(?:TV|DT|HD|LD|CD|LP)\b",
        specimen="EXAM 7.1 and WXYZ 32.1 decode; WXYZ-TV is locked",
    ),
    Rule(
        rid="BSID",
        severity="MEDIUM",
        disposition="CONSENT",
        why="A bsid is the globally unique ID of one ATSC 3.0 broadcast facility. "
            "It is transmitted in the clear and anyone in range reads it -- but it "
            "names the transmitter exactly, so publishing 'we decode bsid 540' says "
            "'we are within about 60 km of that tower'.",
        fix="CONSENT for gate fixtures (the L1 values only mean something against a "
            "named bsid). PARAMETERIZE where it is only a default.",
        pattern=r"\bbsid\b[\"'=\s:(,]{1,4}\d{1,5}\b",
        flags=re.I,
        specimen='<SLT bsid="1234"',
    ),
    Rule(
        rid="MCAST_SLS",
        severity="MEDIUM",
        disposition="CONSENT",
        why="239.255.x.x SLS destinations are chosen by the BROADCASTER and carried "
            "over the air, so they are not the operator's addresses -- but the set "
            "of them is a fingerprint of one multiplex, i.e. one market.",
        fix="CONSENT in decode evidence. In library code these must be DISCOVERED "
            "from the SLT, never hardcoded -- a receiver that only works on one "
            "multiplex is not a receiver. PARAMETERIZE any that are defaults.",
        pattern=r"\b239\.255\.\d{1,3}\.\d{1,3}\b|\b239_255_\d{1,3}_\d{1,3}\b",
        specimen="obj_239_255_32_1_8321_tsi20.bin at 239.255.32.1",
    ),
    Rule(
        rid="RF_MARKET",
        severity="INFO",
        disposition="CONSENT",
        why="An RF channel number or a centre frequency is only a market hint on its "
            "own -- the same channel is in use nationwide. The LIST of channels that "
            "carry ATSC 3.0 near this receiver is the fingerprint, not any one entry.",
        fix="Keep as measured evidence, but never ship a hardcoded --rf DEFAULT: a "
            "stranger's RF 33 is a different station. Require the flag.",
        pattern=r"--rf[\"'\s=]+\d{1,2}\b"
                r"|--rf\b[^\n]{0,40}\bdefault\s*=\s*\d{1,2}\b"
                r"|\brf\s*=\s*\d{1,2}\b"
                r"|\bRF\s?(?:1[4-9]|[2-9]\d)\b"
                r"|\b(?:4[7-9]\d|[5-8]\d\d)(?:\.\d+)?\s*MHz\b",
        specimen='ap.add_argument("--rf", type=int, default=33)',
    ),
    Rule(
        rid="BROADCAST_MEDIA",
        severity="HIGH",
        disposition="SCRUB",
        why="A tracked video/still/guide artifact is COPYRIGHTED BROADCAST CONTENT "
            "that we do not own, and its content (station bug, local news, the "
            "programme grid with times) fingerprints market AND date AND what was "
            "on this receiver's screen. This is the worst ratio of evidence value "
            "to exposure in the repo.",
        fix="Untrack. Replace proof-of-picture with a NUMBER the gate already "
            "computes (flat_frac, decoded frame count, SHA of the elementary "
            "stream), or a heavily cropped/downscaled still with consent. Keep the "
            "originals in the gitignored lab/ tree.",
        path_pattern=r"\.(mp4|m4s|ts|png|jpg|jpeg|gif|wav|sup|srt)$",
        specimen_path="e99_out/station_frame.png",
        specimen="\x00binary-broadcast-frame",
    ),
    Rule(
        rid="BROADCAST_SIGNALLING",
        severity="MEDIUM",
        disposition="CONSENT",
        why="Captured LLS/SLT/ESG XML and tuner lineup dumps are verbatim off-air "
            "signalling. They are the best possible spec evidence AND a complete "
            "market fingerprint in one file.",
        fix="CONSENT for one small SLT as a parser fixture, with the operator "
            "knowing it names the market. SCRUB the full-lineup and full-guide "
            "dumps -- a complete receivable lineup is close to a position fix.",
        path_pattern=r"(lls|slt|esg|sgdu|lineup|guide|grid_now|grid_guide|flex_probe)"
                     r".*\.(xml|json|html)$",
        specimen_path="e99_out/e99_hdhr_lineup.json",
        specimen='{"GuideNumber": "2.1", "GuideName": "EXAMPLE-HD"}',
    ),
    Rule(
        rid="SECRET_MATERIAL",
        severity="FATAL",
        disposition="SCRUB",
        why="Keys, tokens and DRM material must never be published -- and in this "
            "project's case a decryption key would also contradict the standing "
            "scope limit that encrypted services are labeled locked, never attacked.",
        fix="Delete and rotate. If a DRM field must be discussed, describe the "
            "FIELD, never a value.",
        pattern=r"(?i)\b(api[_-]?key|secret[_-]?key|access[_-]?token|password|"
                r"BEGIN [A-Z ]*PRIVATE KEY|widevine[_-]?key|content[_-]?key)\b"
                r"\s*[:=]\s*[\"']?[A-Za-z0-9/+_-]{12,}",
        specimen='api_key = "AbCdEf0123456789xyz"',
    ),
]

RULES_BY_ID = {r.rid: r for r in RULES}

LOCAL_DEFAULT = HERE.parent / "scrub_local.json"


def apply_local(rules: list[Rule], path: Path) -> tuple[int, list[str]]:
    """OR per-installation literals into the shape rules.

    Returns (count, rule_ids) -- the COUNT of literals loaded, never the
    literals themselves. This function is the only place in the tool that ever
    holds a real identifier, and nothing it loads is echoed to the report.
    """
    if not path.exists():
        return 0, []
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except Exception as e:
        print("WARNING: %s could not be parsed (%s); running on SHAPES ONLY"
              % (path, e))
        return 0, []
    n, touched = 0, []
    by_id = {r.rid: r for r in rules}
    for rid, literals in data.items():
        r = by_id.get(rid)
        if r is None or not literals:
            continue
        alts = "|".join(re.escape(str(x)) for x in literals if str(x).strip())
        if not alts:
            continue
        extra = r"(?:%s)" % alts
        r.pattern = (r.pattern + "|" + extra) if r.pattern else extra
        r.rx = re.compile(r.pattern, r.flags)
        n += len(literals)
        touched.append(rid)
    return n, sorted(set(touched))

SKIP_DIRS = {".git", "__pycache__", "node_modules", ".venv"}
# Accepts any comment syntax -- `# scrub-allow: RULE`, `// scrub-allow: RULE`,
# `<!-- scrub-allow: RULE -->`. A marker that only works in Python would push
# people to leave markdown and C unsuppressible, and an allowlist nobody can
# use is an allowlist nobody uses.
ALLOW_RX = re.compile(r"scrub-allow:\s*([A-Z_]+)")


@dataclass
class Finding:
    rid: str
    severity: str
    disposition: str
    path: str
    line: int
    text: str


def _git_ls(root: Path, *args: str) -> list[str]:
    out = subprocess.run(["git", "-C", str(root), "ls-files", *args],
                         capture_output=True, text=True)
    if out.returncode != 0:
        return []
    return [p for p in out.stdout.splitlines() if p.strip()]


def git_tracked(root: Path) -> list[str]:
    return _git_ls(root, "--cached")


def walk_all(root: Path) -> list[str]:
    """Tracked PLUS untracked-but-not-ignored.

    This is exactly the set a careless `git add -A` would sweep into a commit,
    which is the set that matters before a publish.  It deliberately does NOT
    walk the filesystem: lab/ and data/ hold ~200 GB of gitignored captures and
    an os.walk over them would take longer than the audit is worth -- and they
    are, by being ignored, not candidates for publication anyway.  If the
    ignore rules themselves are what you doubt, audit .gitignore, not the disk.
    """
    seen, paths = set(), []
    for p in _git_ls(root, "--cached", "--others", "--exclude-standard"):
        if p not in seen:
            seen.add(p)
            paths.append(p)
    return paths


def is_binary(blob: bytes) -> bool:
    return b"\x00" in blob[:8192]


def scan_file(root: Path, rel: str, rules: list[Rule]) -> tuple[list[Finding], dict]:
    findings: list[Finding] = []
    allowed: dict[str, int] = {}
    full = root / rel

    # path rules first -- they fire on binaries too
    for r in rules:
        if r.prx and r.prx.search(rel):
            findings.append(Finding(r.rid, r.severity, r.disposition, rel, 0,
                                    "<path matches %s>" % r.path_pattern))

    try:
        blob = full.read_bytes()
    except OSError:
        return findings, allowed
    if is_binary(blob):
        return findings, allowed

    text = blob.decode("utf-8", errors="replace")
    for lineno, line in enumerate(text.splitlines(), 1):
        am = ALLOW_RX.search(line)
        allow_id = am.group(1) if am else None
        for r in rules:
            if not r.rx:
                continue
            m = r.rx.search(line)
            if not m:
                continue
            if allow_id == r.rid:
                allowed[r.rid] = allowed.get(r.rid, 0) + 1
                continue
            snippet = line.strip()
            if len(snippet) > 160:
                start = max(0, m.start() - 60)
                snippet = "..." + line[start:m.end() + 60].strip() + "..."
            findings.append(Finding(r.rid, r.severity, r.disposition,
                                    rel, lineno, snippet))
    return findings, allowed


def scan(root: Path, paths: list[str], rules: list[Rule],
         include_self: bool) -> tuple[list[Finding], dict]:
    findings: list[Finding] = []
    allowed: dict[str, int] = {}
    self_rel = None
    try:
        self_rel = str(HERE.relative_to(root)).replace("\\", "/")
    except ValueError:
        pass
    for rel in paths:
        if not include_self and self_rel and rel == self_rel:
            continue
        f, a = scan_file(root, rel, rules)
        findings.extend(f)
        for k, v in a.items():
            allowed[k] = allowed.get(k, 0) + v
    return findings, allowed


def scan_history(root: Path, rules: list[Rule], min_sev: str) -> list[Finding]:
    """Audit every blob that ever existed, not just the checkout.

    Publishing a repo publishes its HISTORY.  A file scrubbed in the working
    tree but still reachable from an old commit is still public the moment the
    repo is, and this is precisely how the fleet's one known email leak
    survived a clean worktree.  So the pre-publish gate has to read the object
    database, and a fix here means rewriting history (filter-repo), not an
    edit-and-commit.
    """
    out = subprocess.run(["git", "-C", str(root), "rev-list", "--objects", "--all"],
                         capture_output=True, text=True)
    if out.returncode != 0:
        return []
    findings: list[Finding] = []
    seen_pairs: set[tuple[str, str, str]] = set()
    active = [r for r in rules
              if r.rx and SEVERITY_ORDER[r.severity] >= SEVERITY_ORDER[min_sev]]

    # Commit METADATA is published too, and it is not a blob: author name and
    # email, committer, and the message body.  This project writes long
    # narrative commit messages, which is exactly where a hostname goes to hide.
    meta = subprocess.run(
        ["git", "-C", str(root), "log", "--all", "--format=%H%x01%an%x01%ae%x01%cn%x01%ce%x01%B%x02"],
        capture_output=True, text=True, errors="replace")
    for rec in meta.stdout.split("\x02"):
        rec = rec.strip("\n")
        if not rec.strip():
            continue
        bits = rec.split("\x01")
        if len(bits) < 6:
            continue
        sha, an, ae, cn, ce, body = bits[0], bits[1], bits[2], bits[3], bits[4], bits[5]
        for label, blob_txt in (("author/committer", "%s <%s> / %s <%s>" % (an, ae, cn, ce)),
                                ("message", body)):
            for r in active:
                m = r.rx.search(blob_txt)
                if not m:
                    continue
                key = (r.rid, "commit-" + label, m.group(0))
                if key in seen_pairs:
                    continue
                seen_pairs.add(key)
                snippet = blob_txt[max(0, m.start() - 40):m.end() + 40].replace("\n", " ")
                findings.append(Finding(r.rid, r.severity, r.disposition,
                                        "<commit %s %s>" % (sha[:8], label), 0,
                                        snippet.strip()))
    for line in out.stdout.splitlines():
        parts = line.split(" ", 1)
        if len(parts) != 2:
            continue          # a commit or a root tree: no path, nothing to read
        sha, path = parts
        blob = subprocess.run(["git", "-C", str(root), "cat-file", "-p", sha],
                              capture_output=True)
        if blob.returncode != 0 or is_binary(blob.stdout):
            continue
        text = blob.stdout.decode("utf-8", errors="replace")
        for r in active:
            m = r.rx.search(text)
            if not m:
                continue
            snippet = text[max(0, m.start() - 40):m.end() + 40].replace("\n", " ")
            key = (r.rid, path, m.group(0))
            if key in seen_pairs:
                continue      # same secret, same path, many revisions: report once
            seen_pairs.add(key)
            findings.append(Finding(r.rid, r.severity, r.disposition,
                                    "%s @%s" % (path, sha[:8]), 0, snippet.strip()))
    return findings


def report(findings: list[Finding], allowed: dict, root: Path, paths: list[str],
           include_self: bool, min_sev: str, show: int) -> str:
    lines = []
    A = lines.append
    A("=" * 78)
    A("SCRUB AUDIT -- %s" % root)
    A("=" * 78)
    A("files scanned : %d %s" % (len(paths),
      "(self INCLUDED)" if include_self else "(this script EXCLUDED from its own "
      "scan -- it holds every specimen by construction; --include-self to override)"))
    A("minimum severity reported: %s" % min_sev)
    A("")

    by_rule: dict[str, list[Finding]] = {}
    for f in findings:
        if SEVERITY_ORDER[f.severity] < SEVERITY_ORDER[min_sev]:
            continue
        by_rule.setdefault(f.rid, []).append(f)

    A("%-24s %-7s %-13s %7s  %s" % ("RULE", "SEV", "DISPOSITION", "HITS", "FILES"))
    A("-" * 78)
    for r in RULES:
        hits = by_rule.get(r.rid, [])
        if not hits:
            continue
        nfiles = len({h.path for h in hits})
        A("%-24s %-7s %-13s %7d  %d" % (r.rid, r.severity, r.disposition,
                                        len(hits), nfiles))
    A("-" * 78)
    total = sum(len(v) for v in by_rule.values())
    A("TOTAL %d findings in %d files" % (total, len({f.path for f in findings})))
    if allowed:
        A("allowlisted (suppressed by an in-line `# scrub-allow:` marker): %s"
          % ", ".join("%s x%d" % kv for kv in sorted(allowed.items())))
    A("")

    for r in RULES:
        hits = by_rule.get(r.rid, [])
        if not hits:
            continue
        A("=" * 78)
        A("%s  [%s]  -> %s" % (r.rid, r.severity, r.disposition))
        A("-" * 78)
        A("WHY : " + r.why)
        A("FIX : " + r.fix)
        A("")
        per_file: dict[str, list[Finding]] = {}
        for h in hits:
            per_file.setdefault(h.path, []).append(h)
        for path in sorted(per_file):
            hs = per_file[path]
            A("  %s  (%d)" % (path, len(hs)))
            for h in hs[:show]:
                A("      %5d | %s" % (h.line, h.text))
            if len(hs) > show:
                A("      ... %d more" % (len(hs) - show))
        A("")
    return "\n".join(lines)


# --------------------------------------------------------------------------
# THE NEGATIVE CONTROL
# --------------------------------------------------------------------------
def selftest(keep: bool) -> int:
    """Plant one specimen per rule in a scratch copy and prove every rule fires.

    Two controls, both required:
      POSITIVE  a planted file per rule -> that rule MUST fire
      NEGATIVE  a clean file with no identifiers -> NO rule may fire

    Plus a THIRD, the one that actually matters for a repo audit: plant a
    specimen into a copy of a REAL tracked file and prove the normal scan path
    (git ls-files -> read -> match) surfaces it.  Regexes that work in isolation
    but never get handed the file are the failure mode this catches.
    """
    scratch = Path(tempfile.mkdtemp(prefix="scrub_selftest_"))
    print("selftest scratch: %s" % scratch)
    ok = True

    # --- NEGATIVE control: a file that must produce nothing -----------------
    clean = scratch / "clean_control.py"
    clean.write_text(
        "# a clean module: spec constants only, no identifiers\n"
        "LDPC_N = 64800\n"
        "FFT_SIZES = (8192, 16384, 32768)\n"
        "def demap(llr):\n"
        "    return llr * 2.0\n", encoding="utf-8")
    f, _ = scan_file(scratch, "clean_control.py", RULES)
    if f:
        ok = False
        print("  FAIL  negative control fired: %s"
              % ", ".join("%s@%d" % (x.rid, x.line) for x in f))
    else:
        print("  pass  negative control: clean file -> 0 findings")

    # --- POSITIVE control: one planted specimen per rule --------------------
    print()
    for r in RULES:
        if not r.specimen:
            print("  ....  %-24s no built-in specimen (local-literals-only "
                  "rule); covered by the local-file leg below" % r.rid)
            continue
        name = r.specimen_path or ("plant_%s.txt" % r.rid.lower())
        target = scratch / name
        target.parent.mkdir(parents=True, exist_ok=True)
        body = ("line one is innocent\n"
                + r.specimen + "\n"
                + "line three is innocent\n")
        target.write_bytes(body.encode("utf-8", errors="replace"))
        rel = str(target.relative_to(scratch)).replace("\\", "/")
        found, _ = scan_file(scratch, rel, RULES)
        fired = [x for x in found if x.rid == r.rid]
        if fired:
            print("  pass  %-24s fired on its specimen at %s:%d"
                  % (r.rid, rel, fired[0].line))
        else:
            ok = False
            print("  FAIL  %-24s did NOT fire on its own specimen (%r in %s)"
                  % (r.rid, r.specimen, rel))

    # --- the local-literals path --------------------------------------------
    # The shape rules cannot catch a project name or a hostname, so the exact
    # ones live in a gitignored file. That path needs its own control: if the
    # loader silently failed, the audit would look identical and be blind to
    # every literal the operator cares most about.
    print()
    import copy
    local_specimens = {
        "PRIVATE_REPO": "zzsecretproj",
        "HOST_USER": "zzworkstation",
        "SITE_ANTENNA": "The ZZ Antenna",
        "CALLSIGN": "ZQZQ",
        "LOCAL_PATH": "zz-sibling-repo",
    }
    local_file = scratch / "scrub_local.json"
    local_file.write_text(json.dumps({k: [v] for k, v in local_specimens.items()}),
                          encoding="utf-8")
    local_rules = copy.deepcopy(RULES)
    n_loaded, touched = apply_local(local_rules, local_file)
    if n_loaded != len(local_specimens):
        ok = False
        print("  FAIL  local file: loaded %d of %d literals"
              % (n_loaded, len(local_specimens)))
    else:
        print("  pass  local file: %d literals loaded into %s"
              % (n_loaded, ", ".join(touched)))
    for rid, literal in local_specimens.items():
        f = scratch / ("local_%s.txt" % rid.lower())
        f.write_text("harmless line\n%s\nharmless line\n" % literal,
                     encoding="utf-8")
        rel = str(f.relative_to(scratch)).replace("\\", "/")
        found, _ = scan_file(scratch, rel, local_rules)
        if any(x.rid == rid for x in found):
            print("  pass  %-24s fired on its LOCAL literal" % rid)
        else:
            ok = False
            print("  FAIL  %-24s did NOT fire on its local literal %r"
                  % (rid, literal))
    # and the negative: without the local file, those literals must be invisible
    quiet = [rid for rid in local_specimens
             if not any(x.rid == rid for x, in
                        [(y,) for y in scan_file(
                            scratch, "local_%s.txt" % rid.lower(), RULES)[0]])]
    if len(quiet) == len(local_specimens):
        print("  pass  without the local file, all %d literals are invisible "
              "(so the published tool carries no dossier)" % len(quiet))
    else:
        print("  note  %d of %d local literals are ALSO caught by a shape rule"
              % (len(local_specimens) - len(quiet), len(local_specimens)))

    # --- the real one: plant into a copy of a real tracked file -------------
    print()
    tracked = [p for p in git_tracked(REPO) if p.endswith(".py")]
    canary_src = None
    for cand in ("atsc3/bootstrap.py", "selftest.py"):
        if cand in tracked:
            canary_src = cand
            break
    if canary_src is None and tracked:
        canary_src = tracked[0]
    if canary_src:
        copy_root = scratch / "repo_copy"
        (copy_root / Path(canary_src).parent).mkdir(parents=True, exist_ok=True)
        shutil.copy2(REPO / canary_src, copy_root / canary_src)
        with open(copy_root / canary_src, "a", encoding="utf-8") as fh:
            fh.write("\n# planted by --selftest: "
                     "HDHR = \"192.168.0.99\"  # device 0123ABCD\n")
        # BEFORE: the untouched original must be quiet for these two rules
        base, _ = scan_file(REPO, canary_src, RULES)
        base_ids = {x.rid for x in base} & {"LAN_ADDR", "DEVICE_ID"}
        after, _ = scan_file(copy_root, canary_src, RULES)
        after_ids = {x.rid for x in after}
        if base_ids:
            print("  note  %s already carries %s before planting"
                  % (canary_src, sorted(base_ids)))
        if {"LAN_ADDR", "DEVICE_ID"} <= after_ids:
            print("  pass  canary: planting into a copy of %s surfaced LAN_ADDR "
                  "+ DEVICE_ID through the normal scan path" % canary_src)
        else:
            ok = False
            print("  FAIL  canary: planted identifiers in %s were NOT surfaced "
                  "(got %s)" % (canary_src, sorted(after_ids)))

    print()
    if not keep:
        shutil.rmtree(scratch, ignore_errors=True)
        print("scratch removed (use --keep-scratch to inspect it)")
    else:
        print("scratch kept: %s" % scratch)

    print()
    print("SELFTEST %s" % ("PASSED -- the instrument is proven to fire" if ok
                           else "FAILED -- a clean scan from this build means NOTHING"))
    return 0 if ok else 2


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", default=str(REPO))
    ap.add_argument("--all-files", action="store_true",
                    help="walk the working tree instead of git ls-files "
                         "(catches untracked files that a careless `git add -A` "
                         "would sweep in)")
    ap.add_argument("--include-self", action="store_true")
    ap.add_argument("--severity", default="INFO", choices=list(SEVERITY_ORDER))
    ap.add_argument("--fail-on", default="FATAL", choices=list(SEVERITY_ORDER))
    ap.add_argument("--rule", action="append", help="only these rule ids")
    ap.add_argument("--show", type=int, default=6, help="example lines per file")
    ap.add_argument("--json", help="write findings as json")
    ap.add_argument("--out", help="write the text report here as well as stdout")
    ap.add_argument("--local", default=str(LOCAL_DEFAULT),
                    help="per-installation literals (gitignored). Default: "
                         "tools/scrub_local.json. Without it the audit runs on "
                         "SHAPES ONLY, which catches the classes but not a "
                         "specific hostname or call sign.")
    ap.add_argument("--no-local", action="store_true",
                    help="ignore scrub_local.json -- shows what a reader with "
                         "only the published tool would find")
    ap.add_argument("--history", action="store_true",
                    help="also audit every blob in git history -- REQUIRED before "
                         "a first publish, because publishing a repo publishes "
                         "everything it ever contained")
    ap.add_argument("--selftest", action="store_true",
                    help="NEGATIVE CONTROL: prove every rule fires on a planted "
                         "specimen before you trust a clean scan")
    ap.add_argument("--keep-scratch", action="store_true")
    args = ap.parse_args()

    if args.selftest:
        return selftest(args.keep_scratch)

    root = Path(args.root).resolve()
    rules = RULES
    n_local, local_rules = (0, [])
    if not args.no_local:
        n_local, local_rules = apply_local(RULES, Path(args.local))
    print("local literals: %s" % (
        "%d loaded from %s, extending %s" % (n_local, args.local,
                                             ", ".join(local_rules))
        if n_local else
        "NONE -- running on shapes only. A specific hostname, antenna nickname "
        "or call sign will NOT be caught by name; see --local."))
    if args.rule:
        want = set(args.rule)
        rules = [r for r in RULES if r.rid in want]
        if not rules:
            print("no such rule; known: %s" % ", ".join(RULES_BY_ID))
            return 2

    paths = walk_all(root) if args.all_files else git_tracked(root)
    if not paths:
        print("no files to scan (not a git repo? try --all-files)")
        return 2

    findings, allowed = scan(root, paths, rules, args.include_self)

    if args.history:
        hist = scan_history(root, rules, args.severity)
        print("=" * 78)
        print("HISTORY AUDIT -- every blob ever committed (%d distinct leaks)"
              % len(hist))
        print("Publishing a repo publishes its history. A finding here is NOT")
        print("fixed by an edit -- it needs git-filter-repo and a force-push to")
        print("a repo nobody has cloned yet.")
        print("=" * 78)
        if not hist:
            print("  clean at severity >= %s" % args.severity)
        for h in sorted(hist, key=lambda x: (-SEVERITY_ORDER[x.severity], x.rid)):
            print("  %-7s %-22s %s" % (h.severity, h.rid, h.path))
            print("          %s" % h.text[:150])
        print()
        findings = findings + hist
    text = report(findings, allowed, root, paths, args.include_self,
                  args.severity, args.show)
    print(text)
    if args.out:
        Path(args.out).write_text(text, encoding="utf-8")
    if args.json:
        Path(args.json).write_text(json.dumps({
            "root": str(root),
            "files_scanned": len(paths),
            "include_self": args.include_self,
            "allowlisted": allowed,
            "rules": [{"id": r.rid, "severity": r.severity,
                       "disposition": r.disposition, "why": r.why, "fix": r.fix}
                      for r in rules],
            "findings": [f.__dict__ for f in findings],
        }, indent=2), encoding="utf-8")

    worst = max([SEVERITY_ORDER[f.severity] for f in findings], default=-1)
    return 1 if worst >= SEVERITY_ORDER[args.fail_on] else 0


if __name__ == "__main__":
    sys.exit(main())

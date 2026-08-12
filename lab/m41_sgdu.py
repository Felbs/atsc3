"""M41 -- the Service Guide: parse SGDUs, build the TV listings.

WHAT WE WERE MISSING
--------------------
The guide arrives on tsi2 as Service Guide Delivery Units: a BINARY container
holding many XML fragments end to end, often gzipped. `m7_objects` treated
each object as one XML blob, so it recovered only the fragments that happened
to be plain text and mangled the container framing on the way out. The SGDD
advertises 377 fragments; we were seeing about twenty, by accident.

THE CONTAINER, REVERSED FROM OUR OWN CAPTURES
---------------------------------------------
    offset  size  meaning
    0       8     reserved / unit id (all zero in every unit we hold)
    8       1     fragment count N
    9       12*N  per-fragment descriptor:
                    +0  4  fragment transport id
                    +4  4  version (zero in every unit we hold)
                    +8  4  OFFSET of this fragment from the start of the data
    9+12N   2     reserved
    9+12N+2 ...   the fragments, concatenated

Verified on two units before being trusted: a 3-fragment unit gives data
start 9+36+2 = 47 and offsets 0/1100/1838, and the fragments are found at
exactly 47/1147/1885 (+39 each for the XML declaration). A 9-fragment unit
gives 9+108+2 = 119 with offsets 0/595/1189, again exact. If the arithmetic
did not land on a '<' we would rather fail than guess -- see GATE 1.

FRAGMENT TYPES (OMA BCAST Service Guide, ATSC A/332)
    Service          a channel
    Schedule         binds a Content to a PresentationWindow (the times)
    Content          a programme: name, description, ratings, components
    PreviewData      artwork / preview
"""
import glob, gzip, os, re, struct, sys, collections
import xml.etree.ElementTree as ET

LAB = os.path.dirname(os.path.abspath(__file__))
OMA = "{urn:oma:xml:bcast:sg:fragments:1.0}"
SA = "{tag:atsc.org,2016:XMLSchemas/ATSC3/SA/1.0/}"

# PresentationWindow times are NTP seconds (epoch 1900), not Unix.
NTP_EPOCH = 2208988800


def unwrap(raw):
    return gzip.decompress(raw) if raw[:2] == b"\x1f\x8b" else raw


def parse_sgdu(d):
    """-> list of (fragment_transport_id, xml_bytes), or [] if not an SGDU."""
    if len(d) < 11:
        return []
    n = d[8]
    if not 1 <= n <= 255:
        return []
    need = 9 + 12 * n + 2
    if len(d) < need:
        return []
    desc = []
    for i in range(n):
        o = 9 + 12 * i
        fid, ver, off = struct.unpack(">III", d[o:o + 12])
        desc.append((fid, ver, off))
    base = need
    out = []
    for i, (fid, ver, off) in enumerate(desc):
        start = base + off
        end = base + desc[i + 1][2] if i + 1 < n else len(d)
        if not (0 <= start < end <= len(d)):
            return []
        frag = d[start:end]
        if not frag.lstrip()[:8].startswith(b"<"):
            return []                       # arithmetic did not land: refuse
        # A fragment does not fill its slot: the next offset is padded, so
        # the slice carries trailing bytes and an XML parser rejects the lot
        # as junk after the document element. Cut at the root's own closing
        # tag. (Symptom before this: exactly one fragment per unit parsed --
        # the LAST, the only one with nothing after it.)
        m = re.match(rb"\s*(?:<\?xml[^>]*\?>\s*)?<([A-Za-z_][\w.:-]*)", frag)
        if m:
            close = b"</" + m.group(1) + b">"
            k = frag.rfind(close)
            if k >= 0:
                frag = frag[:k + len(close)]
        out.append((fid, frag))
    return out


def text_attr(el, tag, attr="text"):
    x = el.find(OMA + tag)
    return x.get(attr) if x is not None else None


def main():
    src = sys.argv[1] if len(sys.argv) > 1 else "stress_0807/live_out2"
    d = os.path.join(LAB, src)
    files = sorted(glob.glob(os.path.join(d, "*")))
    units = frags = 0
    kinds = collections.Counter()
    services, contents, schedules = {}, {}, []
    not_sgdu = 0

    for f in files:
        if f.endswith((".json", ".srt")):
            continue
        raw = open(f, "rb").read()
        parts = parse_sgdu(unwrap(raw))
        if not parts:
            not_sgdu += 1
            continue
        units += 1
        for fid, x in parts:
            frags += 1
            try:
                root = ET.fromstring(x.decode("utf-8", "replace"))
            except ET.ParseError:
                kinds["<unparseable>"] += 1
                continue
            tag = root.tag.replace(OMA, "")
            kinds[tag] += 1
            if tag == "Service":
                services[root.get("id")] = (
                    text_attr(root, "Name") or "?")
            elif tag == "Content":
                # The rating is the element's TEXT -- <sa:RatingDescription>
                # TV-MA</sa:RatingDescription> -- not a `text` attribute and
                # not a child <Name>. Reading it the OMA way returned nothing
                # and reported 0/361 rated, which was wrong, not empty.
                rat = [r.text.strip() for r in root.iter(SA + "RatingDescription")
                       if r.text and r.text.strip()]
                region = [r.text.strip() for r in root.iter(SA + "RegionIdentifier")
                          if r.text and r.text.strip()]
                if region:
                    rat = [f"{v} (region {region[0]})" for v in rat]
                contents[root.get("id")] = dict(
                    name=text_attr(root, "Name") or "?",
                    desc=text_attr(root, "Description") or "",
                    rating="; ".join(x for x in rat if x))
            elif tag == "Schedule":
                svc = root.find(OMA + "ServiceReference")
                for cr in root.iter(OMA + "ContentReference"):
                    pw = cr.find(OMA + "PresentationWindow")
                    if pw is None:
                        continue
                    schedules.append(dict(
                        svc=svc.get("idRef") if svc is not None else "?",
                        cid=cr.get("idRef"),
                        start=int(pw.get("startTime", 0)),
                        end=int(pw.get("endTime", 0))))

    print(f"scanned {len(files)} objects")
    print(f"  SGDUs parsed      {units}")
    print(f"  fragments         {frags}")
    print(f"  not an SGDU       {not_sgdu}")
    print(f"  by type           {dict(kinds)}")

    # ---- gates ------------------------------------------------------
    print()
    g1 = frags > 0 and kinds.get("<unparseable>", 0) == 0
    print(f"GATE 1  every fragment the header located parses as XML"
          f"      {'PASS' if g1 else 'FAIL'}")
    g2 = len(schedules) > 0 and len(contents) > 0
    print(f"GATE 2  schedules {len(schedules)} and contents {len(contents)}"
          f" both non-empty   {'PASS' if g2 else 'FAIL'}")
    linked = sum(1 for s in schedules if s["cid"] in contents)
    g3 = linked > 0
    print(f"GATE 3  schedule->content references resolve: "
          f"{linked}/{len(schedules)}   {'PASS' if g3 else 'FAIL'}")

    # ---- the listings -----------------------------------------------
    import datetime as dt
    print(f"\n=== TV LISTINGS off our own antenna ===")
    rows = []
    for s in schedules:
        c = contents.get(s["cid"])
        if not c:
            continue
        t0 = dt.datetime.fromtimestamp(s["start"] - NTP_EPOCH, dt.timezone.utc)
        t1 = dt.datetime.fromtimestamp(s["end"] - NTP_EPOCH, dt.timezone.utc)
        rows.append((t0, t1, services.get(s["svc"], s["svc"].split("/")[-1]),
                     c["name"], c["rating"], c["desc"]))
    rows.sort()
    for t0, t1, svc, name, rating, desc in rows[:40]:
        mins = int((t1 - t0).total_seconds() / 60)
        print(f"  {t0:%m-%d %H:%M}Z {mins:>4}m  {svc[:22]:<22} {name[:38]:<38}"
              f" {rating[:18]}")
        if desc and desc != name:
            print(f"                             {desc[:88]}")
    print(f"\n  {len(rows)} scheduled programmes resolved")
    if rows:
        span = (rows[-1][1] - rows[0][0]).total_seconds() / 3600
        print(f"  guide spans {rows[0][0]:%m-%d %H:%M}Z .. "
              f"{rows[-1][1]:%m-%d %H:%M}Z  = {span:.1f} hours")
        chans = sorted({r[2] for r in rows})
        print(f"  {len(chans)} services: {', '.join(chans)}")
        rated = [r for r in rows if r[4]]
        print(f"  {len(rated)}/{len(rows)} carry a content advisory rating")
        for r in rated[:4]:
            print(f"      {r[3][:44]:<44} {r[4]}")
        desc = sum(1 for r in rows if r[5] and r[5] != r[3])
        print(f"  {desc}/{len(rows)} carry an episode description")
    return 0


if __name__ == "__main__":
    sys.exit(main())

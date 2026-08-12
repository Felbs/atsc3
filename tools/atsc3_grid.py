#!/usr/bin/env python3
"""atsc3_grid -- the unified EPG grid for the 4K tuner UI  (E70, task #65)

Every service on every carrier in one view, now/next per cell, straight from
the broadcaster's own service guide.

THE PRODUCT DECISION THIS IMPLEMENTS (operator, 8/10)
    "Add what's playing on the encrypted channels to the grid ... just add a
     disclaimer that those channels are encrypted. We should allow our users
     to see whatever decodable data is still accessible."

So a Widevine-protected service gets a full row: its name, its virtual
channel, its resolution and track list, its listings, and its live captions --
because ALL of that is transmitted in the clear (E67 measured it: 98.3% of the
VIDEO payload is encrypted, and none of the signalling, none of the manifest,
and none of the stpp subtitle track).  What it does not get is any suggestion
that you can watch it.  The badge says

    ENCRYPTED -- guide and captions only

and the row carries no play affordance at all.  No circumvention is attempted,
hinted at, or linked.

AND IT NEVER RENDERS A PICTURE.  E67's law: "the decoder ran" is not "the
decoder decoded" -- CENC output passes ffprobe, emits 600 frames with zero
errors and tiles its NAL prefixes 120/120, and is still a flat green field.
A thumbnail would have to pass lab/e67_picture_gate.py's flat_frac test first,
so this page has no thumbnails.  See `PICTURE_POLICY` below and gate leg 5.

Conventions borrowed from the fleet's panels (tv_tuna_panel.py :8642,
prop_history.py :8649): stdlib only, all CSS/JS inline, no external requests,
localhost bind, a BUILD stamp the page checks so a stale tab announces itself,
and the grid DOM is rebuilt only when its content key changes -- never from a
timer, because rebuilding innerHTML under a click destroys the button.

    python tools/atsc3_grid.py                  # serve on 127.0.0.1:8650
    python tools/atsc3_grid.py --port 8651
    python tools/atsc3_grid.py --out grid.html  # one self-contained file
    python tools/atsc3_grid.py --at guide       # anchor on the banked guide
"""
from __future__ import annotations

import argparse
import datetime as dt
import html
import json
import os
import sys
import threading
import time
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if HERE not in sys.path:
    sys.path.insert(0, HERE)

import atsc3_inspect as I   # noqa: E402
import atsc3_vs as V        # noqa: E402

# Bump on every UI change.  An open tab keeps running its old JS across a
# restart of this server; the page compares this against its own copy and
# tells the user to hard-refresh rather than silently showing stale code.
BUILD = "0811-e81grid5"

PORT = 8650                 # 8642-8649 are taken; 8650 was the next free one
SLOT_SECS = 1800            # 30-minute cells, same as stvt_epg.render_grid
DEFAULT_SLOTS = 6
REBUILD_SECS = 300          # the lineup is banked data; it does not churn

PICTURE_POLICY = {
    "renders_thumbnails": False,
    "gate_required_first": "lab/e67_picture_gate.py",
    "why": ("E67: ffprobe, frame count, error count and NAL tiling all PASS "
            "on CENC-encrypted video. flat_frac is the only measurement that "
            "discriminated (0.016 real vs 0.97-0.99 garbage), and temporal "
            "correlation does NOT. So this page claims no picture at all."),
}

BADGE = {
    "PROTECTED": ("enc", "ENCRYPTED — guide and captions only"),
    "CLEAR": ("ok", "clear"),
    "UNKNOWN": ("unk", "protection UNKNOWN"),
}

_lock = threading.Lock()
_cache = {"model": None, "built": 0.0, "sources": None}
_vs_cache = {"sb": None, "built": 0.0}
VS_REBUILD_SECS = 900     # the scoreboard reads banked probes + rolling logs


def get_scoreboard(force=False):
    """The E81 head-to-head, cached.  Never probes the radio or the FLEX --
    `atsc3_vs --probe` does that separately and banks its readings."""
    with _lock:
        now = time.time()
        if (_vs_cache["sb"] is None or force
                or now - _vs_cache["built"] > VS_REBUILD_SECS):
            _vs_cache["sb"] = V.scoreboard()
            _vs_cache["built"] = now
        return _vs_cache["sb"]


# ------------------------------------------------------------------- model

def get_model(sources=None, force=False):
    with _lock:
        now = time.time()
        if (_cache["model"] is None or force
                or now - _cache["built"] > REBUILD_SECS
                or _cache["sources"] != sources):
            _cache["model"] = I.build_lineup(sources, at=now)
            _cache["built"] = now
            _cache["sources"] = sources
        return _cache["model"]


def _anchor(model, at):
    """Where does the grid start?  Honest about a guide that does not reach now.

    Returns (origin_unix, status, message).  `status` is 'live' when the guide
    covers the requested time and 'stale' when it does not -- and when it does
    not we do NOT quietly slide to a time that works.  We say so, and offer the
    guide's own window as an explicit choice.
    """
    span = (model.get("esg") or {}).get("span")
    now = time.time()
    if at == "guide" and span:
        origin = _busiest_slot(model["esg"])
        return origin, "anchored", (
            "Anchored on the banked guide's busiest half-hour (%s) -- "
            "the guide is stitched from captures on different days, and its "
            "earliest slot covers only the services that were on air then. "
            "This is what the broadcaster was transmitting at that moment, "
            "not what is on now." % _ts(origin))
    try:
        t = float(at)
    except (TypeError, ValueError):
        t = now
    origin = (int(t) // SLOT_SECS) * SLOT_SECS
    if not span:
        return origin, "no-guide", "No service-guide fragments are banked."
    if span[0] <= t <= span[1]:
        return origin, "live", ""
    return origin, "stale", (
        "The banked service guide runs %s .. %s and does not reach the "
        "current time, so no now/next is claimed for any service — "
        "clear or encrypted. The rest of the row (protection, transport, "
        "tracks, captions) is read from the air and is current to its "
        "capture." % (_ts(span[0]), _ts(span[1])))


def _busiest_slot(esg):
    """The half-hour where the most distinct services have a programme.

    The banked guide is stitched from captures taken on different days, so its
    first slot is not its representative one -- WRC, for instance, only
    appears in the 8/07 capture.  Picking the busiest slot shows the guide at
    its fullest instead of at its earliest.
    """
    best, best_n = None, -1
    for r in esg["rows"]:
        origin = (r["start"] // SLOT_SECS) * SLOT_SECS
        n = len({x["gsid"] or x["svc_name"] for x in esg["rows"]
                 if x["start"] <= origin < x["end"]})
        if n > best_n:
            best, best_n = origin, n
    return best or 0


def _ts(u):
    return dt.datetime.fromtimestamp(u, dt.timezone.utc).strftime(
        "%m-%d %H:%MZ")


def _slot_labels(origin, n):
    return [dt.datetime.fromtimestamp(origin + i * SLOT_SECS,
                                      dt.timezone.utc).strftime("%H:%M")
            for i in range(n)]


def _clean(s, limit=60):
    """Broadcaster text into a cell.  Strip anything that would mojibake a
    console or smuggle markup; the page escapes on top of this."""
    if not s:
        return ""
    s = "".join(c if 32 <= ord(c) < 0x2500 else "?" for c in str(s))
    return s[:limit]


def grid_json(sources=None, at=None, slots=DEFAULT_SLOTS):
    model = get_model(sources)
    origin, status, message = _anchor(model, at)
    esg = model["esg"]
    rows = []
    for bsid, c in model["carriers"].items():
        rf = c.get("rf") or {}
        for svc in c["services"]:
            prot = svc["protection"]
            cls, label = BADGE[prot["state"]]
            play = svc["playable"]
            cells = []
            for i in range(slots):
                t0 = origin + i * SLOT_SECS
                prog = None
                if status != "stale":
                    for r in esg["rows"]:
                        if r["gsid"] and r["gsid"] != svc.get("globalServiceID"):
                            continue
                        if not r["gsid"] and r["svc_name"] != svc.get("name"):
                            continue
                        if r["start"] <= t0 < r["end"]:
                            prog = r
                            break
                if prog is None:
                    cells.append({"title": "", "empty": True})
                else:
                    cells.append({
                        "title": _clean(prog["name"], 44),
                        "rating": _clean(prog["rating"], 18),
                        "desc": _clean(prog["desc"], 220),
                        "cont": prog["start"] < t0,
                    })
            m = svc.get("media") or {}
            v = m.get("video") or {}
            rows.append({
                "bsid": bsid, "sid": svc["serviceId"],
                "rf": rf.get("rf"), "center_mhz": c.get("center_mhz"),
                "name": _clean(svc.get("name"), 12) or "?",
                "vchannel": svc.get("vchannel") or "-",
                "state": prot["state"], "badge_class": cls,
                "badge": label, "drm": prot["drm"],
                "layers_read": prot["layers_read"],
                "transport": svc.get("transport") or "?",
                "plp": (svc.get("plp") or {}).get("plp"),
                "category": svc.get("category"),
                "playable": play["playable"], "reason": play["reason"],
                "detail": play["detail"],
                "realtime": play.get("realtime"),
                "realtime_note": play.get("realtime_note"),
                "video": (("%s %s" % (v.get("resolution") or "",
                                      (v.get("codec") or "").split(".")[0])
                           ).strip() or None),
                "bitrate": v.get("bitrate_bps"),
                "audio": [a.get("lang") or (a.get("codec") or "?").split(".")[0]
                           for a in (m.get("audio") or [])],
                "captions_n": len(svc.get("captions") or []),
                "captions_status": svc.get("captions_status"),
                "live": bool(svc.get("live")),
                "cells": cells,
            })
    rows.sort(key=lambda r: (-(r["rf"] or 0), _chan_key(r["vchannel"])))
    return {
        "build": BUILD, "generated": time.time(),
        "origin": origin, "slots": _slot_labels(origin, slots),
        "slot_secs": SLOT_SECS, "status": status, "message": message,
        "guide": {"programmes": len(esg["rows"]), "services": esg["services"],
                  "span": esg["span"], "sources": esg["sources"]},
        "sources": model["sources"],
        "picture_policy": PICTURE_POLICY,
        "rows": rows,
    }


def _chan_key(v):
    try:
        a, b = str(v).split(".")
        return (int(a), int(b))
    except Exception:
        return (9999, 0)


def service_json(bsid, sid, sources=None):
    model = get_model(sources)
    c = model["carriers"].get(str(bsid))
    if not c:
        return None
    for svc in c["services"]:
        if str(svc["serviceId"]) != str(sid):
            continue
        return {
            "build": BUILD,
            "name": svc.get("name"), "vchannel": svc.get("vchannel"),
            "bsid": bsid, "sid": sid,
            "rf": (c.get("rf") or {}).get("rf"),
            "rf_source": (c.get("rf") or {}).get("source"),
            "category": svc.get("category"),
            "transport": svc.get("transport"),
            "transport_source": svc.get("transport_source"),
            "plp": svc.get("plp"),
            "protection": svc["protection"],
            "playable": svc["playable"],
            "picture": svc["picture"],
            "media": svc.get("media"),
            "captions": svc.get("captions") or [],
            "captions_kind": svc.get("captions_kind"),
            "captions_source": svc.get("captions_source"),
            "captions_status": svc.get("captions_status"),
            "live": svc.get("live"),
            "globalServiceID": svc.get("globalServiceID"),
            "flow": svc.get("flow"),
        }
    return None


# -------------------------------------------------------------------- page

PAGE = """<!DOCTYPE html>
<html><head><meta charset="utf-8"><title>ATSC 3.0 guide</title>
<style>
body{font-family:Segoe UI,system-ui,sans-serif;background:#05080f;color:#dce6f2;
     margin:0;padding:12px}
h1{font-size:20px;margin:0 0 2px}
.sub{color:#7f96b3;font-size:12px;margin-bottom:10px}
#banner{background:#0d1626;border:1px solid #26436b;border-radius:10px;
        padding:9px 14px;margin-bottom:12px;font-size:13px;line-height:1.5}
#banner.warn{border-color:#6b5a26;background:#161206;color:#e7c96a}
#banner a{color:#7fb0ee}
#wrap{overflow-x:auto;max-width:100%}
table{border-collapse:collapse;width:100%;min-width:900px;font-size:12px}
th{position:sticky;top:0;background:#05080f;text-align:left;padding:6px;
   color:#7f96b3;border-bottom:1px solid #26436b;z-index:2}
td{padding:5px 6px;border-bottom:1px solid #101c30;vertical-align:top}
.chan{white-space:nowrap}
.vch{font-weight:700;color:#e8f1fc;font-size:13px}
.call{color:#9fb4d0}
.meta{color:#5f7591;font-size:10.5px;margin-top:2px}
.badge{display:inline-block;border-radius:4px;padding:1px 6px;font-size:10px;
       font-weight:700;letter-spacing:.3px;white-space:nowrap}
.b-ok{background:#124a2a;color:#67d18a}
.b-enc{background:#4a1a1a;color:#e77}
.b-unk{background:#4a3b12;color:#e7c96a}
.lockrow td{background:#0b0f18}
.show{color:#c7d5e8}
.cont{color:#5f7591}
.rating{color:#7f96b3;font-size:10px}
.empty{color:#2c3b52}
.now{outline:1px solid #2f79d4;border-radius:4px}
.cap{color:#8fb7e6;font-size:10.5px;margin-top:3px;font-style:italic}
tr.svc{cursor:pointer}
tr.svc:hover td{background:#0d1626}
#detail{position:fixed;left:50%;transform:translateX(-50%);top:4vh;
        width:min(560px,94vw);max-height:88vh;overflow:auto;
        background:#0d1626;border:1px solid #26436b;border-radius:10px;
        padding:12px 14px;font-size:12px;display:none;
        box-shadow:0 8px 40px #000d;z-index:9}
#scrim{position:fixed;inset:0;background:#000a;display:none;z-index:8}
#detail h2{font-size:15px;margin:0 0 6px}
#detail .k{color:#7f96b3;font-size:10px;text-transform:uppercase;
           letter-spacing:.5px;margin-top:9px}
#detail table{font-size:10.5px}
#detail td{padding:2px 4px;border-bottom:1px solid #101c30}
#detail .src{color:#5f7591;font-size:9.5px;word-break:break-all}
#close{float:right;background:#152238;color:#9fb4d0;border:1px solid #26436b;
       border-radius:6px;padding:2px 9px;cursor:pointer}
.edu{color:#5f7591;font-size:11px;margin-top:14px;line-height:1.6;
     max-width:1000px}
.edu b{color:#9fb4d0}
.good{color:#67d18a}.warn{color:#e7c96a}.bad{color:#e77}
#tabs{margin:0 0 10px}
.tab{background:#0d1626;color:#9fb4d0;border:1px solid #26436b;
     border-radius:7px;padding:5px 14px;margin-right:6px;cursor:pointer;
     font-size:12px;font-family:inherit}
.tab.on{background:#17335c;color:#e8f1fc;border-color:#3d6ba8}
.vsrow{border-bottom:1px solid #101c30;padding:10px 4px}
.vsq{color:#7f96b3;font-size:11px}
.vsn{font-size:13px;margin:4px 0}
.vsn b{color:#e8f1fc}
.vsc{color:#5f7591;font-size:10.5px;margin-top:3px;line-height:1.5}
.v-WIN{background:#124a2a;color:#67d18a}
.v-LOSE{background:#4a1a1a;color:#e77}
.v-TIE{background:#1e3350;color:#9fb4d0}
.v-NC{background:#4a3b12;color:#e7c96a}
.vbadge{display:inline-block;border-radius:4px;padding:1px 8px;font-size:10px;
        font-weight:700;letter-spacing:.4px;margin-right:8px}
</style></head><body>
<h1>ATSC 3.0 &mdash; unified guide</h1>
<div class="sub" id="sub">loading&hellip;</div>
<div id="tabs">
  <button class="tab on" id="t-guide" onclick="setTab('guide')">guide</button>
  <button class="tab" id="t-vs" onclick="setTab('vs')">vs&nbsp;HDHomeRun</button>
</div>
<div id="banner"></div>
<div id="wrap"></div>
<div id="scrim" onclick="hideDetail()"></div>
<div id="detail"><button id="close" onclick="hideDetail()">close</button>
  <div id="dbody"></div></div>
<div class="edu">
<b>Why encrypted channels are still listed.</b> An ATSC 3.0 broadcaster sends
far more in the clear than the picture. The service list, the virtual channel,
the DASH manifest (resolution, codecs, bitrates, languages), the segment
timing and the <b>subtitle track</b> are all unencrypted even when the video
and audio are locked with Widevine &mdash; measured, not assumed: 98.3&nbsp;% of
the video payload bytes on 5.1/4.1/9.1 are protected, and 0&nbsp;% of the
stpp captions are. So those rows carry everything the station transmits
openly, and a badge that says exactly what you cannot do with them.
<b>Nothing here attempts or assists decryption.</b>
<br><br>
<b>Why there are no thumbnails.</b> Encrypted video decodes <i>cleanly</i>:
it parses, it emits frames, it reports zero errors, and it tiles its NAL
prefixes perfectly &mdash; CENC is designed to leave the structure intact and
encrypt only the payload. The only measurement that told the two apart was the
fraction of pixels sharing one luma value. A picture would have to pass that
gate before it could appear here, so no picture appears here.
</div>
<script>
var AT = new URLSearchParams(location.search).get('at') || '';
var KEY = null, ROWS = [];
function esc(s){return (s===null||s===undefined?'':String(s))
  .replace(/&/g,'&amp;').replace(/</g,'&lt;').replace(/>/g,'&gt;')
  .replace(/"/g,'&quot;');}
function cell(c){
  if(c.empty) return '<td class="empty">&middot;</td>';
  var t = (c.cont?'<span class="cont">&raquo; ':'<span class="show">')
        + esc(c.title) + '</span>';
  if(c.rating) t += '<div class="rating">'+esc(c.rating)+'</div>';
  return '<td>'+t+'</td>';
}
function render(g){
  document.getElementById('sub').textContent =
    g.rows.length+' services on '+
    (new Set(g.rows.map(function(r){return r.bsid;}))).size+
    ' carrier(s)  \\u2022  guide: '+g.guide.programmes+' programmes'+
    '  \\u2022  build '+g.build;
  var b = document.getElementById('banner');
  b.className = (g.status==='live') ? '' : 'warn';
  var msg = g.message || ('Guide window '+
      (g.guide.span? fmt(g.guide.span[0])+' .. '+fmt(g.guide.span[1]) : '?'));
  if(g.status==='stale')
    msg += ' <a href="?at=guide">Show the guide\\u2019s own window &rarr;</a>';
  else if(g.status==='anchored')
    msg += ' <a href="?">Back to now &rarr;</a>';
  b.innerHTML = msg;
  // STABLE DOM: this table holds click targets. Rebuild only when the
  // content actually changes, never on a timer.
  var key = JSON.stringify([g.origin,g.status,g.rows.map(function(r){
      return [r.vchannel,r.state,r.cells.map(function(c){return c.title;})];})]);
  if(key===KEY) return;
  KEY = key; ROWS = g.rows;
  var h = '<table><thead><tr><th style="width:170px">channel</th>';
  for(var i=0;i<g.slots.length;i++) h += '<th>'+esc(g.slots[i])+'</th>';
  h += '</tr></thead><tbody>';
  for(var i=0;i<g.rows.length;i++){
    var r = g.rows[i];
    h += '<tr class="svc'+(r.state==='PROTECTED'?' lockrow':'')+
         '" onclick="showDetail('+i+')">';
    h += '<td class="chan"><span class="vch">'+esc(r.vchannel)+'</span> '+
         '<span class="call">'+esc(r.name)+'</span>'+
         (r.live?' <span class="badge b-ok">LIVE</span>':'')+
         '<div style="margin-top:3px"><span class="badge b-'+r.badge_class+
         '">'+esc(r.badge)+'</span></div>'+
         '<div class="meta">'+esc(r.transport)+
         (r.plp?' &middot; PLP'+esc(r.plp):'')+
         (r.rf?' &middot; RF'+esc(r.rf):'')+
         (r.video?' &middot; '+esc(r.video):'')+
         (r.audio.length?' &middot; '+esc(r.audio.join('/')):'')+'</div>';
    if(r.captions_n)
      h += '<div class="cap">captions available ('+r.captions_n+' cues)</div>';
    h += '</td>';
    for(var j=0;j<r.cells.length;j++) h += cell(r.cells[j]);
    h += '</tr>';
  }
  h += '</tbody></table>';
  document.getElementById('wrap').innerHTML = h;
}
function fmt(u){var d=new Date(u*1000);
  return d.toISOString().slice(5,16).replace('T',' ')+'Z';}
// the key is OUR label, not broadcaster text: it may hold entities.
function kv(k,v){return '<div class="k">'+k+'</div><div>'+v+'</div>';}
function showDetail(i){
  var r = ROWS[i];
  fetch('/api/service?bsid='+encodeURIComponent(r.bsid)+
        '&sid='+encodeURIComponent(r.sid))
    .then(function(x){return x.json();}).then(function(d){
    var h = '<h2>'+esc(d.vchannel)+' '+esc(d.name)+'</h2>';
    h += '<span class="badge b-'+r.badge_class+'">'+esc(r.badge)+'</span>';
    h += kv('playable', (d.playable.playable
        ? '<span class="good">YES</span>'
        : '<span class="bad">NO &mdash; '+esc(d.playable.reason)+'</span>')+
        '<br>'+esc(d.playable.detail)+
        (d.playable.realtime===false
          ? '<br><span class="warn">'+esc(d.playable.realtime_note)+'</span>'
          : ''));
    h += kv('protection &mdash; every field, and the file it was read from',
        '<table>'+d.protection.evidence.map(function(e){
          return '<tr><td>'+esc(e.layer)+'</td><td>'+esc(e.field)+'</td><td>'+
                 esc(e.value)+'</td></tr><tr><td></td><td colspan="2" '+
                 'class="src">'+esc(e.source)+
                 (e.note?' &mdash; '+esc(e.note):'')+'</td></tr>';
        }).join('')+'</table>');
    h += kv('transport', esc(d.transport)+
        '<div class="src">'+esc(d.transport_source)+'</div>'+
        'PLP '+esc(d.plp.plp||'unknown')+
        '<div class="src">'+esc(d.plp.source||d.plp.basis)+'</div>');
    var m = d.media||{};
  function bps(b){return !b?'':(b>=1e6?(b/1e6).toFixed(2)+' Mbps'
                                     :Math.round(b/1e3)+' kbps');}
    var tr = '';
    if(m.video) tr += '<tr><td>video</td><td>'+esc(m.video.resolution)+' '+
        esc(m.video.codec)+'</td><td>'+bps(m.video.bitrate_bps)+'</td></tr>';
    (m.audio||[]).forEach(function(a){tr += '<tr><td>audio</td><td>'+
        esc(a.lang||'-')+' '+esc(a.codec)+'</td><td>'+
        bps(a.bitrate_bps)+'</td></tr>';});
    (m.text||[]).forEach(function(a){tr += '<tr><td>captions</td><td>'+
        esc(a.lang||'-')+' '+esc(a.codec)+'</td><td></td></tr>';});
    h += kv('tracks (from the manifest, which is in the clear)',
        '<table>'+tr+'</table><div class="src">'+esc(m.source||'')+'</div>');
    if(d.captions.length){
      h += kv('captions &mdash; '+esc(d.captions_kind||''),
        d.captions.map(function(c){
          return '<div>'+esc(c.text)+'</div>';}).join('')+
        '<div class="src">'+esc(d.captions_source||'')+'</div>');
    } else if(d.captions_status==='empty'){
      h += kv('captions', 'stpp track present and unencrypted, but the '+
        'banked segments carry an empty document &mdash; no cues in this '+
        'window.');
    }
    h += kv('picture', '<span class="warn">no thumbnail, by policy</span>'+
        '<div class="src">'+esc(d.picture.why)+' &mdash; gate: '+
        esc(d.picture.gate)+'</div>');
    document.getElementById('dbody').innerHTML = h;
    document.getElementById('detail').style.display='block';
    document.getElementById('scrim').style.display='block';
  });
}
function hideDetail(){document.getElementById('detail').style.display='none';
  document.getElementById('scrim').style.display='none';}
// The tab is part of the URL, not just page state: it is deep-linkable,
// survives a reload, and can be verified without simulating a click.
var TAB = (new URLSearchParams(location.search).get('tab')==='vs')
          ? 'vs' : 'guide';
var VSKEY = null;
function setTab(t){
  TAB = t;
  try{
    var u = new URL(location.href);
    if(t==='vs') u.searchParams.set('tab','vs');
    else u.searchParams.delete('tab');
    history.replaceState(null,'',u);
  }catch(e){}
  document.getElementById('t-guide').className = 'tab'+(t==='guide'?' on':'');
  document.getElementById('t-vs').className = 'tab'+(t==='vs'?' on':'');
  KEY = null; VSKEY = null;              // force one rebuild on a real click
  load();
}
function vsbadge(v){
  var cls = v==='NOT-COMPARABLE' ? 'NC' : v;
  return '<span class="vbadge v-'+cls+'">'+esc(v)+'</span>';
}
function renderVS(sb){
  document.getElementById('sub').textContent =
    'head-to-head \\u2022 win '+sb.tally.WIN+'  lose '+sb.tally.LOSE+
    '  tie '+sb.tally.TIE+'  n/c '+sb.tally['NOT-COMPARABLE']+
    '  \\u2022 '+sb.headline+'  \\u2022 build '+sb.build;
  var b = document.getElementById('banner');
  b.className = sb.tally.LOSE > sb.tally.WIN ? 'warn' : '';
  b.innerHTML = 'One antenna, passive splitter \\u2014 every paired reading '+
    'is same-air and same-instant. A favourable number we cannot defend is '+
    'worse than an honest loss, so a dimension the method cannot support is '+
    'marked NOT-COMPARABLE rather than scored.';
  var key = JSON.stringify(sb.tally)+sb.build+(sb.probe_when||'');
  if(key===VSKEY) return;
  VSKEY = key;
  var h = '';
  for(var i=0;i<sb.dimensions.length;i++){
    var d = sb.dimensions[i];
    h += '<div class="vsrow">'+vsbadge(d.verdict)+
         '<b>'+esc(d.dimension)+'</b>'+
         '<div class="vsq">'+esc(d.question)+'</div>'+
         '<div class="vsn">ours <b>'+esc(d.ours)+'</b>'+
         ' &nbsp;&middot;&nbsp; theirs <b>'+esc(d.theirs)+'</b></div>'+
         '<div class="vsc">method: '+esc(d.method)+'</div>';
    for(var j=0;j<d.caveats.length;j++)
      h += '<div class="vsc">&bull; '+esc(d.caveats[j])+'</div>';
    h += '</div>';
  }
  var cov = null, cap = null;
  for(var i=0;i<sb.dimensions.length;i++){
    if(sb.dimensions[i].dimension==='coverage') cov = sb.dimensions[i];
    if(sb.dimensions[i].dimension==='capability parity') cap = sb.dimensions[i];
  }
  if(cov){
    h += '<h1 style="font-size:15px;margin:16px 0 6px">per service</h1>'+
         '<table><thead><tr><th>RF</th><th>service</th><th>vchan</th>'+
         '<th>ours</th><th>theirs</th><th></th></tr></thead><tbody>';
    cov.detail.rows.forEach(function(r){
      h += '<tr><td>'+esc(r.rf||'-')+'</td><td>'+esc(r.name)+'</td><td>'+
           esc(r.vchannel||'-')+'</td><td>'+esc(r.ours)+'</td><td>'+
           esc(r.theirs)+'</td><td>'+vsbadge(r.verdict)+'</td></tr>';
    });
    h += '</tbody></table>';
  }
  if(cap){
    h += '<h1 style="font-size:15px;margin:16px 0 6px">capability</h1>'+
         '<table><thead><tr><th></th><th>capability</th><th>ours</th>'+
         '<th>theirs</th><th>evidence</th></tr></thead><tbody>';
    cap.detail.rows.forEach(function(r){
      h += '<tr><td>'+vsbadge(r.verdict)+'</td><td>'+esc(r.capability)+
           '</td><td>'+esc(r.ours)+'</td><td>'+esc(r.theirs)+'</td>'+
           '<td class="vsc">'+esc(r.evidence)+'</td></tr>';
    });
    h += '</tbody></table>';
  }
  if(sb.corroboration && sb.corroboration.length){
    h += '<h1 style="font-size:15px;margin:16px 0 6px">independent '+
         'corroboration <span class="vsq">(their L1 parser checking our '+
         'claims &mdash; unscored)</span></h1><table><tbody>';
    sb.corroboration.forEach(function(c){
      h += '<tr><td>RF'+esc(c.rf)+'</td><td>'+esc(c.check)+'</td><td>ours '+
           esc(c.ours_claim)+'</td><td>flex '+esc(c.flex_reports)+'</td><td>'+
           (c.agrees===null?'<span class="cont">n/a</span>'
             :(c.agrees?'<span class="good">agrees</span>'
                       :'<span class="bad">DISAGREES</span>'))+
           '</td><td class="vsc">'+esc(c.our_source)+'</td></tr>';
    });
    h += '</tbody></table>';
  }
  document.getElementById('wrap').innerHTML = h;
}
function load(){
  if(TAB==='vs'){
    fetch('/api/vs').then(function(x){return x.json();}).then(renderVS);
    return;
  }
  fetch('/api/grid'+(AT?('?at='+encodeURIComponent(AT)):''))
    .then(function(x){return x.json();}).then(function(g){
      if(window._b===undefined) window._b = g.build;
      else if(g.build!==window._b)
        document.getElementById('sub').textContent =
          '\\u21bb this page was updated behind your tab \\u2014 press '+
          'Ctrl+Shift+R to load the new version';
      render(g);
    });
}
setTab(TAB);
setInterval(load, 300000);
</script></body></html>
"""


# ------------------------------------------------------------------ server

class H(BaseHTTPRequestHandler):
    sources = None

    def log_message(self, *a):
        pass

    def _send(self, body, ctype="application/json", code=200):
        data = body.encode("utf-8") if isinstance(body, str) else body
        self.send_response(code)
        self.send_header("Content-Type", ctype)
        self.send_header("Content-Length", str(len(data)))
        self.send_header("Cache-Control", "no-store")
        self.end_headers()
        try:
            self.wfile.write(data)
        except OSError:
            pass

    def do_GET(self):
        path, _, q = self.path.partition("?")
        args = {}
        for kv in q.split("&"):
            if "=" in kv:
                k, v = kv.split("=", 1)
                args[k] = _unquote(v)
        try:
            if path == "/":
                return self._send(PAGE, "text/html; charset=utf-8")
            if path == "/api/grid":
                g = grid_json(self.sources, at=args.get("at") or None,
                              slots=int(args.get("slots", DEFAULT_SLOTS)))
                return self._send(json.dumps(g, default=str))
            if path == "/api/service":
                d = service_json(args.get("bsid"), args.get("sid"),
                                 self.sources)
                if d is None:
                    return self._send(json.dumps({"error": "no such service"}),
                                      code=404)
                return self._send(json.dumps(I.strip_private(d), default=str))
            if path == "/api/vs":
                return self._send(json.dumps(get_scoreboard(), default=str))
            if path == "/api/health":
                return self._send(json.dumps(
                    {"build": BUILD, "picture_policy": PICTURE_POLICY}))
        except Exception as exc:                       # never drop the socket
            return self._send(json.dumps({"error": repr(exc)}), code=500)
        self.send_error(404)


def _unquote(s):
    s = s.replace("+", " ")
    out, i = [], 0
    while i < len(s):
        if s[i] == "%" and i + 2 < len(s):
            try:
                out.append(chr(int(s[i + 1:i + 3], 16)))
                i += 3
                continue
            except ValueError:
                pass
        out.append(s[i])
        i += 1
    return "".join(out)


class Grid(ThreadingHTTPServer):
    request_queue_size = 128
    daemon_threads = True
    allow_reuse_address = True


# ------------------------------------------------------------ static export

def render_static(g, out_path):
    """One self-contained HTML file: same grid, no server, no fetch.

    For handing someone a snapshot, and for the gate -- a static render is
    exactly what an automated check can assert against.
    """
    rows = []
    for r in g["rows"]:
        cells = []
        for c in r["cells"]:
            if c.get("empty"):
                cells.append('<td class="empty">&middot;</td>')
            else:
                t = ('<span class="cont">&raquo; ' if c.get("cont")
                     else '<span class="show">') + html.escape(c["title"]) \
                    + "</span>"
                if c.get("rating"):
                    t += '<div class="rating">%s</div>' % html.escape(
                        c["rating"])
                cells.append("<td>%s</td>" % t)
        meta = " &middot; ".join(x for x in [
            html.escape(r["transport"]),
            ("PLP%s" % r["plp"]) if r["plp"] else "",
            ("RF%s" % r["rf"]) if r["rf"] else "",
            html.escape(r["video"] or ""),
            html.escape("/".join(r["audio"])) if r["audio"] else ""] if x)
        cap = ('<div class="cap">captions available (%d cues)</div>'
               % r["captions_n"]) if r["captions_n"] else ""
        rows.append(
            '<tr class="svc%s"><td class="chan"><span class="vch">%s</span> '
            '<span class="call">%s</span>%s'
            '<div style="margin-top:3px"><span class="badge b-%s">%s</span>'
            '</div><div class="meta">%s</div>%s</td>%s</tr>' % (
                " lockrow" if r["state"] == "PROTECTED" else "",
                html.escape(r["vchannel"]), html.escape(r["name"]),
                ' <span class="badge b-ok">LIVE</span>' if r["live"] else "",
                r["badge_class"], html.escape(r["badge"]), meta, cap,
                "".join(cells)))
    head = "".join("<th>%s</th>" % html.escape(s) for s in g["slots"])
    banner = html.escape(g["message"]) if g["message"] else (
        "Guide window %s .. %s" % (_ts(g["guide"]["span"][0]),
                                   _ts(g["guide"]["span"][1]))
        if g["guide"]["span"] else "no guide")
    page = PAGE.split("<script>")[0]
    page = page.replace('<div id="wrap"></div>',
                        '<div id="wrap"><table><thead><tr>'
                        '<th style="width:170px">channel</th>%s</tr></thead>'
                        '<tbody>%s</tbody></table></div>'
                        % (head, "".join(rows)))
    page = page.replace('<div class="sub" id="sub">loading&hellip;</div>',
                        '<div class="sub" id="sub">%d services &bull; guide '
                        '%d programmes &bull; build %s &bull; static export '
                        '%s</div>' % (len(g["rows"]), g["guide"]["programmes"],
                                      BUILD, _ts(time.time())))
    page = page.replace('<div id="banner"></div>',
                        '<div id="banner" class="%s">%s</div>'
                        % ("warn" if g["status"] != "live" else "", banner))
    page += "</body></html>\n"
    with open(out_path, "w", encoding="utf-8") as fh:
        fh.write(page)
    return out_path


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--port", type=int, default=PORT)
    ap.add_argument("--bind", default=os.environ.get("ATSC3_GRID_BIND",
                                                     "127.0.0.1"))
    ap.add_argument("--source", action="append", default=None)
    ap.add_argument("--at", default=None,
                    help="unix time, or 'guide' to anchor on the banked guide")
    ap.add_argument("--slots", type=int, default=DEFAULT_SLOTS)
    ap.add_argument("--out", default=None,
                    help="write one self-contained HTML file and exit")
    ap.add_argument("--json", action="store_true")
    a = ap.parse_args(argv)
    I.be_polite()
    if a.out or a.json:
        g = grid_json(a.source, at=a.at, slots=a.slots)
        if a.json:
            json.dump(g, sys.stdout, indent=2, default=str)
            sys.stdout.write("\n")
        if a.out:
            print("wrote %s" % render_static(g, a.out), flush=True)
        return 0
    H.sources = a.source
    print("ATSC 3.0 guide grid: http://%s:%d/   build %s"
          % (a.bind, a.port, BUILD), flush=True)
    Grid((a.bind, a.port), H).serve_forever()
    return 0


if __name__ == "__main__":
    sys.exit(main())

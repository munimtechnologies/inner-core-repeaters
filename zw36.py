#!/usr/bin/env python3
# Copyright 2026 MV2, Munim, Inc.
# Licensed under the Apache License, Version 2.0 (Apache-2.0). See LICENSE.
# Made with AI agents working in MV2's name.
"""zw36.py - pre-registration 36: forward rotation closure on multiplet F, 2024-2026.

Pre-registered in findings/2026-11-24-prereg36-multiplet-f-forward-closure-plan.md, committed
(d12433b) before any waveform of the 11 new events was fetched, decoded or looked at. This code
was written on 2026-11-25, still before any such waveform was read and before zw33-v1 or
zw35g2-v1 had returned anything.

Steps (all fixed in the plan):
  C  core phase: the zw33.py 0.1.0 measurement (pinned commit 7f92852), unchanged, on short-
     period vertical records of YKA (CN YK*) and ILAR (IM IL*) fetched from EarthScope FDSN for
     the 3 old + 11 new (+2 report-only) multiplet-F events. Graded window D_c310, 1-2 Hz.
  S  source: mantlep2.py 0.13.0 (pinned b4712f0e) on the 8 M03-24 events + 3 old F + 13 new,
     scored per pair with the G2 rule of pre-registration 35 (q >= 0.93 hit; network vote).
  P36a  forward closure, closure35.py 0.1.0 (pinned 39b2f1c), qualifying triples only.
  P36b  forward quiescence at YKA: same(2022-05-19, n) for new n with source N.

Readings of the plan made in code, before any data (stated here so they are public):
  * Event keys are UTC dates (YYYY-MM-DD); all 16 dates differ.
  * ILAR hardware flag: raised when at least one measured (same/changed) ILAR pair spans
    2018-09-25, every measured spanning pair is 'changed', and no measured within-epoch pair is
    'changed' (an empty within-epoch set does not exclude hardware). The M03-24 part of the flag
    (six events, from zw33-v1) is added by hand at grading, because zw33-v1's output is not
    visible inside this job.
  * The zw35g2-v1 decision on (2005, 2012) is likewise compared by hand at grading
    (combine_external); a disagreement makes that pair 'inconclusive'.
  * Old F members are fetched from EarthScope with the same code as the new ones; zw33-v1 used
    the authors' Zenodo files. Both are reported.

Top-level imports are stdlib only (the tests run without numpy/obspy).
"""

import argparse
import csv
import json
import math
import os
import sys
import time
import traceback
import urllib.request

VERSION = "0.1.0"
NAN = float("nan")

REPO = "https://raw.githubusercontent.com/munimtechnologies/inner-core-repeaters/"
ZW33_COMMIT = "7f928523dbdb5c1e5006110cbea1cdc3c2334e20"       # zw33.py 0.1.0
CLOSURE35_COMMIT = "39b2f1cc05a231fd3b9a4cd71b133ced6285d9c1"  # closure35.py 0.1.0
MANTLEP2_COMMIT = "b4712f0ebbc70b896231620c905ff53a4b8270fe"   # mantlep2.py 0.13.0
PREREG_URL = REPO + "main/findings/2026-11-24-prereg36-multiplet-f-forward-closure-plan.md"
UA = "inner-core-repeaters (MV2, Munim, Inc.)"
FDSN = "https://service.earthscope.org"

OLD8 = ("usp000bt5x", "usp000fz4d", "usp000hute", "usc000t4wj", "us7000am0e", "us6000h2wz",
        "us6000mk65", "us6000runq")
# USGS ComCat: (id, origin UTC, lat, lon, depth km, mag)
F_OLD = [
    ("usp000dsnt", "2005-06-10T13:08:12.920", -58.706, -25.207, 9.8, 4.7),
    ("usp000jwcf", "2012-12-03T03:11:10.790", -58.864, -25.565, 41.2, 4.9),
    ("us6000hm6g", "2022-05-19T02:24:53.021", -58.796, -25.3186, 35.0, 5.1),
]
# The 11 graded new events: the table of pre-registration 36 (ComCat query of 2026-11-24).
F_NEW = [
    ("us7000md2s", "2024-04-19T04:22:13.351", -58.7859, -25.1511, 31.6, 4.7),
    ("us6000mt07", "2024-04-22T17:21:22.541", -58.7238, -25.3084, 28.7, 4.9),
    ("us7000myv7", "2024-07-12T11:48:07.064", -58.8516, -25.2828, 33.2, 4.9),
    ("us7000plsg", "2025-03-22T04:24:23.383", -58.7555, -25.2674, 30.7, 4.9),
    ("us6000rc8b", "2025-09-23T05:36:33.685", -58.7722, -25.3301, 33.9, 5.1),
    ("us6000rk5n", "2025-10-29T13:28:34.928", -58.7537, -25.1169, 35.0, 4.8),
    ("us6000rtg8", "2025-12-08T18:41:15.496", -58.8007, -25.4875, 30.0, 5.5),
    ("us7000rmeh", "2026-01-04T10:35:54.685", -58.6920, -25.3448, 35.0, 5.0),
    ("us6000sg9g", "2026-03-14T08:34:34.421", -58.7514, -25.4670, 35.0, 5.1),
    ("us7000shsi", "2026-05-02T10:31:45.258", -58.7598, -25.2720, 26.4, 5.0),
    ("us6000tmzw", "2026-08-23T17:23:29.488", -58.7568, -25.2117, 35.0, 4.9),
]
# Report-only (just outside 15 km): ComCat values read 2026-11-25.
F_REPORT = [
    ("us6000rvhe", "2025-12-17T17:05:57.315", -58.6643, -25.2234, 35.0, 5.6),
    ("us7000t2ws", "2026-07-24T00:37:03.351", -58.6939, -25.1649, 35.0, 4.9),
]
F2022_DATE = "2022-05-19"
F2005_DATE = "2005-06-10"
F2012_DATE = "2012-12-03"
ILAR_HW_DATE = "2018-09-25"

HIT_Q = 0.93
MIN_REF_PAIRS = 10
MIN_STATIONS = 3
ARRAY_REQ = {"YKA": ("CN", "YK*"), "ILAR": ("IM", "IL*")}
FETCH_MARGIN_S = 90.0


# ----------------------------------------------------------------------------- pure helpers
def isf(x):
    try:
        return math.isfinite(float(x))
    except (TypeError, ValueError):
        return False


def key(a, b):
    return (a, b) if a <= b else (b, a)


def date_of(ev):
    return ev[1][:10]


def graded_dates():
    return [date_of(e) for e in F_OLD + F_NEW]


def new_dates():
    return set(date_of(e) for e in F_NEW)


def station_pair_score(D, gated, ref_ids, a, c, hit_q=HIT_Q, min_ref=MIN_REF_PAIRS):
    """As zw35g2.station_pair_score 0.1.0. D: {(e1, e2) sorted: D_late}; gated: {event: bool}."""
    out = dict(counted=False, q=None, hit=None, d_pair=None, n_ref=0)
    ref = [e for e in ref_ids if gated.get(e)]
    vals = []
    for i in range(len(ref)):
        for j in range(i + 1, len(ref)):
            v = D.get(key(ref[i], ref[j]))
            if isf(v):
                vals.append(float(v))
    out["n_ref"] = len(vals)
    if not (gated.get(a) and gated.get(c)):
        out["reason"] = "event_not_gated"
        return out
    d = D.get(key(a, c))
    if not isf(d):
        out["reason"] = "no_pair_D"
        return out
    out["d_pair"] = float(d)
    if len(vals) < min_ref:
        out["reason"] = "too_few_reference_pairs"
        return out
    q = sum(1 for v in vals if v < float(d)) / float(len(vals))
    out.update(counted=True, q=q, hit=bool(q >= hit_q - 1e-12))
    return out


def decide(hits, n_min=MIN_STATIONS):
    """Pre-registration 11's network vote on per-station hit flags."""
    n = len(hits)
    h = sum(1 for x in hits if x)
    if n < n_min:
        return "untestable", h
    if h >= int(math.ceil(n / 2.0)):
        return "S", h
    if h <= max(1, int(math.floor(0.2 * n))):
        return "N", h
    return "inconclusive", h


def source_decisions(scores):
    """scores: {station: {(date_a, date_b): score dict}} -> {(date_a, date_b): dict}."""
    pairs = sorted({p for sc in scores.values() for p in sc})
    out = {}
    for p in pairs:
        rows = [(st, sc[p]) for st, sc in sorted(scores.items()) if sc.get(p, {}).get("counted")]
        dec, h = decide([r["hit"] for _, r in rows])
        out[p] = dict(decision=dec, n_stations=len(rows), hits=h,
                      q=[round(r["q"], 4) for _, r in rows])
    return out


def combine_external(own, external):
    """The (2005, 2012) rule: if zw35g2-v1 graded the pair too, a disagreement -> inconclusive."""
    if external in (None, "", "untestable") or own == external:
        return own
    if own == "untestable":
        return own
    return "inconclusive"


def src_of(src, a, b):
    d = src.get(key(a, b))
    if isinstance(d, dict):
        return d.get("decision", "untestable")
    return d or "untestable"


def p36a_grade(triple_rows, src, new):
    """triple_rows: closure35.grade rows (recovery_channel, a, b, c, other_channel, grade).
    src: {(date, date): decision or dict}. new: set of new-event dates.
    Returns (qualifying rows, summary)."""
    qual = []
    for r in triple_rows:
        a, b, c = r["a"], r["b"], r["c"]
        if not ({a, b, c} & new):
            continue
        if src_of(src, a, b) != "N" or src_of(src, a, c) != "N":
            continue
        qual.append(dict(r))
    measured = [r for r in qual if r["grade"] in ("holds", "violated")]
    ev = {}
    for r in measured:
        for e in (r["a"], r["b"], r["c"]):
            if e in new:
                ev.setdefault(e, set()).add(r["grade"])
    v_ev = sorted(e for e, g in ev.items() if g == {"violated"})
    h_ev = sorted(e for e, g in ev.items() if g == {"holds"})
    m_ev = sorted(e for e, g in ev.items() if len(g) == 2)
    V, H = len(v_ev), len(h_ev)
    if not measured:
        grade = "untestable"
    elif V >= 2 and V > H:
        grade = "broken forward"
    elif H >= 2 and V == 0:
        grade = "holds forward"
    else:
        grade = "mixed"
    return qual, dict(grade=grade, n_qualifying=len(qual), n_measured=len(measured),
                      V_events=v_ev, H_events=h_ev, M_events=m_ev, V=V, H=H, M=len(m_ev))


def hardware_flag(ilar_labels, hw_date=ILAR_HW_DATE):
    """ilar_labels: {(date_a, date_b): label}. See module docstring for the reading."""
    span, within = [], []
    for (a, b), lab in ilar_labels.items():
        if lab not in ("same", "changed"):
            continue
        a2, b2 = key(a, b)
        (span if (a2 < hw_date <= b2) else within).append(lab)
    flag = bool(span) and all(l == "changed" for l in span) and all(l == "same" for l in within)
    return flag, dict(n_span=len(span), span_changed=sum(1 for l in span if l == "changed"),
                      n_within=len(within), within_same=sum(1 for l in within if l == "same"))


def apply_hw(grade, flag):
    if grade == "broken forward" and flag:
        return "broken forward, ILAR hardware not excluded"
    return grade


def p36b_grade(yka_labels, src, new, ref=F2022_DATE):
    """YKA same(2022, n) for new n with source N on (2022, n)."""
    items = []
    for n in sorted(new):
        if src_of(src, ref, n) != "N":
            continue
        items.append(dict(n=n, label=yka_labels.get(key(ref, n), "missing")))
    lab = [i["label"] for i in items if i["label"] in ("same", "changed")]
    ns, nc = lab.count("same"), lab.count("changed")
    if len(lab) < 2:
        grade = "untestable"
    elif nc == 0:
        grade = "holds"
    elif nc >= 2 and nc >= ns:
        grade = "broken"
    else:
        grade = "mixed"
    return items, dict(grade=grade, n_source_N=len(items), n_measured=len(lab), same=ns,
                       changed=nc)


def labels_from_pairs(pair_rows, array, band, window, dates=None):
    out = {}
    for r in pair_rows:
        if r["array"] != array or r["band"] != band:
            continue
        if dates is not None and (r["a"] not in dates or r["b"] not in dates):
            continue
        out[key(r["a"], r["b"])] = r.get("label_" + window, "nd")
    return out


def closure_input(pair_rows, band, window, dates):
    return [dict(channel=r["array"], event_i=r["a"], event_j=r["b"],
                 label=r.get("label_" + window, "nd"))
            for r in pair_rows
            if r["band"] == band and r["a"] in dates and r["b"] in dates]


def jsonable(o):
    if isinstance(o, dict):
        return {(str(k) if not isinstance(k, tuple) else "|".join(map(str, k))): jsonable(v)
                for k, v in o.items()}
    if isinstance(o, (list, tuple, set)):
        return [jsonable(v) for v in (sorted(o) if isinstance(o, set) else o)]
    if isinstance(o, float):
        return o if math.isfinite(o) else None
    if hasattr(o, "item") and not isinstance(o, (str, bytes)):
        try:
            return jsonable(o.item())
        except Exception:  # noqa: BLE001
            return str(o)
    return o


def write_csv(path, rows):
    keys = []
    for r in rows:
        for k in r:
            if k not in keys:
                keys.append(k)
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=keys or ["empty"])
        w.writeheader()
        for r in rows:
            w.writerow({k: (json.dumps(jsonable(v)) if isinstance(v, (list, dict, set)) else
                            jsonable(v)) for k, v in r.items()})


# ----------------------------------------------------------------------------- I/O
def fetch(url, tries=4, wait=10):
    for k in range(tries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": UA})
            with urllib.request.urlopen(req, timeout=60) as r:
                if r.status == 200:
                    return r.read().decode("utf-8", "replace")
        except Exception as ex:  # noqa: BLE001
            print("fetch error", url, repr(ex)[:200], flush=True)
        time.sleep(wait * (k + 1))
    return None


def load_pinned(dest, name, commit, version):
    os.makedirs(dest, exist_ok=True)
    p = os.path.join(dest, name + ".py")
    if not os.path.exists(p):
        body = fetch(REPO + commit + "/" + name + ".py")
        if not body or ('VERSION = "%s"' % version) not in body:
            return None
        compile(body, name + ".py", "exec")
        with open(p, "w", encoding="utf-8") as f:
            f.write(body)
    if dest not in sys.path:
        sys.path.insert(0, dest)
    return __import__(name)


def event_dict(e):
    from obspy import UTCDateTime
    return dict(id=e[0], date=date_of(e), origin=UTCDateTime(e[1]), lat=e[2], lon=e[3],
                depth=e[4], mag=e[5])


def fetch_array(client, zw, model, array, ev, cache):
    """Raw SHZ stream and channel inventory for one array and event (cached)."""
    from obspy import read, read_inventory
    from obspy.geodetics import locations2degrees
    mpath = os.path.join(cache, "%s_%s.mseed" % (array, ev["id"]))
    ipath = os.path.join(cache, "%s_%s.xml" % (array, ev["id"]))
    if os.path.exists(mpath) and os.path.exists(ipath):
        return read(mpath), read_inventory(ipath), "cache"
    if os.path.exists(mpath + ".none"):
        return None, None, "no data (cached)"
    clat, clon = zw.ARRAY_CENTER[array]
    dist = locations2degrees(ev["lat"], ev["lon"], clat, clon)
    tt, _ph = zw.pkikp_time(model, ev["depth"], dist)
    if tt is None:
        return None, None, "no PKIKP"
    t1 = ev["origin"] + tt - zw.PRE - FETCH_MARGIN_S
    t2 = ev["origin"] + tt + zw.POST + FETCH_MARGIN_S
    net, sta = ARRAY_REQ[array]
    last = ""
    for k in range(3):
        try:
            st = client.get_waveforms(net, sta, "*", "SHZ", t1, t2)
            inv = client.get_stations(network=net, station=sta, channel="SHZ", starttime=t1,
                                      endtime=t2, level="channel")
            st.write(mpath, format="MSEED")
            inv.write(ipath, format="STATIONXML")
            return st, inv, "fetched"
        except Exception as ex:  # noqa: BLE001
            last = repr(ex)[:160]
            if "No data" in last or "204" in last:
                open(mpath + ".none", "w").close()
                return None, None, "no data"
            time.sleep(20 * (k + 1))
    return None, None, "fetch failed " + last


def process_stream(zw, st, inv, ev, array, band, model):
    """zw33.process_folder, applied to an EarthScope stream instead of SAC files."""
    from obspy import Stream
    from obspy.geodetics import locations2degrees
    import numpy as np
    traces, meta, by_el = {}, {}, {}
    for tr in st:
        el = zw.canonical_element(tr.stats.station)
        if el is None or zw.element_array(el) != array or not tr.stats.channel.upper().startswith("SH"):
            continue
        by_el.setdefault(el, []).append(tr)
    for el in sorted(by_el):
        locs = sorted({t.stats.location for t in by_el[el]})
        sub = Stream([t.copy() for t in by_el[el] if t.stats.location == locs[0]])
        try:
            sub.merge(method=0)
        except Exception as ex:  # noqa: BLE001
            meta[el] = "merge_error " + repr(ex)[:60]
            continue
        if len(sub) != 1:
            meta[el] = "merge left %d traces" % len(sub)
            continue
        tr = sub[0]
        if np.ma.isMaskedArray(tr.data):
            if np.ma.getmaskarray(tr.data).any():
                meta[el] = "gap"
                continue
            tr.data = np.asarray(tr.data)
        try:
            c = inv.get_coordinates(tr.id, ev["origin"])
            stla, stlo = c["latitude"], c["longitude"]
        except Exception:  # noqa: BLE001
            stla, stlo = zw.ARRAY_CENTER[array]
        dist = locations2degrees(ev["lat"], ev["lon"], stla, stlo)
        tt, ph = zw.pkikp_time(model, ev["depth"], dist)
        if tt is None:
            meta[el] = "no PKIKP"
            continue
        t0 = ev["origin"] + tt - zw.PRE
        if tr.stats.starttime > t0 - 2.0 or tr.stats.endtime < t0 + zw.PRE + zw.POST + 2.0:
            meta[el] = "record does not cover the window"
            continue
        try:
            tr.data = tr.data.astype("float64")
            tr.detrend("demean")
            tr.detrend("linear")
            tr.taper(0.05)
            tr.filter("bandpass", freqmin=band[0], freqmax=band[1], corners=4, zerophase=False)
            tr.interpolate(sampling_rate=zw.FS, method="lanczos", a=12, starttime=t0, npts=zw.N)
        except Exception as ex:  # noqa: BLE001
            meta[el] = "processing_error " + repr(ex)[:80]
            continue
        x = np.asarray(tr.data, dtype=float)
        if len(x) != zw.N or not np.all(np.isfinite(x)) or float(np.std(x)) == 0.0:
            meta[el] = "bad samples"
            continue
        traces[el] = x
        meta[el] = "ok %s loc=%s %.3f deg %s %.2f s fs=%.1f" % (tr.stats.station, locs[0], dist,
                                                                ph, tt, tr.stats.sampling_rate)
    return traces, meta


# ----------------------------------------------------------------------------- steps
def core_step(zw, model, cache, s):
    from obspy.clients.fdsn import Client
    client = Client(base_url=FDSN, _discover_services=False, user_agent=UA, timeout=120)
    events = [event_dict(e) for e in F_OLD + F_NEW + F_REPORT]
    raw = {}
    fetch_rows = []
    for arr in ("YKA", "ILAR"):
        for ev in events:
            st, inv, how = fetch_array(client, zw, model, arr, ev, cache)
            raw[(arr, ev["date"])] = (st, inv)
            fetch_rows.append(dict(array=arr, event=ev["id"], date=ev["date"], how=how,
                                   n_traces=len(st) if st is not None else 0))
            print("FETCH", arr, ev["date"], how, len(st) if st is not None else 0, flush=True)
    pair_rows, event_rows = [], []
    for band in (zw.BAND,) + tuple(zw.BANDS_REPORT):
        bl = "%.1f-%.1f" % band
        for arr in ("YKA", "ILAR"):
            stacks = {}
            for ev in events:
                st, inv = raw[(arr, ev["date"])]
                row = dict(band=bl, array=arr, date=ev["date"], event=ev["id"],
                           new=ev["date"] in new_dates(),
                           report_only=ev["id"] in [e[0] for e in F_REPORT])
                if st is None:
                    row.update(usable=False, n_read=0)
                    event_rows.append(row)
                    continue
                tr, meta = process_stream(zw, st, inv, ev, arr, band, model)
                stk, info = zw.event_stack(tr)
                row.update(info)
                row["usable"] = len(stk) >= zw.MIN_COMMON
                row["element_notes"] = " | ".join("%s: %s" % kv for kv in sorted(meta.items()))
                event_rows.append(row)
                if stk:
                    stacks[ev["date"]] = stk
            dates = sorted(stacks)
            for i in range(len(dates)):
                for j in range(i + 1, len(dates)):
                    a, b = dates[i], dates[j]
                    try:
                        r = zw.measure_pair(stacks[a], stacks[b])
                    except Exception as ex:  # noqa: BLE001
                        traceback.print_exc()
                        r = dict(ok_direct=False, reason="exception " + repr(ex)[:120])
                    row = dict(band=bl, array=arr, a=a, b=b)
                    for k2, _w in zw.WINDOWS:
                        r.setdefault("label_" + k2, "nd")
                    row.update(r)
                    pair_rows.append(row)
                    print("PAIR", bl, arr, a, b, r.get("label_" + zw.PRIMARY), flush=True)
    s["fetch"] = fetch_rows
    return pair_rows, event_rows


def source_step(mp, cache, plots, budget_min, t0, s):
    import requests
    from obspy.taup import TauPyModel
    fam = F_OLD + F_NEW + F_REPORT
    mp.EVENTS = list(mp.EVENTS) + list(fam)
    mp.IDS = [e[0] for e in mp.EVENTS]
    mp.DATES = {e[0]: e[1][:10] for e in mp.EVENTS}
    id2date = {e[0]: date_of(e) for e in fam}
    session = requests.Session()
    session.headers["User-Agent"] = UA
    found = None
    for _ in range(3):
        found = mp.discover(session)
        if found:
            break
        time.sleep(30)
    source = "station-service"
    if not found:
        found, source = mp.FALLBACK_STATIONS, "fallback"
    cand = []
    for net, sta, lat, lon in found:
        d = mp.gc_deg(mp.SRC_LAT, mp.SRC_LON, lat, lon)
        if mp.DIST_RANGE[0] <= d <= mp.DIST_RANGE[1]:
            cand.append((d, net, sta, lat, lon))
    cand.sort()
    cand = cand[:mp.MAX_STATIONS]
    s.update(station_source=source, n_candidate_stations=len(cand))
    model = TauPyModel("ak135")
    scores, stations, skipped, ev_all, pr_all = {}, [], [], [], []
    ids = [e[0] for e in fam]
    for d, net, sta, lat, lon in cand:
        name = "%s.%s" % (net, sta)
        if (time.time() - t0) / 60.0 > budget_min:
            skipped.append(name)
            continue
        ev_rows, pair_rows = [], []
        try:
            res = mp.run_station(session, model, net, sta, lat, lon, cache, plots, ev_rows,
                                 pair_rows, {})
        except Exception as ex:  # noqa: BLE001
            traceback.print_exc()
            res = dict(status="exception", error=repr(ex)[:300])
        gated = {r["event"]: str(r.get("gated")) == "True" for r in ev_rows}
        D = {key(r["ev1"], r["ev2"]): float(r["D_late"]) for r in pair_rows if isf(r.get("D_late"))}
        sc = {}
        for i in range(len(ids)):
            for j in range(i + 1, len(ids)):
                pa = key(id2date[ids[i]], id2date[ids[j]])
                sc[pa] = station_pair_score(D, gated, OLD8, ids[i], ids[j])
        scores[name] = sc
        stations.append(dict(station=name, dist_deg=round(d, 2), status=res.get("status"),
                             n_gated_old=sum(1 for e in OLD8 if gated.get(e)),
                             gated_family=sorted(id2date[e] for e in ids if gated.get(e))))
        ev_all.extend(ev_rows)
        pr_all.extend(pair_rows)
        print("STATION", name, res.get("status"), flush=True)
    s.update(source_stations=stations, source_skipped_time_budget=skipped)
    return scores, ev_all, pr_all


# ----------------------------------------------------------------------------- main
def run(args):
    t0 = time.time()
    os.makedirs(args.out_dir, exist_ok=True)
    s = dict(version=VERSION, started=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
             pins=dict(zw33=ZW33_COMMIT, closure35=CLOSURE35_COMMIT, mantlep2=MANTLEP2_COMMIT),
             events=dict(old=F_OLD, new=F_NEW, report_only=F_REPORT))
    s["prereg_public_at_start"] = bool(fetch(PREREG_URL, tries=2, wait=5))
    pin = "/tmp/zw36_pinned"
    zw = load_pinned(pin, "zw33", ZW33_COMMIT, "0.1.0")
    cl = load_pinned(pin, "closure35", CLOSURE35_COMMIT, "0.1.0")
    mp = load_pinned(pin, "mantlep2", MANTLEP2_COMMIT, "0.13.0")
    if zw is None or cl is None or mp is None or [e[0] for e in mp.EVENTS] != list(OLD8):
        s["error"] = "pinned code could not be loaded or differs"
        with open(os.path.join(args.out_dir, "zw36_summary.json"), "w") as f:
            json.dump(jsonable(s), f, indent=1)
        return 2
    if (zw.SAME_MAX_RATIO, zw.CHANGED_MIN_RATIO, zw.FLOOR_MIN) != (
            cl.SAME_MAX_RATIO, cl.CHANGED_MIN_RATIO, cl.FLOOR_MIN):
        s["error"] = "zw33 and closure35 labelling constants differ"
        with open(os.path.join(args.out_dir, "zw36_summary.json"), "w") as f:
            json.dump(jsonable(s), f, indent=1)
        return 2

    from obspy.taup import TauPyModel
    model = TauPyModel("ak135")
    wcache = os.path.join(args.data_dir, "zw36", "core")
    os.makedirs(wcache, exist_ok=True)
    pair_rows, event_rows = [], []
    try:
        pair_rows, event_rows = core_step(zw, model, wcache, s)
    except Exception as ex:  # noqa: BLE001
        traceback.print_exc()
        s["core_error"] = repr(ex)[:400]
    write_csv(os.path.join(args.out_dir, "zw36_core_pairs.csv"), pair_rows)
    write_csv(os.path.join(args.out_dir, "zw36_core_events.csv"), event_rows)

    scache = os.path.join(args.data_dir, "zw36", "mantlep")
    plots = os.path.join(args.out_dir, "plots")
    os.makedirs(scache, exist_ok=True)
    os.makedirs(plots, exist_ok=True)
    scores, sev, spr = {}, [], []
    try:
        scores, sev, spr = source_step(mp, scache, plots, args.time_budget_min, t0, s)
    except Exception as ex:  # noqa: BLE001
        traceback.print_exc()
        s["source_error"] = repr(ex)[:400]
    write_csv(os.path.join(args.out_dir, "zw36_source_events.csv"), sev)
    write_csv(os.path.join(args.out_dir, "zw36_source_pairs.csv"), spr)
    rows = []
    for st, sc in sorted(scores.items()):
        for (a, b), v in sorted(sc.items()):
            r = dict(station=st, a=a, b=b)
            r.update({k: v.get(k) for k in ("counted", "q", "hit", "d_pair", "n_ref", "reason")})
            rows.append(r)
    write_csv(os.path.join(args.out_dir, "zw36_source_station_scores.csv"), rows)
    src = source_decisions(scores)
    write_csv(os.path.join(args.out_dir, "zw36_source_decisions.csv"),
              [dict(a=a, b=b, **v) for (a, b), v in sorted(src.items())])

    prim = "%.1f-%.1f" % zw.BAND
    gd = set(graded_dates())
    new = new_dates()
    results = {}
    for window, _w in zw.WINDOWS:
        trip, clos = cl.grade(cl.build_labels(closure_input(pair_rows, prim, window, gd)))
        qual, p36a = p36a_grade(trip, src, new)
        ilar = labels_from_pairs(pair_rows, "ILAR", prim, window, gd)
        flag, flag_info = hardware_flag(ilar)
        p36a["grade_with_hw_flag"] = apply_hw(p36a["grade"], flag)
        p36a.update(hw_flag_F=flag, hw_info_F=flag_info)
        yka = labels_from_pairs(pair_rows, "YKA", prim, window, gd)
        items, p36b = p36b_grade(yka, src, new)
        ilar_items, ilar_rep = p36b_grade(ilar, src, new)
        results[window] = dict(graded=(window == zw.PRIMARY), closure35=clos, P36a=p36a,
                               P36b=p36b, P36b_items=items, report_ILAR_2022_n=ilar_rep,
                               report_ILAR_items=ilar_items)
        if window == zw.PRIMARY:
            write_csv(os.path.join(args.out_dir, "zw36_closure_triples.csv"), trip)
            write_csv(os.path.join(args.out_dir, "zw36_p36a_qualifying.csv"), qual)
    # report-only: all 16 events, and the 0.5-2 Hz band
    alld = gd | set(date_of(e) for e in F_REPORT)
    rep = {}
    trip16, clos16 = cl.grade(cl.build_labels(closure_input(pair_rows, prim, zw.PRIMARY, alld)))
    rep["closure_16_events"] = dict(closure35=clos16, P36a=p36a_grade(trip16, src, new | set(
        date_of(e) for e in F_REPORT))[1])
    for band in zw.BANDS_REPORT:
        bl = "%.1f-%.1f" % band
        tripb, closb = cl.grade(cl.build_labels(closure_input(pair_rows, bl, zw.PRIMARY, gd)))
        rep[bl] = dict(closure35=closb, P36a=p36a_grade(tripb, src, new)[1],
                       P36b=p36b_grade(labels_from_pairs(pair_rows, "YKA", bl, zw.PRIMARY, gd),
                                       src, new)[1])
    pf = src.get(key(F2005_DATE, F2012_DATE), {})
    s.update(results=results, report_only=rep,
             source_2005_2012=dict(this_run=pf.get("decision", "untestable"),
                                   note="compare with zw35g2-v1 by hand (combine_external)"),
             headline=dict(P36a=results.get(zw.PRIMARY, {}).get("P36a", {}).get("grade_with_hw_flag"),
                           P36b=results.get(zw.PRIMARY, {}).get("P36b", {}).get("grade")),
             n_core_pairs=len(pair_rows), n_source_stations=len(scores),
             minutes=round((time.time() - t0) / 60.0, 1),
             finished=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()))
    with open(os.path.join(args.out_dir, "zw36_summary.json"), "w", encoding="utf-8") as f:
        json.dump(jsonable(s), f, indent=1)
    print(json.dumps(jsonable(s["headline"])), flush=True)
    return 0


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--data-dir", default="/data")
    ap.add_argument("--out-dir", default="/out")
    ap.add_argument("--time-budget-min", type=float, default=200.0)
    return run(ap.parse_args(argv))


if __name__ == "__main__":
    sys.exit(main())

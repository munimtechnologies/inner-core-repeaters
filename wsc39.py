#!/usr/bin/env python3
# Copyright 2026 MV2, Munim, Inc.
# Licensed under the Apache License, Version 2.0 (Apache-2.0). See LICENSE.
# Made with AI agents working in MV2's name.
"""wsc39.py - pre-registration 39: are the members of Wang et al.'s 16 South Sandwich multiplets one source?

Pre-registered in findings/2026-11-27-prereg39-wang-source-census.md, written and published before
any mantle-P record of a non-M03-24 Wang multiplet member was fetched (the three multiplet-F events
2005/2012/2022 were measured in zw35g2-v1 and are graded again here as a replication check).

Wang, Vidale et al. 2024 (Nature 631:340, Supplementary Data Table 3) list 16 multiplets of
repeating South Sandwich Islands earthquakes; both camps of the inner-core debate read waveform
changes between their members. Gate G2 of pre-registration 35 found that the 2012 member of
multiplet F is not a repeat in mantle P. This census applies the same rule to every member pair:

  * Measurement: mantlep2.py v0.13.0 unchanged (pinned commit), at every IU/II/GT BHZ station
    30-85 deg from 59.03 S 25.68 W (up to 60, nearest first), as in zw35g2. Per station:
    one REFERENCE pass on the eight M03-24 events of pre-registration 11, then one pass per
    multiplet on that multiplet's ComCat-located members only (so each multiplet is aligned on
    its own members, not on M03-24).
  * Station score of a member pair (a, c): q = share of the station's gated reference pairs with
    D_late < D_late(a, c). Counted if a and c are gated in the multiplet pass and the reference
    pass has >= 10 gated pairs. Hit if q >= 0.93.
  * Pair decision (pre-registration 11 rule): n counted stations, h hits; n < 3 untestable;
    h >= ceil(n/2) S; h <= max(1, floor(0.2 n)) N; otherwise inconclusive.
  * Pairs whose two members are both M03-24 reference events are not graded (they are the yardstick).
  * Multiplet verdict: odd (any pair S) / untestable (no testable pair) / clean (N in >= 2/3 of
    testable pairs) / unresolved (otherwise). Odd member: in >= 1 S pair, S in >= half of its
    testable pairs, and some N pair exists that does not involve it.
  * Census headline over the 14 graded multiplets (all but I and M): fewer than 5 testable ->
    untestable; odd <= 1 -> source-clean; odd >= 3 and odd >= 25% of testable -> source problems
    common; otherwise some source problems.
  * Positive control: multiplet I, which Wang et al. label "not good repeaters": passed if odd,
    failed if clean, else inconclusive. A failed control is attached to the headline.
  * Replication: Wang's multiplet F contains Zhang & Wen's multiplet F; its pairs
    (2005, 2012) and (2005, 2022) are compared with zw35g2-v1's S and N.

Top-level imports are standard library only; numpy/obspy/requests are imported inside run().
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

REPO = "https://raw.githubusercontent.com/munimtechnologies/inner-core-repeaters/"
MANTLEP2_COMMIT = "b4712f0ebbc70b896231620c905ff53a4b8270fe"   # mantlep2.py 0.13.0 (as zw35g2.py)
PREREG_URL = REPO + "main/findings/2026-11-27-prereg39-wang-source-census.md"

OLD8 = ("usp000bt5x", "usp000fz4d", "usp000hute", "usc000t4wj", "us7000am0e", "us6000h2wz",
        "us6000mk65", "us6000runq")

# Wang et al. 2024 Supplementary Data Table 1 event number -> USGS ComCat match
# (id, origin UTC, lat, lon, depth km, mag), from icr-wang-v1 (matched by origin time; all dt <= 0.1 s).
WANG_EVENTS = {
    701: ("usp0004qcq", "1991-04-15T20:39:47.810", -56.192, -26.54, 33.0, 5.2),
    702: ("usp0005hc3", "1992-11-25T12:55:28.250", -56.177, -26.778, 33.0, 4.9),
    704: ("usp00065uk", "1993-12-30T17:45:00.390", -58.956, -25.356, 10.0, 5.0),
    705: ("usp0007zxb", "1997-04-04T02:35:44.800", -57.893, -25.599, 33.0, 4.8),
    706: ("usp0008ct4", "1997-12-17T15:09:36.940", -58.962, -25.677, 50.0, 4.6),
    707: ("usp0008ds5", "1998-01-05T04:01:59.890", -56.315, -25.845, 37.7, 5.4),
    708: ("usp0008k0w", "1998-04-12T21:33:47.420", -56.116, -26.768, 100.0, 4.7),
    709: ("usp0008t10", "1998-08-09T20:03:30.790", -58.767, -25.349, 33.0, 4.9),
    711: None,   # 1999-05-14T05:05:08.20: no ComCat event within 90 s
    712: ("usp0009aw9", "1999-07-01T18:56:04.000", -59.05, -25.718, 33.0, 4.8),
    713: ("usp000a8mp", "2001-01-29T09:14:29.820", -58.909, -25.473, 33.0, 4.8),
    715: ("usp000b3ys", "2002-05-08T14:59:53.580", -56.129, -26.625, 33.0, 4.8),
    716: ("usp000b8b4", "2002-07-18T14:26:51.630", -58.688, -25.275, 33.0, 5.0),
    718: ("usp000bt5x", "2003-03-21T10:47:47.620", -59.071, -25.47, 33.0, 4.8),
    719: ("usp000by62", "2003-05-26T08:45:21.770", -57.535, -25.989, 73.3, 4.7),
    720: ("usp000c3yx", "2003-08-03T05:38:57.550", -58.915, -25.561, 33.0, 5.0),
    723: ("usp000cqp2", "2004-03-23T06:20:00.030", -56.19, -26.998, 77.5, 4.8),
    724: ("usp000cucb", "2004-05-07T09:50:31.630", -57.907, -25.535, 71.5, 5.5),
    726: ("usp000djau", "2005-03-16T12:27:06.940", -56.248, -26.778, 70.5, 4.9),
    727: ("usp000djxx", "2005-03-24T06:08:38.490", -55.985, -27.078, 69.9, 5.0),
    729: ("usp000dr5u", "2005-05-19T15:16:40.890", -57.921, -25.549, 52.3, 4.7),
    730: ("usp000dsnt", "2005-06-10T13:08:12.920", -58.706, -25.207, 9.8, 4.7),
    733: ("usp000f787", "2007-03-23T17:12:51.620", -58.653, -25.181, 20.5, 5.0),
    736: ("usp000gwse", "2009-05-02T22:27:07.920", -58.966, -25.591, 48.6, 4.7),
    737: ("usp000gwuv", "2009-05-05T06:16:07.350", -58.704, -25.115, 36.2, 5.1),
    739: ("usp000h4ph", "2009-12-02T23:14:00.830", -56.156, -26.772, 59.7, 5.4),
    742: ("usp000hute", "2011-02-20T17:37:37.320", -58.979, -25.857, 10.0, 5.0),
    744: ("usp000jwcf", "2012-12-03T03:11:10.790", -58.864, -25.565, 41.2, 4.9),
    745: ("usc000hjka", "2013-06-09T00:21:34.640", -58.724, -25.362, 33.0, 5.2),
    752: ("us10004xsg", "2016-03-15T13:35:17.740", -56.0147, -27.0579, 41.0, 5.2),
    753: ("us20006ugf", "2016-08-27T04:51:05.460", -56.2915, -26.8872, 58.55, 5.4),
    754: ("us10007tyy", "2017-01-18T14:52:56.800", -57.5213, -25.9955, 56.98, 5.0),
    756: ("us10008rd4", "2017-05-11T01:16:59.340", -56.3574, -25.7766, 10.0, 5.4),
    757: ("us20009npn", "2017-06-20T12:54:35.910", -56.2587, -27.0145, 56.19, 5.3),
    758: ("us2000c537", "2017-12-15T21:53:19.860", -56.3586, -26.7386, 51.06, 5.1),
    801: ("usp00071k7", "1995-08-08T04:37:55.830", -55.897, -26.995, 59.9, 4.8),
    802: ("usp00075p4", "1995-10-24T09:25:33.080", -56.297, -26.558, 33.0, 4.7),
    807: ("usp000bmza", "2003-01-03T06:04:51.640", -57.432, -25.828, 33.0, 4.5),
    808: ("usp000bz6u", "2003-06-06T09:01:08.080", -57.154, -25.725, 33.0, 4.9),
    809: ("usp000d3z7", "2004-09-09T12:50:09.580", -56.317, -25.734, 10.0, 4.9),
    810: ("usp000d9p9", "2004-11-30T20:16:12.300", -57.392, -25.764, 83.7, 4.7),
    811: ("usp000dkyk", "2005-04-01T22:43:22.770", -57.363, -25.801, 59.3, 4.5),
    812: ("usp000dxt2", "2005-08-21T18:06:56.620", -55.977, -27.932, 144.3, 4.9),
    822: ("usp000hhbu", "2010-08-02T04:18:04.260", -57.416, -25.772, 52.5, 4.5),
    823: ("usp000ht84", "2011-01-16T21:22:20.450", -55.925, -27.755, 108.8, 4.6),
    824: ("usp000jc8h", "2011-12-16T03:49:46.250", -57.391, -26.043, 69.3, 4.5),
    825: ("usp000jg1n", "2012-03-10T14:08:27.630", -57.479, -25.799, 68.5, 4.6),
    827: ("usp000jm37", "2012-05-29T15:14:58.940", -58.881, -25.563, 83.7, 4.9),
    830: ("usb000kh8x", "2013-10-19T13:43:46.310", -56.2908, -26.6858, 62.8, 4.9),
    833: ("us10003609", "2015-08-29T14:07:37.400", -57.4271, -25.6848, 56.34, 4.6),
    837: ("us20006uah", "2016-08-26T14:46:02.890", -57.4592, -25.9054, 75.87, 4.8),
    842: ("us1000efdr", "2018-05-27T15:05:52.060", -55.9513, -27.782, 111.82, 4.8),
    845: ("us70007pxg", "2020-02-14T09:55:28.590", -58.905, -25.5422, 35.0, 5.1),
    848: ("us7000am0e", "2020-07-11T18:58:06.501", -59.0899, -25.7164, 57.13, 5.3),
    849: ("us6000c6sf", "2020-10-08T15:48:26.691", -58.9366, -25.5815, 77.2, 5.2),
    851: ("us7000f8lp", "2021-08-27T12:27:49.950", -57.3942, -25.8894, 46.92, 4.7),
    852: ("us7000g8fw", "2022-01-02T08:30:46.266", -58.7064, -25.2481, 35.0, 4.8),
    853: ("us7000ga77", "2022-01-08T15:56:06.741", -58.8471, -25.1434, 10.0, 4.9),
    854: ("us6000h2wz", "2022-03-07T23:25:40.208", -59.0811, -25.7858, 35.0, 5.1),
    856: ("us6000hm6g", "2022-05-19T02:24:53.021", -58.796, -25.3186, 35.0, 5.1),
    857: ("us7000hem5", "2022-05-30T06:40:20.729", -57.6282, -25.8732, 53.29, 4.6),
    859: ("us6000ittv", "2022-10-15T00:12:35.059", -59.1057, -25.595, 38.897, 5.2),
    860: ("us6000j8es", "2022-12-10T12:01:47.583", -56.3417, -26.8518, 78.794, 5.0),
    861: ("us6000j9pd", "2022-12-16T10:47:50.125", -58.7849, -25.2723, 35.0, 5.0),
    862: ("us7000j2tv", "2023-01-07T16:18:55.214", -56.2146, -27.033, 62.042, 4.9),
}

# Wang et al. 2024 Supplementary Data Table 3 (event numbers as printed)
MULTIPLETS = {
    "A": (712, 718, 724, 742, 848, 854, 859),
    "B": (702, 708, 723, 739, 757, 862),
    "C": (704, 706, 713, 827, 849),
    "D": (807, 810, 825, 833),
    "E": (719, 754, 811, 822, 857),
    "F": (709, 730, 744, 856),
    "G": (802, 715, 830, 758, 860),
    "H": (808, 824, 837, 851),
    "I": (733, 745, 853, 861),      # "not good repeaters" (Wang et al.): positive control
    "J": (716, 737, 852),
    "K": (705, 711, 729),
    "L": (701, 726, 753),
    "M": (812, 823, 842),           # "too noisy" (Wang et al.): report only
    "N": (809, 707, 756),
    "O": (720, 736, 845),
    "P": (727, 752, 801),
}
POSITIVE_CONTROL = "I"
REPORT_ONLY = ("M",)
GRADED = tuple(k for k in sorted(MULTIPLETS) if k != POSITIVE_CONTROL and k not in REPORT_ONLY)
# Multiplets whose catalogue locations lie in the M03-24 cluster (C02 in icr-wang-v1): the yardstick
# is closest to them. Results are also reported split near / far (report only).
NEAR = ("A", "C", "F", "I", "J", "O")
# zw35g2-v1 (graded 2026-11-27): Wang 730 = 2005-06-10, 744 = 2012-12-03, 856 = 2022-05-19
REPLICATION = {(730, 744): "S", (730, 856): "N", (744, 856): "S"}

HIT_Q = 0.93
MIN_REF_PAIRS = 10
MIN_STATIONS = 3
NAN = float("nan")


# ----------------------------------------------------------------------------- pure helpers
def isf(x):
    try:
        return math.isfinite(float(x))
    except (TypeError, ValueError):
        return False


def key(a, b):
    return (a, b) if a <= b else (b, a)


def cid(w):
    e = WANG_EVENTS.get(w)
    return e[0] if e else None


def located(letter):
    """ComCat-located members of a multiplet, in origin-time order (Wang numbers)."""
    ms = [w for w in MULTIPLETS[letter] if WANG_EVENTS.get(w)]
    return sorted(ms, key=lambda w: WANG_EVENTS[w][1])


def member_pairs(letter):
    """Graded member pairs (Wang numbers, earlier first). Pairs of two M03-24 reference events are
    the yardstick and are not graded."""
    ms = located(letter)
    out = []
    for i in range(len(ms)):
        for j in range(i + 1, len(ms)):
            if cid(ms[i]) in OLD8 and cid(ms[j]) in OLD8:
                continue
            out.append((ms[i], ms[j]))
    return out


def gc_km(lat1, lon1, lat2, lon2):
    p1, p2 = math.radians(lat1), math.radians(lat2)
    dl = math.radians(lon2 - lon1)
    c = math.sin(p1) * math.sin(p2) + math.cos(p1) * math.cos(p2) * math.cos(dl)
    return 6371.0 * math.acos(max(-1.0, min(1.0, c)))


def catalogue_sep(wa, wb):
    a, b = WANG_EVENTS[wa], WANG_EVENTS[wb]
    return round(gc_km(a[2], a[3], b[2], b[3]), 1), round(abs(a[4] - b[4]), 1)


def reference_values(D_ref, gated_ref, ref_ids=OLD8):
    ref = [e for e in ref_ids if gated_ref.get(e)]
    vals = []
    for i in range(len(ref)):
        for j in range(i + 1, len(ref)):
            v = D_ref.get(key(ref[i], ref[j]))
            if isf(v):
                vals.append(float(v))
    return vals


def station_pair_score(ref_vals, D, gated, a, c, hit_q=HIT_Q, min_ref=MIN_REF_PAIRS):
    """ref_vals: reference D_late list at the station; D: {(id1, id2) sorted: D_late} from the
    multiplet pass; gated: {id: bool} from the multiplet pass; a, c: ComCat ids."""
    out = dict(counted=False, q=None, hit=None, d_pair=None, n_ref=len(ref_vals), reason="")
    if not (gated.get(a) and gated.get(c)):
        out["reason"] = "member_not_gated"
        return out
    d = D.get(key(a, c))
    if not isf(d):
        out["reason"] = "no_pair_D"
        return out
    out["d_pair"] = float(d)
    if len(ref_vals) < min_ref:
        out["reason"] = "too_few_reference_pairs"
        return out
    q = sum(1 for v in ref_vals if v < float(d)) / float(len(ref_vals))
    out.update(counted=True, q=q, hit=bool(q >= hit_q - 1e-12))
    return out


def decide(hits, n_min=MIN_STATIONS):
    n = len(hits)
    h = sum(1 for x in hits if x)
    if n < n_min:
        return "untestable", h
    if h >= int(math.ceil(n / 2.0)):
        return "S", h
    if h <= max(1, int(math.floor(0.2 * n))):
        return "N", h
    return "inconclusive", h


def multiplet_verdict(decisions):
    """decisions: {pair: S|N|inconclusive|untestable}."""
    vals = list(decisions.values())
    testable = [d for d in vals if d != "untestable"]
    if any(d == "S" for d in vals):
        return "odd"
    if not testable:
        return "untestable"
    n_n = sum(1 for d in testable if d == "N")
    if 3 * n_n >= 2 * len(testable):
        return "clean"
    return "unresolved"


def odd_members(decisions):
    """Members in >= 1 S pair with S in >= half of their testable pairs, provided an N pair exists
    that does not involve them."""
    members = sorted({m for p in decisions for m in p})
    out = []
    for m in members:
        mine = [d for p, d in decisions.items() if m in p and d != "untestable"]
        s = sum(1 for d in mine if d == "S")
        other_n = any(d == "N" for p, d in decisions.items() if m not in p)
        if s >= 1 and 2 * s >= len(mine) and other_n:
            out.append(m)
    return out


def census(verdicts, graded=GRADED):
    v = [verdicts.get(k, "untestable") for k in graded]
    counts = {x: sum(1 for y in v if y == x) for x in ("odd", "clean", "unresolved", "untestable")}
    testable = counts["odd"] + counts["clean"] + counts["unresolved"]
    if testable < 5:
        head = "untestable"
    elif counts["odd"] <= 1:
        head = "source-clean"
    elif counts["odd"] >= 3 and 4 * counts["odd"] >= testable:
        head = "source problems common"
    else:
        head = "some source problems"
    return head, counts, testable


def control_status(verdict):
    return {"odd": "passed", "clean": "failed"}.get(verdict, "inconclusive")


def final_headline(head, control):
    return head + (" [positive control failed]" if control == "failed" else "")


def grade_multiplet(letter, scores):
    """scores: {station: {(wa, wb): score dict}} -> dict with pair decisions and the verdict."""
    pairs = {}
    for p in member_pairs(letter):
        rows = [(st, sc[p]) for st, sc in sorted(scores.items()) if sc.get(p, {}).get("counted")]
        dec, h = decide([r["hit"] for _, r in rows])
        km, ddep = catalogue_sep(*p)
        pairs[p] = dict(decision=dec, n_stations=len(rows), hits=h, stations=[st for st, _ in rows],
                        q=[round(r["q"], 4) for _, r in rows], cat_sep_km=km, cat_ddepth_km=ddep)
    dec = {p: v["decision"] for p, v in pairs.items()}
    return dict(pairs=pairs, verdict=multiplet_verdict(dec), odd_members=odd_members(dec))


def replication(f_result):
    out = {}
    for p, want in REPLICATION.items():
        got = f_result["pairs"].get(p, {}).get("decision")
        out["%d-%d" % p] = dict(zw35g2=want, here=got, agree=(got == want))
    return out


def summarise(graded_by_letter):
    verdicts = {k: v["verdict"] for k, v in graded_by_letter.items()}
    head, counts, testable = census(verdicts)
    ctrl = control_status(verdicts.get(POSITIVE_CONTROL, "untestable"))
    near = census(verdicts, tuple(k for k in GRADED if k in NEAR))
    far = census(verdicts, tuple(k for k in GRADED if k not in NEAR))
    return dict(headline=final_headline(head, ctrl), headline_raw=head, counts=counts,
                n_testable=testable, positive_control=dict(multiplet=POSITIVE_CONTROL,
                verdict=verdicts.get(POSITIVE_CONTROL), status=ctrl),
                near_report_only=dict(headline=near[0], counts=near[1]),
                far_report_only=dict(headline=far[0], counts=far[1]), verdicts=verdicts)


# ----------------------------------------------------------------------------- I/O
def fetch(url, tries=4, wait=10):
    for k in range(tries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "inner-core-repeaters (MV2, Munim, Inc.)"})
            with urllib.request.urlopen(req, timeout=60) as r:
                if r.status == 200:
                    return r.read().decode("utf-8", "replace")
        except Exception as ex:  # noqa: BLE001
            print("fetch error", url, repr(ex)[:200], flush=True)
        time.sleep(wait * (k + 1))
    return None


def load_mantlep2(dest):
    os.makedirs(dest, exist_ok=True)
    p = os.path.join(dest, "mantlep2.py")
    if not os.path.exists(p):
        body = fetch(REPO + MANTLEP2_COMMIT + "/mantlep2.py")
        if not body or 'VERSION = "0.13.0"' not in body:
            return None
        compile(body, "mantlep2.py", "exec")
        with open(p, "w") as f:
            f.write(body)
    if dest not in sys.path:
        sys.path.insert(0, dest)
    import mantlep2  # noqa: E402
    return mantlep2


def write_csv(path, rows):
    keys = []
    for r in rows:
        for k in r:
            if k not in keys:
                keys.append(k)
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=keys or ["empty"])
        w.writeheader()
        w.writerows(rows)


def set_events(mp, events):
    mp.EVENTS = list(events)
    mp.IDS = [e[0] for e in events]
    mp.DATES = {e[0]: e[1][:10] for e in events}


def run_pass(mp, session, model, events, net, sta, lat, lon, cache, plotdir, tag, ev_all, pair_all):
    set_events(mp, events)
    os.makedirs(plotdir, exist_ok=True)
    ev_rows, pair_rows = [], []
    try:
        res = mp.run_station(session, model, net, sta, lat, lon, cache, plotdir, ev_rows, pair_rows, {})
    except Exception as ex:  # noqa: BLE001
        traceback.print_exc()
        res = dict(status="exception", error=repr(ex)[:300])
    gated = {r["event"]: str(r.get("gated")) == "True" for r in ev_rows}
    D = {key(r["ev1"], r["ev2"]): float(r["D_late"]) for r in pair_rows if isf(r.get("D_late"))}
    for r in ev_rows:
        r["pass"] = tag
    for r in pair_rows:
        r["pass"] = tag
    ev_all.extend(ev_rows)
    pair_all.extend(pair_rows)
    return res, gated, D


def write_outputs(out_dir, s, scores, ev_all, pair_all, mp):
    graded = {k: grade_multiplet(k, scores.get(k, {})) for k in sorted(MULTIPLETS)}
    s["multiplets"] = {k: dict(verdict=v["verdict"], odd_members=v["odd_members"],
                               odd_member_ids=[cid(w) for w in v["odd_members"]],
                               pairs={"%d-%d" % p: d for p, d in v["pairs"].items()})
                       for k, v in graded.items()}
    s["census"] = summarise(graded)
    s["replication_F"] = replication(graded["F"])
    rows = []
    for k, v in graded.items():
        for p, d in v["pairs"].items():
            rows.append(dict(multiplet=k, wang_1=p[0], wang_2=p[1], id_1=cid(p[0]), id_2=cid(p[1]),
                             date_1=WANG_EVENTS[p[0]][1][:10], date_2=WANG_EVENTS[p[1]][1][:10],
                             decision=d["decision"], n_stations=d["n_stations"], hits=d["hits"],
                             cat_sep_km=d["cat_sep_km"], cat_ddepth_km=d["cat_ddepth_km"],
                             graded=(k in GRADED), stations=" ".join(d["stations"]),
                             q=" ".join(str(x) for x in d["q"])))
    write_csv(os.path.join(out_dir, "wsc39_pairs.csv"), rows)
    srows = []
    for k, by_st in sorted(scores.items()):
        for st, sc in sorted(by_st.items()):
            for p, v in sc.items():
                r = dict(multiplet=k, station=st, wang_1=p[0], wang_2=p[1])
                r.update({x: v.get(x) for x in ("counted", "q", "hit", "d_pair", "n_ref", "reason")})
                srows.append(r)
    write_csv(os.path.join(out_dir, "wsc39_station_scores.csv"), srows)
    write_csv(os.path.join(out_dir, "wsc39_events.csv"), ev_all)
    write_csv(os.path.join(out_dir, "wsc39_mantlep_pairs.csv"), pair_all)
    with open(os.path.join(out_dir, "wsc39_summary.json"), "w") as f:
        json.dump(mp.jsonable(s) if mp else s, f, indent=1)


def run(args):
    t0 = time.time()
    os.makedirs(args.out_dir, exist_ok=True)
    s = dict(version=VERSION, started=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
             mantlep2_commit=MANTLEP2_COMMIT, hit_q=HIT_Q, min_ref_pairs=MIN_REF_PAIRS,
             min_stations=MIN_STATIONS, graded=GRADED, positive_control=POSITIVE_CONTROL,
             report_only=REPORT_ONLY, near=NEAR, finished=None)
    s["prereg_public_at_start"] = bool(fetch(PREREG_URL, tries=2, wait=5))
    mp = load_mantlep2("/tmp/wsc39_pinned")
    if mp is None or [e[0] for e in mp.EVENTS] != list(OLD8):
        s["error"] = "mantlep2.py v0.13.0 download failed or its EVENTS differ from M03-24"
        with open(os.path.join(args.out_dir, "wsc39_summary.json"), "w") as f:
            json.dump(s, f, indent=1)
        return 2
    ref_events = list(mp.EVENTS)

    import requests
    from obspy.taup import TauPyModel
    session = requests.Session()
    session.headers["User-Agent"] = "inner-core-repeaters (MV2, Munim, Inc.)"
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

    cache = os.path.join(args.data_dir, "mantlep")
    plots = os.path.join(args.out_dir, "plots")
    os.makedirs(cache, exist_ok=True)
    model = TauPyModel("ak135")
    scores = {k: {} for k in MULTIPLETS}
    ev_all, pair_all, stations, skipped = [], [], [], []
    for d, net, sta, lat, lon in cand:
        name = "%s.%s" % (net, sta)
        if (time.time() - t0) / 60.0 > args.time_budget_min:
            skipped.append(name)
            continue
        res, g_ref, D_ref = run_pass(mp, session, model, ref_events, net, sta, lat, lon, cache,
                                     os.path.join(plots, "ref"), "ref", ev_all, pair_all)
        ref_vals = reference_values(D_ref, g_ref)
        st = dict(station=name, dist_deg=round(d, 2), ref_status=res.get("status"),
                  n_ref_pairs=len(ref_vals), multiplets={})
        for k in sorted(MULTIPLETS):
            evs = [WANG_EVENTS[w] for w in located(k)]
            r2, g, D = run_pass(mp, session, model, evs, net, sta, lat, lon, cache,
                                os.path.join(plots, k), k, ev_all, pair_all)
            sc = {}
            for p in member_pairs(k):
                sc[p] = station_pair_score(ref_vals, D, g, cid(p[0]), cid(p[1]))
            scores[k][name] = sc
            st["multiplets"][k] = dict(status=r2.get("status"), n_gated=sum(1 for e in evs if g.get(e[0])),
                                       n_counted=sum(1 for v in sc.values() if v["counted"]))
        stations.append(st)
        print(name, "ref pairs", len(ref_vals), json.dumps(st["multiplets"])[:600], flush=True)
        s.update(stations=stations, skipped_time_budget=skipped,
                 minutes=round((time.time() - t0) / 60.0, 1))
        write_outputs(args.out_dir, s, scores, ev_all, pair_all, mp)
    set_events(mp, ref_events)
    s.update(stations=stations, skipped_time_budget=skipped,
             finished=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
             minutes=round((time.time() - t0) / 60.0, 1))
    write_outputs(args.out_dir, s, scores, ev_all, pair_all, mp)
    print(json.dumps(mp.jsonable(dict(census=s["census"], replication_F=s["replication_F"])))[:3000], flush=True)
    return 0


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--data-dir", default="/data")
    ap.add_argument("--out-dir", default="/out")
    ap.add_argument("--time-budget-min", type=float, default=170.0)
    return run(ap.parse_args(argv))


if __name__ == "__main__":
    sys.exit(main())

#!/usr/bin/env python3
# Copyright 2026 MV2, Munim, Inc.
# Licensed under the Apache License, Version 2.0 (Apache-2.0). See LICENSE.
# Made with AI agents working in MV2's name; no human expert has checked it.
"""m2measure.py: measurement front end of pre-registration 38, arm M (one clock, two arrays).

Part of inner-core-repeaters by MV2, Munim, Inc. (mv2@munimtech.com). Licence: Apache-2.0.
Decision rules: findings/2026-11-27-prereg38-armM-plan.md and m2core/m2core.py v0.1.0 (frozen
2026-11-27 before any M2 waveform was fetched). This file is frozen separately, also before any
M2 waveform is fetched: findings/2026-11-27-prereg38-armM-m2measure-frozen.md. It writes m2core's
input file (array, mult, ev1, ev2, t1, t2, D, ok) and changes nothing in the plan.

What it measures, per Wang, Vidale et al. 2024 (Nature 631:340) multiplet (Supplementary Data
Table 3, with multiplet A's misprinted 724 read as 734: pre-registration 39 amendment 1; I and M
excluded), one multiplet at a time, on its ComCat-matched members only:

  YKA   D = D_c310 (1 - cc in the 3-10 s late-PKIKP-coda window) and ok from yka.py v0.2.0 via
        m0324_check v0.8.0's run_yka, with yka_repick 0.9.2's constrained (+/-2 s) picker,
        exactly as m34.py did for M03-24. Only the event list changes.
  ILAR  each event is beamed with the ilar-v1 code (run.py v0.1.0 fetch_event / process_event,
        pinned commit 7067d4f); repeater edges among the multiplet's members come from run.py's
        detect_pairs with ilar-v1's thresholds (cc >= 0.90, max lag 8 s, snr_group >= 3.0) and no
        catalogue-distance limit (membership is Wang's, not ours); then D comes from
        m0324_check v0.8.0's ilar_frame + all_D, unchanged. ok = D finite, both events beamed, and
        both pass ilar-v1's SNR gates in the multiplet's frame (snr_group >= 3.0, reference-window
        SNR >= 5.0, DF-window SNR >= 3.0, computed as run.analyze_multiplet does).

Gate (before any grade is read): the same code run on the nine M03-24 (Z004) events must
reproduce our published rows, results/m0324chk-v2/yka_z004_pairs.csv (D_c310) and
ilar_z004_matrix.csv (D_reproduced): max |difference| <= 0.001 over the pairs common to both, and
at least 90% of the published finite rows present, at each array. Only if the gate passes is
m2core v0.1.0 run on the output (its results go to <out>/m2core/). If it fails, m2core is not
run and nothing is graded until a recount note explains why.

Top-level imports are standard library only (numpy, obspy and requests are used by the pinned
modules, imported inside run()).
"""
import argparse
import calendar
import csv
import io
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
YKA_PUB = REPO + "main/results/m0324chk-v2/yka_z004_pairs.csv"
ILAR_PUB = REPO + "main/results/m0324chk-v2/ilar_z004_matrix.csv"
M2CORE_COMMIT = "063ca49589f4f99e798bfac81eaf31b32036be23"   # m2core v0.1.0 as frozen 2026-11-27
M2CORE_URL = REPO + M2CORE_COMMIT + "/m2core/m2core.py"
YKA_REPICK_URL = REPO + "main/yka_repick.py"                   # checked: VERSION = "0.9.2"
PLAN_URL = REPO + "main/findings/2026-11-27-prereg38-armM-plan.md"
FROZEN_URL = REPO + "main/findings/2026-11-27-prereg38-armM-m2measure-frozen.md"

ILAR_PROC_TAG = "proc_m2_v1"
ILAR_MIN_SNR_GROUP = 3.0     # ilar-v1 (run.py defaults, results/ilar-v1/run_summary.json)
ILAR_MIN_SNR_REF = 5.0
ILAR_MIN_SNR_DF = 3.0
ILAR_DETECT_THR = 0.90
ILAR_MAX_LAG_S = 8.0
ILAR_MAX_KM = 20000.0        # no catalogue-distance limit inside a Wang multiplet
GATE_TOL = 0.001
GATE_MIN_FRAC = 0.9
EXCLUDED = ("I", "M")

# Wang et al. 2024 Supplementary Data Table 1 event number -> USGS ComCat match
# (id, origin UTC, lat, lon, depth km, mag); identical to wsc39.py 0.2.0 (icr-wang-v1 matches).
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
    726: ("usp000djau", "2005-03-16T12:27:06.940", -56.248, -26.778, 70.5, 4.9),
    727: ("usp000djxx", "2005-03-24T06:08:38.490", -55.985, -27.078, 69.9, 5.0),
    729: ("usp000dr5u", "2005-05-19T15:16:40.890", -57.921, -25.549, 52.3, 4.7),
    730: ("usp000dsnt", "2005-06-10T13:08:12.920", -58.706, -25.207, 9.8, 4.7),
    733: ("usp000f787", "2007-03-23T17:12:51.620", -58.653, -25.181, 20.5, 5.0),
    734: ("usp000fz4d", "2008-02-17T19:32:09.100", -58.986, -25.697, 66.5, 4.9),
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

# Wang et al. 2024 Supplementary Data Table 3 (A: 734 for the misprinted 724, amendment 1)
MULTIPLETS = {
    "A": (712, 718, 734, 742, 848, 854, 859),
    "B": (702, 708, 723, 739, 757, 862),
    "C": (704, 706, 713, 827, 849),
    "D": (807, 810, 825, 833),
    "E": (719, 754, 811, 822, 857),
    "F": (709, 730, 744, 856),
    "G": (802, 715, 830, 758, 860),
    "H": (808, 824, 837, 851),
    "I": (733, 745, 853, 861),
    "J": (716, 737, 852),
    "K": (705, 711, 729),
    "L": (701, 726, 753),
    "M": (812, 823, 842),
    "N": (809, 707, 756),
    "O": (720, 736, 845),
    "P": (727, 752, 801),
}
MEASURED = tuple(k for k in sorted(MULTIPLETS) if k not in EXCLUDED)


# ----------------------------------------------------------------------------- pure helpers
def isf(x):
    try:
        return math.isfinite(float(x))
    except (TypeError, ValueError):
        return False


def fnum(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return NAN


def truthy(x):
    return x is True or str(x).strip().lower() in ("true", "1", "yes")


def decimal_year(iso):
    """z004_events.csv convention (m34.py): year + (day_of_year - 1 + fraction of day) / days."""
    s = iso.strip().rstrip("Z")
    date, hms = s.split("T")
    y, mo, d = (int(v) for v in date.split("-"))
    h, mi, sec = hms.split(":")
    t0 = calendar.timegm((y, 1, 1, 0, 0, 0, 0, 0, 0))
    t = calendar.timegm((y, mo, d, int(h), int(mi), 0, 0, 0, 0)) + float(sec)
    days = 366 if calendar.isleap(y) else 365
    return y + (t - t0) / 86400.0 / days


def event_dict(w):
    """Wang number -> event dict in the format the pinned pipelines use (None if not matched)."""
    e = WANG_EVENTS.get(w)
    if not e:
        return None
    cid, t, lat, lon, dep, mag = e
    iso = t + "Z"
    return dict(id=cid, time=iso, lat=float(lat), lon=float(lon), depth=float(dep), mag=float(mag),
                t=decimal_year(iso), wang=w)


def multiplet_events(letter):
    """{ComCat id: event dict} for the matched members of one multiplet."""
    out = {}
    for w in MULTIPLETS[letter]:
        e = event_dict(w)
        if e is not None:
            out[e["id"]] = e
    return out


def member_pairs(evs):
    """All member pairs, earlier event first (ties by id)."""
    order = sorted(evs, key=lambda i: (evs[i]["t"], i))
    return [(order[i], order[j]) for i in range(len(order)) for j in range(i + 1, len(order))]


def ilar_ok(D, g1, g2, r1, r2, f1, f2):
    """ilar-v1 gates: snr_group >= 3.0, reference SNR >= 5.0, DF SNR >= 3.0 for both events."""
    vals = (D, g1, g2, r1, r2, f1, f2)
    if not all(isf(v) for v in vals):
        return False
    return (min(g1, g2) >= ILAR_MIN_SNR_GROUP and min(r1, r2) >= ILAR_MIN_SNR_REF
            and min(f1, f2) >= ILAR_MIN_SNR_DF)


def ilar_rows(letter, evs, D, snr):
    """m2core rows for one multiplet at ILAR. D: {(a, b): D} from all_D (a earlier);
    snr: {id: (snr_group, snr_ref, snr_df)}. Pairs without a value get D nan and ok False."""
    rows = []
    for a, b in member_pairs(evs):
        d = D.get((a, b), D.get((b, a), NAN))
        ga, ra, fa = snr.get(a, (NAN, NAN, NAN))
        gb, rb, fb = snr.get(b, (NAN, NAN, NAN))
        rows.append(dict(array="ILAR", mult=letter, ev1=a, ev2=b, t1=evs[a]["t"], t2=evs[b]["t"],
                         D=fnum(d), ok=ilar_ok(fnum(d), ga, gb, ra, rb, fa, fb),
                         snr_group1=ga, snr_group2=gb, snr_ref1=ra, snr_ref2=rb, snr_df1=fa, snr_df2=fb))
    return rows


def yka_rows(letter, evs, pipeline_rows):
    """m2core rows for one multiplet at YKA from run_yka's pairs csv rows (D = D_c310)."""
    by = {}
    for r in pipeline_rows:
        k = tuple(sorted((str(r.get("ev1")), str(r.get("ev2")))))
        by.setdefault(k, r)
    rows = []
    for a, b in member_pairs(evs):
        r = by.get(tuple(sorted((a, b))))
        d = fnum(r.get("D_c310")) if r else NAN
        ok = bool(r is not None and truthy(r.get("ok")) and isf(d))
        rows.append(dict(array="YKA", mult=letter, ev1=a, ev2=b, t1=evs[a]["t"], t2=evs[b]["t"],
                         D=d, ok=ok, D_direct=fnum(r.get("D")) if r else NAN,
                         D_c13=fnum(r.get("D_c13")) if r else NAN))
    return rows


def compare(ours, published, key_ours, key_pub):
    """ours/published: lists of dict rows with ev1, ev2. Max |diff| over pairs common to both
    (finite on both sides) and the number of published finite rows."""
    def k(r):
        return tuple(sorted((str(r.get("ev1")), str(r.get("ev2")))))
    pub = {k(r): fnum(r.get(key_pub)) for r in published if isf(r.get(key_pub))}
    diffs = []
    for r in ours:
        kk = k(r)
        if kk in pub and isf(r.get(key_ours)):
            diffs.append(abs(fnum(r.get(key_ours)) - pub[kk]))
    return dict(n_published=len(pub), n_common=len(diffs), max_abs_diff=max(diffs) if diffs else None)


def gate_ok(c, tol=GATE_TOL, min_frac=GATE_MIN_FRAC):
    if not c or not c.get("n_published") or c.get("max_abs_diff") is None:
        return False
    return c["max_abs_diff"] <= tol and c["n_common"] >= min_frac * c["n_published"]


def frame_snr(beam, off, w, n0, n1, n):
    """SNR of a beam in window w after shifting by off samples (run.analyze_multiplet's snr)."""
    s = int(round(off))
    a, c = w[0] + s, w[1] + s
    if a < 0 or c > n or c <= a:
        return NAN
    nb = math.sqrt(sum(float(x) ** 2 for x in beam[n0:n1]) / float(max(1, n1 - n0)))
    if not nb > 0:
        return NAN
    seg = beam[a:c]
    return math.sqrt(sum(float(x) ** 2 for x in seg) / float(len(seg))) / nb


def read_csv_text(text):
    return list(csv.DictReader(io.StringIO(text))) if text else []


def read_rows(path):
    if os.path.exists(path) and os.path.getsize(path) > 0:
        with open(path, newline="") as f:
            return list(csv.DictReader(f))
    return []


def jsonable(o):
    if isinstance(o, dict):
        return {str(k): jsonable(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [jsonable(v) for v in o]
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
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys or ["empty"])
        w.writeheader()
        w.writerows(rows)


def fetch(url, tries=5, wait=15, timeout=90):
    for k in range(tries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "inner-core-repeaters (MV2, Munim, Inc.)"})
            with urllib.request.urlopen(req, timeout=timeout) as r:
                if r.status == 200:
                    return r.read().decode("utf-8")
        except Exception as ex:  # noqa: BLE001
            print("fetch error", url, repr(ex)[:160], flush=True)
        time.sleep(wait * (k + 1))
    return None


def load_module(name, url, marker, here, dest="/tmp/m2_pinned"):
    """Import `name` from next to this file if present, else download it once (must contain
    `marker`, e.g. its VERSION line) into dest and import it from there."""
    for d in (here, dest):
        if os.path.exists(os.path.join(d, name + ".py")):
            if d not in sys.path:
                sys.path.insert(0, d)
            break
    else:
        body = fetch(url, tries=4, wait=10)
        if not body or marker not in body:
            raise RuntimeError("could not download %s with %s" % (name, marker))
        compile(body, name + ".py", "exec")
        os.makedirs(dest, exist_ok=True)
        with open(os.path.join(dest, name + ".py"), "w") as f:
            f.write(body)
        if dest not in sys.path:
            sys.path.insert(0, dest)
    return __import__(name)


# ----------------------------------------------------------------------------- ILAR
class Ilar:
    def __init__(self, mods, data_dir):
        import m0324_check as mc  # noqa: E402  (pinned v0.8.0)
        import run as ilrun  # noqa: E402  (pinned run.py v0.1.0, the ilar-v1 code)
        from phases import Phases, load_or_build  # noqa: E402
        self.mc, self.run = mc, ilrun
        self.data_dir = data_dir
        self.ph = Phases(load_or_build(os.path.join(data_dir, "ak135_pkp_grid_v1.npz")))
        self.proc = os.path.join(data_dir, ILAR_PROC_TAG)
        os.makedirs(self.proc, exist_ok=True)
        self.client = None
        self.coords = None
        self.data = {}
        self.status = {}

    def ensure(self, e, deadline):
        i = e["id"]
        if i in self.data or i in self.status:
            return
        pp = os.path.join(self.proc, i + ".npz")
        try:
            if not os.path.exists(pp):
                if time.time() > deadline:
                    self.status[i] = "skipped_time_budget"
                    return
                if self.client is None:
                    self.client = self.run.get_client()
                    self.coords = self.run.get_coords(self.client, self.data_dir)
                st, why = self.run.fetch_event(self.client, e, self.ph, self.data_dir)
                if st is None:
                    self.status[i] = why
                    return
                d, why2 = self.run.process_event(st, e, self.ph, self.coords)
                if d is None:
                    self.status[i] = why2
                    return
                self.run.save_proc(pp, d)
                time.sleep(0.2)
            self.data[i] = self.mc.load_ilar(pp)
            self.status[i] = "ok"
        except Exception as ex:  # noqa: BLE001
            self.status[i] = "exception:" + repr(ex)[:150]

    def measure(self, evs):
        """-> ({(a, b): D}, {id: (snr_group, snr_ref, snr_df)}, info) for one event set."""
        mc, run = self.mc, self.run
        ids = sorted([i for i in evs if i in self.data], key=lambda i: (evs[i]["t"], i))
        info = dict(n_events=len(evs), n_beamed=len(ids))
        snr = {i: (float(self.data[i]["snr_group"]), NAN, NAN) for i in ids}
        use = [i for i in ids if self.data[i]["snr_group"] >= ILAR_MIN_SNR_GROUP]
        info["n_snr_group_ok"] = len(use)
        if len(use) < 2:
            info["status"] = "fewer than 2 beamed events with snr_group >= 3"
            return {}, snr, info
        edges, ncomp, nuse = run.detect_pairs(use, evs, self.data, ILAR_MAX_KM, ILAR_DETECT_THR,
                                              ILAR_MAX_LAG_S, ILAR_MIN_SNR_GROUP)
        info.update(n_compared=ncomp, n_edges=len(edges),
                    edges=[dict(ev1=a, ev2=b, cc=round(c, 4)) for (a, b), (c, l) in sorted(edges.items())])
        frame = mc.ilar_frame(use, self.data, edges)
        if frame is None:
            info["status"] = "no DF pick on the multiplet stack"
            return {}, snr, info
        info.update(status="ok", root=frame["root"], frame_members=list(frame["members"]),
                    pick_minus_pred_s=frame["pick_minus_pred_s"], dfw=list(frame["dfw"]), refw=list(frame["refw"]))
        for m in frame["members"]:
            b = self.data[m]["beam"]
            snr[m] = (float(self.data[m]["snr_group"]),
                      frame_snr(b, frame["off"][m], frame["refw"], mc.N0, mc.N1, mc.N),
                      frame_snr(b, frame["off"][m], frame["dfw"], mc.N0, mc.N1, mc.N))
        beams = {m: self.data[m]["beam"] for m in frame["members"]}
        D = mc.all_D(beams, frame, evs)
        return {k: float(v) for k, v in D.items()}, snr, info


# ----------------------------------------------------------------------------- YKA
class Yka:
    def __init__(self, data_dir, out_dir, here):
        import m0324_check as mc  # noqa: E402  (pinned v0.8.0)
        import xcorr  # noqa: E402
        yka_repick = load_module("yka_repick", YKA_REPICK_URL, 'VERSION = "0.9.2"', here)
        self.mc, self.xcorr, self.rp = mc, xcorr, yka_repick
        self.data_dir, self.out_dir = data_dir, out_dir
        self.yka = mc.load_yka_module("/tmp/icr_mods")

    def measure(self, name, evs):
        if self.yka is None:
            return [], dict(status="yka.py download failed")
        sub = os.path.join(self.out_dir, "yka_runs", name)
        plots = os.path.join(sub, "plots")
        os.makedirs(plots, exist_ok=True)
        rec = []
        self.yka.pick_onset = self.rp.constrained_picker(self.xcorr.pick_onset, self.rp.REPICK_S, record=rec)
        sr = {}
        try:
            self.mc.run_yka(evs, list(evs), self.data_dir, sub, plots, sr)
        except Exception as ex:  # noqa: BLE001
            traceback.print_exc()
            sr["yka"] = dict(status="exception", error=repr(ex)[:300])
        rows = read_rows(os.path.join(sub, "yka_z004_pairs.csv"))
        y = sr.get("yka", {}) if isinstance(sr.get("yka"), dict) else {}
        mi = y.get("multiplet", {}) if isinstance(y.get("multiplet"), dict) else {}
        info = dict(status=y.get("status"), n_processed=y.get("n_processed"), n_pairs=len(rows),
                    n_ok=sum(1 for r in rows if truthy(r.get("ok"))),
                    pick_minus_pred_s=mi.get("pick_minus_pred_s"),
                    pick_valid=self.rp.pick_valid(mi.get("pick_minus_pred_s")),
                    event_status=read_rows(os.path.join(sub, "yka_z004_status.csv")))
        return rows, info


# ----------------------------------------------------------------------------- m2core
def run_core(here, pairs_path, out_dir):
    m2core = load_module("m2core", M2CORE_URL, 'VERSION = "0.1.0"', here)   # frozen v0.1.0
    if m2core.VERSION != "0.1.0":
        raise RuntimeError("m2core version %s, expected 0.1.0" % m2core.VERSION)
    os.makedirs(out_dir, exist_ok=True)
    pairs = m2core.read_pairs(pairs_path)
    res, units, null, drops, closure_rows = m2core.analyse(pairs)
    with open(os.path.join(out_dir, "m2_summary.json"), "w") as f:
        json.dump(res, f, indent=1)
    m2core.write_csv(os.path.join(out_dir, "m2_units.csv"), units)
    m2core.write_csv(os.path.join(out_dir, "m2_null.csv"), null)
    m2core.write_csv(os.path.join(out_dir, "m2_drop_one.csv"), drops)
    m2core.write_csv(os.path.join(out_dir, "m2_closure_pairs.csv"), closure_rows)
    return dict(version=m2core.VERSION, checks_ok=res["checks_ok"], M2a_grade=res["M2a"]["grade"],
                M2b_grade=res["M2b"]["grade"])


# ----------------------------------------------------------------------------- main
def run_all(args):
    t0 = time.time()
    deadline = t0 + 60.0 * args.time_budget_min
    here = os.path.dirname(os.path.abspath(__file__))
    mods = os.path.abspath(args.mods)
    for p in (here, mods):
        if p not in sys.path:
            sys.path.insert(0, p)
    os.makedirs(args.out_dir, exist_ok=True)
    s = dict(version=VERSION, started=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()), measured=MEASURED,
             multiplet_A=list(MULTIPLETS["A"]), ilar_proc_tag=ILAR_PROC_TAG, gate_tol=GATE_TOL,
             gate_min_frac=GATE_MIN_FRAC, finished=None)
    s["plan_public_at_start"] = bool(fetch(PLAN_URL, tries=2, wait=5))
    s["frozen_note_public_at_start"] = bool(fetch(FROZEN_URL, tries=2, wait=5))

    import m0324_check as mc  # noqa: E402
    gate_evs = mc.read_events(os.path.join(mods, "z004_events.csv"))
    ilar = Ilar(mods, args.data_dir)
    yka = Yka(args.data_dir, args.out_dir, here)

    # ---- gate: M03-24 (Z004), both arrays
    gate = {}
    for e in sorted(gate_evs.values(), key=lambda e: e["t"]):
        ilar.ensure(e, deadline)
    Dg, snrg, infog = ilar.measure(gate_evs)
    ours_i = [dict(ev1=a, ev2=b, D=v) for (a, b), v in Dg.items()]
    pub_i = read_csv_text(fetch(ILAR_PUB))
    gate["ILAR"] = compare(ours_i, pub_i, "D", "D_reproduced")
    gate["ILAR"]["info"] = {k: v for k, v in infog.items() if k != "edges"}
    yrows_g, yinfo_g = yka.measure("gate_Z004", gate_evs)
    pub_y = read_csv_text(fetch(YKA_PUB))
    gate["YKA"] = compare(yrows_g, pub_y, "D_c310", "D_c310")
    gate["YKA"]["info"] = {k: v for k, v in yinfo_g.items() if k != "event_status"}
    gate["passed"] = gate_ok(gate["ILAR"]) and gate_ok(gate["YKA"])
    s["gate"] = gate
    print("GATE", json.dumps(jsonable(gate))[:1500], flush=True)

    # ---- the Wang multiplets, one at a time
    rows, per = [], {}
    for letter in MEASURED:
        evs = multiplet_events(letter)
        rec = dict(n_matched=len(evs), members=sorted(evs, key=lambda i: evs[i]["t"]))
        if time.time() > deadline:
            rec["status"] = "skipped_time_budget"
            per[letter] = rec
            rows += ilar_rows(letter, evs, {}, {}) + yka_rows(letter, evs, [])
            continue
        for e in sorted(evs.values(), key=lambda e: e["t"]):
            ilar.ensure(e, deadline)
        try:
            D, snr, info = ilar.measure(evs)
        except Exception as ex:  # noqa: BLE001
            traceback.print_exc()
            D, snr, info = {}, {}, dict(status="exception:" + repr(ex)[:200])
        ir = ilar_rows(letter, evs, D, snr)
        yr_raw, yinfo = yka.measure(letter, evs)
        yr = yka_rows(letter, evs, yr_raw)
        rec.update(ilar=info, yka={k: v for k, v in yinfo.items() if k != "event_status"},
                   yka_event_status=yinfo.get("event_status"),
                   n_ok_ilar=sum(1 for r in ir if r["ok"]), n_ok_yka=sum(1 for r in yr if r["ok"]))
        per[letter] = rec
        rows += ir + yr
        print("MULT", letter, json.dumps(jsonable(dict(n_ok_ilar=rec["n_ok_ilar"], n_ok_yka=rec["n_ok_yka"],
                                                       ilar_status=info.get("status"))))[:400], flush=True)
        s["multiplets"] = per
        write_csv(os.path.join(args.out_dir, "m2_pairs.csv"), rows)
        with open(os.path.join(args.out_dir, "m2measure_summary.json"), "w") as f:
            json.dump(jsonable(s), f, indent=1)

    s["multiplets"] = per
    s["ilar_event_status"] = ilar.status
    s["n_rows"] = len(rows)
    s["n_ok"] = dict(ILAR=sum(1 for r in rows if r["array"] == "ILAR" and r["ok"]),
                     YKA=sum(1 for r in rows if r["array"] == "YKA" and r["ok"]))
    pairs_path = os.path.join(args.out_dir, "m2_pairs.csv")
    write_csv(pairs_path, rows)
    if gate["passed"]:
        try:
            s["m2core"] = run_core(here, pairs_path, os.path.join(args.out_dir, "m2core"))
        except Exception as ex:  # noqa: BLE001
            traceback.print_exc()
            s["m2core"] = dict(status="exception", error=repr(ex)[:300])
    else:
        s["m2core"] = dict(status="not run: the M03-24 reproduction gate failed (nothing is graded)")
    s["finished"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    s["minutes"] = round((time.time() - t0) / 60.0, 1)
    with open(os.path.join(args.out_dir, "m2measure_summary.json"), "w") as f:
        json.dump(jsonable(s), f, indent=1)
    print(json.dumps(jsonable(dict(gate_passed=gate["passed"], n_ok=s["n_ok"], m2core=s["m2core"])))[:2000],
          flush=True)
    return 0


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--mods", required=True, help="folder with the pinned v0.8.0 modules (fetch_pinned.py)")
    ap.add_argument("--data-dir", default="/data")
    ap.add_argument("--out-dir", default="/out")
    ap.add_argument("--time-budget-min", type=float, default=200.0)
    return run_all(ap.parse_args(argv))


if __name__ == "__main__":
    sys.exit(main())

#!/usr/bin/env python3
"""mantlep3.py: source-verified clock test across South Sandwich Islands families.

Part of inner-core-repeaters by MV2, Munim, Inc. (made with AI agents). Licence: Apache-2.0.
Pre-registered in findings/2026-11-18-h13-prereg.md (pre-registration 13), written and
published before this script fetched any new waveform.

Why: at Yellowknife (YKA) the inner-core coda difference of catalogue "repeaters" tracks how far
apart the sources are, not time (diag-v1), so the published retrospective tests are limited by
the source. For one family, M03-24, mantle P showed the sources are alike while the core-path
coda changed and changed back (mantlep2-v1, pre-registration 11, N). This script applies that
source control to every published pair and asks which clock the core change follows.

Step 1  For every ok pair of the published YKA pair table (results/ykab-v2/yka_pairs.csv) and
        ILAR pair table (results/ilar-v1/pairs.csv), measure mantle-P late-coda dissimilarity
        D_late (3-10 s, 1-2 Hz) at the eight stations usable in mantlep2-v1, one family at a
        time, with mantlep2's alignment and gates.
Step 2  Source score S of a pair = median over stations of the percentile of its D_late among
        all gated pairs at that station (stations with >= 20 gated pairs; pairs with >= 2).
Step 3  Calibration: does S predict the core D (Spearman >= 0.15, family-bootstrap 5% > 0)?
        Clock: among source-verified pairs (S <= 1/3), does core D follow the backtracking
        clock dB (Wang et al. 2024: turn 2008.5, back at 1/2.5 the rate), the monotone clock dM,
        or neither? Decision rule: decide13().

Writes to --out-dir: h13_summary.json, mantlep3_pairs.csv, mantlep3_families.csv,
mantlep3_events.csv, h13_scores.csv, h13_pairs_scored.csv, plots/.
"""
import argparse
import bisect
import csv
import datetime
import io
import itertools
import json
import math
import os
import random
import sys
import time
import traceback

import numpy as np

import mantlep2 as m2

VERSION = "0.15.0"
RAW = "https://raw.githubusercontent.com/munimtechnologies/inner-core-repeaters/"
TABLES = {   # array: (commit, path) of the published pair tables; only rows with ok == True are used
    "yka": ("f65c966adf2cce3a22ffd9e4629c842959121407", "results/ykab-v2/yka_pairs.csv"),
    "ilar": ("7067d4fc2b0a2fdb1378a1c406dd680c8d0b69b1", "results/ilar-v1/pairs.csv"),
}
PREREG = "findings/2026-11-18-h13-prereg.md"
USGS = "https://earthquake.usgs.gov/fdsnws/event/1/query"
USGS_PARAMS = dict(format="csv", starttime="1990-01-01", endtime="2026-11-19", minlatitude=-61.5,
                   maxlatitude=-54, minlongitude=-31, maxlongitude=-22, minmagnitude=4.3,
                   orderby="time-asc", limit=20000)

# The eight stations usable in mantlep2-v1, quietest first; coordinates are fallbacks only.
STATIONS = [
    ("IU", "QSPA", -89.9289, 144.4382), ("GT", "VNDA", -77.5172, 161.8528),
    ("IU", "LCO", -29.011, -70.7005), ("GT", "DBIC", 6.6702, -4.8566),
    ("GT", "PLCA", -40.7328, -70.5508), ("GT", "LPAZ", -16.2879, -68.1307),
    ("GT", "LBTB", -25.015, 25.597), ("IU", "LSZ", -15.2779, 28.1882),
]

# M03-24 (merged family Y007): already used in pre-registrations 7-11. Pairs with both events
# in this set are left out of every primary statistic (sensitivity only).
M0324 = frozenset(("usp000bt5x", "usp000fz4d", "usp000hute", "usc000t4wj", "us7000am0e",
                   "us6000h2wz", "us6000ittv", "us6000mk65", "us6000runq"))

TURN, RATIO = 2008.5, 2.5
EDGE_S = 0.25                 # reference onset this close to the search start = noise, discarded
MIN_STATION_PAIRS = 20
MIN_ST = 2
S_MAX = 1.0 / 3.0
CAL_RHO = 0.15
RHO_MIN = 0.10
RHO_NONE_P95 = 0.20
FRAC = 0.90
MIN_V, MIN_FAM, MIN_CROSS = 40, 4, 10
NBOOT = 1000
NAN = float("nan")


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


def equiv(t):
    """Backtracking-clock time: the inner core is back where it was at equiv(t)."""
    return t if t <= TURN else TURN - (t - TURN) / RATIO


def d_B(t1, t2):
    return abs(equiv(t2) - equiv(t1))


def d_M(t1, t2):
    return abs(t2 - t1)


def is_cross(t1, t2):
    return (t1 < TURN) != (t2 < TURN)


def median(xs):
    xs = sorted(float(x) for x in xs if isf(x))
    n = len(xs)
    if n == 0:
        return NAN
    m = n // 2
    return xs[m] if n % 2 else 0.5 * (xs[m - 1] + xs[m])


def rank_avg(v):
    """Ranks 1..n, ties averaged."""
    r = m2.rankdata(np.asarray(v, dtype=float))
    return [float(x) + 1.0 for x in r]


def spearman(a, b):
    return m2.spearman(a, b)


def sep_km(lat1, lon1, lat2, lon2):
    p1, p2 = math.radians(lat1), math.radians(lat2)
    dp, dl = p2 - p1, math.radians(lon2 - lon1)
    h = math.sin(dp / 2) ** 2 + math.cos(p1) * math.cos(p2) * math.sin(dl / 2) ** 2
    return 2 * 6371.0 * math.asin(min(1.0, math.sqrt(h)))


def iso_to_epoch(s):
    s = s.strip().rstrip("Z")
    fmt = "%Y-%m-%dT%H:%M:%S.%f" if "." in s else "%Y-%m-%dT%H:%M:%S"
    return datetime.datetime.strptime(s, fmt).replace(tzinfo=datetime.timezone.utc).timestamp()


def read_catalog(text):
    by_id, items = {}, []
    for r in csv.DictReader(io.StringIO(text)):
        try:
            e = dict(lat=float(r["latitude"]), lon=float(r["longitude"]), depth=fnum(r.get("depth")),
                     epoch=iso_to_epoch(r["time"]))
        except (KeyError, ValueError):
            continue
        by_id[r.get("id", "")] = e
        items.append(e)
    items.sort(key=lambda e: e["epoch"])
    return dict(by_id=by_id, items=items, keys=[e["epoch"] for e in items])


def locate(ev_id, iso, cat, tol=5.0):
    if ev_id in cat["by_id"]:
        return cat["by_id"][ev_id], "id"
    if not iso:
        return None, "none"
    try:
        t = iso_to_epoch(iso)
    except ValueError:
        return None, "none"
    i = bisect.bisect_left(cat["keys"], t)
    best = None
    for j in (i - 1, i):
        if 0 <= j < len(cat["items"]):
            d = abs(cat["keys"][j] - t)
            if d <= tol and (best is None or d < best[0]):
                best = (d, cat["items"][j])
    return (best[1], "time") if best else (None, "none")


def read_pairs(text, array):
    """Rows with ok == True of a published pair table, with both clocks recomputed."""
    rows = []
    for r in csv.DictReader(io.StringIO(text)):
        if str(r.get("ok", "")).strip() != "True":
            continue
        t1, t2 = fnum(r.get("t1")), fnum(r.get("t2"))
        a, b = (r.get("ev1") or "").strip(), (r.get("ev2") or "").strip()
        if not (isf(t1) and isf(t2) and a and b and a != b):
            continue
        rows.append(dict(array=array, mult=("I:" if array == "ilar" else "") + r.get("mult", ""),
                         ev1=a, ev2=b, key=tuple(sorted((a, b))),
                         time1=r.get("time1", ""), time2=r.get("time2", ""), t1=t1, t2=t2,
                         mag1=fnum(r.get("mag1")), mag2=fnum(r.get("mag2")),
                         D=fnum(r.get("D")), D_c13=fnum(r.get("D_c13")), D_c310=fnum(r.get("D_c310")),
                         dB=d_B(t1, t2), dM=d_M(t1, t2), cross=is_cross(t1, t2),
                         m0324=(a in M0324 and b in M0324)))
    return rows


def build_events(rows, cat):
    ev, miss = {}, []
    for r in rows:
        for e, iso in ((r["ev1"], r["time1"]), (r["ev2"], r["time2"])):
            if e in ev or e in miss:
                continue
            loc, how = locate(e, iso, cat)
            if loc is None or not iso:
                miss.append(e)
                continue
            dep = loc["depth"] if isf(loc["depth"]) else 33.0
            ev[e] = dict(id=e, time=iso, lat=loc["lat"], lon=loc["lon"], depth=dep, how=how)
    return ev, miss


def families(rows, events):
    fam = {}
    for r in rows:
        if r["ev1"] not in events or r["ev2"] not in events:
            continue
        f = fam.setdefault(r["mult"], dict(members=[], need=set()))
        for e in (r["ev1"], r["ev2"]):
            if e not in f["members"]:
                f["members"].append(e)
        f["need"].add(r["key"])
    return fam


# ----------------------------------------------------------------------------- per family
def measure_family(X, raw_snr, ids, need, fs=None, center=None):
    """mantlep2's alignment and gates for one family at one station.
    X: {event: trace on the common grid, ak135 P at center}. Returns (res, gates, pairs) with
    pairs {(ev1, ev2) sorted: pair_D dict} for the needed pairs whose two events pass the gates."""
    fs = m2.FS if fs is None else fs
    center = int(round(m2.PRE * fs)) if center is None else center
    ids = [e for e in ids if e in X]
    res = dict(n_events=len(ids))
    cands = [e for e in ids if np.isfinite(raw_snr.get(e, NAN))]
    if len(cands) < 2:
        res["status"] = "too_few"
        return res, {}, {}
    ref = max(cands, key=lambda e: raw_snr[e])
    anchor = m2.onset_in(X[ref], fs, m2.ONSET_SEARCH[0], m2.ONSET_SEARCH[1], center)
    if anchor is None:
        res["status"] = "no_ref_onset"
        return res, {}, {}
    res.update(ref=ref, ref_onset_minus_pred_s=round((anchor - center) / fs, 2))
    if anchor - (center + int(round(m2.ONSET_SEARCH[0] * fs))) <= int(round(EDGE_S * fs)):
        res["status"] = "edge"
        return res, {}, {}
    A = {}
    for e in ids:
        if e == ref:
            A[e] = X[e]
            continue
        lag, c = m2.align2(X[ref], X[e], fs, anchor)
        if lag is None:
            continue
        A[e] = m2.shift(X[e], lag)
    norm = []
    for a in A.values():
        seg = a[anchor + int(round(m2.ALIGN_WIN[0] * fs)):anchor + int(round(m2.ALIGN_WIN[1] * fs))]
        mx = float(np.max(np.abs(seg))) if len(seg) else 0.0
        norm.append(a / mx if mx > 0 else a)
    stack = np.mean(norm, axis=0)
    onset = m2.onset_in(stack, fs, -m2.STACK_HALF_S, m2.STACK_HALF_S, anchor)
    if onset is None:
        res["status"] = "no_onset"
        return res, {}, {}
    gates = {}
    for e in ids:
        if e not in A:
            gates[e] = False
            continue
        s = m2.window_ratio(A[e], fs, onset, m2.SIG, m2.NOISE, "peak")
        l = m2.window_ratio(A[e], fs, onset, m2.WINDOWS["D_late"], m2.NOISE, "rms")
        gates[e] = bool(np.isfinite(s) and s >= m2.SNR_DIRECT_MIN and np.isfinite(l) and l >= m2.SNR_LATE_MIN)
    pairs = {}
    for k in sorted(need):
        a, b = k
        if not (gates.get(a) and gates.get(b)):
            continue
        r = m2.pair_D(A[a], A[b], onset, fs)
        if r is None or not np.isfinite(r.get("D_late", NAN)):
            continue
        pairs[k] = r
    res.update(status="ok", onset_minus_pred_s=round((onset - center) / fs, 2),
               n_gated=int(sum(1 for v in gates.values() if v)), n_pairs=len(pairs))
    return res, gates, pairs


# ----------------------------------------------------------------------------- scores and test
def source_scores(pair_rows, min_station_pairs=MIN_STATION_PAIRS):
    """{(ev1, ev2): dict(S, n_st, D_late_median)}, list of stations used."""
    by, seen = {}, set()
    for r in pair_rows:
        key = tuple(sorted((r["ev1"], r["ev2"])))
        k = (r["station"],) + key
        if k in seen or not isf(r.get("D_late")):
            continue
        seen.add(k)
        by.setdefault(r["station"], []).append((key, float(r["D_late"])))
    pct, raw, used = {}, {}, []
    for st, lst in sorted(by.items()):
        if len(lst) < min_station_pairs:
            continue
        used.append(st)
        rk = rank_avg([d for _, d in lst])
        n = float(len(lst))
        for (key, d), q in zip(lst, rk):
            pct.setdefault(key, []).append(q / n)
            raw.setdefault(key, []).append(d)
    out = {k: dict(S=median(v), n_st=len(v), D_late_median=median(raw[k])) for k, v in pct.items()}
    return out, used


def decide13(rB, rM, frac_B, frac_M, p95_B, p95_M):
    """Pre-registered clock decision among source-verified pairs."""
    if not all(isf(x) for x in (rB, rM, frac_B, frac_M, p95_B, p95_M)):
        return "untestable"
    if rB >= RHO_MIN and frac_B >= FRAC:
        return "B"
    if rM >= RHO_MIN and frac_M >= FRAC:
        return "M"
    if max(rB, rM) < RHO_MIN and max(p95_B, p95_M) < RHO_NONE_P95:
        return "none"
    return "inconclusive"


def _pct(xs, q):
    xs = sorted(x for x in xs if isf(x))
    if not xs:
        return NAN
    return xs[int(round(q * (len(xs) - 1)))]


def family_boot(rows, stat, nboot=NBOOT, seed=13):
    by = {}
    for r in rows:
        by.setdefault(r["mult"], []).append(r)
    fams = sorted(by)
    rng = random.Random(seed)
    out = []
    for _ in range(nboot if len(fams) >= 2 else 0):
        sub = [r for f in (rng.choice(fams) for _ in fams) for r in by[f]]
        out.append(stat(sub))
    return out


def clock(V, key, nboot=NBOOT, seed=13):
    fams = sorted({r["mult"] for r in V})
    n_cross = sum(1 for r in V if r["cross"])
    out = dict(n_pairs=len(V), n_families=len(fams), n_cross=n_cross,
               testable=bool(len(V) >= MIN_V and len(fams) >= MIN_FAM and n_cross >= MIN_CROSS))
    rB = spearman([r["dB"] for r in V], [r[key] for r in V])
    rM = spearman([r["dM"] for r in V], [r[key] for r in V])

    def st(sub):
        return (spearman([r["dB"] for r in sub], [r[key] for r in sub]),
                spearman([r["dM"] for r in sub], [r[key] for r in sub]))

    boots = [b for b in family_boot(V, st, nboot, seed) if isf(b[0]) and isf(b[1])]
    nb = len(boots)
    frac_B = sum(1 for b, m in boots if b > m) / float(nb) if nb else NAN
    frac_M = sum(1 for b, m in boots if m > b) / float(nb) if nb else NAN
    p95_B = _pct([b for b, _ in boots], 0.95)
    p95_M = _pct([m for _, m in boots], 0.95)
    out.update(rho_B=rB, rho_M=rM, n_boot=nb, frac_B_gt_M=frac_B, frac_M_gt_B=frac_M,
               p05_B=_pct([b for b, _ in boots], 0.05), p95_B=p95_B,
               p05_M=_pct([m for _, m in boots], 0.05), p95_M=p95_M)
    out["decision"] = decide13(rB, rM, frac_B, frac_M, p95_B, p95_M) if out["testable"] else "untestable"
    return out


def clock_test(rows, scores, key="D", s_max=S_MAX, include_m0324=False, nboot=NBOOT, seed=13):
    base = [r for r in rows if isf(r.get(key)) and (include_m0324 or not r["m0324"])]
    scored = [r for r in base if r["key"] in scores and scores[r["key"]]["n_st"] >= MIN_ST]
    out = dict(metric=key, s_max=s_max, include_m0324=include_m0324, n_pairs=len(base),
               n_scored=len(scored), n_families_scored=len({r["mult"] for r in scored}))
    S = [scores[r["key"]]["S"] for r in scored]
    rho = spearman(S, [r[key] for r in scored])
    cb = family_boot(scored, lambda sub: spearman([scores[r["key"]]["S"] for r in sub],
                                                  [r[key] for r in sub]), nboot, seed + 1)
    p05 = _pct(cb, 0.05)
    out["calibration"] = dict(rho=rho, p05=p05, p95=_pct(cb, 0.95),
                              passed=bool(isf(rho) and rho >= CAL_RHO and isf(p05) and p05 > 0))
    V = [r for r in scored if scores[r["key"]]["S"] <= s_max + 1e-12]
    out["clock"] = clock(V, key, nboot, seed)
    if not out["calibration"]["passed"]:
        out["decision"] = "untestable-calibration"
    else:
        out["decision"] = out["clock"]["decision"]
    sep_all = spearman([r.get("sep_km", NAN) for r in scored], [r[key] for r in scored])
    sep_v = spearman([r.get("sep_km", NAN) for r in V], [r[key] for r in V])
    out["diagnostic_sep_rho"] = dict(scored=sep_all, verified=sep_v)
    return out


def grade(rows_by_array, scores, nboot=NBOOT):
    yka, ilar = rows_by_array.get("yka", []), rows_by_array.get("ilar", [])
    g = dict(primary=clock_test(yka, scores, "D", nboot=nboot))
    g["decision"] = g["primary"]["decision"]
    g["secondary_ilar_D"] = clock_test(ilar, scores, "D", nboot=nboot)
    g["secondary_yka_D_c310"] = clock_test(yka, scores, "D_c310", nboot=nboot)
    g["sensitivity"] = {
        "yka_D_S_le_0.25": clock_test(yka, scores, "D", s_max=0.25, nboot=nboot),
        "yka_D_S_le_0.5": clock_test(yka, scores, "D", s_max=0.5, nboot=nboot),
        "yka_D_with_m0324": clock_test(yka, scores, "D", include_m0324=True, nboot=nboot),
        "yka_D_c13": clock_test(yka, scores, "D_c13", nboot=nboot),
    }
    return g


# ----------------------------------------------------------------------------- data access
def fetch_text(session, url, params=None, tries=10, wait=60, must=None):
    for k in range(tries):
        body = m2.get(session, url, params or {})
        if body:
            text = body.decode("utf-8", "replace")
            if must is None or must in text:
                return text
        if k + 1 < tries:
            time.sleep(wait)
    return None


def prereg_public(session, wait_min=20.0):
    t0 = time.time()
    while True:
        try:
            r = session.get(RAW + "main/" + PREREG, timeout=30)
            if r.status_code == 200 and b"pre-registration 13" in r.content.lower():
                return True, round((time.time() - t0) / 60.0, 1)
        except Exception as ex:  # noqa: BLE001
            print("prereg check", repr(ex)[:120], flush=True)
        if time.time() - t0 > wait_min * 60:
            return False, round((time.time() - t0) / 60.0, 1)
        time.sleep(60)


def station_coords(session):
    body = m2.get(session, m2.STATION, dict(net="IU,GT", sta=",".join(s for _, s, _, _ in STATIONS),
                                            cha="BHZ", level="station", format="text"))
    found = {}
    if body:
        try:
            for net, sta, lat, lon in m2.parse_station_text(body.decode("utf-8", "replace")):
                found[(net, sta)] = (lat, lon)
        except Exception:  # noqa: BLE001
            found = {}
    return [(net, sta) + found.get((net, sta), (lat, lon)) + ((net, sta) in found,)
            for net, sta, lat, lon in STATIONS]


def run_station(session, model, net, sta, slat, slon, events, fams, cache, deadline):
    from obspy import UTCDateTime
    from obspy.geodetics import locations2degrees
    name = "%s.%s" % (net, sta)
    center = int(round(m2.PRE * m2.FS))
    X, raw, ev_rows = {}, {}, []
    needed = sorted({e for f in fams.values() for e in f["members"]}, key=lambda e: events[e]["time"])
    stopped = False
    for eid in needed:
        if time.time() > deadline:
            stopped = True
            break
        E = events[eid]
        row = dict(station=name, event=eid, time=E["time"], status="", raw_snr="")
        ev_rows.append(row)
        try:
            o = UTCDateTime(E["time"])
            dist = locations2degrees(E["lat"], E["lon"], slat, slon)
            tp = m2.p_time(model, E["depth"], dist)
            if tp is None:
                row["status"] = "no_P"
                continue
            tr, inv, loc = m2.load_trace(session, net, sta, o + tp - m2.FETCH_PRE, o + tp + m2.FETCH_POST, cache)
            if tr is None:
                row["status"] = "no_data"
                continue
            x, resp = m2.process(tr, inv, o + tp - m2.PRE)
        except Exception as ex:  # noqa: BLE001
            row["status"] = "error:" + repr(ex)[:80]
            continue
        if not np.all(np.isfinite(x)) or not np.any(x):
            row["status"] = "bad_samples"
            continue
        X[eid] = x
        raw[eid] = m2.window_ratio(x, m2.FS, center, m2.RAW_SIG, m2.RAW_NOISE, "peak")
        row.update(status="ok", raw_snr=round(raw[eid], 2) if np.isfinite(raw[eid]) else "")
    fam_rows, pair_rows = [], []
    for fname, f in sorted(fams.items()):
        try:
            res, gates, pairs = measure_family(X, raw, f["members"], f["need"])
        except Exception as ex:  # noqa: BLE001
            res, gates, pairs = dict(status="exception:" + repr(ex)[:80]), {}, {}
        fam_rows.append(dict(station=name, mult=fname, n_members=len(f["members"]), n_need=len(f["need"]),
                             **{k: v for k, v in res.items()}))
        for (a, b), r in sorted(pairs.items()):
            pair_rows.append(dict(station=name, mult=fname, ev1=a, ev2=b, cc_direct=round(r["cc_direct"], 4),
                                  D_early=round(r["D_early"], 5) if isf(r["D_early"]) else "",
                                  D_late=round(r["D_late"], 5)))
    st = dict(station=name, n_events_needed=len(needed), n_traces=len(X), stopped_time_budget=stopped,
              n_family_ok=sum(1 for r in fam_rows if r.get("status") == "ok"), n_pairs=len(pair_rows))
    return st, ev_rows, fam_rows, pair_rows


def plot(out_dir, rows_by_array, scores):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    os.makedirs(os.path.join(out_dir, "plots"), exist_ok=True)
    yka = [r for r in rows_by_array.get("yka", []) if not r["m0324"] and isf(r["D"]) and r["key"] in scores
           and scores[r["key"]]["n_st"] >= MIN_ST]
    V = [r for r in yka if scores[r["key"]]["S"] <= S_MAX + 1e-12]
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.5))
    ax[0].scatter([scores[r["key"]]["S"] for r in yka], [r["D"] for r in yka], s=6, alpha=0.5)
    ax[0].axvline(S_MAX, color="k", lw=0.8)
    ax[0].set_xlabel("mantle-P source score S (low = sources alike)")
    ax[0].set_ylabel("YKA core D")
    for j, (k, lab) in enumerate((("dB", "backtracking-clock distance dB (yr)"), ("dM", "time apart dM (yr)"))):
        c = ["tab:red" if r["cross"] else "tab:blue" for r in V]
        ax[j + 1].scatter([r[k] for r in V], [r["D"] for r in V], s=8, c=c, alpha=0.7)
        ax[j + 1].set_xlabel(lab)
        ax[j + 1].set_ylabel("YKA core D (verified pairs; red = across 2008.5)")
    fig.suptitle("Pre-registration 13: source-verified pairs (mantlep3 %s)" % VERSION)
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, "plots", "h13_overview.png"), dpi=100)
    plt.close(fig)


def write_csv(path, rows):
    m2.write_csv(path, rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", default="/data")
    ap.add_argument("--out-dir", default="/out")
    ap.add_argument("--time-budget-min", type=float, default=200.0)
    args = ap.parse_args()
    t0 = time.time()
    deadline = t0 + 60.0 * args.time_budget_min
    import requests
    from obspy.taup import TauPyModel
    cache = os.path.join(args.data_dir, "mantlep")       # shared with mantlep v1/v2
    os.makedirs(cache, exist_ok=True)
    os.makedirs(args.out_dir, exist_ok=True)
    session = requests.Session()
    session.headers["User-Agent"] = "inner-core-repeaters (MV2, Munim, Inc.)"
    s = dict(version=VERSION, started=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
             tables=TABLES, stations=[("%s.%s" % (n, t)) for n, t, _, _ in STATIONS],
             rules=dict(S_MAX=S_MAX, MIN_ST=MIN_ST, MIN_STATION_PAIRS=MIN_STATION_PAIRS, CAL_RHO=CAL_RHO,
                        RHO_MIN=RHO_MIN, RHO_NONE_P95=RHO_NONE_P95, FRAC=FRAC, MIN_V=MIN_V,
                        MIN_FAM=MIN_FAM, MIN_CROSS=MIN_CROSS, NBOOT=NBOOT, EDGE_S=EDGE_S,
                        TURN=TURN, RATIO=RATIO))
    pub, waited = prereg_public(session)
    s["prereg_public_at_start"] = pub
    s["prereg_wait_min"] = waited
    print("prereg public:", pub, waited, flush=True)
    rows_by_array, all_rows = {}, []
    for arr, (commit, path) in TABLES.items():
        text = fetch_text(session, RAW + commit + "/" + path, must="ev1")
        rows = read_pairs(text, arr) if text else []
        rows_by_array[arr] = rows
        all_rows += rows
        s.setdefault("n_ok_rows", {})[arr] = len(rows)
    cat_text = fetch_text(session, USGS, USGS_PARAMS, must="latitude")
    cat = read_catalog(cat_text) if cat_text else dict(by_id={}, items=[], keys=[])
    s["catalog_events"] = len(cat["items"])
    events, miss = build_events(all_rows, cat)
    s["n_events"], s["events_unlocated"] = len(events), miss
    for r in all_rows:
        if r["ev1"] in events and r["ev2"] in events:
            a, b = events[r["ev1"]], events[r["ev2"]]
            r["sep_km"] = sep_km(a["lat"], a["lon"], b["lat"], b["lon"])
    fams = families(all_rows, events)
    s["n_families"] = len(fams)
    print("events", len(events), "families", len(fams), "rows", {k: len(v) for k, v in rows_by_array.items()},
          flush=True)
    model = TauPyModel("ak135")
    st_rows, ev_rows, fam_rows, pair_rows = [], [], [], []
    for net, sta, lat, lon, from_service in station_coords(session):
        if time.time() > deadline:
            st_rows.append(dict(station="%s.%s" % (net, sta), skipped_time_budget=True))
            continue
        try:
            st, er, fr, pr = run_station(session, model, net, sta, lat, lon, events, fams, cache, deadline)
        except Exception as ex:  # noqa: BLE001
            traceback.print_exc()
            st, er, fr, pr = dict(station="%s.%s" % (net, sta), error=repr(ex)[:300]), [], [], []
        st["coords_from_service"] = from_service
        st_rows.append(st)
        ev_rows += er
        fam_rows += fr
        pair_rows += pr
        print(json.dumps(m2.jsonable(st)), "%.1f min" % ((time.time() - t0) / 60.0), flush=True)
    s["stations_run"] = st_rows
    scores, used = source_scores(pair_rows)
    s["stations_scored"] = used
    s["n_pairs_scored"] = sum(1 for v in scores.values() if v["n_st"] >= MIN_ST)
    s["grade"] = grade(rows_by_array, scores)
    s["finished"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    s["minutes"] = round((time.time() - t0) / 60.0, 1)
    write_csv(os.path.join(args.out_dir, "mantlep3_pairs.csv"), pair_rows)
    write_csv(os.path.join(args.out_dir, "mantlep3_families.csv"), fam_rows)
    write_csv(os.path.join(args.out_dir, "mantlep3_events.csv"), ev_rows)
    write_csv(os.path.join(args.out_dir, "h13_scores.csv"),
              [dict(ev1=k[0], ev2=k[1], S=round(v["S"], 4), n_st=v["n_st"],
                    D_late_median=round(v["D_late_median"], 5)) for k, v in sorted(scores.items())])
    out_rows = []
    for r in all_rows:
        sc = scores.get(r["key"], {})
        out_rows.append(dict(array=r["array"], mult=r["mult"], ev1=r["ev1"], ev2=r["ev2"], t1=round(r["t1"], 4),
                             t2=round(r["t2"], 4), dB=round(r["dB"], 4), dM=round(r["dM"], 4), cross=r["cross"],
                             m0324=r["m0324"], D=r["D"], D_c310=r["D_c310"], sep_km=r.get("sep_km", ""),
                             S=sc.get("S", ""), n_st=sc.get("n_st", 0)))
    write_csv(os.path.join(args.out_dir, "h13_pairs_scored.csv"), out_rows)
    with open(os.path.join(args.out_dir, "h13_summary.json"), "w") as f:
        json.dump(m2.jsonable(s), f, indent=1)
    try:
        plot(args.out_dir, rows_by_array, scores)
    except Exception as ex:  # noqa: BLE001
        print("plot failed", repr(ex)[:200], flush=True)
    print(json.dumps(m2.jsonable(dict(decision=s["grade"]["decision"],
                                      primary=s["grade"]["primary"])))[:4000], flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())

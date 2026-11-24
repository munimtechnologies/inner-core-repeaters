#!/usr/bin/env python3
# Copyright 2026 MV2, Munim, Inc.
# Licensed under the Apache License, Version 2.0 (Apache-2.0). See LICENSE.
# Made with AI agents working in MV2's name.
"""ks32.py - pre-registration 32: is the reported 1 s SKS advance of the Kermadec doublet
in the direct SKS pulse (outer core), or only in a long low-frequency window that also holds
other phases (Zhang & Wen's SPdKS/SKPdS reading)?

Plan: findings/2026-11-24-prereg32-kermadec-sks-plan.md (published before this code
downloads any waveform). All decision rules are in PARAMS and the grade_* functions.

Claim under test: Zhou (2022), Commun. Earth Environ. 3:97 (doi 10.1038/s43247-022-00432-7):
SKS of the Kermadec pair 1997-05-03 (Mw 6.9, 108 km) / 2018-09-10 (Mw 6.9, 115 km) arrives
about 1 s earlier in 2018 at COLA, ULN, INK and YAK "when any of the reference waves (Pdiff,
Sdiff, PP, or SS) from the two earthquakes are aligned"; band 20-100 s; SKS windows of
Table S2 (s after origin): COLA 1411-1463, ULN 1439-1482, INK 1447-1494, YAK 1429-1471;
values -1.05, -1.0, -0.8, -1.1 s. Zhang & Wen (ESSOAr 2023, doi
10.22541/essoar.167276417.78056130/v1) see no SKS change at higher frequency and attribute
the shift to SPdKS/SKPdS. Zhou replied (doi 10.22541/essoar.167458062.24504379/v1).

Measured: for each station and event pair (A earlier, B later), each phase's arrival change
relative to its own iasp91 prediction (each event at its own USGS location, depth, origin),
d = (T_B,obs - T_B,pred) - (T_A,obs - T_A,pred), by normalised cross-correlation of A and B.
Station value x = median over passing references of d(SKS) - d(ref). Origin-time errors
cancel in x.

Top level imports are stdlib + numpy only; ObsPy is imported inside functions. Data access,
response removal and rotation reuse sdv23.py (pre-registration 23), downloaded unchanged by
the job from tag v0.32.0.
"""

import argparse
import csv
import hashlib
import json
import math
import os
import sys
import time

import numpy as np

VERSION = "0.1.0"
PLAN_PATH = "findings/2026-11-24-prereg32-kermadec-sks-plan.md"

# USGS ComCat, M >= 5.8, depth 70-160 km, within 40 km of (31.77 S, 179.38 W), 1990-2026,
# queried 2026-11-23. Columns: id, origin (UTC), lat, lon, depth km, mag, magType, usgs id.
EVENTS = [
    ("E1997", "1997-05-03T16:46:02.02", -31.791, -179.382, 108.3, 6.9, "mwc", "usp00081pj"),
    ("E1999", "1999-04-20T19:04:08.32", -31.888, -179.040, 95.7, 6.5, "mwb", "usp00096mm"),
    ("E2006", "2006-05-16T10:39:23.34", -31.810, -179.307, 152.0, 7.4, "mwb", "usp000eh87"),
    ("E2018", "2018-09-10T04:19:02.63", -31.7447, -179.3728, 115.0, 6.9, "mww", "us2000hbpv"),
]
PRIMARY_PAIR = ("E1997", "E2018")
STATIONS = [("IU", "COLA"), ("IU", "ULN"), ("CN", "INK"), ("IU", "YAK")]

# Zhou (2022) Table S2: SKS measurement windows (s after origin) and the reported smallest
# differential traveltimes (s), as copied from the supplement (hands report 2026-11-23).
ZHOU_WIN = {"COLA": (1411.0, 1463.0), "ULN": (1439.0, 1482.0),
            "INK": (1447.0, 1494.0), "YAK": (1429.0, 1471.0)}
ZHOU_VAL = {"COLA": -1.05, "ULN": -1.0, "INK": -0.8, "YAK": -1.1}

# phase key -> (TauP names (first arrival among them), component)
PHASES = {
    "SKS": (["SKS"], "R"),
    "SKKS": (["SKKS"], "R"),
    "P1": (["P", "Pdiff"], "Z"),
    "PP": (["PP"], "Z"),
    "S1": (["S", "Sdiff"], "T"),
    "SS": (["SS"], "T"),
}
REFS = ["P1", "PP", "S1", "SS"]

PARAMS = dict(
    band_long=(0.01, 0.05),    # Hz, Zhou's 20-100 s band (T1)
    band_short=(0.1, 0.5),     # Hz, direct-pulse band (T2, graded)
    band_mid=(0.05, 0.2),      # Hz, reported only
    sr=20.0,                   # Hz, every trace resampled to this rate first
    tb=600.0, te=2600.0,       # s after origin
    pre_filt=(0.002, 0.004, 2.0, 4.0),
    output="VEL",
    ref_pre=10.0, ref_post=40.0,        # long-band reference windows around prediction
    fallback_pre=10.0, fallback_post=40.0,
    maxlag_long=8.0,
    peak_search=5.0, short_pre=3.0, short_post=4.0, maxlag_short=1.5,
    noise_rel=(-60.0, -10.0),
    cc_long=0.90, cc_short=0.80, snr_min=3.0,
    advance_s=-0.5, min_stations=3,
)

T1_LABELS = ("reproduced", "partial", "not-reproduced")
T2_LABELS = ("present", "mixed", "absent")


def log(*a):
    print(time.strftime("%H:%M:%S"), "[ks32]", *a, flush=True)


# ----------------------------------------------------------------- pure numerics

def window_lag(a, b, i1, i2, maxlag):
    """Same algorithm as sdv23.window_lag: normalised cross-correlation of a[i1:i2] with
    b[i1+k:i2+k], |k| <= maxlag samples. Returns (lag, cc, edge); lag > 0 means b is
    delayed relative to a; parabolic refinement."""
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    i1, i2, maxlag = int(i1), int(i2), int(maxlag)
    if i2 - i1 < 3 or i1 < 0 or i2 > len(a) or i1 - maxlag < 0 or i2 + maxlag > len(b):
        return float("nan"), float("nan"), True
    sa = a[i1:i2] - a[i1:i2].mean()
    na = math.sqrt(float(np.dot(sa, sa)))
    ccs = np.full(2 * maxlag + 1, np.nan)
    for j, k in enumerate(range(-maxlag, maxlag + 1)):
        sb = b[i1 + k:i2 + k]
        sb = sb - sb.mean()
        nb = math.sqrt(float(np.dot(sb, sb)))
        if na > 0 and nb > 0:
            ccs[j] = float(np.dot(sa, sb)) / (na * nb)
    if not np.any(np.isfinite(ccs)):
        return float("nan"), float("nan"), True
    j = int(np.nanargmax(ccs))
    if j == 0 or j == len(ccs) - 1:
        return float(j - maxlag), float(ccs[j]), True
    y0, y1, y2 = ccs[j - 1], ccs[j], ccs[j + 1]
    den = y0 - 2.0 * y1 + y2
    d = 0.5 * (y0 - y2) / den if den != 0 else 0.0
    return float(j - maxlag + d), float(y1 - 0.25 * (y0 - y2) * d), False


def shift_left(b, s):
    """c[i] = b[i + s] (zeros outside b). s may be negative."""
    b = np.asarray(b, float)
    s = int(s)
    c = np.zeros_like(b)
    n = len(b)
    if s >= 0:
        if s < n:
            c[:n - s] = b[s:]
    elif -s < n:
        c[-s:] = b[:n + s]
    return c


def lag_in_window(a, b, sr, tb, w1, w2, shift_s, maxlag_s):
    """Delay (s) of trace b relative to trace a, both sampled at sr from tb s after their own
    origins. The window [w1, w2] (s after origin) is set on a; b is first moved earlier by
    shift_s (rounded to samples), then searched within +/- maxlag_s.
    Returns (total_delay_s, cc, edge)."""
    s = int(round(shift_s * sr))
    c = shift_left(b, s)
    i1 = int(round((w1 - tb) * sr))
    i2 = int(round((w2 - tb) * sr))
    k, cc, edge = window_lag(a, c, i1, i2, int(round(maxlag_s * sr)))
    if not np.isfinite(k):
        return float("nan"), float("nan"), True
    return float((k + s) / sr), float(cc), bool(edge)


def phase_change(a, b, sr, tb, t_pred_a, t_pred_b, pre, post, maxlag_s):
    """Arrival change of one phase, each event relative to its own prediction; window
    [t_pred_a - pre, t_pred_a + post] on a. Returns (d_s, cc, edge)."""
    tot, cc, edge = lag_in_window(a, b, sr, tb, t_pred_a - pre, t_pred_a + post,
                                  t_pred_b - t_pred_a, maxlag_s)
    if not np.isfinite(tot):
        return float("nan"), float("nan"), True
    return float(tot - (t_pred_b - t_pred_a)), cc, edge


def snr_win(x, sr, tb, n1, n2, s1, s2):
    """Peak |x| in [s1, s2] over rms in [n1, n2] (s after origin)."""
    x = np.asarray(x, float)
    i = [int(round((t - tb) * sr)) for t in (n1, n2, s1, s2)]
    noise = x[max(0, i[0]):max(0, i[1])]
    sig = x[max(0, i[2]):max(0, i[3])]
    if len(noise) < 3 or len(sig) < 1:
        return float("nan")
    rms = math.sqrt(float(np.mean(noise ** 2)))
    return float(np.max(np.abs(sig)) / rms) if rms > 0 else float("inf")


def snr_at(x, sr, tb, t_phase, noise_rel, sig_pre, sig_post):
    """Peak |x| in [t_phase - sig_pre, t_phase + sig_post] over rms in the noise window
    t_phase + noise_rel."""
    return snr_win(x, sr, tb, t_phase + noise_rel[0], t_phase + noise_rel[1],
                   t_phase - sig_pre, t_phase + sig_post)


def peak_time(x, sr, tb, t0, half):
    """Time (s after origin) of the largest |x| within [t0 - half, t0 + half]."""
    x = np.asarray(x, float)
    i1 = max(0, int(round((t0 - half - tb) * sr)))
    i2 = min(len(x), int(round((t0 + half - tb) * sr)) + 1)
    if i2 - i1 < 1:
        return float("nan")
    j = int(np.argmax(np.abs(x[i1:i2])))
    return float(tb + (i1 + j) / sr)


def zhou_window(sta, tsks_1997, tsks_a, P=None):
    """SKS window (s after A's origin) for the long-band measurement. Zhou's Table S2
    window is used if it contains the iasp91 SKS time of E1997 at that station; it is kept
    at the same offsets from the predicted SKS of event A. Otherwise the fallback
    [tSKS_A - fallback_pre, tSKS_A + fallback_post]. Returns (w1, w2, source)."""
    P = P or PARAMS
    w = ZHOU_WIN.get(sta)
    if w is not None and w[0] <= tsks_1997 <= w[1]:
        return (float(tsks_a + (w[0] - tsks_1997)), float(tsks_a + (w[1] - tsks_1997)),
                "table-s2")
    return (float(tsks_a - P["fallback_pre"]), float(tsks_a + P["fallback_post"]), "fallback")


def passes(m, cc_min, snr_min):
    """m: dict(d, cc, edge, snr_a, snr_b)."""
    try:
        return bool(np.isfinite(m["d"]) and not m["edge"] and m["cc"] >= cc_min
                    and m["snr_a"] >= snr_min and m["snr_b"] >= snr_min)
    except (KeyError, TypeError):
        return False


def station_value(sks, refs, cc_min, snr_min):
    """x = median over passing references of d(SKS) - d(ref); c = median d(ref).
    Returns dict(x, c, n_ref, refs, status)."""
    nan = float("nan")
    if not passes(sks, cc_min, snr_min):
        return dict(x=nan, c=nan, n_ref=0, refs=[], status="sks-gate")
    ok = sorted(k for k, m in refs.items() if passes(m, cc_min, snr_min))
    if not ok:
        return dict(x=nan, c=nan, n_ref=0, refs=[], status="no-ref")
    diffs = [sks["d"] - refs[k]["d"] for k in ok]
    return dict(x=float(np.median(diffs)), c=float(np.median([refs[k]["d"] for k in ok])),
                n_ref=len(ok), refs=ok, status="ok")


def differentials(changes, target="SKS", refs=REFS):
    """changes: phase -> d (s). Returns ref -> d(target) - d(ref) for finite pairs."""
    out = {}
    dt = changes.get(target)
    if dt is None or not np.isfinite(dt):
        return out
    for r in refs:
        dr = changes.get(r)
        if dr is not None and np.isfinite(dr):
            out[r] = float(dt - dr)
    return out


def grade_count(values, thr, min_n, labels):
    """values: station -> x (nan = untestable). labels = (positive, middle, negative).
    N finite < min_n -> untestable; k = #(x <= thr): k >= 3 positive, k <= 1 negative,
    else middle."""
    xs = [v for v in values.values() if v is not None and np.isfinite(v)]
    if len(xs) < min_n:
        return "untestable"
    k = sum(1 for v in xs if v <= thr)
    if k >= 3:
        return labels[0]
    if k <= 1:
        return labels[2]
    return labels[1]


def grade_joint(t1, t2):
    """Pre-registered reading of T1 (long window, Zhou's band) with T2 (direct pulse)."""
    if t1 == "untestable":
        return "untestable"
    if t1 == "not-reproduced":
        return "not-reproduced"
    if t1 == "reproduced":
        if t2 == "present":
            return "direct-SKS-advance"
        if t2 == "absent":
            return "window-content"
        return "mixed"
    return "mixed"


def jsonable(o):
    if isinstance(o, np.bool_):
        return bool(o)
    if isinstance(o, np.integer):
        return int(o)
    if isinstance(o, np.floating):
        return float(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    return str(o)


# ----------------------------------------------------------------- data side (ObsPy)

def predict(model, depth, dist_deg):
    """Phase key -> first predicted time (s) among its TauP names."""
    names = sorted({n for v in PHASES.values() for n in v[0]})
    arr = model.get_travel_times(source_depth_in_km=depth, distance_in_degree=dist_deg,
                                 phase_list=names)
    out = {}
    for key, (cands, _c) in PHASES.items():
        ts = [float(a.time) for a in arr if a.name in cands]
        if ts:
            out[key] = min(ts)
    return out


def station_event(base, D, model, ev, net, sta, P, errors):
    """One event at one station -> dict(dist, baz, preds, loc, ids, bands{band: R/T/Z})."""
    from obspy import UTCDateTime
    from obspy.geodetics import gps2dist_azimuth, locations2degrees
    t0 = UTCDateTime(ev[1])
    inv = D.stations("ks32_%s_%s_%s" % (ev[0], net, sta), network=net, station=sta,
                     channel="BH?", level="response", starttime=t0, endtime=t0 + 3600)
    if inv is None:
        errors.append("%s %s.%s: no station metadata" % (ev[0], net, sta))
        return None
    st_lat = inv[0][0].latitude
    st_lon = inv[0][0].longitude
    dist = locations2degrees(ev[2], ev[3], st_lat, st_lon)
    _, _, baz = gps2dist_azimuth(st_lat, st_lon, ev[2], ev[3])
    st = D.waveforms(net, sta, "*", "BH?", t0 + P["tb"] - 300, t0 + P["te"] + 300)
    if st is None or len(st) == 0:
        errors.append("%s %s.%s: no waveforms" % (ev[0], net, sta))
        return None
    locs = sorted({tr.stats.location for tr in st})
    loc = "00" if "00" in locs else ("" if "" in locs else locs[0])
    st = st.select(location=loc)
    for tr in st:
        if abs(tr.stats.sampling_rate - P["sr"]) > 1e-6:
            tr.resample(P["sr"])
    preds = predict(model, ev[4], dist)
    bands = {}
    ids = sorted({tr.id for tr in st})
    for band in (P["band_long"], P["band_short"], P["band_mid"]):
        try:
            comps = base.process(st, inv, t0, P["tb"], P["te"], P["output"], P["pre_filt"], band)
        except Exception as ex:  # noqa: BLE001
            errors.append("%s %s.%s band %s: process failed: %r" % (ev[0], net, sta, band, ex))
            continue
        if "Z" not in comps:
            continue
        rec = {"Z": comps["Z"]["x"], "sr": comps["Z"]["sr"]}
        if "1" in comps and "2" in comps:
            a1, a2 = comps["1"]["az"], comps["2"]["az"]
            if np.isfinite(a1) and np.isfinite(a2):
                n = min(len(comps["1"]["x"]), len(comps["2"]["x"]))
                r, t = base.to_rt(comps["1"]["x"][:n], comps["2"]["x"][:n], a1, a2, baz)
                rec["R"], rec["T"] = r, t
                rec["az"] = (a1, a2)
        bands[band] = rec
    return dict(dist=float(dist), baz=float(baz), preds=preds, loc=loc, ids=ids, bands=bands)


def measure_pair(A, B, sta, tsks97, P):
    """All measurements for one station and one pair. Returns (rows, values) where values
    = {"long": station_value, "short": ..., "mid": ..., plus raw (uncorrected) x}."""
    rows = []
    values = {}
    tb, sr = P["tb"], P["sr"]
    bl = P["band_long"]
    ra, rb = A["bands"].get(bl), B["bands"].get(bl)
    if ra is None or rb is None or "SKS" not in A["preds"] or "SKS" not in B["preds"]:
        return rows, {"long": dict(x=float("nan"), status="no-data")}

    def one(band, key, w1, w2, shift_extra, maxlag, rec_a, rec_b):
        comp = PHASES[key][1]
        ta, tbp = A["preds"][key], B["preds"][key]
        n = min(len(rec_a[comp]), len(rec_b[comp]))
        tot, cc, edge = lag_in_window(rec_a[comp][:n], rec_b[comp][:n], sr, tb, w1, w2,
                                      (tbp - ta) + shift_extra, maxlag)
        d = tot - (tbp - ta) if np.isfinite(tot) else float("nan")
        na1, na2 = ta + P["noise_rel"][0], ta + P["noise_rel"][1]
        nb1, nb2 = tbp + P["noise_rel"][0], tbp + P["noise_rel"][1]
        sh = (tbp - ta) + shift_extra
        sa = snr_win(rec_a[comp], sr, tb, na1, na2, w1, w2)
        sb = snr_win(rec_b[comp], sr, tb, nb1 + shift_extra, nb2 + shift_extra, w1 + sh, w2 + sh)
        m = dict(d=float(d), cc=float(cc), edge=bool(edge), snr_a=sa, snr_b=sb,
                 dp=float(tbp - ta))
        rows.append(dict(band="%.2f-%.2f" % band, phase=key, comp=comp, w1=round(w1, 2),
                         w2=round(w2, 2), pred_a=round(ta, 3), pred_b=round(tbp, 3), **m))
        return m

    # ---- long band (Zhou's band): T1
    have = [k for k in PHASES if k in A["preds"] and k in B["preds"] and PHASES[k][1] in ra
            and PHASES[k][1] in rb]
    if "SKS" not in have:
        return rows, {"long": dict(x=float("nan"), status="no-R")}
    w1, w2, wsrc = zhou_window(sta, tsks97, A["preds"]["SKS"], P)
    long_m = {"SKS": one(bl, "SKS", w1, w2, 0.0, P["maxlag_long"], ra, rb)}
    for key in REFS + ["SKKS"]:
        if key in have:
            t = A["preds"][key]
            long_m[key] = one(bl, key, t - P["ref_pre"], t + P["ref_post"], 0.0,
                              P["maxlag_long"], ra, rb)
    refs_l = {k: long_m[k] for k in REFS if k in long_m}
    v_long = station_value(long_m["SKS"], refs_l, P["cc_long"], P["snr_min"])
    v_long["window"] = [w1, w2, wsrc]
    if v_long["status"] == "ok":
        raw = [(long_m["SKS"]["d"] + long_m["SKS"]["dp"]) - (refs_l[k]["d"] + refs_l[k]["dp"])
               for k in v_long["refs"]]
        v_long["x_raw"] = float(np.median(raw))
        v_long["per_ref"] = {k: long_m["SKS"]["d"] - refs_l[k]["d"] for k in v_long["refs"]}
        if "SKKS" in long_m and passes(long_m["SKKS"], P["cc_long"], P["snr_min"]):
            v_long["x_skks"] = float(np.median([long_m["SKKS"]["d"] - refs_l[k]["d"]
                                                for k in v_long["refs"]]))
    values["long"] = v_long
    c = v_long.get("c", float("nan"))

    # ---- short (graded T2) and mid (reported) bands: direct pulses, centred on c
    for tag, band in (("short", P["band_short"]), ("mid", P["band_mid"])):
        rsa, rsb = A["bands"].get(band), B["bands"].get(band)
        if rsa is None or rsb is None or not np.isfinite(c):
            values[tag] = dict(x=float("nan"), status="no-data" if np.isfinite(c) else "no-c")
            continue
        ms = {}
        for key in ["SKS"] + REFS:
            comp = PHASES[key][1]
            if key not in have or comp not in rsa or comp not in rsb:
                continue
            tpk = peak_time(rsa[comp], sr, tb, A["preds"][key], P["peak_search"])
            if not np.isfinite(tpk):
                continue
            ms[key] = one(band, key, tpk - P["short_pre"], tpk + P["short_post"], c,
                          P["maxlag_short"], rsa, rsb)
        if "SKS" not in ms:
            values[tag] = dict(x=float("nan"), status="no-SKS")
            continue
        v = station_value(ms["SKS"], {k: ms[k] for k in REFS if k in ms}, P["cc_short"],
                          P["snr_min"])
        if v["status"] == "ok":
            v["per_ref"] = {k: ms["SKS"]["d"] - ms[k]["d"] for k in v["refs"]}
        values[tag] = v
    return rows, values


def run(args):
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    import sdv23 as base  # pinned copy downloaded by the job's setup (tag v0.32.0)
    from obspy.taup import TauPyModel
    P = PARAMS
    os.makedirs(args.out_dir, exist_ok=True)
    here = os.path.dirname(os.path.abspath(__file__))
    with open(os.path.join(here, "sdv23.py"), "rb") as f:
        sdv_sha = hashlib.sha256(f.read()).hexdigest()
    with open(os.path.abspath(__file__), "rb") as f:
        code_sha = hashlib.sha256(f.read()).hexdigest()
    plan_ok, plan_sha = base.check_plan(args.wait_min, PLAN_PATH)
    log("plan public:", plan_ok, plan_sha[:16])
    errors = []
    D = base.Data(args.data_dir)
    model = TauPyModel("iasp91")
    ev = {e[0]: e for e in EVENTS}
    rec = {}
    for k in ev:
        for net, sta in STATIONS:
            try:
                r = station_event(base, D, model, ev[k], net, sta, P, errors)
            except Exception as ex:  # noqa: BLE001
                errors.append("%s %s.%s: %r" % (k, net, sta, ex))
                r = None
            if r is not None:
                rec[(k, sta)] = r
                log("event", k, net, sta, "dist %.2f" % r["dist"], "loc '%s'" % r["loc"],
                    "bands", len(r["bands"]))
    rows_all, st_rows = [], []
    values = {}
    keys = sorted(ev)
    for i, a in enumerate(keys):
        for b in keys[i + 1:]:
            for net, sta in STATIONS:
                if (a, sta) not in rec or (b, sta) not in rec or ("E1997", sta) not in rec:
                    continue
                tsks97 = rec[("E1997", sta)]["preds"].get("SKS", float("nan"))
                try:
                    rows, vals = measure_pair(rec[(a, sta)], rec[(b, sta)], sta, tsks97, P)
                except Exception as ex:  # noqa: BLE001
                    errors.append("measure %s-%s %s: %r" % (a, b, sta, ex))
                    continue
                for r in rows:
                    r.update(A=a, B=b, net=net, sta=sta)
                rows_all.extend(rows)
                values[(a, b, sta)] = vals
                for tag, v in vals.items():
                    st_rows.append(dict(A=a, B=b, net=net, sta=sta,
                                        dist=round(rec[(a, sta)]["dist"], 3), band=tag,
                                        x=v.get("x"), x_raw=v.get("x_raw", ""),
                                        x_skks=v.get("x_skks", ""), c=v.get("c", ""),
                                        n_ref=v.get("n_ref", ""), refs="|".join(v.get("refs", [])),
                                        status=v.get("status", ""),
                                        window="|".join(str(w) for w in v.get("window", []))))
                    log(a, b, sta, tag, v.get("status"), v.get("x"))
    pa, pb = PRIMARY_PAIR
    t1_vals = {s: values.get((pa, pb, s), {}).get("long", {}).get("x", float("nan"))
               for _n, s in STATIONS}
    t2_vals = {s: values.get((pa, pb, s), {}).get("short", {}).get("x", float("nan"))
               for _n, s in STATIONS}
    t1 = grade_count(t1_vals, P["advance_s"], P["min_stations"], T1_LABELS)
    t2 = grade_count(t2_vals, P["advance_s"], P["min_stations"], T2_LABELS)
    joint = grade_joint(t1, t2)
    log("GRADED T1", t1, "T2", t2, "joint", joint)
    cols = ["A", "B", "net", "sta", "band", "phase", "comp", "w1", "w2", "pred_a", "pred_b",
            "dp", "d", "cc", "edge", "snr_a", "snr_b"]
    with open(os.path.join(args.out_dir, "ks32_phases.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        for r in rows_all:
            w.writerow({c: r.get(c, "") for c in cols})
    scols = ["A", "B", "net", "sta", "dist", "band", "x", "x_raw", "x_skks", "c", "n_ref", "refs",
             "status", "window"]
    with open(os.path.join(args.out_dir, "ks32_stations.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=scols)
        w.writeheader()
        for r in st_rows:
            w.writerow({c: r.get(c, "") for c in scols})
    meta = {"%s|%s" % k: dict(dist=v["dist"], baz=v["baz"], loc=v["loc"], ids=v["ids"],
                              preds=v["preds"], bands=["%.2f-%.2f" % b for b in v["bands"]])
            for k, v in rec.items()}
    summary = dict(
        version=VERSION, code_sha256=code_sha, sdv23_sha256=sdv_sha,
        plan_public_at_start=plan_ok, plan_sha256=plan_sha,
        primary_pair=list(PRIMARY_PAIR),
        graded=dict(T1=t1, T2=t2, joint=joint),
        T1_values=t1_vals, T2_values=t2_vals, zhou_values=ZHOU_VAL,
        primary_detail={s: values.get((pa, pb, s), {}) for _n, s in STATIONS},
        series_long={"%s-%s|%s" % k: v.get("long", {}).get("x") for k, v in values.items()},
        station_meta=meta, errors=errors, params=P,
        finished=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()))
    with open(os.path.join(args.out_dir, "ks32_summary.json"), "w") as f:
        json.dump(summary, f, indent=1, default=jsonable)
    log("done", len(rows_all), "phase rows;", len(errors), "errors")


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", default="/data")
    ap.add_argument("--out-dir", default="/out")
    ap.add_argument("--wait-min", type=float, default=120.0)
    run(ap.parse_args(argv))


if __name__ == "__main__":
    main()

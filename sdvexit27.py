#!/usr/bin/env python3
# Copyright 2026 MV2, Munim, Inc.
# Licensed under the Apache License, Version 2.0 (Apache-2.0). See LICENSE.
# Made with AI agents working in MV2's name.
"""sdvexit27.py - pre-registration 27: does SKS at IU.SDV change over a decade on OTHER
repeating-earthquake pairs whose SKS leaves the core near the same exit point as the
Aleutian doublet 2000-07-07 / 2009-12-17 (Zhang & Wen 2025, Nat. Commun. 16:8754)?

Plan: findings/2026-11-21-prereg27-sdv-exit-pairs-plan.md (published before this code runs;
the job waits for it and records its sha256).

Zhang & Wen model the decadal ULVZ change near the SKS entry point (source side) but say the
exit point (receiver side) "could explain the data equally well". Pre-registration 26 tests
the entry side with neighbour stations. This job tests the exit side with neighbour SOURCES:
pairs of nearly co-located earthquakes (one before the 2009-09-18 metadata step, one after
2009-12-17, >= 5 years apart) whose SKS to IU.SDV exits the core within 100 km of the
doublet's exit point but enters it far away. Controls are pairs at other azimuths.

Processing, relocation and membership are those of sdv23.py 0.1.1, imported unchanged from
its pinned commit (fetched at run time, sha256 recorded). Top-level imports: stdlib + numpy.
"""

import argparse
import csv
import hashlib
import io
import json
import math
import os
import sys
import tempfile
import time
import urllib.parse
import urllib.request

import numpy as np

VERSION = "0.1.0"
PLAN_PATH = "findings/2026-11-21-prereg27-sdv-exit-pairs-plan.md"
REPO = "https://raw.githubusercontent.com/munimtechnologies/inner-core-repeaters/"
SDV23_COMMIT = "dbbd91d3f26f9dc2165cbf59a5c05f7d82d8521d"
SDV23_MUST = b'VERSION = "0.1.1"'
P23_SUMMARY = "results/sdv23-v2/sdv23_summary.json"
USGS_CSV = "https://earthquake.usgs.gov/fdsnws/event/1/query"
R_CMB = 3480.0
R_EARTH = 6371.0
YEAR_S = 365.25 * 86400.0
WINDOWS = ("SKS_R", "SKS_Z", "S_T", "P_Z")
# IU.SDV location 00 nominal frame before the 2009-09-18 metadata step (BHN 0, BHE 90).
NOMINAL = (0.0, 90.0)

PARAMS = dict(
    cat_start="1999-11-06", pre_end="2009-09-17T23:59:59", post_start="2009-12-17",
    cat_end="2026-11-01", post_split="2018-01-01",
    min_mag=5.3, max_depth=100.0, ring=(85.0, 125.0),
    pair_sep_km=15.0, pair_dz_km=25.0, pair_dmag=0.8, min_span_yr=5.0,
    wedge_daz=20.0, wedge_dist=(85.0, 110.0), ctrl_daz=25.0, ctrl_dist=(88.0, 110.0),
    max_geom_events=400, max_ctrl_cands=40,
    entry_nb_km=100.0, exit_nb_km=100.0, control_min_km=300.0,
    caps=dict(exit=(25, 12), control=(15, 8), near=(10, 6)),   # (max checks, max verified)
    p_max_stations=12,
    min_s_sks_sep_s=10.0, min_skks_sep_s=6.0,
    peak_search_s=5.0,
    maxlag_s=0.6, sks_cc_min=0.80, sks_snr_min=5.0, ref_cc_min=0.90,
    changed_s=0.08, stable_s=0.04, z_agree_s=0.04,
    stable_frac=2.0 / 3.0,
    g0_tol_s=0.010, ctrl_min_n=3, ctrl_false_max=0.20,
)


def log(*a):
    print(time.strftime("%H:%M:%S"), "[sdvexit27]", *a, flush=True)


# ----------------------------------------------------------------- pure numerics

def _fin(x):
    try:
        return x is not None and math.isfinite(float(x))
    except (TypeError, ValueError):
        return False


def wrap180(x):
    return (float(x) + 180.0) % 360.0 - 180.0


def dest_point(lat, lon, az, ang):
    """Point at angular distance ang (deg) from (lat, lon) along initial azimuth az (deg)."""
    p1, l1, a, d = (math.radians(v) for v in (lat, lon, az, ang))
    sp2 = math.sin(p1) * math.cos(d) + math.cos(p1) * math.sin(d) * math.cos(a)
    p2 = math.asin(max(-1.0, min(1.0, sp2)))
    l2 = l1 + math.atan2(math.sin(a) * math.sin(d) * math.cos(p1),
                         math.cos(d) - math.sin(p1) * math.sin(p2))
    return math.degrees(p2), (math.degrees(l2) + 180.0) % 360.0 - 180.0


def ang_dist(lat1, lon1, lat2, lon2):
    """Great-circle angular distance (deg), haversine."""
    p1, p2 = math.radians(lat1), math.radians(lat2)
    dp, dl = p2 - p1, math.radians(lon2 - lon1)
    h = math.sin(dp / 2) ** 2 + math.cos(p1) * math.cos(p2) * math.sin(dl / 2) ** 2
    return math.degrees(2.0 * math.asin(min(1.0, math.sqrt(h))))


def azimuth(lat1, lon1, lat2, lon2):
    """Initial azimuth (deg, 0..360) from point 1 to point 2."""
    p1, p2 = math.radians(lat1), math.radians(lat2)
    dl = math.radians(lon2 - lon1)
    y = math.sin(dl) * math.cos(p2)
    x = math.cos(p1) * math.sin(p2) - math.sin(p1) * math.cos(p2) * math.cos(dl)
    return math.degrees(math.atan2(y, x)) % 360.0


def arc_km(lat1, lon1, lat2, lon2, r=R_CMB):
    return math.radians(ang_dist(lat1, lon1, lat2, lon2)) * r


def hav_km_vec(lat, lon, lats, lons):
    """Surface distance (km) from one point to arrays of points."""
    p1 = math.radians(lat)
    p2 = np.radians(np.asarray(lats, float))
    dp = p2 - p1
    dl = np.radians(np.asarray(lons, float) - lon)
    h = np.sin(dp / 2) ** 2 + math.cos(p1) * np.cos(p2) * np.sin(dl / 2) ** 2
    return 2.0 * R_EARTH * np.arcsin(np.minimum(1.0, np.sqrt(h)))


def cmb_indices(depths, cmb, tol=0.5):
    idx = [i for i, d in enumerate(depths) if abs(float(d) - cmb) <= tol]
    if not idx:
        return None, None
    return idx[0], idx[-1]


def iso_epoch(iso):
    """'YYYY-MM-DDTHH:MM:SS(.fff)(Z)' or 'YYYY-MM-DD' -> POSIX seconds (UTC)."""
    import calendar
    s = iso.strip().rstrip("Z")
    date, _, clock = s.partition("T")
    y, mo, d = (int(v) for v in date.split("-"))
    hh, mm, ss = 0, 0, 0.0
    if clock:
        parts = clock.split(":")
        hh, mm = int(parts[0]), int(parts[1])
        ss = float(parts[2]) if len(parts) > 2 else 0.0
    return calendar.timegm((y, mo, d, hh, mm, 0, 0, 0, 0)) + ss


def parse_usgs_csv(text):
    """USGS fdsnws-event CSV -> list of dict(id, iso, lat, lon, depth, mag, magtype, ts)."""
    out = []
    for r in csv.DictReader(io.StringIO(text)):
        if (r.get("type") or "earthquake") != "earthquake":
            continue
        try:
            ev = dict(id=r["id"], iso=r["time"].rstrip("Z"), lat=float(r["latitude"]),
                      lon=float(r["longitude"]), depth=float(r["depth"] or "nan"),
                      mag=float(r["mag"]), magtype=r.get("magType", ""))
        except (KeyError, ValueError, TypeError):
            continue
        if not all(_fin(ev[k]) for k in ("lat", "lon", "depth", "mag")):
            continue
        ev["ts"] = iso_epoch(ev["iso"])
        out.append(ev)
    return out


def find_pairs(pre, post, P=PARAMS):
    """Candidate pairs (A from pre, B from post) by separation, depth, magnitude and span."""
    if not pre or not post:
        return []
    la = np.array([e["lat"] for e in post])
    lo = np.array([e["lon"] for e in post])
    dz = np.array([e["depth"] for e in post])
    mg = np.array([e["mag"] for e in post])
    ts = np.array([e["ts"] for e in post])
    out = []
    for A in pre:
        d = hav_km_vec(A["lat"], A["lon"], la, lo)
        m = ((d <= P["pair_sep_km"]) & (np.abs(dz - A["depth"]) <= P["pair_dz_km"])
             & (np.abs(mg - A["mag"]) <= P["pair_dmag"])
             & (ts - A["ts"] >= P["min_span_yr"] * YEAR_S))
        for j in np.nonzero(m)[0]:
            B = post[int(j)]
            out.append(dict(A=A, B=B, sep_km=float(d[j]),
                            span_yr=float((B["ts"] - A["ts"]) / YEAR_S),
                            mag_min=float(min(A["mag"], B["mag"]))))
    return out


def rank_key(p):
    return (round(p["sep_km"], 3), -p["mag_min"], p["A"]["iso"], p["B"]["iso"])


def greedy_unique(pairs):
    """Pairs in rank order, each event used at most once."""
    used, out = set(), []
    for p in sorted(pairs, key=rank_key):
        a, b = p["A"]["id"], p["B"]["id"]
        if a in used or b in used:
            continue
        used.add(a)
        used.add(b)
        out.append(p)
    return out


def classify_pair(d_entry, d_exit, P=PARAMS):
    e_ok, x_ok = _fin(d_entry), _fin(d_exit)
    if e_ok and x_ok and d_entry <= P["entry_nb_km"] and d_exit <= P["exit_nb_km"]:
        return "near"
    if x_ok and d_exit <= P["exit_nb_km"] and e_ok and d_entry > P["entry_nb_km"]:
        return "exit"
    if e_ok and x_ok and d_entry > P["control_min_km"] and d_exit > P["control_min_km"]:
        return "control"
    return "other"


def is_doublet(p, t_a, t_b, tol=60.0):
    return abs(p["A"]["ts"] - t_a) <= tol and abs(p["B"]["ts"] - t_b) <= tol


def shift_left(b, s):
    b = np.asarray(b, float)
    s = int(s)
    c = np.zeros_like(b)
    n = len(b)
    if s >= 0:
        if s < n:
            c[:n - s] = b[s:]
    else:
        if -s < n:
            c[-s:] = b[:n + s]
    return c


def window_lag(a, b, i1, i2, maxlag):
    """Normalised cross-correlation of a[i1:i2] with b[i1+k:i2+k], |k| <= maxlag samples.
    Returns (lag, cc, edge); lag > 0 means b is delayed (as sdv23.window_lag)."""
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


def snr(x, n1, n2, s1, s2):
    x = np.asarray(x, float)
    noise = x[max(0, int(n1)):max(0, int(n2))]
    sig = x[max(0, int(s1)):max(0, int(s2))]
    if len(noise) < 3 or len(sig) < 1:
        return float("nan")
    rms = math.sqrt(float(np.mean(noise ** 2)))
    return float(np.max(np.abs(sig)) / rms) if rms > 0 else float("inf")


def pair_lags(compsA, compsB, wins_idx, sh, maxlag):
    """B shifted left by sh samples, then fine lag per window. Returns window ->
    (total lag in samples = k + sh, cc, edge)."""
    out = {}
    shifted = {c: shift_left(x, sh) for c, x in compsB.items()}
    for w, (c, i1, i2) in wins_idx.items():
        k, cc, edge = window_lag(compsA[c], shifted[c], i1, i2, maxlag)
        out[w] = (k + sh if _fin(k) else float("nan"), cc, edge)
    return out


def loc_pred(dxN, dxE, az_deg, dp_skm):
    a = math.radians(az_deg)
    return -dp_skm * (dxN * math.cos(a) + dxE * math.sin(a))


def station_x(lags, ccs, edges, dt0, locp, P=PARAMS):
    """As sksnbr26.station_x: x = SKS_R lag - reference lag - loc_pred difference."""
    out = {}
    for sks in ("SKS_R", "SKS_Z"):
        for ref in ("S_T", "P_Z"):
            out["%s_vs_%s" % (sks, ref)] = lags[sks] - lags[ref] - (locp[sks] - locp[ref])
        out["%s_abs" % sks] = lags[sks] - dt0 - locp[sks]

    def ok(w):
        return _fin(ccs.get(w)) and ccs[w] >= P["ref_cc_min"] and not edges.get(w, True)
    ref = "S_T" if ok("S_T") else ("P_Z" if ok("P_Z") else "abs")
    if ref == "abs":
        out["x"], out["xz"] = out["SKS_R_abs"], out["SKS_Z_abs"]
    else:
        out["x"], out["xz"] = out["SKS_R_vs_" + ref], out["SKS_Z_vs_" + ref]
    out["ref"] = ref
    return out


def grade_pair(x, xz, usable, P=PARAMS):
    if not usable or not _fin(x):
        return "unusable"
    if abs(x) <= P["stable_s"]:
        return "stable"
    if (abs(x) >= P["changed_s"] and _fin(xz) and xz * x > 0
            and abs(xz) >= P["z_agree_s"]):
        return "changed"
    return "unclear"


def grade_class(gates, group, controls, P=PARAMS):
    """group / controls: lists of (pair grade, x). Pre-registered class grade."""
    ctrl = [x for g, x in controls if g != "unusable" and _fin(x)]
    n_ctrl = len(ctrl)
    f = (sum(abs(x) >= P["changed_s"] for x in ctrl) / n_ctrl) if n_ctrl else float("nan")
    quiet = bool(n_ctrl >= P["ctrl_min_n"] and f <= P["ctrl_false_max"])
    use = [(g, x) for g, x in group if g != "unusable" and _fin(x)]
    n = len(use)
    changed = [x for g, x in use if g == "changed"]
    n_stable = sum(g == "stable" for g, _ in use)
    res = dict(n_usable=n, n_changed=len(changed), n_stable=n_stable, n_ctrl=n_ctrl,
               ctrl_false_rate=f, controls_quiet=quiet,
               median_x=float(np.median([x for _, x in use])) if n else float("nan"))
    same_sign = bool(changed) and (all(x > 0 for x in changed) or all(x < 0 for x in changed))
    if not gates or not all(bool(v) for v in gates.values()) or n == 0:
        res["grade"] = "untestable"
    elif len(changed) >= max(1, math.ceil(n / 2.0)) and same_sign and quiet:
        res["grade"] = "change"
    elif not changed and n_stable >= P["stable_frac"] * n - 1e-9:
        res["grade"] = "stable"
    else:
        res["grade"] = "unclear"
    return res


# ----------------------------------------------------------------- helpers

def sha256(data):
    return hashlib.sha256(data).hexdigest()


def fetch(url, tries=4, timeout=120):
    err = None
    for k in range(tries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "inner-core-repeaters (MV2, Munim, Inc.)"})
            with urllib.request.urlopen(req, timeout=timeout) as r:
                return r.read()
        except Exception as ex:  # noqa: BLE001
            err = ex
            time.sleep(15 * (k + 1))
    raise RuntimeError("fetch failed %s: %r" % (url, err))


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


def wcsv(path, rows):
    if not rows:
        return
    keys = []
    for r in rows:
        for k in r:
            if k not in keys:
                keys.append(k)
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        for r in rows:
            w.writerow(r)


def load_sdv23():
    body = fetch(REPO + SDV23_COMMIT + "/sdv23.py")
    if SDV23_MUST not in body:
        raise RuntimeError("unexpected sdv23.py at " + SDV23_COMMIT)
    d = tempfile.mkdtemp(prefix="sdv23code_")
    with open(os.path.join(d, "sdv23.py"), "wb") as f:
        f.write(body)
    if d not in sys.path:
        sys.path.insert(0, d)
    import sdv23  # noqa: E402
    return sdv23, sha256(body)


def catalogue(data_dir, sdv, t1, t2, P=PARAMS):
    """USGS CSV ring query around IU.SDV, cached in data_dir."""
    os.makedirs(os.path.join(data_dir, "sdvexit27"), exist_ok=True)
    q = dict(format="csv", starttime=t1, endtime=t2, latitude=sdv["lat"],
             longitude=sdv["lon"], minradius=P["ring"][0], maxradius=P["ring"][1],
             minmagnitude=P["min_mag"], maxdepth=P["max_depth"], orderby="time-asc")
    f = os.path.join(data_dir, "sdvexit27", "cat_%s_%s.csv" % (t1[:10], t2[:10]))
    if os.path.exists(f):
        with open(f, "rb") as fh:
            body = fh.read()
    else:
        body = fetch(USGS_CSV + "?" + urllib.parse.urlencode(q))
        with open(f, "wb") as fh:
            fh.write(body)
    return parse_usgs_csv(body.decode("utf-8", "replace")), sha256(body)


def to_event(S, e):
    depth = min(max(float(e["depth"]), 0.0), 700.0)
    return S.ev_dict((e["id"], e["iso"], e["lat"], e["lon"], depth, e["mag"], e["magtype"]))


# ----------------------------------------------------------------- analysis

def geometry(S, model, cmb, M, lat, lon):
    """As sksnbr26.geometry: distance, azimuths, phase times/slownesses, SKS entry/exit."""
    from obspy.geodetics import gps2dist_azimuth, locations2degrees
    dist = float(locations2degrees(M["lat"], M["lon"], lat, lon))
    dkm, az, baz = gps2dist_azimuth(M["lat"], M["lon"], lat, lon)
    g = dict(dist=dist, dist_km=dkm / 1000.0, az=float(az), baz=float(baz),
             tP=float("nan"), tS=float("nan"), tSKS=float("nan"), tSKKS=float("nan"),
             pP=float("nan"), pS=float("nan"), pSKS=float("nan"), entry=None, exit=None)
    ph = S.phase_info(model, M["depth"], dist)
    for k in ("P", "S", "SKS"):
        if k in ph:
            g["t" + k], g["p" + k] = ph[k]
    sk = model.get_travel_times(source_depth_in_km=M["depth"], distance_in_degree=dist,
                                phase_list=["SKKS"])
    if sk:
        g["tSKKS"] = float(sk[0].time)
    pp = model.get_pierce_points(source_depth_in_km=M["depth"], distance_in_degree=dist,
                                 phase_list=["SKS"])
    if pp:
        pc = pp[0].pierce
        dd = np.asarray(pc["dist"], float)
        conv = math.degrees if abs(math.degrees(float(dd[-1])) - dist) < 1.0 else (lambda v: v)
        ie, ix = cmb_indices(pc["depth"], cmb)
        if ie is not None:
            g["entry"] = dest_point(M["lat"], M["lon"], az, conv(float(dd[ie])))
            g["exit"] = dest_point(M["lat"], M["lon"], az, conv(float(dd[ix])))
    return g


def sep_ok(g, P=PARAMS):
    return bool(_fin(g["tS"]) and _fin(g["tSKS"])
                and abs(g["tS"] - g["tSKS"]) >= P["min_s_sks_sep_s"]
                and (not _fin(g["tSKKS"]) or g["tSKKS"] - g["tSKS"] >= P["min_skks_sep_s"]))


def measure_sdv(S, D, SP, g, A, B, reloc, P=PARAMS):
    """SKS change at IU.SDV for pair (A earlier, B later), as the SDV branch of
    sksnbr26.run, in two frames: 'meta' (each record's own metadata) and 'nom' (both in
    the pre-step nominal frame 0/90 deg). Returns a row dict."""
    tbS, teS = g["tP"] - 300.0, g["dist_km"] / 2.6 + 150.0
    body, locs = {}, {}
    for E in (A, B):
        st, inv, loc = S.sdv_event(D, E, tbS, teS)
        if st is None:
            raise RuntimeError("nodata " + E["id"])
        b = S.process(st, inv, E["t"], g["tP"] - 150.0, g["tS"] + 120.0, "DISP",
                      SP["pre_filt_body"], SP["band_body"], SP["sr_up"])
        if not all(c in b for c in ("1", "2", "Z")):
            raise RuntimeError("missing components " + E["id"])
        if not all(_fin(b[c]["az"]) for c in ("1", "2")):
            raise RuntimeError("no orientation " + E["id"])
        body[E["id"]], locs[E["id"]] = b, loc
    a, b = body[A["id"]], body[B["id"]]
    sr = a["Z"]["sr"]
    if abs(b["Z"]["sr"] - sr) > 1e-6:
        raise RuntimeError("sampling rates differ")
    tbb = g["tP"] - 150.0
    offA = float(np.mean([a[c]["off"] for c in ("1", "2", "Z")]))
    offB = float(np.mean([b[c]["off"] for c in ("1", "2", "Z")]))
    offZ = a["Z"]["off"]
    idx = lambda t: int(round((t - tbb - offA) * sr))  # noqa: E731
    idxz = lambda t: int(round((t - tbb - offZ) * sr))  # noqa: E731
    tt = lambda i: tbb + offZ + i / sr  # noqa: E731
    rots = {"meta": (a["1"]["az"], a["2"]["az"], b["1"]["az"], b["2"]["az"]),
            "nom": (NOMINAL[0], NOMINAL[1], NOMINAL[0], NOMINAL[1])}
    RA, TA = S.to_rt(a["1"]["x"], a["2"]["x"], rots["meta"][0], rots["meta"][1], g["baz"])

    def peak(x, tc, lo=-P["peak_search_s"], hi=P["peak_search_s"]):
        i1, i2 = max(0, idxz(tc + lo)), min(len(x), idxz(tc + hi))
        return tt(i1 + int(np.argmax(np.abs(x[i1:i2]))))
    kR, kZ = peak(RA, g["tSKS"]), peak(a["Z"]["x"], g["tSKS"])
    kT, kP = peak(TA, g["tS"], hi=8.0), peak(a["Z"]["x"], g["tP"])
    wins = {"SKS_R": ("R", kR - 2.0, kR + 3.0), "SKS_Z": ("Z", kZ - 2.0, kZ + 3.0),
            "S_T": ("T", kT - 3.0, kT + 10.0), "P_Z": ("Z", kP - 1.0, kP + 4.0)}
    wi = {w: (c, idx(t1), idx(t2)) for w, (c, t1, t2) in wins.items()}
    sh = int(round((reloc["dt0"] - (offB - offA)) * sr))
    slow = {"SKS_R": g["pSKS"], "SKS_Z": g["pSKS"], "S_T": g["pS"], "P_Z": g["pP"]}
    locp = {w: loc_pred(reloc["dxN"], reloc["dxE"], g["az"], slow[w]) for w in WINDOWS}
    row = dict(loc_a=locs[A["id"]], loc_b=locs[B["id"]], az1_a=a["1"]["az"],
               az1_b=b["1"]["az"], coarse_s=sh / sr)
    for rn, (a1, a2, b1, b2) in rots.items():
        RAr, TAr = S.to_rt(a["1"]["x"], a["2"]["x"], a1, a2, g["baz"])
        RBr, TBr = S.to_rt(b["1"]["x"], b["2"]["x"], b1, b2, g["baz"])
        cA = {"R": RAr, "T": TAr, "Z": a["Z"]["x"]}
        cB = {"R": RBr, "T": TBr, "Z": b["Z"]["x"]}
        res = pair_lags(cA, cB, wi, sh, int(P["maxlag_s"] * sr))
        lags = {w: (v[0] / sr + offB - offA) if _fin(v[0]) else float("nan")
                for w, v in res.items()}
        ccs = {w: v[1] for w, v in res.items()}
        edges = {w: bool(v[2]) for w, v in res.items()}
        xs = station_x(lags, ccs, edges, reloc["dt0"], locp)
        n1, n2 = idx(g["tP"] - 120.0), idx(g["tP"] - 20.0)
        _, s1, s2 = wi["SKS_R"]
        snrA = snr(RAr, n1, n2, s1, s2)
        snrB = snr(shift_left(RBr, sh), n1, n2, s1, s2)
        usable = bool(_fin(ccs["SKS_R"]) and ccs["SKS_R"] >= P["sks_cc_min"]
                      and not edges["SKS_R"] and _fin(snrA) and _fin(snrB)
                      and snrA >= P["sks_snr_min"] and snrB >= P["sks_snr_min"]
                      and _fin(xs["x"]))
        for w in WINDOWS:
            row["%s_lag_%s" % (rn, w)] = lags[w]
            row["%s_cc_%s" % (rn, w)] = ccs[w]
            row["%s_edge_%s" % (rn, w)] = edges[w]
        for k2, v in xs.items():
            row["%s_%s" % (rn, k2)] = v
        row["%s_snr_a" % rn], row["%s_snr_b" % rn] = snrA, snrB
        row["%s_usable" % rn] = usable
        row["%s_grade" % rn] = grade_pair(xs["x"], xs["xz"], usable)
    for w in WINDOWS:
        row["t1_" + w] = round(wins[w][1], 3)
    return row


def run(args):
    P = PARAMS
    t_start = time.time()
    budget = args.time_budget_min * 60
    os.makedirs(args.out_dir, exist_ok=True)
    summary = dict(version=VERSION, started=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
                   params=P, errors=[])
    with open(os.path.abspath(__file__), "rb") as f:
        summary["code_sha256"] = sha256(f.read())
    S, s23 = load_sdv23()
    summary["sdv23"] = dict(commit=SDV23_COMMIT, sha256=s23, version=S.VERSION)
    ok, psha = S.check_plan(args.wait_min, PLAN_PATH)
    summary["plan_public_at_start"] = ok
    summary["plan_sha256"] = psha
    log("plan public:", ok, psha[:16])
    p23 = {}
    try:
        p23 = json.loads(fetch(REPO + "main/" + P23_SUMMARY).decode("utf-8"))
    except Exception as ex:  # noqa: BLE001
        summary["errors"].append("p23 summary: " + repr(ex)[:160])
    x23 = (p23.get("primary") or {}).get("x_meta")
    summary["p23_x_meta"] = x23

    from obspy.taup import TauPyModel
    SP = S.PARAMS
    SPm = dict(SP, p_max_stations=P["p_max_stations"])
    D = S.Data(args.data_dir)
    model = TauPyModel("iasp91")
    cmb = float(getattr(model.model, "cmb_depth", 2889.0))
    summary["cmb_depth_km"] = cmb
    EV = {r[0]: S.ev_dict(r) for r in S.EVENTS}
    M0, C0 = EV[S.MASTER], EV[S.PARTNER]
    sdv = S.SDV
    g0 = geometry(S, model, cmb, M0, sdv["lat"], sdv["lon"])
    if g0["entry"] is None:
        raise RuntimeError("no SKS CMB geometry for the doublet")
    summary["doublet_geometry"] = g0
    az0 = azimuth(sdv["lat"], sdv["lon"], M0["lat"], M0["lon"])
    summary["doublet_az_from_sdv"] = az0
    rows, mrows = [], []

    # ---- G0/G1: the doublet itself, sdv23 membership with its own full station set
    mres, mr = S.membership(D, model, M0, C0, SP, {})
    for r in mr:
        r["pair"] = "doublet"
    mrows += mr
    reloc0 = mres.get("reloc")
    summary["doublet_membership"] = mres
    drow = dict(pair="doublet", cls="doublet", a_iso=M0["iso"], b_iso=C0["iso"],
                d_entry_km=0.0, d_exit_km=0.0, status="")
    if reloc0:
        try:
            drow.update(measure_sdv(S, D, SP, g0, M0, C0, reloc0))
            drow["status"] = "ok"
        except Exception as ex:  # noqa: BLE001
            drow["status"] = "error: " + repr(ex)[:160]
    rows.append(drow)
    x_d = drow.get("meta_SKS_R_vs_S_T", float("nan"))
    gates = dict(
        G0_p23_reproduced=bool(_fin(x23) and _fin(x_d)
                               and abs(float(x_d) - float(x23)) <= P["g0_tol_s"]),
        G1_member_relocated=bool(mres.get("member") and reloc0),
    )
    log("gates", gates, "x_doublet_meta", x_d)

    # ---- catalogue and candidate pairs
    pre, post = [], []
    try:
        pre, h1 = catalogue(args.data_dir, sdv, P["cat_start"], P["pre_end"])
        p1, h2 = catalogue(args.data_dir, sdv, P["post_start"], P["post_split"])
        p2, h3 = catalogue(args.data_dir, sdv, P["post_split"], P["cat_end"])
        post = p1 + p2
        summary["catalogue"] = dict(n_pre=len(pre), n_post=len(post), sha256=[h1, h2, h3])
    except Exception as ex:  # noqa: BLE001
        summary["errors"].append("catalogue: " + repr(ex)[:200])
    for e in pre + post:
        e["dist_sdv"] = ang_dist(sdv["lat"], sdv["lon"], e["lat"], e["lon"])
        e["daz"] = wrap180(azimuth(sdv["lat"], sdv["lon"], e["lat"], e["lon"]) - az0)
    cands = find_pairs(pre, post)
    t_a0, t_b0 = iso_epoch(M0["iso"]), iso_epoch(C0["iso"])
    cands = [p for p in cands if not is_doublet(p, t_a0, t_b0)]
    wedge = [p for p in cands if abs(p["A"]["daz"]) <= P["wedge_daz"]
             and P["wedge_dist"][0] <= p["A"]["dist_sdv"] <= P["wedge_dist"][1]]
    ctrl_c = [p for p in cands if abs(p["A"]["daz"]) >= P["ctrl_daz"]
              and P["ctrl_dist"][0] <= p["A"]["dist_sdv"] <= P["ctrl_dist"][1]]
    summary["n_candidate_pairs"] = len(cands)
    summary["n_wedge_pairs"] = len(wedge)
    summary["n_ctrl_candidate_pairs"] = len(ctrl_c)
    log("pairs", len(cands), "wedge", len(wedge), "ctrl", len(ctrl_c))

    # ---- geometry (event A to SDV) and classes
    geo = {}
    order_a = []
    for p in sorted(wedge, key=rank_key):
        if p["A"]["id"] not in order_a:
            order_a.append(p["A"]["id"])
    order_a = order_a[:P["max_geom_events"]]
    ctrl_rank = greedy_unique(ctrl_c)[:P["max_ctrl_cands"]]
    ctrl_ids = {id(p) for p in ctrl_rank}
    need = set(order_a) | {p["A"]["id"] for p in ctrl_rank}
    evA = {p["A"]["id"]: p["A"] for p in wedge + ctrl_rank}
    for eid in sorted(need):
        try:
            geo[eid] = geometry(S, model, cmb, to_event(S, evA[eid]), sdv["lat"], sdv["lon"])
        except Exception as ex:  # noqa: BLE001
            summary["errors"].append("geometry %s: %s" % (eid, repr(ex)[:100]))
    cand_rows, by_cls = [], {"near": [], "exit": [], "control": []}
    for p in wedge + ctrl_rank:
        g = geo.get(p["A"]["id"])
        if g is None or g["entry"] is None:
            continue
        de = arc_km(g["entry"][0], g["entry"][1], g0["entry"][0], g0["entry"][1])
        dx = arc_km(g["exit"][0], g["exit"][1], g0["exit"][0], g0["exit"][1])
        c = classify_pair(de, dx)
        is_ctrl_cand = id(p) in ctrl_ids
        if is_ctrl_cand and c != "control":
            c = "other"
        if not is_ctrl_cand and c == "control":
            c = "other"
        p.update(d_entry_km=de, d_exit_km=dx, cls=c, sep_ok=sep_ok(g))
        cand_rows.append(dict(a_id=p["A"]["id"], a_iso=p["A"]["iso"], b_id=p["B"]["id"],
                              b_iso=p["B"]["iso"], a_mag=p["A"]["mag"], b_mag=p["B"]["mag"],
                              a_lat=p["A"]["lat"], a_lon=p["A"]["lon"],
                              a_depth=p["A"]["depth"], sep_km=round(p["sep_km"], 2),
                              span_yr=round(p["span_yr"], 2),
                              dist_sdv=round(p["A"]["dist_sdv"], 3),
                              daz=round(p["A"]["daz"], 3), d_entry_km=round(de, 1),
                              d_exit_km=round(dx, 1), cls=c, sep_ok=p["sep_ok"]))
        if c in by_cls and p["sep_ok"]:
            by_cls[c].append(p)
    for c in by_cls:
        by_cls[c] = greedy_unique(by_cls[c])
    summary["n_class_candidates"] = {c: len(v) for c, v in by_cls.items()}
    log("class candidates", summary["n_class_candidates"])

    # ---- membership + measurement, class by class (exit first: primary)
    checks = {}
    for c in ("exit", "control", "near"):
        max_chk, max_ver = P["caps"][c]
        n_chk = n_ver = 0
        for p in by_cls[c]:
            if n_chk >= max_chk or n_ver >= max_ver:
                break
            if time.time() - t_start > budget:
                summary["errors"].append("time budget reached in class " + c)
                break
            n_chk += 1
            A, B = to_event(S, p["A"]), to_event(S, p["B"])
            tag = "%s|%s" % (A["id"], B["id"])
            row = dict(pair=tag, cls=c, a_iso=A["iso"], b_iso=B["iso"], a_mag=A["mag"],
                       b_mag=B["mag"], sep_km=round(p["sep_km"], 2),
                       span_yr=round(p["span_yr"], 2), d_entry_km=p["d_entry_km"],
                       d_exit_km=p["d_exit_km"], dist_sdv=p["A"]["dist_sdv"], status="")
            try:
                mres_p, mr_p = S.membership(D, model, A, B, SPm, {})
                for r in mr_p:
                    r["pair"] = tag
                mrows += mr_p
                row.update(member=mres_p.get("member"), n_cc090=mres_p.get("n_cc090"),
                           n_cc095=mres_p.get("n_cc095"))
                rl = mres_p.get("reloc")
                if not (mres_p.get("member") and rl):
                    row["status"] = "not a verified member"
                else:
                    n_ver += 1
                    row.update(dt0=rl["dt0"], dxN=rl["dxN"], dxE=rl["dxE"], reloc_rms=rl["rms"])
                    row.update(measure_sdv(S, D, SP, geo[p["A"]["id"]], A, B, rl))
                    row["status"] = "ok"
                    log(c, tag, "x_nom %.4f grade %s" % (row["nom_x"], row["nom_grade"]))
            except Exception as ex:  # noqa: BLE001
                row["status"] = "error: " + repr(ex)[:160]
                log(c, tag, row["status"])
            rows.append(row)
        checks[c] = dict(checked=n_chk, verified=n_ver)
    summary["checks"] = checks

    # ---- grades (primary frame: nom)
    def grp(c, fr="nom"):
        return [(r.get(fr + "_grade", "unusable"), r.get(fr + "_x")) for r in rows
                if r["cls"] == c and r["status"] == "ok"]
    ctrl = grp("control")
    prim = grade_class(gates, grp("exit"), ctrl)
    prim.update(gates=gates, x_doublet_meta=x_d, x_doublet_nom=drow.get("nom_x"), x_p23=x23)
    sec = grade_class(gates, grp("near"), ctrl)
    rep_meta = grade_class(gates, grp("exit", "meta"), grp("control", "meta"))
    allx = [r.get("nom_x") for r in rows if r["cls"] != "doublet" and r["status"] == "ok"
            and r.get("nom_usable")]
    summary["primary_exit"] = prim
    summary["secondary_near"] = sec
    summary["reported_exit_meta_frame"] = rep_meta
    summary["station_side"] = dict(n=len(allx), median_x=float(np.median(allx)) if allx
                                   else float("nan"))
    log("PRIMARY", prim)

    wcsv(os.path.join(args.out_dir, "sdvexit27_candidates.csv"), cand_rows)
    wcsv(os.path.join(args.out_dir, "sdvexit27_rows.csv"), rows)
    wcsv(os.path.join(args.out_dir, "sdvexit27_members.csv"), mrows)
    summary["finished"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    summary["runtime_min"] = round((time.time() - t_start) / 60.0, 1)
    with open(os.path.join(args.out_dir, "sdvexit27_summary.json"), "w") as f:
        json.dump(summary, f, indent=1, default=jsonable)
    try:
        plot(args.out_dir, rows)
    except Exception as ex:  # noqa: BLE001
        log("plot failed", repr(ex)[:200])
    return summary


def plot(out_dir, rows):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    os.makedirs(os.path.join(out_dir, "plots"), exist_ok=True)
    fig, ax = plt.subplots(figsize=(7, 4.5))
    col = {"doublet": "k", "near": "C1", "exit": "C0", "control": "0.5"}
    for r in rows:
        if r.get("status") != "ok" or not r.get("nom_usable"):
            continue
        ax.scatter(max(r["d_exit_km"], 1.0), 1000 * r["nom_x"], c=col.get(r["cls"], "C3"),
                   label=r["cls"])
    h, lab = ax.get_legend_handles_labels()
    uniq = dict(zip(lab, h))
    ax.legend(uniq.values(), uniq.keys(), fontsize=8)
    ax.axhline(-160, color="C3", ls="--", lw=0.8)
    ax.axhline(0, color="k", lw=0.5)
    ax.set_xscale("log")
    ax.set_xlabel("distance of SKS exit point from the doublet's (km)")
    ax.set_ylabel("x = SKS change at IU.SDV, later - earlier (ms)")
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, "plots", "sdvexit27.png"), dpi=110)


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", default="/data")
    ap.add_argument("--out-dir", default="/out")
    ap.add_argument("--wait-min", type=float, default=120.0)
    ap.add_argument("--time-budget-min", type=float, default=240.0)
    a = ap.parse_args(argv)
    run(a)
    return 0


if __name__ == "__main__":
    sys.exit(main())

#!/usr/bin/env python3
# Copyright 2026 MV2, Munim, Inc.
# Licensed under the Apache License, Version 2.0 (Apache-2.0). See LICENSE.
# Made with AI agents working in MV2's name.
"""sksnbr26.py - pre-registration 26: do other SKS paths that enter the core near IU.SDV's
SKS entry point change between the Aleutian doublet 2000-07-07 / 2009-12-17?

Plan: findings/2026-11-21-prereg26-sks-entry-neighbours-plan.md (published before this code
runs; the job waits for it and records its sha256).

Zhang & Wen 2025 (Nat. Commun. 16:8754) model a decadal ULVZ change near the SKS entry point
(source side) of the Aleutian -> IU.SDV path and say the exit point fits equally well. Stations
at nearly the same azimuth from the source (e.g. VE.FUNV, IU.PTGA) enter the core within tens of
km of SDV's entry point but exit far from SDV's exit point. This job measures the same SKS
differential-time change at every station that recorded both events at 85-125 deg, classifies
stations by CMB entry/exit distance from SDV's, and grades the source-side hypothesis.

Processing, relocation and membership are those of sdv23.py 0.1.1, imported unchanged from its
pinned commit (fetched at run time, sha256 recorded). Top-level imports: stdlib + numpy only.
"""

import argparse
import csv
import hashlib
import json
import math
import os
import sys
import tempfile
import time
import urllib.request

import numpy as np

VERSION = "0.1.0"
PLAN_PATH = "findings/2026-11-21-prereg26-sks-entry-neighbours-plan.md"
REPO = "https://raw.githubusercontent.com/munimtechnologies/inner-core-repeaters/"
SDV23_COMMIT = "dbbd91d3f26f9dc2165cbf59a5c05f7d82d8521d"
SDV23_MUST = b'VERSION = "0.1.1"'
P23_SUMMARY = "results/sdv23-v2/sdv23_summary.json"
R_CMB = 3480.0
CONTROL_NETS = ("IU", "II", "GE", "G", "GT", "MN")
WINDOWS = ("SKS_R", "SKS_Z", "S_T", "P_Z")

PARAMS = dict(
    dist_range=(85.0, 125.0),
    entry_nb_km=100.0, exit_nb_km=100.0, control_min_km=300.0, report_km=300.0,
    max_controls=40,
    min_s_sks_sep_s=10.0, min_skks_sep_s=6.0,
    peak_search_s=5.0,
    maxlag_s=0.6, sks_cc_min=0.80, sks_snr_min=5.0, ref_cc_min=0.90,
    changed_s=0.08, stable_s=0.04, z_agree_s=0.04,
    g0_tol_s=0.010, ctrl_min_n=5, ctrl_false_max=0.20,
)


def log(*a):
    print(time.strftime("%H:%M:%S"), "[sksnbr26]", *a, flush=True)


# ----------------------------------------------------------------- pure numerics

def _fin(x):
    try:
        return x is not None and math.isfinite(float(x))
    except (TypeError, ValueError):
        return False


def dest_point(lat, lon, az, ang):
    """Point at angular distance ang (deg) from (lat, lon) along initial azimuth az (deg) on a
    sphere. Returns (lat, lon) in degrees, lon in [-180, 180)."""
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


def arc_km(lat1, lon1, lat2, lon2, r=R_CMB):
    return math.radians(ang_dist(lat1, lon1, lat2, lon2)) * r


def cmb_indices(depths, cmb, tol=0.5):
    """First and last pierce-point index at the CMB depth (SKS entry and exit)."""
    idx = [i for i, d in enumerate(depths) if abs(float(d) - cmb) <= tol]
    if not idx:
        return None, None
    return idx[0], idx[-1]


def select_loc(chans_by_loc):
    """chans_by_loc: location code -> set of channel codes. Preference 00, '', 10, others;
    the location needs a Z channel and a horizontal pair (N/E or 1/2)."""
    def key(l):
        return ({"00": 0, "": 1, "10": 2}.get(l, 3), l)
    for loc in sorted(chans_by_loc, key=key):
        ends = {c[-1] for c in chans_by_loc[loc]}
        if "Z" in ends and ({"N", "E"} <= ends or {"1", "2"} <= ends):
            return loc
    return None


def shift_left(b, s):
    """c[i] = b[i + s] (zeros outside). Same as sdv23.shift_left."""
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
    """Same as sdv23.window_lag: normalised cross-correlation of a[i1:i2] with b[i1+k:i2+k],
    |k| <= maxlag samples. Returns (lag, cc, edge); lag > 0 means b is delayed."""
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
    """compsA/compsB: comp -> array; wins_idx: window -> (comp, i1, i2) on A. B is shifted
    left by sh samples, then the fine lag (|k| <= maxlag) is measured. Returns
    window -> (total lag in samples = k + sh, cc, edge)."""
    out = {}
    shifted = {c: shift_left(x, sh) for c, x in compsB.items()}
    for w, (c, i1, i2) in wins_idx.items():
        k, cc, edge = window_lag(compsA[c], shifted[c], i1, i2, maxlag)
        out[w] = (k + sh if _fin(k) else float("nan"), cc, edge)
    return out


def loc_pred(dxN, dxE, az_deg, dp_skm):
    """Same as sdv23.loc_pred: predicted change of a (differential) time for a source shift."""
    a = math.radians(az_deg)
    return -dp_skm * (dxN * math.cos(a) + dxE * math.sin(a))


def station_x(lags, ccs, edges, dt0, locp, P=PARAMS):
    """lags: window -> lag (s, 2009 minus 2000, relative to catalogue origins); ccs, edges:
    window -> cc / edge flag; dt0: relocation origin shift (s); locp: window -> loc_pred term
    for that window's slowness (s). Returns dict(ref, x, xz, and every version)."""
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


def classify(d_entry, d_exit, net, P=PARAMS):
    if _fin(d_entry) and d_entry <= P["entry_nb_km"]:
        return "entry-neighbour"
    if _fin(d_exit) and d_exit <= P["exit_nb_km"]:
        return "exit-neighbour"
    if (_fin(d_entry) and _fin(d_exit) and d_entry > P["control_min_km"]
            and d_exit > P["control_min_km"] and net in CONTROL_NETS):
        return "control"
    return "other"


def grade_station(x, xz, usable, P=PARAMS):
    if not usable or not _fin(x):
        return "unusable"
    if abs(x) <= P["stable_s"]:
        return "stable"
    if (abs(x) >= P["changed_s"] and _fin(xz) and xz * x > 0
            and abs(xz) >= P["z_agree_s"]):
        return "changed"
    return "unclear"


def grade_group(gates, group, controls, P=PARAMS):
    """group / controls: lists of (station grade, x). Returns dict(grade, n_usable, n_ctrl,
    ctrl_false_rate, controls_quiet)."""
    ctrl = [x for g, x in controls if g != "unusable" and _fin(x)]
    n_ctrl = len(ctrl)
    f = (sum(abs(x) >= P["changed_s"] for x in ctrl) / n_ctrl) if n_ctrl else float("nan")
    quiet = bool(n_ctrl >= P["ctrl_min_n"] and f <= P["ctrl_false_max"])
    use = [g for g, _ in group if g != "unusable"]
    res = dict(n_usable=len(use), n_ctrl=n_ctrl, ctrl_false_rate=f, controls_quiet=quiet,
               grades=use)
    if not gates or not all(bool(v) for v in gates.values()) or not use:
        res["grade"] = "untestable"
    elif any(g == "changed" for g in use) and quiet:
        res["grade"] = "shared-change"
    elif all(g == "stable" for g in use):
        res["grade"] = "no-shared-change"
    else:
        res["grade"] = "unclear"
    return res


# ----------------------------------------------------------------- helpers

def sha256(data):
    return hashlib.sha256(data).hexdigest()


def fetch(url, tries=4, timeout=90):
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


# ----------------------------------------------------------------- analysis

def geometry(S, model, cmb, M, lat, lon):
    """Distance, azimuths, phase times/slownesses and SKS CMB entry/exit points."""
    from obspy.geodetics import gps2dist_azimuth, locations2degrees
    dist = float(locations2degrees(M["lat"], M["lon"], lat, lon))
    _, az, baz = gps2dist_azimuth(M["lat"], M["lon"], lat, lon)
    g = dict(dist=dist, az=float(az), baz=float(baz), tP=float("nan"), tS=float("nan"),
             tSKS=float("nan"), tSKKS=float("nan"), pP=float("nan"), pS=float("nan"),
             pSKS=float("nan"), entry=None, exit=None)
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
            g["entry_ang"] = conv(float(dd[ie]))
            g["exit_ang"] = conv(float(dd[ix]))
    return g


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
    D = S.Data(args.data_dir)
    model = TauPyModel("iasp91")
    cmb = float(getattr(model.model, "cmb_depth", 2889.0))
    summary["cmb_depth_km"] = cmb
    EV = {r[0]: S.ev_dict(r) for r in S.EVENTS}
    M, C = EV[S.MASTER], EV[S.PARTNER]

    # ---- membership + relocation (sdv23 rules, unchanged)
    mres, mrows = S.membership(D, model, M, C, SP, {})
    reloc = mres.get("reloc")
    summary["membership"] = mres
    log("membership", mres.get("n_cc090"), mres.get("n_cc095"), reloc)

    # ---- candidate stations: BHZ at both event times, 85-125 deg
    lists = {}
    for E in (M, C):
        inv = D.stations("nb26_list_" + E["id"], latitude=M["lat"], longitude=M["lon"],
                         minradius=P["dist_range"][0], maxradius=P["dist_range"][1],
                         channel="BHZ", starttime=E["t"], endtime=E["t"] + 3600,
                         level="station")
        d = {}
        if inv is not None:
            for net in inv:
                for s in net:
                    d[(net.code, s.code)] = (float(s.latitude), float(s.longitude))
        lists[E["id"]] = d
    keys = sorted(set(lists[M["id"]]) & set(lists[C["id"]]))
    coords = {k: lists[M["id"]][k] for k in keys}
    coords[("IU", "SDV")] = (S.SDV["lat"], S.SDV["lon"])
    summary["n_candidates"] = len(coords)

    geo = {}
    for k, (lat, lon) in coords.items():
        try:
            geo[k] = geometry(S, model, cmb, M, lat, lon)
        except Exception as ex:  # noqa: BLE001
            summary["errors"].append("geometry %s.%s: %s" % (k[0], k[1], repr(ex)[:120]))
    gs = geo.get(("IU", "SDV"))
    if not gs or gs["entry"] is None:
        raise RuntimeError("no SKS CMB geometry for IU.SDV")
    summary["sdv_geometry"] = gs
    sta_rows = []
    for k, g in geo.items():
        de = arc_km(g["entry"][0], g["entry"][1], gs["entry"][0], gs["entry"][1]) \
            if g["entry"] else float("nan")
        dx = arc_km(g["exit"][0], g["exit"][1], gs["exit"][0], gs["exit"][1]) \
            if g["exit"] else float("nan")
        sep_ok = (_fin(g["tS"]) and _fin(g["tSKS"]) and abs(g["tS"] - g["tSKS"])
                  >= P["min_s_sks_sep_s"] and (not _fin(g["tSKKS"])
                                               or g["tSKKS"] - g["tSKS"] >= P["min_skks_sep_s"]))
        cls = "sdv" if k == ("IU", "SDV") else classify(de, dx, k[0])
        sta_rows.append(dict(net=k[0], sta=k[1], lat=coords[k][0], lon=coords[k][1],
                             dist=round(g["dist"], 3), az=round(g["az"], 3),
                             baz=round(g["baz"], 3), d_entry_km=de, d_exit_km=dx,
                             entry_lat=g["entry"][0] if g["entry"] else float("nan"),
                             entry_lon=g["entry"][1] if g["entry"] else float("nan"),
                             exit_lat=g["exit"][0] if g["exit"] else float("nan"),
                             exit_lon=g["exit"][1] if g["exit"] else float("nan"),
                             sep_ok=bool(sep_ok), cls=cls))
    order = {"sdv": 0, "entry-neighbour": 1, "exit-neighbour": 2, "other": 3, "control": 4}
    sta_rows.sort(key=lambda r: (order[r["cls"]], r["net"], r["sta"]))
    todo, n_ctrl = [], 0
    for r in sta_rows:
        if not r["sep_ok"] and r["cls"] != "sdv":
            r["plan"] = "skip: phase separation"
            continue
        if r["cls"] == "other" and not (min(r["d_entry_km"], r["d_exit_km"]) <= P["report_km"]):
            r["plan"] = "skip: other"
            continue
        if r["cls"] == "control":
            if n_ctrl >= P["max_controls"]:
                r["plan"] = "skip: control cap"
                continue
            n_ctrl += 1
        r["plan"] = "measure"
        todo.append(r)
    log("candidates", len(sta_rows), "to measure", len(todo),
        {c: sum(r["cls"] == c for r in todo) for c in order})

    rows = []
    if not reloc:
        summary["errors"].append("no relocation: E2009a not relocated; nothing measured")
    tbS, teS = gs["tP"] - 300.0, None
    for r in (todo if reloc else []):
        if time.time() - t_start > budget:
            summary["errors"].append("time budget reached before %s.%s" % (r["net"], r["sta"]))
            break
        k = (r["net"], r["sta"])
        g = geo[k]
        row = dict(net=r["net"], sta=r["sta"], cls=r["cls"], dist=r["dist"], az=r["az"],
                   d_entry_km=r["d_entry_km"], d_exit_km=r["d_exit_km"], status="")
        try:
            body, locs = {}, {}
            for E in (M, C):
                if k == ("IU", "SDV"):
                    from obspy.geodetics import gps2dist_azimuth
                    dist_km = gps2dist_azimuth(M["lat"], M["lon"], S.SDV["lat"],
                                               S.SDV["lon"])[0] / 1000.0
                    teS = dist_km / 2.6 + 150.0
                    st, inv, loc = S.sdv_event(D, E, tbS, teS)
                    if st is None:
                        raise RuntimeError("nodata " + E["id"])
                else:
                    t1, t2 = E["t"] + g["tP"] - 300.0, E["t"] + g["tS"] + 300.0
                    st = D.waveforms(r["net"], r["sta"], "*", "BH?", t1, t2)
                    if st is None:
                        raise RuntimeError("nodata " + E["id"])
                    cbl = {}
                    for tr in st:
                        cbl.setdefault(tr.stats.location, set()).add(tr.stats.channel)
                    loc = select_loc(cbl)
                    if loc is None:
                        raise RuntimeError("no 3-component location " + E["id"])
                    st = st.select(location=loc)
                    inv = D.stations("nb26_%s_%s_%s" % (r["net"], r["sta"], E["id"]),
                                     network=r["net"], station=r["sta"], location=loc or "--",
                                     channel="BH?", starttime=t1, endtime=t2, level="response")
                b = S.process(st, inv, E["t"], g["tP"] - 150.0, g["tS"] + 120.0, "DISP",
                              SP["pre_filt_body"], SP["band_body"], SP["sr_up"])
                if not all(c in b for c in ("1", "2", "Z")):
                    raise RuntimeError("missing components " + E["id"])
                if not all(_fin(b[c]["az"]) for c in ("1", "2")):
                    raise RuntimeError("no orientation " + E["id"])
                body[E["id"]], locs[E["id"]] = b, loc
            A, B = body[M["id"]], body[C["id"]]
            sr = A["Z"]["sr"]
            if abs(B["Z"]["sr"] - sr) > 1e-6:
                raise RuntimeError("sampling rates differ")
            tbb = g["tP"] - 150.0
            offA = float(np.mean([A[c]["off"] for c in ("1", "2", "Z")]))
            offB = float(np.mean([B[c]["off"] for c in ("1", "2", "Z")]))
            offZ = A["Z"]["off"]
            idx = lambda t: int(round((t - tbb - offA) * sr))  # noqa: E731
            idxz = lambda t: int(round((t - tbb - offZ) * sr))  # noqa: E731
            tt = lambda i: tbb + offZ + i / sr  # noqa: E731
            RA, TA = S.to_rt(A["1"]["x"], A["2"]["x"], A["1"]["az"], A["2"]["az"], g["baz"])
            RB, TB = S.to_rt(B["1"]["x"], B["2"]["x"], B["1"]["az"], B["2"]["az"], g["baz"])
            cA = {"R": RA, "T": TA, "Z": A["Z"]["x"]}
            cB = {"R": RB, "T": TB, "Z": B["Z"]["x"]}

            def peak(x, tc, lo=-P["peak_search_s"], hi=P["peak_search_s"]):
                i1, i2 = max(0, idxz(tc + lo)), min(len(x), idxz(tc + hi))
                return tt(i1 + int(np.argmax(np.abs(x[i1:i2]))))
            kR, kZ = peak(RA, g["tSKS"]), peak(A["Z"]["x"], g["tSKS"])
            kT, kP = peak(TA, g["tS"], hi=8.0), peak(A["Z"]["x"], g["tP"])
            wins = {"SKS_R": ("R", kR - 2.0, kR + 3.0), "SKS_Z": ("Z", kZ - 2.0, kZ + 3.0),
                    "S_T": ("T", kT - 3.0, kT + 10.0), "P_Z": ("Z", kP - 1.0, kP + 4.0)}
            wi = {w: (c, idx(t1), idx(t2)) for w, (c, t1, t2) in wins.items()}
            sh = int(round((reloc["dt0"] - (offB - offA)) * sr))
            res = pair_lags(cA, cB, wi, sh, int(P["maxlag_s"] * sr))
            lags = {w: (v[0] / sr + offB - offA) if _fin(v[0]) else float("nan")
                    for w, v in res.items()}
            ccs = {w: v[1] for w, v in res.items()}
            edges = {w: bool(v[2]) for w, v in res.items()}
            slow = {"SKS_R": g["pSKS"], "SKS_Z": g["pSKS"], "S_T": g["pS"], "P_Z": g["pP"]}
            locp = {w: loc_pred(reloc["dxN"], reloc["dxE"], g["az"], slow[w]) for w in WINDOWS}
            xs = station_x(lags, ccs, edges, reloc["dt0"], locp)
            n1, n2 = idx(g["tP"] - 120.0), idx(g["tP"] - 20.0)
            _, s1, s2 = wi["SKS_R"]
            snrA = snr(RA, n1, n2, s1, s2)
            snrB = snr(shift_left(RB, sh), n1, n2, s1, s2)
            usable = bool(_fin(ccs["SKS_R"]) and ccs["SKS_R"] >= P["sks_cc_min"]
                          and not edges["SKS_R"] and _fin(snrA) and _fin(snrB)
                          and snrA >= P["sks_snr_min"] and snrB >= P["sks_snr_min"]
                          and _fin(xs["x"]))
            row.update(loc2000=locs[M["id"]], loc2009=locs[C["id"]],
                       az1_2000=A["1"]["az"], az1_2009=B["1"]["az"],
                       snr2000=snrA, snr2009=snrB, coarse_s=sh / sr)
            for w in WINDOWS:
                row["lag_" + w] = lags[w]
                row["cc_" + w] = ccs[w]
                row["edge_" + w] = edges[w]
                row["t1_" + w] = round(wins[w][1], 3)
            row.update({k2: v for k2, v in xs.items()})
            row["usable"] = usable
            row["grade"] = grade_station(xs["x"], xs["xz"], usable)
            row["status"] = "ok"
            log(r["net"], r["sta"], r["cls"], "x %.4f ref %s cc %.3f snr %.1f/%.1f -> %s"
                % (xs["x"], xs["ref"], ccs["SKS_R"], snrA, snrB, row["grade"]))
        except Exception as ex:  # noqa: BLE001
            row["status"] = "error: " + repr(ex)[:160]
            row["grade"] = "unusable"
            log(r["net"], r["sta"], row["status"])
        rows.append(row)

    # ---- grades
    sdv = next((x for x in rows if x["net"] == "IU" and x["sta"] == "SDV"), {})
    x_sdv_st = sdv.get("SKS_R_vs_S_T", float("nan"))
    gates = dict(
        G0_p23_reproduced=bool(_fin(x23) and _fin(x_sdv_st)
                               and abs(float(x_sdv_st) - float(x23)) <= P["g0_tol_s"]),
        G1_member_relocated=bool(mres.get("member") and reloc),
    )
    ctrl = [(x["grade"], x.get("x")) for x in rows if x["cls"] == "control"]
    ent = [(x["grade"], x.get("x")) for x in rows if x["cls"] == "entry-neighbour"]
    ext = [(x["grade"], x.get("x")) for x in rows if x["cls"] == "exit-neighbour"]
    prim = grade_group(gates, ent, ctrl)
    sec = grade_group(gates, ext, ctrl)
    prim.update(gates=gates, x_sdv_vs_ST=x_sdv_st, x_p23=x23,
                entry_neighbours=[{k2: x.get(k2) for k2 in ("net", "sta", "d_entry_km",
                                                             "d_exit_km", "x", "xz", "ref",
                                                             "grade", "status")}
                                  for x in rows if x["cls"] == "entry-neighbour"])
    summary["primary"] = prim
    summary["secondary_exit"] = sec
    summary["n_measured"] = len(rows)
    log("PRIMARY", prim)

    wcsv(os.path.join(args.out_dir, "sksnbr26_stations.csv"), sta_rows)
    wcsv(os.path.join(args.out_dir, "sksnbr26_rows.csv"), rows)
    wcsv(os.path.join(args.out_dir, "sksnbr26_members.csv"), mrows)
    summary["finished"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    summary["runtime_min"] = round((time.time() - t_start) / 60.0, 1)
    with open(os.path.join(args.out_dir, "sksnbr26_summary.json"), "w") as f:
        json.dump(summary, f, indent=1, default=jsonable)
    try:
        plot(args.out_dir, sta_rows, rows, gs)
    except Exception as ex:  # noqa: BLE001
        log("plot failed", repr(ex)[:200])
    return summary


def plot(out_dir, sta_rows, rows, gs):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    os.makedirs(os.path.join(out_dir, "plots"), exist_ok=True)
    fig, ax = plt.subplots(1, 2, figsize=(12, 5))
    xs = {(r["net"], r["sta"]): r for r in rows if r.get("usable")}
    for r in sta_rows:
        if not _fin(r["entry_lat"]):
            continue
        u = xs.get((r["net"], r["sta"]))
        col = "0.7" if u is None else ("C3" if u["x"] < 0 else "C0")
        size = 10 if u is None else 10 + 400 * min(abs(u["x"]), 0.3)
        ax[0].scatter(r["entry_lon"], r["entry_lat"], s=size, c=col)
        if u is not None or r["cls"] in ("entry-neighbour", "sdv"):
            ax[0].annotate("%s.%s" % (r["net"], r["sta"]), (r["entry_lon"], r["entry_lat"]),
                           fontsize=6)
    ax[0].plot(gs["entry"][1], gs["entry"][0], "k*", ms=12)
    ax[0].set_title("SKS CMB entry points (size ~ |x|, red = earlier in 2009)")
    ax[0].set_xlabel("lon")
    ax[0].set_ylabel("lat")
    u = [r for r in rows if r.get("usable")]
    ax[1].scatter([r["d_entry_km"] for r in u], [1000 * r["x"] for r in u], c="C0")
    for r in u:
        ax[1].annotate(r["sta"], (r["d_entry_km"], 1000 * r["x"]), fontsize=6)
    ax[1].axhline(-160, color="C3", ls="--", lw=0.8)
    ax[1].axhline(0, color="k", lw=0.5)
    ax[1].set_xscale("symlog", linthresh=100)
    ax[1].set_xlabel("distance of SKS entry point from SDV's (km)")
    ax[1].set_ylabel("x = SKS change 2009 - 2000 (ms)")
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, "plots", "sksnbr26.png"), dpi=110)


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", default="/data")
    ap.add_argument("--out-dir", default="/out")
    ap.add_argument("--wait-min", type=float, default=120.0)
    ap.add_argument("--time-budget-min", type=float, default=180.0)
    a = ap.parse_args(argv)
    run(a)
    return 0


if __name__ == "__main__":
    sys.exit(main())

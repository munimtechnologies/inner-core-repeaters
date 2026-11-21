#!/usr/bin/env python3
# Copyright 2026 MV2, Munim, Inc.
# Licensed under the Apache License, Version 2.0 (Apache-2.0). See LICENSE.
# Made with AI agents working in MV2's name.
"""sdvwidth25.py - pre-registration 25: is the reported 160 ms SKS advance at IU.SDV an
onset reading of a broader 2009 pulse, and if the 2009 pulse is broader, is that the
earthquake source (seen in every phase) or the path (SKS only)?

Plan: findings/2026-11-21-prereg25-sdv-onset-width-plan.md (published before this code
runs; the job waits for it and records its sha256).

Pre-registration 23 (sdv23.py 0.1.1) found that the whole-pulse cross-correlation change of
radial SKS minus transverse S at IU.SDV between the Aleutian doublet 2000-07-07 / 2009-12-17
is +28 ms, not the -160 ms reported by Zhang & Wen 2025 (Nat. Commun. 16:8754), who describe
the 2009 SKS as earlier and "broader". This test reads the same records in two further ways:
 1. onset times: where the envelope of each pulse first rises to 10/20/30% of its peak;
 2. pulse width: the stretch factor s that best maps the 2000 pulse onto the 2009 pulse
    (s > 1: 2009 broader), for radial and vertical SKS, the direct transverse S pulse and
    vertical P at SDV (same azimuth as SKS), and vertical P at the global stations that verify
    the doublet (azimuth-averaged source reference).
Processing, windows, rotations and relocation are those of sdv23.py 0.1.1, imported unchanged
from its pinned commit (fetched at run time, sha256 recorded). Waveforms are read from the
sdv23 data cache (same data_key). Top-level imports: stdlib + numpy only.
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
PLAN_PATH = "findings/2026-11-21-prereg25-sdv-onset-width-plan.md"
REPO = "https://raw.githubusercontent.com/munimtechnologies/inner-core-repeaters/"
SDV23_COMMIT = "dbbd91d3f26f9dc2165cbf59a5c05f7d82d8521d"
SDV23_MUST = b'VERSION = "0.1.1"'
P23_SUMMARY = "results/sdv23-v2/sdv23_summary.json"
PULSES = ("SKS_R", "SKS_Z", "S_T_pulse", "P_Z")

PARAMS = dict(
    onset_fracs=(0.1, 0.2, 0.3),   # envelope fractions of the peak; 0.2 is primary
    onset_back_s=3.0,              # scan at most this far back from the envelope peak
    s_min=0.70, s_max=1.40, s_step=0.005,
    stretch_maxlag_s=0.3,          # free lag while fitting the stretch (after coarse alignment)
    coarse_p_s=3.0,                # coarse alignment of global P pairs
    p23_tol_s=0.005,               # G0: our SKS_R x (meta) must match P23's within 5 ms
    sks_cc_min=0.80,
    onset_adv_s=0.10, onset_none_s=0.05,
    w_broad=1.05, w_narrow=0.95, w_same=0.03,
    ref_cc_min=0.90, glob_cc_min=0.90, glob_min_n=5,
    src_tol_ref=0.03, src_tol_glob=0.05, path_excess=0.05,
    nboot=2000, seed=25,
)


def log(*a):
    print(time.strftime("%H:%M:%S"), "[sdvwidth25]", *a, flush=True)


# ----------------------------------------------------------------- pure numerics

def hilbert(x):
    """Hilbert transform (imaginary part of the analytic signal): cos -> sin."""
    x = np.asarray(x, float)
    n = len(x)
    X = np.fft.fft(x)
    h = np.zeros(n)
    if n % 2 == 0:
        h[0] = h[n // 2] = 1.0
        h[1:n // 2] = 2.0
    else:
        h[0] = 1.0
        h[1:(n + 1) // 2] = 2.0
    return np.imag(np.fft.ifft(X * h))


def envelope(x):
    x = np.asarray(x, float)
    return np.sqrt(x ** 2 + hilbert(x) ** 2)


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


def peak_index(env, i1, i2):
    env = np.asarray(env, float)
    i1, i2 = max(0, int(i1)), min(len(env), int(i2))
    if i2 - i1 < 1:
        return -1
    return i1 + int(np.argmax(env[i1:i2]))


def onset_index(env, i1, i2, frac, back):
    """Onset of a pulse: take the envelope peak in [i1, i2), scan backwards at most `back`
    samples to the first sample below frac * peak, and interpolate linearly to the crossing.
    Returns (onset as a float sample index, peak index); onset is nan if not found."""
    env = np.asarray(env, float)
    ip = peak_index(env, i1, i2)
    if ip < 0:
        return float("nan"), ip
    thr = frac * env[ip]
    lo = max(0, ip - int(back))
    j = ip
    while j > lo and env[j] >= thr:
        j -= 1
    if env[j] >= thr:
        return float("nan"), ip
    e0, e1 = float(env[j]), float(env[j + 1])
    pos = (thr - e0) / (e1 - e0) if e1 != e0 else 0.0
    return float(j + pos), ip


def stretch_fit(a, b, i1, i2, ic, smin, smax, step, maxlag):
    """Stretch factor s such that b compressed about sample ic best matches a in [i1, i2):
    bs[i] = b(ic + (i - ic) * s), with a free lag |k| <= maxlag. s > 1 means b is broader
    than a. Returns dict(s, lag, cc, cc_s1, edge)."""
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    n = len(b)
    pos_all = np.arange(n, dtype=float)
    lo, hi = int(i1) - int(maxlag), int(i2) + int(maxlag)
    _, cc1, _ = window_lag(a, b, i1, i2, maxlag)
    best = dict(s=float("nan"), lag=float("nan"), cc=float("nan"), cc_s1=cc1, edge=True)
    if lo < 0 or hi > n:
        return best
    seg = np.arange(lo, hi, dtype=float)
    for s in np.arange(smin, smax + step / 2.0, step):
        bs = np.zeros(n)
        bs[lo:hi] = np.interp(ic + (seg - ic) * s, pos_all, b, left=0.0, right=0.0)
        k, cc, edge = window_lag(a, bs, i1, i2, maxlag)
        if np.isfinite(cc) and (not np.isfinite(best["cc"]) or cc > best["cc"]):
            best.update(s=float(round(s, 6)), lag=k, cc=cc, edge=edge)
    return best


def boot_median(vals, nboot, seed):
    v = np.asarray([x for x in vals if np.isfinite(x)], float)
    if len(v) == 0:
        return float("nan"), float("nan"), float("nan")
    rng = np.random.default_rng(seed)
    meds = np.array([np.median(rng.choice(v, len(v))) for _ in range(nboot)])
    lo, hi = np.percentile(meds, [5, 95])
    return float(np.median(v)), float(lo), float(hi)


# ----------------------------------------------------------------- decision rules

def _fin(x):
    try:
        return x is not None and math.isfinite(float(x))
    except (TypeError, ValueError):
        return False


def grade_onset(x, ok, P=PARAMS):
    """x: onset change (2009 minus 2000) of radial SKS at 20% of peak, relative to the
    transverse-S lag, location-corrected, metadata rotation (s; negative = earlier)."""
    if not ok or not _fin(x):
        return "untestable"
    if x <= -P["onset_adv_s"]:
        return "onset-advance"
    if x >= -P["onset_none_s"]:
        return "no-onset-advance"
    return "partial"


def grade_width(s, ok, P=PARAMS):
    if not ok or not _fin(s):
        return "untestable"
    if s >= P["w_broad"]:
        return "broader"
    if s <= P["w_narrow"]:
        return "narrower"
    if abs(s - 1.0) <= P["w_same"]:
        return "same"
    return "slight"


def grade_source(w_grade, s_r, s_ref, s_glob, P=PARAMS):
    """Is a width change of radial SKS shared by the same-azimuth references at SDV (s_ref)
    and the azimuth-averaged global P (s_glob)?"""
    if w_grade == "untestable":
        return "untestable"
    if w_grade == "same":
        return "not-needed"
    if not (_fin(s_r) and _fin(s_ref) and _fin(s_glob)):
        return "untestable"
    if abs(s_r - s_ref) <= P["src_tol_ref"] and abs(s_r - s_glob) <= P["src_tol_glob"]:
        return "source"
    if abs(s_r - 1.0) - max(abs(s_ref - 1.0), abs(s_glob - 1.0)) >= P["path_excess"]:
        return "path-specific"
    return "unclear"


def headline(o_grade, s_grade):
    if o_grade == "untestable":
        return "untestable"
    if o_grade == "no-onset-advance":
        return "no-advance-in-either-reading"
    if o_grade == "partial":
        return "partial-onset-advance"
    if s_grade == "source":
        return "onset-advance-from-source-width"
    if s_grade == "path-specific":
        return "onset-advance-path-specific"
    return "onset-advance-unexplained"


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
    delta23 = ((p23.get("pairs") or {}).get(S.PARTNER) or {}).get("delta")
    summary["p23"] = dict(x_meta=x23, delta=delta23)

    from obspy.geodetics import gps2dist_azimuth, locations2degrees
    from obspy.taup import TauPyModel
    SP = S.PARAMS
    D = S.Data(args.data_dir)
    model = TauPyModel("iasp91")
    EV = {r[0]: S.ev_dict(r) for r in S.EVENTS}
    M, C = EV[S.MASTER], EV[S.PARTNER]
    dist_deg = locations2degrees(M["lat"], M["lon"], S.SDV["lat"], S.SDV["lon"])
    dist_km, az_src, baz = gps2dist_azimuth(M["lat"], M["lon"], S.SDV["lat"], S.SDV["lon"])
    dist_km /= 1000.0
    ph = S.phase_info(model, M["depth"], dist_deg)
    tP, tS, tSKS = ph["P"][0], ph["S"][0], ph["SKS"][0]
    pP, pS, pSKS = ph["P"][1], ph["S"][1], ph["SKS"][1]
    summary["geometry"] = dict(dist_deg=dist_deg, dist_km=dist_km, az_src=az_src, baz=baz,
                               tP=tP, tS=tS, tSKS=tSKS)
    tb, te = tP - 300.0, dist_km / 2.6 + 150.0      # identical to sdv23 (cache keys)
    tbb, teb = tP - 150.0, tS + 120.0
    body, meta = {}, {}
    for eid, E in ((S.MASTER, M), (S.PARTNER, C)):
        st, inv, _loc = S.sdv_event(D, E, tb, te)
        if st is None:
            raise RuntimeError("no IU.SDV data for " + eid)
        b = S.process(st, inv, E["t"], tbb, teb, "DISP", SP["pre_filt_body"], SP["band_body"],
                      SP["sr_up"])
        body[eid] = b
        meta[eid] = (b["1"]["az"], b["2"]["az"])
    A, B = body[S.MASTER], body[S.PARTNER]
    sr = A["Z"]["sr"]
    off = A["Z"]["off"]
    idx = lambda t: int(round((t - tbb - off) * sr))  # noqa: E731
    tt = lambda i: tbb + off + i / sr  # noqa: E731
    RA, TA = S.to_rt(A["1"]["x"], A["2"]["x"], meta[S.MASTER][0], meta[S.MASTER][1], baz)

    def peak(x, tc, lo=-SP["peak_search_s"], hi=SP["peak_search_s"]):
        i1, i2 = idx(tc + lo), idx(tc + hi)
        return tt(i1 + int(np.argmax(np.abs(x[i1:i2]))))
    kR, kZ = peak(RA, tSKS), peak(A["Z"]["x"], tSKS)
    kT, kRs, kP = peak(TA, tS, hi=8.0), peak(RA, tS, hi=8.0), peak(A["Z"]["x"], tP)
    wins = {"SKS_R": ("R", kR - 2.0, kR + 3.0), "postSKS_R": ("R", kR + 5.0, kR + 11.0),
            "SKS_Z": ("Z", kZ - 2.0, kZ + 3.0), "S_T": ("T", kT - 3.0, kT + 10.0),
            "S_R_early": ("R", kRs - 2.0, kRs + 1.5), "S_R_mid": ("R", kRs + 1.5, kRs + 5.0),
            "P_Z": ("Z", kP - 1.0, kP + 4.0), "S_T_pulse": ("T", kT - 2.0, kT + 3.0)}
    summary["windows"] = {k: [v[0], round(v[1], 3), round(v[2], 3)] for k, v in wins.items()}
    slow = {"SKS_R": pSKS, "SKS_Z": pSKS, "S_T_pulse": pS, "P_Z": pP}

    mres, mrows = S.membership(D, model, M, C, SP, {})
    reloc = mres.get("reloc")
    summary["membership"] = mres
    log("membership", mres.get("n_cc090"), mres.get("n_cc095"), reloc)

    aA1, aA2 = meta[S.MASTER]
    aB1, aB2 = meta[S.PARTNER]
    rots = {"meta": (aA1, aA2, aB1, aB2), "naive": (0.0, 90.0, 0.0, 90.0)}
    if _fin(delta23):
        rots["cons"] = (aA1, aA2, aA1 + float(delta23), aA2 + float(delta23))
    lagrows = S.measure_sdv_pair(A, B, rots, baz, wins, tbb, SP)
    lag = {(r["rotation"], r["window"]): r for r in lagrows}
    offA = float(np.mean([A[c]["off"] for c in ("1", "2", "Z")]))
    offB = float(np.mean([B[c]["off"] for c in ("1", "2", "Z")]))
    back = int(P["onset_back_s"] * sr)
    sml = int(P["stretch_maxlag_s"] * sr)

    def lp(w, ref_slow):
        if not reloc:
            return 0.0
        return S.loc_pred(reloc["dxN"], reloc["dxE"], az_src, slow[w] - ref_slow)

    rows, traces = [], {}
    for rname, (a1, a2, b1, b2) in rots.items():
        RA_, TA_ = S.to_rt(A["1"]["x"], A["2"]["x"], a1, a2, baz)
        RB_, TB_ = S.to_rt(B["1"]["x"], B["2"]["x"], b1, b2, baz)
        sh = int(round(lag[(rname, "S_T")]["coarse_s"] * sr))
        comps = {"R": (RA_, S.shift_left(RB_, sh)), "T": (TA_, S.shift_left(TB_, sh)),
                 "Z": (A["Z"]["x"], S.shift_left(B["Z"]["x"], sh))}
        envs = {c: (envelope(a), envelope(b)) for c, (a, b) in comps.items()}
        lagST = lag[(rname, "S_T")]["lag_s"]
        lagPZ = lag[(rname, "P_Z")]["lag_s"]
        for w in PULSES:
            comp, t1, t2 = wins[w]
            a, b = comps[comp]
            ea, eb = envs[comp]
            i1, i2 = idx(t1), idx(t2)
            ic = peak_index(ea, i1, i2)
            fit = stretch_fit(a, b, i1, i2, ic, P["s_min"], P["s_max"], P["s_step"], sml)
            L = lag[(rname, w)]
            row = dict(rotation=rname, window=w, comp=comp, t1=round(t1, 3), t2=round(t2, 3),
                       lag_s=L["lag_s"], cc=L["cc"],
                       dtc_vs_ST=L["lag_s"] - lagST - lp(w, pS),
                       dtc_vs_PZ=L["lag_s"] - lagPZ - lp(w, pP),
                       s=fit["s"], s_cc=fit["cc"], cc_s1=fit["cc_s1"], s_edge=fit["edge"])
            for f in P["onset_fracs"]:
                tag = "on%02d" % int(round(f * 100))
                onA, ipA = onset_index(ea, i1, i2, f, back)
                onB, _ = onset_index(eb, i1, i2, f, back)
                d_on = ((onB + sh - onA) / sr + offB - offA) if (_fin(onA) and _fin(onB)) \
                    else float("nan")
                row[tag + "_vs_ST"] = d_on - lagST - lp(w, pS)
                row[tag + "_vs_PZ"] = d_on - lagPZ - lp(w, pP)
                row[tag + "_vs_own"] = d_on - L["lag_s"]
                row[tag + "_rise_A_s"] = (ipA - onA) / sr if _fin(onA) else float("nan")
            row["pred_on20_from_s"] = (-(fit["s"] - 1.0) * row["on20_rise_A_s"]
                                       if _fin(fit["s"]) and _fin(row["on20_rise_A_s"])
                                       else float("nan"))
            rows.append(row)
            if rname == "meta" and w == "SKS_R":
                traces = dict(a=a[i1 - back:i2], b=b[i1 - back:i2], t0=tt(i1 - back), sr=sr)
        log("rotation", rname, "done")

    # ---- global P width reference (stations that verified the doublet)
    glob = []
    for r in mrows:
        if r.get("status") != "ok" or not (_fin(r.get("cc")) and r["cc"] >= SP["p_member_cc"]):
            continue
        if time.time() - t_start > budget:
            summary["errors"].append("global P stopped by time budget")
            break
        g = dict(net=r["net"], sta=r["sta"], dist=r["dist"], az=r["az"], status="")
        try:
            tPs = S.phase_info(model, M["depth"], r["dist"])["P"][0]
            arrs = []
            for E in (M, C):
                st = D.waveforms(r["net"], r["sta"], "*", "BHZ", E["t"] + tPs - 120,
                                 E["t"] + tPs + 120)
                if st is None:
                    raise RuntimeError("nodata " + E["id"])
                locs = sorted({tr.stats.location for tr in st},
                              key=lambda l: {"00": 0, "": 1, "10": 2}.get(l, 3))
                st = st.select(location=locs[0])
                inv = D.stations("p_%s_%s_%s" % (r["net"], r["sta"], E["id"]), network=r["net"],
                                 station=r["sta"], location=locs[0] or "--", channel="BHZ",
                                 starttime=E["t"] + tPs - 120, endtime=E["t"] + tPs + 120,
                                 level="response")
                pr = S.process(st, inv, E["t"], tPs - 60, tPs + 60, "DISP", SP["pre_filt_body"],
                               SP["band_body"], SP["sr_up"])
                arrs.append(pr["Z"])
            a, b = arrs
            srp = a["sr"]
            s1 = int(round((tPs - 3 - (tPs - 60) - a["off"]) * srp))
            s2 = int(round((tPs + 5 - (tPs - 60) - a["off"]) * srp))
            pk = s1 + int(np.argmax(np.abs(a["x"][s1:s2])))
            i1, i2 = pk - int(1 * srp), pk + int(4 * srp)
            kc, _, edc = window_lag(a["x"], b["x"], i1, i2, int(P["coarse_p_s"] * srp))
            if not _fin(kc) or edc:
                g["status"] = "coarse-fail"
                glob.append(g)
                continue
            bsh = S.shift_left(b["x"], int(round(kc)))
            ea, eb = envelope(a["x"]), envelope(bsh)
            ic = peak_index(ea, i1, i2)
            fit = stretch_fit(a["x"], bsh, i1, i2, ic, P["s_min"], P["s_max"], P["s_step"],
                              int(P["stretch_maxlag_s"] * srp))
            k1, cc1, _ = window_lag(a["x"], bsh, i1, i2, int(P["stretch_maxlag_s"] * srp))
            onA, ipA = onset_index(ea, i1, i2, 0.2, int(P["onset_back_s"] * srp))
            onB, _ = onset_index(eb, i1, i2, 0.2, int(P["onset_back_s"] * srp))
            g.update(s=fit["s"], s_cc=fit["cc"], cc_s1=cc1,
                     on20_vs_own=((onB - onA - k1) / srp) if (_fin(onA) and _fin(onB))
                     else float("nan"),
                     on20_rise_A_s=((ipA - onA) / srp) if _fin(onA) else float("nan"),
                     status="ok")
        except Exception as ex:  # noqa: BLE001
            g["status"] = "error: " + repr(ex)[:120]
        glob.append(g)
    log("global P stations", sum(g["status"] == "ok" for g in glob))

    # ---- grades
    def get(rot, w, k):
        return next((x[k] for x in rows if x["rotation"] == rot and x["window"] == w),
                    float("nan"))
    x_meta_ours = get("meta", "SKS_R", "dtc_vs_ST")
    x_on = get("meta", "SKS_R", "on20_vs_ST")
    gates = dict(
        G0_p23_reproduced=bool(_fin(x23) and _fin(x_meta_ours)
                               and abs(x_meta_ours - float(x23)) <= P["p23_tol_s"]),
        G1_sks_cc=bool(_fin(get("meta", "SKS_R", "cc"))
                       and get("meta", "SKS_R", "cc") >= P["sks_cc_min"]),
        G2_onsets=bool(_fin(x_on)),
    )
    all_ok = all(gates.values())
    og = grade_onset(x_on, all_ok)
    s_R = get("meta", "SKS_R", "s")
    wg = grade_width(s_R, all_ok and _fin(get("meta", "SKS_R", "s_cc"))
                     and get("meta", "SKS_R", "s_cc") >= P["sks_cc_min"])
    refs = [get("meta", w, "s") for w in ("P_Z", "S_T_pulse")
            if _fin(get("meta", w, "s_cc")) and get("meta", w, "s_cc") >= P["ref_cc_min"]]
    s_ref = float(np.median(refs)) if refs else float("nan")
    gl = [g["s"] for g in glob if g.get("status") == "ok" and _fin(g.get("s_cc"))
          and g["s_cc"] >= P["glob_cc_min"]]
    med, lo, hi = boot_median(gl, P["nboot"], P["seed"])
    s_glob = med if len(gl) >= P["glob_min_n"] else float("nan")
    sg = grade_source(wg, s_R, s_ref, s_glob)
    prim = dict(gates=gates, x_meta_ours=x_meta_ours, x_meta_p23=x23, x_onset20=x_on,
                x_onset10=get("meta", "SKS_R", "on10_vs_ST"),
                x_onset30=get("meta", "SKS_R", "on30_vs_ST"),
                x_onset20_naive=get("naive", "SKS_R", "on20_vs_ST"),
                x_onset20_vs_PZ=get("meta", "SKS_R", "on20_vs_PZ"),
                onset_grade=og, s_R=s_R, width_grade=wg, s_ref=s_ref, n_ref=len(refs),
                s_glob=s_glob, s_glob_5_95=[lo, hi], n_glob=len(gl), source_grade=sg,
                headline=headline(og, sg))
    summary["primary"] = prim
    log("PRIMARY", prim)

    wcsv(os.path.join(args.out_dir, "sdvwidth25_rows.csv"), rows)
    wcsv(os.path.join(args.out_dir, "sdvwidth25_globalP.csv"), glob)
    wcsv(os.path.join(args.out_dir, "sdvwidth25_lags.csv"), lagrows)
    summary["finished"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    summary["runtime_min"] = round((time.time() - t_start) / 60.0, 1)
    with open(os.path.join(args.out_dir, "sdvwidth25_summary.json"), "w") as f:
        json.dump(summary, f, indent=1, default=jsonable)
    try:
        plot(args.out_dir, traces, rows, glob)
    except Exception as ex:  # noqa: BLE001
        log("plot failed", repr(ex)[:200])
    return summary


def plot(out_dir, traces, rows, glob):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    os.makedirs(os.path.join(out_dir, "plots"), exist_ok=True)
    fig, ax = plt.subplots(2, 1, figsize=(9, 8))
    if traces:
        for k, lab in (("a", "2000-07-07"), ("b", "2009-12-17 (aligned on transverse S)")):
            x = np.asarray(traces[k], float)
            t = traces["t0"] + np.arange(len(x)) / traces["sr"]
            m = np.max(np.abs(x)) or 1.0
            ax[0].plot(t, x / m, label=lab)
            ax[0].plot(t, envelope(x) / m, lw=0.6, color=ax[0].lines[-1].get_color())
        ax[0].set_title("IU.SDV radial SKS, 0.4-2 Hz displacement, metadata rotation")
        ax[0].set_xlabel("time after the 2000 origin (s)")
        ax[0].legend(fontsize=8)
    gs = [g["s"] for g in glob if g.get("status") == "ok" and _fin(g.get("s"))]
    if gs:
        ax[1].hist(gs, bins=20, color="0.6", label="vertical P, global stations")
    for w, col in (("SKS_R", "C3"), ("SKS_Z", "C1"), ("P_Z", "C0"), ("S_T_pulse", "C2")):
        v = next((r["s"] for r in rows if r["rotation"] == "meta" and r["window"] == w), None)
        if _fin(v):
            ax[1].axvline(v, color=col, label="IU.SDV " + w)
    ax[1].set_xlabel("stretch factor s (2009 pulse width / 2000 pulse width)")
    ax[1].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, "plots", "sdvwidth25.png"), dpi=110)


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", default="/data")
    ap.add_argument("--out-dir", default="/out")
    ap.add_argument("--wait-min", type=float, default=120.0)
    ap.add_argument("--time-budget-min", type=float, default=90.0)
    a = ap.parse_args(argv)
    run(a)
    return 0


if __name__ == "__main__":
    sys.exit(main())

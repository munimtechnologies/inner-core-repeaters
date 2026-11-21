#!/usr/bin/env python3
# Copyright 2026 MV2, Munim, Inc.
# Licensed under the Apache License, Version 2.0 (Apache-2.0). See LICENSE.
# Made with AI agents working in MV2's name.
"""sdvz24.py - pre-registration 24: which sensor orientation do the published IU.SDV
waveforms of Zhang & Wen (2025) assume for each event of doublet AI-0009, and how much
does the claimed 160 ms SKS advance depend on it?

Plan: findings/2026-11-21-prereg24-sdv-zenodo-frame-plan.md (written and published
before this code downloads the authors' data or any waveform).

Data: Zhang & Wen, "Preprocessed Doublet Data of AI-0009 ...", Zenodo
doi 10.5281/zenodo.16421554 (CC BY 4.0), file data_IU.SDV.zip (MD5
9ba0853da7a537a7a8f1117407675112). Their horizontals are rotated to the great-circle
path (R/T). We compare them with raw IU.SDV location-00 records from the EarthScope
FDSN service, converted to geographic N/E with the station metadata of each event time.

Parts:
 F. frame fit: for each event, the rotation a_fit (deg) of the metadata N/E frame and
    the time lag that make our R/T match theirs best (closed-form angle at each lag, 0.4-2
    Hz displacement, window SKS-30 s .. S+60 s). The channel-1 azimuth their processing
    assumed is A = metadata azimuth + a_fit. rho = A(2009) - A(2000).
 M. on their own traces: x = lag(SKS_R) - lag(S_T) (2009 relative to 2000), after a
    coarse alignment on S_T, and x as a function of an extra rotation da of the 2009 R/T.
 RZ. reported: Rayleigh-wave (Rayleigh-Z) azimuth of each event's location-00 sensor.

Top-level imports are stdlib + numpy (+ sdv23, same repository); ObsPy inside functions.
"""

import argparse
import csv
import hashlib
import json
import math
import os
import sys
import time
import zipfile

import numpy as np

import sdv23 as base

VERSION = "0.1.0"
PLAN_PATH = "findings/2026-11-21-prereg24-sdv-zenodo-frame-plan.md"
ZENODO_URL = "https://zenodo.org/api/records/16421554/files/data_IU.SDV.zip/content"
ZENODO_MD5 = "9ba0853da7a537a7a8f1117407675112"
# USGS ComCat origins (as in sdv23.EVENTS): id -> (time, lat, lon, depth km)
EVENTS = {"E2000a": ("2000-07-07T15:46:44.56", 51.411, 179.978, 31.0),
          "E2009a": ("2009-12-17T20:01:21.61", 51.398, 179.963, 35.0)}
MASTER, PARTNER = "E2000a", "E2009a"

PARAMS = dict(
    band=(0.4, 2.0), pre_filt=(0.1, 0.2, 7.0, 9.0),
    fit_sr=50.0, fit_maxlag_s=5.0, fit_pre_s=30.0, fit_post_s=60.0, fit_min_s=40.0,
    fit_c_min=0.95, frame_tol=1.0, meta_step=-3.0,
    sr_up=200.0, grid_pre_s=40.0, grid_post_s=60.0,
    coarse_maxlag_s=5.0, maxlag_s=0.6, peak_search_s=5.0,
    st_cc_min=0.90, sks_cc_min=0.80, reproduce_s=0.10,
    sens_min=-20.0, sens_max=20.0, sens_step=0.5,
    consistent_tol=1.0, misrot_min=2.0,
    rz_band=(0.02, 0.05), rz_pre_filt=(0.003, 0.005, 0.5, 1.0), rz_c_min=0.80,
    rz_half=45.0, rz_step=0.1,
    max_member_mb=300,
)


def log(*a):
    print(time.strftime("%H:%M:%S"), "[sdvz24]", *a, flush=True)


# ----------------------------------------------------------------- pure numerics

def reframe(r, t, da):
    """R/T recomputed as if channel 1 pointed da degrees further clockwise than the
    azimuth used to make (r, t). reframe(to_rt(c1, c2, a, a+90, b), da) ==
    to_rt(c1, c2, a+da, a+da+90, b)."""
    a = math.radians(da)
    r = np.asarray(r, float)
    t = np.asarray(t, float)
    return math.cos(a) * r - math.sin(a) * t, math.cos(a) * t + math.sin(a) * r


def fit_frame(c1, c2, rp, tp, baz, maxlag):
    """Find azimuth a (deg; channel 2 at a + 90) and lag k (samples) such that
    to_rt(c1[m+k : m+k+n], c2[...], a, a + 90, baz) best matches (rp, tp), with n = len(rp)
    and m = maxlag; c1, c2 must have n + 2m samples. Two conventions are tried: 'proper'
    (ObsPy/SAC: T = R + 90 deg clockwise) and 'reflected' (their T negated).
    Returns dict(a, c, lag, conv) with c the normalised vector correlation."""
    c1 = np.asarray(c1, float)
    c2 = np.asarray(c2, float)
    rp = np.asarray(rp, float)
    tp = np.asarray(tp, float)
    n, m = len(rp), int(maxlag)
    if n < 3 or len(tp) != n or len(c1) != n + 2 * m or len(c2) != n + 2 * m:
        raise ValueError("fit_frame: bad lengths")
    nt = math.sqrt(float(np.dot(rp, rp) + np.dot(tp, tp)))
    best = dict(a=float("nan"), c=-1.0, lag=0, conv="")
    if nt <= 0:
        return best
    for conv, sgn in (("proper", 1.0), ("reflected", -1.0)):
        tq = sgn * tp
        for k in range(-m, m + 1):
            s1 = c1[m + k:m + k + n]
            s2 = c2[m + k:m + k + n]
            ns = math.sqrt(float(np.dot(s1, s1) + np.dot(s2, s2)))
            if ns <= 0:
                continue
            X = -(float(np.dot(s1, rp)) + float(np.dot(s2, tq)))
            Y = float(np.dot(s1, tq)) - float(np.dot(s2, rp))
            c = math.hypot(X, Y) / (ns * nt)
            if c > best["c"]:
                th = math.degrees(math.atan2(Y, X))
                best = dict(a=float((baz - th) % 360.0), c=float(c), lag=int(k), conv=conv)
    return best


def on_grid(x, sr, t0, t1, t2, sr_out):
    """Samples of x (rate sr, first sample at time t0) at t1 + k / sr_out, k = 0 ..
    round((t2 - t1) * sr_out) - 1. FFT-upsampled first if sr < sr_out; linear
    interpolation; NaN outside the trace."""
    x = np.asarray(x, float)
    if sr < sr_out:
        f = int(math.ceil(sr_out / sr - 1e-9))
        x = base.upsample(x, f)
        sr = sr * f
    n = int(round((t2 - t1) * sr_out))
    tg = t1 + np.arange(n) / sr_out
    tx = t0 + np.arange(len(x)) / sr
    return np.interp(tg, tx, x, left=np.nan, right=np.nan)


def classify_frame(rho, P=PARAMS):
    """rho = assumed channel-1 azimuth in 2009 minus in 2000 (deg)."""
    if rho is None or not np.isfinite(rho):
        return "untestable"
    if abs(float(base.wrap180(rho))) <= P["frame_tol"]:
        return "same-frame"
    if abs(float(base.wrap180(rho - P["meta_step"]))) <= P["frame_tol"]:
        return "metadata-frame"
    return "other"


def grade_combined(rho, delta, P=PARAMS):
    """eps = rho - delta: rotation error of their 2009 R/T relative to 2000, given the
    physical sensor rotation delta between the events. Returns (grade, eps)."""
    if rho is None or delta is None or not np.isfinite(rho) or not np.isfinite(delta):
        return "unresolved", float("nan")
    eps = float(base.wrap180(rho - delta))
    if abs(eps) <= P["consistent_tol"]:
        return "consistent", eps
    if abs(eps) >= P["misrot_min"]:
        return "misrotated", eps
    return "unclear", eps


def grade_reproduce(gates, x0, P=PARAMS):
    if not gates or not all(bool(v) for v in gates.values()):
        return "untestable"
    if x0 is None or not np.isfinite(x0):
        return "untestable"
    return "reproduced" if x0 <= -P["reproduce_s"] else "not-reproduced"


def pick_windows(ra, ta, sr, g0, tsks, ts, P=PARAMS):
    """Measurement windows (times after origin) from the master's R and T on a grid
    starting at g0 with rate sr."""
    idx = lambda t: int(round((t - g0) * sr))  # noqa: E731
    tt = lambda i: g0 + i / sr  # noqa: E731

    def peak(x, tc, lo, hi):
        i1, i2 = max(0, idx(tc + lo)), min(len(x), idx(tc + hi))
        if i2 - i1 < 2:
            raise ValueError("peak search outside the grid")
        seg = np.nan_to_num(np.asarray(x[i1:i2], float))
        return tt(i1 + int(np.argmax(np.abs(seg))))
    ps = P["peak_search_s"]
    kR = peak(ra, tsks, -ps, ps)
    kT = peak(ta, ts, -ps, 8.0)
    kRs = peak(ra, ts, -ps, 8.0)
    return {"SKS_R": ("R", kR - 2.0, kR + 3.0), "postSKS_R": ("R", kR + 5.0, kR + 11.0),
            "S_T": ("T", kT - 3.0, kT + 10.0), "S_R_early": ("R", kRs - 2.0, kRs + 1.5),
            "S_R_mid": ("R", kRs + 1.5, kRs + 5.0)}


def measure(ra, ta, rb, tb2, sr, g0, wins, P=PARAMS, shift=None):
    """Lags (s; > 0 = later in b) of b against a in each window, after a coarse integer
    alignment on S_T (+/- coarse_maxlag_s) unless shift (samples) is given.
    Returns (dict window -> dict(lag, cc, edge), shift)."""
    idx = lambda t: int(round((t - g0) * sr))  # noqa: E731
    ra, ta, rb, tb2 = (np.nan_to_num(np.asarray(v, float)) for v in (ra, ta, rb, tb2))
    if shift is None:
        _, t1, t2 = wins["S_T"]
        kc, _, edge = base.window_lag(ta, tb2, idx(t1), idx(t2),
                                      int(round(P["coarse_maxlag_s"] * sr)))
        shift = int(round(kc)) if (np.isfinite(kc) and not edge) else 0
    comps = {"R": (ra, base.shift_left(rb, shift)), "T": (ta, base.shift_left(tb2, shift))}
    ml = int(round(P["maxlag_s"] * sr))
    out = {}
    for name, (c, t1, t2) in wins.items():
        a, b = comps[c]
        k, cc, edge = base.window_lag(a, b, idx(t1), idx(t2), ml)
        out[name] = dict(lag=(k + shift) / sr if np.isfinite(k) else float("nan"),
                         cc=cc, edge=bool(edge))
    return out, int(shift)


def needed_range(wins, P=PARAMS):
    """Time span (after origin) that both events must cover for every window."""
    pad = P["coarse_maxlag_s"] + P["maxlag_s"] + 0.5
    return (min(v[1] for v in wins.values()) - pad, max(v[2] for v in wins.values()) + pad)


def covered(x, g0, sr, t1, t2):
    i1, i2 = int(math.floor((t1 - g0) * sr)), int(math.ceil((t2 - g0) * sr))
    if i1 < 0 or i2 > len(x):
        return False
    return bool(np.all(np.isfinite(np.asarray(x[i1:i2], float))))


def set_key(rel_path, cha, sr, start_iso):
    """Group key for traces of one processing set (same file name apart from the
    component code)."""
    base_name = rel_path.replace("\\", "/")
    if cha and cha in os.path.basename(base_name):
        d, f = os.path.split(base_name)
        return d + "/" + f.replace(cha, cha[:-1] + "?")
    return "%s|%.3f|%s" % (os.path.dirname(base_name), float(sr), start_iso[:21])


def comp_of(cha, fname):
    c = (cha or "").strip().upper()
    if c:
        return c[-1]
    toks = os.path.basename(fname).upper().replace("_", ".").split(".")
    for tok in reversed(toks):
        if tok in ("R", "T", "Z", "N", "E", "BHR", "BHT", "BHZ", "BHN", "BHE", "HHR", "HHT"):
            return tok[-1]
    return ""


def safe_member(name):
    """Normalised relative path of a zip member, or None if unsafe."""
    n = os.path.normpath(name.replace("\\", "/"))
    if n.startswith("..") or os.path.isabs(n) or ":" in n or n in (".", ""):
        return None
    return n


def md5(data):
    return hashlib.md5(data).hexdigest()


# ----------------------------------------------------------------- data steps

def extract_zip(path, dest, max_mb, depth=0):
    """Extract a zip safely (and zips inside it, one level). Returns list of files."""
    out = []
    with zipfile.ZipFile(path) as z:
        for info in z.infolist():
            if info.is_dir():
                continue
            rel = safe_member(info.filename)
            if rel is None or info.file_size > max_mb * 1024 * 1024:
                continue
            target = os.path.join(dest, rel)
            os.makedirs(os.path.dirname(target) or dest, exist_ok=True)
            with z.open(info) as src, open(target, "wb") as dst:
                dst.write(src.read())
            if rel.lower().endswith(".zip") and depth == 0:
                out += extract_zip(target, target + "_x", max_mb, depth + 1)
            else:
                out.append(target)
    return out


SAC_KEYS = ["kstnm", "khole", "kcmpnm", "cmpaz", "cmpinc", "baz", "az", "gcarc", "dist",
            "evla", "evlo", "evdp", "stla", "stlo", "o", "b", "e", "delta", "npts",
            "kevnm", "t0", "t1", "t2", "t3", "kt0", "kt1", "kt2", "kt3",
            "user0", "user1", "user2", "kuser0", "iztype", "idep"]


def manifest(files, root, origins):
    """Header-only read of every file. Returns rows (one per trace)."""
    from obspy import read
    rows = []
    for f in sorted(files):
        rel = os.path.relpath(f, root).replace("\\", "/")
        try:
            st = read(f, headonly=True)
        except Exception as ex:  # noqa: BLE001
            rows.append(dict(path=rel, status="unreadable: " + repr(ex)[:80]))
            continue
        for i, tr in enumerate(st):
            s = tr.stats
            sac = s.get("sac", {}) or {}
            ev = ""
            for eid, t0 in origins.items():
                if t0 - 3600 <= s.starttime <= t0 + 7200:
                    ev = eid
            row = dict(path=rel, index=i, status="ok", format=s.get("_format", ""),
                       net=s.network, sta=s.station, loc=s.location, cha=s.channel,
                       comp=comp_of(s.channel, rel), start=str(s.starttime),
                       end=str(s.endtime), sr=float(s.sampling_rate), npts=int(s.npts),
                       event=ev, start_rel=(float(s.starttime - origins[ev]) if ev
                                            else float("nan")))
            for k in SAC_KEYS:
                v = sac.get(k, "")
                row["sac_" + k] = v.strip() if isinstance(v, str) else v
            row["set"] = set_key(rel, s.channel, s.sampling_rate, str(s.starttime))
            rows.append(row)
    return rows


def their_trace(root, row, origin, P):
    """Full read + demean/detrend/taper + the common band-pass. Returns (x, sr, t0)."""
    from obspy import read
    st = read(os.path.join(root, row["path"]))
    tr = st[int(row["index"])].copy()
    tr.data = np.asarray(tr.data, dtype=np.float64)
    tr.detrend("demean")
    tr.detrend("linear")
    tr.taper(0.05)
    tr.filter("bandpass", freqmin=P["band"][0], freqmax=P["band"][1], corners=2,
              zerophase=True)
    return np.asarray(tr.data, float), float(tr.stats.sampling_rate), float(tr.stats.starttime
                                                                            - origin)


def run(args):
    from obspy import UTCDateTime
    from obspy.geodetics import gps2dist_azimuth, locations2degrees
    from obspy.taup import TauPyModel

    P = PARAMS
    t_start = time.time()
    os.makedirs(args.out_dir, exist_ok=True)
    S = dict(version=VERSION, base_version=base.VERSION, params=P, errors=[],
             started=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()))
    with open(os.path.abspath(__file__), "rb") as f:
        S["code_sha256"] = base.sha256(f.read())
    ok, psha = base.check_plan(args.wait_min, PLAN_PATH)
    S["plan_public_at_start"], S["plan_sha256"] = ok, psha
    log("plan public:", ok, psha[:16])

    # ---- authors' data
    zdir = os.path.join(args.data_dir, "sdvz24")
    os.makedirs(zdir, exist_ok=True)
    zpath = os.path.join(zdir, "data_IU.SDV.zip")
    if not os.path.exists(zpath):
        data = base.fetch_url(ZENODO_URL, tries=3, timeout=600)
        with open(zpath, "wb") as f:
            f.write(data)
    with open(zpath, "rb") as f:
        zmd5 = md5(f.read())
    S["zip_md5"] = zmd5
    S["G0_zip_md5_ok"] = zmd5 == ZENODO_MD5
    log("zip md5", zmd5, S["G0_zip_md5_ok"])
    root = os.path.join(zdir, "x")
    files = extract_zip(zpath, root, P["max_member_mb"])
    S["n_files"] = len(files)
    origins = {eid: UTCDateTime(v[0]) for eid, v in EVENTS.items()}
    man = manifest(files, root, origins)
    S["n_traces"] = sum(1 for r in man if r.get("status") == "ok")
    log("files", len(files), "traces", S["n_traces"])

    # ---- geometry (master)
    model = TauPyModel("iasp91")
    t_m, la_m, lo_m, dp_m = EVENTS[MASTER]
    dist_deg = locations2degrees(la_m, lo_m, base.SDV["lat"], base.SDV["lon"])
    dist_km, az_src, baz = gps2dist_azimuth(la_m, lo_m, base.SDV["lat"], base.SDV["lon"])
    dist_km /= 1000.0
    ph = base.phase_info(model, dp_m, dist_deg)
    tP, tS, tSKS = ph["P"][0], ph["S"][0], ph["SKS"][0]
    S["geometry"] = dict(dist_deg=dist_deg, dist_km=dist_km, baz=baz, tP=tP, tS=tS,
                         tSKS=tSKS)
    log("dist %.3f baz %.3f tSKS %.1f tS %.1f" % (dist_deg, baz, tSKS, tS))

    # ---- their sets per event
    sets = {}
    for r in man:
        if r.get("status") != "ok" or not r.get("event") or r.get("sta") not in ("SDV", ""):
            continue
        sets.setdefault((r["event"], r["set"]), {})[r["comp"]] = r
    fw1, fw2 = tSKS - P["fit_pre_s"], tS + P["fit_post_s"]
    g0, g1 = tSKS - P["grid_pre_s"], tS + P["grid_post_s"]

    def covers(rr, t1, t2):
        return all(rr[c]["start_rel"] <= t1 and
                   rr[c]["start_rel"] + (rr[c]["npts"] - 1) / rr[c]["sr"] >= t2
                   for c in ("R", "T"))
    primary = {}
    for eid in EVENTS:
        keys = sorted(k for (e, k), rr in sets.items() if e == eid and "R" in rr and "T" in rr)
        full = [k for k in keys if covers(sets[(eid, k)], g0, g1)]
        primary[eid] = (full or keys or [None])[0]
    S["sets"] = {"%s|%s" % k: sorted(v) for k, v in sets.items()}
    S["primary_set"] = primary
    log("primary sets", primary)

    # ---- raw IU.SDV data (location 00 primary; 10 reported)
    D = base.Data(args.data_dir)
    raw, meta, rz = {}, {}, {}
    tbb, teb = tP - 150.0, tS + 120.0
    for eid, (tiso, la, lo, dp) in EVENTS.items():
        t0 = origins[eid]
        for loc in ("00", "10"):
            try:
                st = D.waveforms("IU", "SDV", loc, "BH?", t0 - 600, t0 + dist_km / 2.6 + 300)
                if st is None:
                    S["errors"].append("raw nodata %s %s" % (eid, loc))
                    continue
                inv = D.stations("z24_%s_%s" % (eid, loc), network="IU", station="SDV",
                                 location=loc, channel="BH?", starttime=t0 - 600,
                                 endtime=t0 + dist_km / 2.6 + 300, level="response")
                b = base.process(st, inv, t0, tbb, teb, "DISP", P["pre_filt"], P["band"])
                if not all(c in b for c in ("1", "2")):
                    S["errors"].append("raw components %s %s" % (eid, loc))
                    continue
                raw[(eid, loc)] = b
                meta[(eid, loc)] = (b["1"]["az"], b["2"]["az"])
                if loc == "00" and "Z" in b:
                    w1, w2 = dist_km / 4.2 - 150.0, dist_km / 3.0 + 150.0
                    L = base.process(st, inv, t0, w1, w2, "VEL", P["rz_pre_filt"], P["rz_band"])
                    srl = L["Z"]["sr"]
                    i1 = int(round((dist_km / 4.2 - w1 - L["Z"]["off"]) * srl))
                    i2 = int(round((dist_km / 3.0 - w1 - L["Z"]["off"]) * srl))
                    a_hat, c = base.orient_rz(L["1"]["x"][i1:i2], L["2"]["x"][i1:i2],
                                              L["Z"]["x"][i1:i2], L["2"]["az"] - L["1"]["az"],
                                              baz, L["1"]["az"], P["rz_half"], P["rz_step"])
                    rz[eid] = dict(a_hat=a_hat, c=c, az1_meta=L["1"]["az"])
                log("raw", eid, loc, meta[(eid, loc)], rz.get(eid))
            except Exception as ex:  # noqa: BLE001
                S["errors"].append("raw %s %s: %s" % (eid, loc, repr(ex)[:140]))
    S["meta_azimuths"] = {"%s|%s" % k: v for k, v in meta.items()}
    S["rz"] = rz

    # ---- F. frame fits
    fs = P["fit_sr"]
    m = int(round(P["fit_maxlag_s"] * fs))
    fit_rows, fits = [], {}
    loaded = {}

    def load(eid, key):
        if (eid, key) not in loaded:
            rr = sets[(eid, key)]
            loaded[(eid, key)] = {c: their_trace(root, rr[c], origins[eid], P)
                                  for c in ("R", "T") if c in rr}
        return loaded[(eid, key)]
    for (eid, key), rr in sorted(sets.items()):
        if "R" not in rr or "T" not in rr:
            continue
        try:
            th = load(eid, key)
            Rg = on_grid(*th["R"], fw1, fw2, fs)
            Tg = on_grid(*th["T"], fw1, fw2, fs)
            ok_ = np.isfinite(Rg) & np.isfinite(Tg)
            if not ok_.any():
                raise ValueError("no overlap with the fit window")
            i_ok = np.flatnonzero(ok_)
            j1, j2 = int(i_ok[0]), int(i_ok[-1]) + 1
            if not ok_[j1:j2].all() or (j2 - j1) / fs < P["fit_min_s"]:
                raise ValueError("fit window coverage %.1f s" % ((j2 - j1) / fs))
            w1 = fw1 + j1 / fs
            n = j2 - j1
            for loc in ("00", "10"):
                if (eid, loc) not in raw:
                    continue
                b = raw[(eid, loc)]
                az1, az2 = meta[(eid, loc)]
                nn = min(len(b["1"]["x"]), len(b["2"]["x"]))
                N, E = base.rotate_to_ne(b["1"]["x"][:nn], b["2"]["x"][:nn], az1, az2)
                t0r = tbb + b["1"]["off"]
                Ng = on_grid(N, b["1"]["sr"], t0r, w1 - m / fs, w1 + (n + m) / fs, fs)
                Eg = on_grid(E, b["2"]["sr"], t0r, w1 - m / fs, w1 + (n + m) / fs, fs)
                if len(Ng) != n + 2 * m or not (np.all(np.isfinite(Ng))
                                                and np.all(np.isfinite(Eg))):
                    fit_rows.append(dict(event=eid, set=key, loc=loc, status="raw coverage"))
                    continue
                fr = fit_frame(Ng, Eg, Rg[j1:j2], Tg[j1:j2], baz, m)
                a_fit = float(base.wrap180(fr["a"]))
                row = dict(event=eid, set=key, loc=loc, a_fit=a_fit, az1_meta=az1,
                           assumed_az1=float((az1 + a_fit) % 360.0), c=fr["c"],
                           lag_s=fr["lag"] / fs, conv=fr["conv"], fit_s=n / fs,
                           sac_baz=rr["R"].get("sac_baz", ""), sac_cmpaz_R=rr["R"].get(
                               "sac_cmpaz", ""), sac_cmpaz_T=rr["T"].get("sac_cmpaz", ""))
                fit_rows.append(row)
                fits[(eid, key, loc)] = row
                log("fit", eid, key, loc, round(a_fit, 2), round(fr["c"], 4), fr["conv"])
        except Exception as ex:  # noqa: BLE001
            fit_rows.append(dict(event=eid, set=key, loc="", status="error: " + repr(ex)[:120]))

    # ---- frame grade (primary sets, location 00)
    fA = fits.get((MASTER, primary[MASTER], "00"))
    fB = fits.get((PARTNER, primary[PARTNER], "00"))
    gates_f = dict(G0_zip=S["G0_zip_md5_ok"], G1_sets=bool(primary[MASTER] and primary[PARTNER]),
                   G2_raw00=bool(fA is not None and fB is not None),
                   G3_fit_c=bool(fA and fB and fA["c"] >= P["fit_c_min"]
                                 and fB["c"] >= P["fit_c_min"]),
                   G3_same_conv=bool(fA and fB and fA["conv"] == fB["conv"]))
    rho = (float(base.wrap180(fB["assumed_az1"] - fA["assumed_az1"]))
           if (fA and fB) else float("nan"))
    fgrade = classify_frame(rho) if all(gates_f.values()) else "untestable"
    rz_ok = all(e in rz and rz[e]["c"] >= P["rz_c_min"] for e in EVENTS)
    d_rz = float(base.wrap180(rz[PARTNER]["a_hat"] - rz[MASTER]["a_hat"])) if rz_ok \
        else float("nan")
    comb, eps = grade_combined(rho if fgrade != "untestable" else float("nan"), d_rz)
    S["frame"] = dict(grade=fgrade, rho=rho, gates=gates_f, fit_master=fA, fit_partner=fB,
                      delta_rz=d_rz, combined_with_rz=comb, eps_rel_rz=eps)
    log("FRAME", fgrade, "rho", rho, "rz delta", d_rz, comb, eps)

    # ---- M. measurement on their primary traces
    meas = dict(grade="untestable")
    sens_rows = []
    try:
        if not (primary[MASTER] and primary[PARTNER]):
            raise ValueError("no R/T set for both events")
        sr = P["sr_up"]
        g = {}
        for eid in EVENTS:
            th = load(eid, primary[eid])
            conv = (fits.get((eid, primary[eid], "00")) or {}).get("conv") or "proper"
            sgn = -1.0 if conv == "reflected" else 1.0
            g[eid] = (on_grid(*th["R"], g0, g1, sr), sgn * on_grid(*th["T"], g0, g1, sr))
        RA, TA = g[MASTER]
        RB, TB = g[PARTNER]
        wins = pick_windows(RA, TA, sr, g0, tSKS, tS, P)
        n1, n2 = needed_range(wins, P)
        gates_m = dict(G1_coverage=all(covered(x, g0, sr, n1, n2) for x in (RA, TA, RB, TB)))
        res0, shift0 = measure(RA, TA, RB, TB, sr, g0, wins, P)
        gates_m["G4_ST_cc"] = bool(np.isfinite(res0["S_T"]["cc"])
                                   and res0["S_T"]["cc"] >= P["st_cc_min"]
                                   and not res0["S_T"]["edge"])
        gates_m["G4_SKS_cc"] = bool(np.isfinite(res0["SKS_R"]["cc"])
                                    and res0["SKS_R"]["cc"] >= P["sks_cc_min"]
                                    and not res0["SKS_R"]["edge"])
        x0 = res0["SKS_R"]["lag"] - res0["S_T"]["lag"]
        idx = lambda t: int(round((t - g0) * sr))  # noqa: E731

        def rms(x, w):
            s = np.nan_to_num(x[idx(w[1]):idx(w[2])])
            return float(np.sqrt(np.mean(s ** 2))) if len(s) else float("nan")

        def ratio(a, b):
            return a / b if (np.isfinite(a) and np.isfinite(b) and b > 0) else float("nan")
        ratios = {e: dict(T_over_R_SKS=ratio(rms(g[e][1], wins["SKS_R"]),
                                             rms(g[e][0], wins["SKS_R"])),
                          T_over_R_S=ratio(rms(g[e][1], wins["S_T"]), rms(g[e][0], wins["S_T"])))
                  for e in EVENTS}
        xs = {}
        for da in np.arange(P["sens_min"], P["sens_max"] + P["sens_step"] / 2, P["sens_step"]):
            da = round(float(da), 3)
            rb2, tb2 = reframe(RB, TB, da)
            res, _ = measure(RA, TA, rb2, tb2, sr, g0, wins, P, shift=shift0)
            xs[da] = res["SKS_R"]["lag"] - res["S_T"]["lag"]
            row = dict(da=da, x=xs[da])
            for wn, v in res.items():
                row[wn + "_lag"], row[wn + "_cc"], row[wn + "_edge"] = v["lag"], v["cc"], v["edge"]
            sens_rows.append(row)
        slope = (xs.get(1.0, float("nan")) - xs.get(-1.0, float("nan"))) / 2.0
        early = {d: next(r["S_R_early_lag"] for r in sens_rows if r["da"] == d)
                 for d in (-3.0, 0.0, 3.0)}
        null = [d for d in sorted(xs) if np.isfinite(xs[d])]
        da_null = min(null, key=lambda d: abs(xs[d])) if null else float("nan")
        meas = dict(grade=grade_reproduce(gates_m, x0), gates=gates_m, x0=x0,
                    lags0=res0, coarse_shift_s=shift0 / sr, windows=wins,
                    slope_s_per_deg=slope, x_at={str(d): xs.get(d) for d in
                                                  (-17.0, -3.0, 0.0, 3.0, 17.0)},
                    da_min_abs_x=da_null, x_min_abs=xs.get(da_null, float("nan")),
                    S_R_early_change_pm3=(early[3.0] - early[-3.0]), ratios=ratios)
    except Exception as ex:  # noqa: BLE001
        S["errors"].append("measurement: " + repr(ex)[:200])
    S["measurement"] = meas
    log("MEASUREMENT", meas.get("grade"), meas.get("x0"), meas.get("slope_s_per_deg"))

    # ---- outputs
    def wcsv(name, rows):
        if not rows:
            return
        keys = []
        for r in rows:
            for k in r:
                if k not in keys:
                    keys.append(k)
        with open(os.path.join(args.out_dir, name), "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=keys)
            w.writeheader()
            for r in rows:
                w.writerow(r)
    wcsv("sdvz24_manifest.csv", man)
    wcsv("sdvz24_fits.csv", fit_rows)
    wcsv("sdvz24_sens.csv", sens_rows)
    S["finished"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    S["runtime_min"] = round((time.time() - t_start) / 60.0, 1)
    with open(os.path.join(args.out_dir, "sdvz24_summary.json"), "w") as f:
        json.dump(S, f, indent=1, default=base.jsonable)
    try:
        plot(args.out_dir, sens_rows, S)
    except Exception as ex:  # noqa: BLE001
        log("plot failed", repr(ex)[:200])
    return S


def plot(out_dir, sens_rows, S):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    os.makedirs(os.path.join(out_dir, "plots"), exist_ok=True)
    fig, ax = plt.subplots(1, 1, figsize=(8, 4.5))
    if sens_rows:
        d = [r["da"] for r in sens_rows]
        ax.plot(d, [1000.0 * r["x"] for r in sens_rows], "-o", ms=2,
                label="SKS_R minus S_T change (2009 - 2000)")
        ax.plot(d, [1000.0 * (r["S_R_early_lag"] - r["S_T_lag"]) for r in sens_rows], "-",
                label="radial S early minus S_T")
    for v in (-3.0, 0.0):
        ax.axvline(v, color="k", lw=0.6)
    ax.axhline(-160.0, color="r", lw=0.6, ls="--", label="paper: SKS 160 ms earlier")
    ax.set_xlabel("extra rotation applied to the published 2009 R/T (deg)")
    ax.set_ylabel("ms")
    ax.set_title("IU.SDV, Zhang & Wen (2025) published traces; frame: %s" %
                 S.get("frame", {}).get("grade", "?"))
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, "plots", "sdvz24_rotation.png"), dpi=110)


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", default="/data")
    ap.add_argument("--out-dir", default="/out")
    ap.add_argument("--wait-min", type=float, default=120.0)
    a = ap.parse_args(argv)
    run(a)
    return 0


if __name__ == "__main__":
    sys.exit(main())

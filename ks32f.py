#!/usr/bin/env python3
# Copyright 2026 MV2, Munim, Inc.
# Licensed under the Apache License, Version 2.0 (Apache-2.0). See LICENSE.
# Made with AI agents working in MV2's name.
"""ks32f.py - pre-registration 32, amendment 1: the source-fingerprint reading (F).

Plan: findings/2026-11-24-prereg32-amendment-1.md (published before any waveform of
pre-registration 32 was downloaded). The graded T1/T2/joint headline of ks32.py is
unchanged; F is a second, separately labelled reading.

Idea. In one record, the SKS-minus-reference time change between two earthquakes
(B later, A earlier) cancels origin-time and clock errors. What remains is
  r(s, k) = X + [dT_SKS - dT_k](relative source move u = (north, east) km, depth dh km)
where X is a change that only SKS feels (the outer-core term) and the bracket is the
relative-location term, which depends on each phase's horizontal slowness, the
station azimuth and each phase's depth derivative. P-type references (Pdiff, PP)
differ from SKS by about 0.1 s per km of depth; S-type references (Sdiff, SS) by
about 0.03 s per km or less. Fitting X, u and dh together to all passing
(station, reference) pairs separates a core term from a relative mislocation without
trusting any catalogue location or depth.

Input: ks32_phases.csv of job ks-sks32-v1 (columns A, B, sta, band, phase, d, dp, cc,
edge, snr_a, snr_b), committed under results/ks32-v1/. Raw change of a phase =
d + dp (B relative to A, each after its own origin time).

Top-level imports are stdlib + numpy only; ObsPy (TauP) is imported inside run().
"""

import argparse
import csv
import hashlib
import io
import json
import math
import os
import sys
import time
import urllib.request

import numpy as np

VERSION = "0.1.0"
AMEND_PATH = "findings/2026-11-24-prereg32-amendment-1.md"
REPO_RAW = "https://raw.githubusercontent.com/munimtechnologies/inner-core-repeaters/main/"
PHASES_PATH = "results/ks32-v1/ks32_phases.csv"
R_EARTH = 6371.0
KM_PER_DEG = 111.195

# Same USGS values as ks32.EVENTS: (id, origin, lat, lon, depth km)
EV_A = ("E1997", "1997-05-03T16:46:02.02", -31.791, -179.382, 108.3)
EV_B = ("E2018", "2018-09-10T04:19:02.63", -31.7447, -179.3728, 115.0)
# EarthScope station service, 2026-11-24 (station level)
STATIONS = {"COLA": (64.873566, -147.861772), "ULN": (47.8651, 107.0532),
            "INK": (68.3065, -133.5254), "YAK": (62.031, 129.6805)}
# Zhang & Wen (ESSOAr 2023) relocation of the 2018 event relative to 1997, as read
# from their preprint (direction wording ambiguous: "204.33 deg" and "NW").
ZW_RELOC = dict(km=14.20, az_text="204.33 deg (text also says NW)", dh_km=1.45)

TAUP = {"SKS": ["SKS"], "P1": ["P", "Pdiff"], "PP": ["PP"], "S1": ["S", "Sdiff"],
        "SS": ["SS"]}
REFS = ["P1", "PP", "S1", "SS"]
BANDS = {"long": "0.01-0.05", "short": "0.10-0.50", "mid": "0.05-0.20"}
GRADED_BAND = "long"

FPARAMS = dict(
    cc={"long": 0.90, "short": 0.80, "mid": 0.80},  # same gates as ks32
    snr_min=3.0,
    min_rows=7,          # >= 3 degrees of freedom for 4 unknowns
    min_stations=3,
    cond_max=100.0,      # condition number of the column-normalised design
    sigma_max=0.5,       # s, formal 1-sigma of X
    core_thr=-0.5,       # s, same threshold as T1/T2
    null_thr=-0.25,      # s
    ddeg=0.05,           # deg, numerical d/dDelta step
    dh=2.0,              # km, numerical d/dh step
)
F_LABELS = ("core-term", "intermediate", "no-core-term", "untestable")


def log(*a):
    print(time.strftime("%H:%M:%S"), "[ks32f]", *a, flush=True)


# ----------------------------------------------------------------- pure numerics

def fnum(s):
    try:
        return float(s)
    except (TypeError, ValueError):
        return float("nan")


def truthy(s):
    return str(s).strip().lower() in ("true", "1", "yes")


def dist_az(lat1, lon1, lat2, lon2):
    """Spherical epicentral distance (deg) and azimuth (deg, clockwise from north) from
    point 1 to point 2."""
    p1, p2 = math.radians(lat1), math.radians(lat2)
    dl = math.radians(lon2 - lon1)
    c = math.sin(p1) * math.sin(p2) + math.cos(p1) * math.cos(p2) * math.cos(dl)
    d = math.degrees(math.acos(max(-1.0, min(1.0, c))))
    az = math.degrees(math.atan2(math.sin(dl) * math.cos(p2),
                                 math.cos(p1) * math.sin(p2)
                                 - math.sin(p1) * math.cos(p2) * math.cos(dl))) % 360.0
    return d, az


def usgs_vector(a=EV_A, b=EV_B):
    """Relative position of B from A in the USGS catalogue: (north km, east km, dh km)."""
    latm = math.radians(0.5 * (a[2] + b[2]))
    dlon = ((b[3] - a[3] + 180.0) % 360.0) - 180.0
    return ((b[2] - a[2]) * KM_PER_DEG, dlon * KM_PER_DEG * math.cos(latm), b[4] - a[4])


def parse_rows(rows, a, b, band):
    """ks32_phases.csv rows -> sta -> phase -> dict(tot, cc, edge, snr_a, snr_b)."""
    out = {}
    for r in rows:
        if r.get("A") != a or r.get("B") != b or r.get("band") != band:
            continue
        tot = fnum(r.get("d")) + fnum(r.get("dp"))
        out.setdefault(r.get("sta"), {})[r.get("phase")] = dict(
            tot=tot, cc=fnum(r.get("cc")), edge=truthy(r.get("edge")),
            snr_a=fnum(r.get("snr_a")), snr_b=fnum(r.get("snr_b")))
    return out


def ok(m, cc_min, snr_min):
    """Same gate as ks32.passes; NaN never passes."""
    try:
        return bool(np.isfinite(m["tot"]) and np.isfinite(m["cc"]) and not m["edge"]
                    and m["cc"] >= cc_min and m["snr_a"] >= snr_min
                    and m["snr_b"] >= snr_min)
    except (KeyError, TypeError):
        return False


def raw_diffs(parsed, cc_min, snr_min, refs=REFS):
    """[(sta, ref, r)] with r = raw SKS change - raw ref change, both passing."""
    out = []
    for sta in sorted(parsed):
        ph = parsed[sta]
        if "SKS" not in ph or not ok(ph["SKS"], cc_min, snr_min):
            continue
        for k in refs:
            if k in ph and ok(ph[k], cc_min, snr_min):
                out.append((sta, k, float(ph["SKS"]["tot"] - ph[k]["tot"])))
    return out


def design_row(part_sks, part_k, az_deg):
    """part = (dT/dDelta s/deg, dT/dh s/km). Row for unknowns (X, uN, uE, dh):
    dT_phase = -(dT/dDelta in s/rad)/R * (uN cos az + uE sin az) + dT/dh * dh."""
    dp = (part_sks[0] - part_k[0]) * (180.0 / math.pi) / R_EARTH   # s/km
    az = math.radians(az_deg)
    return [1.0, -dp * math.cos(az), -dp * math.sin(az), part_sks[1] - part_k[1]]


def fit(G, y):
    """Least squares y = G m. Returns dict(m, sig, rms, cond, dof, res)."""
    G = np.asarray(G, float).reshape(-1, 4)
    y = np.asarray(y, float).reshape(-1)
    n, k = G.shape
    nan = float("nan")
    if n < k:
        return dict(m=[nan] * k, sig=[nan] * k, rms=nan, cond=float("inf"), dof=n - k,
                    res=[])
    norms = np.sqrt((G ** 2).sum(axis=0))
    if np.any(norms == 0):
        cond = float("inf")
    else:
        cond = float(np.linalg.cond(G / norms))
    m = np.linalg.lstsq(G, y, rcond=None)[0]
    res = y - G @ m
    dof = n - k
    rms = float(math.sqrt(float(np.mean(res ** 2))))
    if dof > 0 and np.isfinite(cond):
        s2 = float(res @ res) / dof
        cov = s2 * np.linalg.pinv(G.T @ G)
        sig = np.sqrt(np.clip(np.diag(cov), 0.0, None))
    else:
        sig = np.full(k, nan)
    return dict(m=[float(v) for v in m], sig=[float(v) for v in sig], rms=rms, cond=cond,
                dof=int(dof), res=[float(v) for v in res])


def grade_f(n, nsta, cond, X, sX, P=FPARAMS):
    """Pre-registered reading F (amendment 1)."""
    for v in (cond, X, sX):
        if v is None or not np.isfinite(v):
            return "untestable"
    if n < P["min_rows"] or nsta < P["min_stations"]:
        return "untestable"
    if cond > P["cond_max"] or sX > P["sigma_max"]:
        return "untestable"
    if X <= P["core_thr"] and X + 2.0 * sX < 0.0:
        return "core-term"
    if X > P["null_thr"]:
        return "no-core-term"
    return "intermediate"


def analyse(rows, parts, azs, band_label, cc_min, P=FPARAMS, a="E1997", b="E2018"):
    """Fit X, uN, uE, dh for one band. parts: sta -> phase -> (dTdD, dTdh); azs: sta -> az."""
    parsed = parse_rows(rows, a, b, band_label)
    diffs = raw_diffs(parsed, cc_min, P["snr_min"])
    G, y, used = [], [], []
    for sta, k, r in diffs:
        if sta not in parts or sta not in azs or "SKS" not in parts[sta] or k not in parts[sta]:
            continue
        row = design_row(parts[sta]["SKS"], parts[sta][k], azs[sta])
        if not all(np.isfinite(v) for v in row) or not np.isfinite(r):
            continue
        G.append(row)
        y.append(r)
        used.append(dict(sta=sta, ref=k, r=r))
    nsta = len({u["sta"] for u in used})
    f = fit(np.asarray(G, float).reshape(-1, 4), np.asarray(y, float))
    X, sX = f["m"][0], f["sig"][0]
    grade = grade_f(len(y), nsta, f["cond"], X, sX, P)
    uN, uE, dh = usgs_vector()
    x_usgs = float("nan")
    if G:
        Ga = np.asarray(G, float)
        pred = Ga[:, 1:] @ np.array([uN, uE, dh])
        corr = np.asarray(y, float) - pred
        x_usgs = float(np.mean(corr))
        for u, p, c in zip(used, pred, corr):
            u["loc_term_usgs"] = float(p)
            u["r_minus_usgs"] = float(c)
    for u, rr in zip(used, f["res"]):
        u["residual"] = rr
    return dict(band=band_label, n=len(y), n_stations=nsta, grade=grade,
                X=X, sigma_X=sX, uN_km=f["m"][1], uE_km=f["m"][2], dh_km=f["m"][3],
                sigmas=f["sig"], rms=f["rms"], cond=f["cond"], dof=f["dof"],
                x_fixed_usgs=x_usgs, used=used)


def first_time(model, depth, dist, names):
    arr = model.get_travel_times(source_depth_in_km=depth, distance_in_degree=dist,
                                 phase_list=list(names))
    ts = [float(x.time) for x in arr if x.name in names]
    return min(ts) if ts else float("nan")


def partials(model, depth, dist, names, ddeg, dh):
    """Central differences of the first arrival among names: (dT/dDelta s/deg, dT/dh s/km)."""
    t1 = first_time(model, depth, dist + ddeg, names)
    t0 = first_time(model, depth, dist - ddeg, names)
    h1 = first_time(model, depth + dh, dist, names)
    h0 = first_time(model, depth - dh, dist, names)
    return ((t1 - t0) / (2.0 * ddeg), (h1 - h0) / (2.0 * dh))


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


# ----------------------------------------------------------------- job side

def fetch(url, tries=3, timeout=60):
    err = ""
    for k in range(tries):
        try:
            with urllib.request.urlopen(url, timeout=timeout) as r:
                return r.read()
        except Exception as ex:  # noqa: BLE001
            err = repr(ex)[:160]
            time.sleep(10 * (k + 1))
    raise RuntimeError("fetch failed %s: %s" % (url, err))


def fetch_wait(url, wait_min):
    deadline = time.time() + wait_min * 60
    while True:
        try:
            return fetch(url, tries=1)
        except Exception:  # noqa: BLE001
            if time.time() > deadline:
                return None
            time.sleep(120)


def run(args):
    from obspy.taup import TauPyModel
    P = FPARAMS
    os.makedirs(args.out_dir, exist_ok=True)
    amend = fetch_wait(REPO_RAW + AMEND_PATH, args.wait_min)
    amend_sha = hashlib.sha256(amend).hexdigest() if amend else ""
    log("amendment public:", bool(amend), amend_sha[:16])
    if args.phases and os.path.exists(args.phases):
        with open(args.phases, "rb") as f:
            raw = f.read()
    else:
        raw = fetch_wait(REPO_RAW + PHASES_PATH, args.wait_min)
    if raw is None:
        raise SystemExit("ks32_phases.csv not available")
    phases_sha = hashlib.sha256(raw).hexdigest()
    rows = list(csv.DictReader(io.StringIO(raw.decode("utf-8"))))
    log("rows", len(rows), "sha", phases_sha[:16])
    model = TauPyModel("iasp91")
    parts, azs, dists = {}, {}, {}
    for sta, (la, lo) in STATIONS.items():
        d, az = dist_az(EV_A[2], EV_A[3], la, lo)
        dists[sta], azs[sta] = d, az
        parts[sta] = {k: partials(model, EV_A[4], d, names, P["ddeg"], P["dh"])
                      for k, names in TAUP.items()}
        log(sta, "dist %.3f az %.2f" % (d, az), parts[sta])
    res = {tag: analyse(rows, parts, azs, BANDS[tag], P["cc"][tag], P) for tag in BANDS}
    for tag, v in res.items():
        log(tag, v["grade"], "X", v["X"], "+/-", v["sigma_X"], "n", v["n"], "cond", v["cond"])
    uN, uE, dh = usgs_vector()
    summary = dict(
        version=VERSION, amendment_public_at_start=bool(amend), amendment_sha256=amend_sha,
        phases_sha256=phases_sha, graded=dict(F=res[GRADED_BAND]["grade"]),
        bands=res, partials=parts, azimuth_deg=azs, distance_deg=dists,
        usgs_relative=dict(uN_km=uN, uE_km=uE, dh_km=dh), zhang_wen_relocation=ZW_RELOC,
        params=P, finished=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()))
    with open(os.path.join(args.out_dir, "ks32f_summary.json"), "w") as f:
        json.dump(summary, f, indent=1, default=jsonable)
    cols = ["band", "sta", "ref", "r", "loc_term_usgs", "r_minus_usgs", "residual"]
    with open(os.path.join(args.out_dir, "ks32f_rows.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        for tag, v in res.items():
            for u in v["used"]:
                w.writerow({c: (tag if c == "band" else u.get(c, "")) for c in cols})
    log("done; F =", res[GRADED_BAND]["grade"])


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-dir", default="/out")
    ap.add_argument("--phases", default="")
    ap.add_argument("--wait-min", type=float, default=30.0)
    run(ap.parse_args(argv))


if __name__ == "__main__":
    sys.exit(main())

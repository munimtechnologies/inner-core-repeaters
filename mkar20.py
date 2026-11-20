#!/usr/bin/env python3
# Copyright 2026 MV2, Munim, Inc.
# Licensed under the Apache License, Version 2.0 (Apache-2.0). See LICENSE.
# Made with AI agents working in MV2's name.
"""mkar20.py - pre-registration 20: does the Chile -> Makanchi (MKAR) inner-core path follow
the backtracking clock?

Plan (written and published before this script fetched any waveform):
findings/2026-11-20-prereg20-mkar-clock-plan.md.

1. Candidate pairs among eligible Chilean events (chilepower.eligible: MKAR era 2000-08-28 to
   2025-08-15, 143-165 deg; here M >= 5.0): epicentres <= 10 km, |dM| <= 0.3, |dz| <= 20 km.
   Families = connected components of the straddle + before pairs.
2. MKAR 1-2 Hz beams (chilepilot.py functions, unchanged) for every family event; core pair
   measurement D (1 - correlation in onset +1..+6 s) after alignment on the direct PKP pulse.
3. Mantle-P source scores at seven IU stations (mantlep3.run_station / source_scores, unchanged)
   for the pairs that pass the MKAR gates.
4. Decision: mantlep3.clock_test(rows, scores, "D") unchanged (pre-registration 13's rule).

Top level needs only the standard library and numpy. The pinned helper modules (chilepower,
chilepilot, mantlep2, mantlep3) and obspy are imported inside functions.

Usage: python mkar20.py --data-dir /data --out-dir /out --wait-min 120 --time-budget-min 300
"""

import argparse
import csv
import hashlib
import json
import math
import os
import random
import sys
import time
import traceback
import urllib.request

import numpy as np

VERSION = "0.1.0"
NAN = float("nan")

RAW = "https://raw.githubusercontent.com/munimtechnologies/inner-core-repeaters/"
PLAN = "findings/2026-11-20-prereg20-mkar-clock-plan.md"
PLAN_TOKEN = b"pre-registration 20"
CATALOGUE = RAW + "fe970e1fbbf8d55bc179cb3d78fa33a29271267b/results/chilepilot-v2/chile_catalogue.csv"
DEPS = (  # (file, commit, line the file must contain)
    ("chilepower.py", "2f803d276539edeb4cbdd6ed9bc2a15914ab8e38", 'VERSION = "0.1.0"'),
    ("chilepilot.py", "3df9ac8dec3a46d75b9be75265f7a9a7938acd7f", 'VERSION = "0.1.0"'),
    ("mantlep2.py", "b4712f0ebbc70b896231620c905ff53a4b8270fe", 'VERSION = "0.13.0"'),
    ("mantlep3.py", "6787707e046d68ebe42f9e1a322d4250b868ed05", 'VERSION = "0.15.0"'),
)

TURN, RATIO = 2008.5, 2.5
PAIR = dict(km=10.0, dm=0.3, dz=20.0, min_mag=5.0)
MKAR = tuple(("KZ", "MK%02d" % i) for i in range(1, 10))
BAND = (1.0, 2.0)

# core pair measurement (seconds relative to ak135 PKIKP, or to the picked onset)
SIG = (-3.0, 7.0)
ONSET_SEARCH = (-3.0, 4.0)
ONSET_FRAC = 0.3
COARSE_LAG = 3.0
DIRECT = (-0.5, 1.0)
FINE_LAG = 0.4
CODA = (1.0, 6.0)
CODA_LATE = (3.0, 10.0)
DT_LAG = 0.3
SNR_MIN = 3.0
CC_DIRECT_MIN = 0.5

# mantle-P stations 45-80 deg from Chile (fallback coordinates; the core is never touched)
MANTLE_STATIONS = (
    ("IU", "ANMO", 34.9459, -106.4572), ("IU", "CCM", 38.0557, -91.2446),
    ("IU", "HRV", 42.5064, -71.5583), ("IU", "SSPA", 40.6358, -77.8876),
    ("IU", "SJG", 18.1091, -66.1500), ("IU", "DWPF", 28.1103, -81.4327),
    ("IU", "TSUM", -19.2022, 17.5838),
)

NEAR, FAR = 2.0, 4.0
PC_MIN, PC_RHO = 15, 0.10
NBOOT = 1000
EPS = 1e-9


def log(*a):
    print(time.strftime("%H:%M:%S"), *a, flush=True)


# ----------------------------------------------------------------------------- pure helpers
def isf(x):
    try:
        return math.isfinite(float(x))
    except (TypeError, ValueError):
        return False


def median(xs):
    v = sorted(float(x) for x in xs if isf(x))
    n = len(v)
    if n == 0:
        return NAN
    return v[n // 2] if n % 2 else 0.5 * (v[n // 2 - 1] + v[n // 2])


def ranks(v):
    """Ranks 1..n, ties averaged."""
    order = sorted(range(len(v)), key=lambda i: v[i])
    r = [0.0] * len(v)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and v[order[j + 1]] == v[order[i]]:
            j += 1
        for k in range(i, j + 1):
            r[order[k]] = 0.5 * (i + j) + 1.0
        i = j + 1
    return r


def spearman(a, b):
    pts = [(float(x), float(y)) for x, y in zip(a, b) if isf(x) and isf(y)]
    if len(pts) < 3:
        return NAN
    ra = ranks([p[0] for p in pts])
    rb = ranks([p[1] for p in pts])
    ma, mb = sum(ra) / len(ra), sum(rb) / len(rb)
    sab = sum((x - ma) * (y - mb) for x, y in zip(ra, rb))
    saa = sum((x - ma) ** 2 for x in ra)
    sbb = sum((y - mb) ** 2 for y in rb)
    if saa <= 0 or sbb <= 0:
        return NAN
    return sab / math.sqrt(saa * sbb)


def pct(xs, q):
    v = sorted(x for x in xs if isf(x))
    if not v:
        return NAN
    return v[int(round(q * (len(v) - 1)))]


def equiv(t):
    """Backtracking-clock time (as mantlep3.equiv)."""
    return t if t <= TURN else TURN - (t - TURN) / RATIO


def d_B(t1, t2):
    return abs(equiv(t2) - equiv(t1))


def d_M(t1, t2):
    return abs(t2 - t1)


def kind_of(t1, t2):
    a, b = t1 < TURN, t2 < TURN
    if a != b:
        return "straddle"
    return "before" if a else "after"


def km_between(lat1, lon1, lat2, lon2):
    p1, p2 = math.radians(lat1), math.radians(lat2)
    dl = math.radians(lon2 - lon1)
    c = math.sin(p1) * math.sin(p2) + math.cos(p1) * math.cos(p2) * math.cos(dl)
    return 6371.0 * math.acos(max(-1.0, min(1.0, c)))


def candidate_pairs(ev, km=PAIR["km"], dm=PAIR["dm"], dz=PAIR["dz"], min_mag=PAIR["min_mag"]):
    """All pairs (earlier event first) of events with mag >= min_mag, epicentres <= km,
    |dM| <= dm, |dz| <= dz (a missing depth passes). ev: dicts with id, t, lat, lon, depth, mag."""
    sel = sorted((r for r in ev if isf(r.get("mag")) and r["mag"] >= min_mag - EPS), key=lambda r: r["t"])
    out = []
    dlat = km / 111.0 + 0.01
    for i in range(len(sel)):
        a = sel[i]
        for j in range(i + 1, len(sel)):
            b = sel[j]
            if abs(a["lat"] - b["lat"]) > dlat or abs(a["mag"] - b["mag"]) > dm + EPS:
                continue
            za, zb = a.get("depth"), b.get("depth")
            if isf(za) and isf(zb) and abs(za - zb) > dz + EPS:
                continue
            d = km_between(a["lat"], a["lon"], b["lat"], b["lon"])
            if d > km + EPS:
                continue
            out.append(dict(ev1=a["id"], ev2=b["id"], t1=a["t"], t2=b["t"], key=tuple(sorted((a["id"], b["id"]))),
                            sep_km=d, kind=kind_of(a["t"], b["t"]), dB=d_B(a["t"], b["t"]), dM=d_M(a["t"], b["t"])))
    return out


def families(pairs, times):
    """Connected components of the straddle + before pairs. Returns ({event: family}, rows) where
    rows are all pairs with both events in the same family, each with 'mult' set. Families are
    named F001, F002, ... in order of their earliest member."""
    parent = {}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for p in pairs:
        if p["kind"] not in ("straddle", "before"):
            continue
        for e in (p["ev1"], p["ev2"]):
            parent.setdefault(e, e)
        ra, rb = find(p["ev1"]), find(p["ev2"])
        if ra != rb:
            parent[rb] = ra
    groups = {}
    for e in parent:
        groups.setdefault(find(e), []).append(e)
    ordered = sorted(groups.values(), key=lambda m: (min(times[e] for e in m), min(m)))
    fam = {}
    for k, members in enumerate(ordered):
        for e in members:
            fam[e] = "F%03d" % (k + 1)
    rows = []
    for p in pairs:
        f1, f2 = fam.get(p["ev1"]), fam.get(p["ev2"])
        if f1 is not None and f1 == f2:
            rows.append(dict(p, mult=f1, m0324=False, cross=(p["kind"] == "straddle")))
    return fam, rows


# ----------------------------------------------------------------------------- signal helpers
def envelope(x, fs, smooth_s=0.5):
    """Smoothed amplitude envelope (analytic signal via FFT, then a moving average)."""
    x = np.asarray(x, dtype=float)
    n = len(x)
    if n == 0:
        return x.copy()
    X = np.fft.fft(x)
    h = np.zeros(n)
    if n % 2 == 0:
        h[0] = h[n // 2] = 1.0
        h[1:n // 2] = 2.0
    else:
        h[0] = 1.0
        h[1:(n + 1) // 2] = 2.0
    a = np.abs(np.fft.ifft(X * h))
    w = max(1, int(round(smooth_s * fs)))
    return np.convolve(a, np.ones(w) / w, mode="same")


def ncc(a, b):
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    if len(a) != len(b) or len(a) < 2:
        return NAN
    a = a - a.mean()
    b = b - b.mean()
    d = math.sqrt(float(np.sum(a * a)) * float(np.sum(b * b)))
    return float(np.sum(a * b)) / d if d > 0 else NAN


def lag_scan(a, b, i0, i1, maxlag, center=0):
    """Best integer lag L in center +/- maxlag maximising ncc(a[i0:i1], b[i0+L:i1+L]).
    Returns (L, cc, cc at L-1, cc at L+1) or (None, nan, nan, nan)."""
    n = len(b)
    ccs = {}
    for L in range(center - maxlag - 1, center + maxlag + 2):
        j0, j1 = i0 + L, i1 + L
        if i0 < 0 or i1 > len(a) or j0 < 0 or j1 > n:
            continue
        ccs[L] = ncc(a[i0:i1], b[j0:j1])
    inside = [L for L in range(center - maxlag, center + maxlag + 1) if isf(ccs.get(L))]
    if not inside:
        return None, NAN, NAN, NAN
    best = max(inside, key=lambda L: ccs[L])
    return best, ccs[best], ccs.get(best - 1, NAN), ccs.get(best + 1, NAN)


def parabolic(cm, c0, cp):
    """Sub-sample offset of a peak from three samples (-1, 0, +1); 0 if not a proper peak."""
    if not (isf(cm) and isf(c0) and isf(cp)):
        return 0.0
    den = cm - 2.0 * c0 + cp
    if den >= 0:
        return 0.0
    return max(-0.5, min(0.5, 0.5 * (cm - cp) / den))


def fshift(x, s):
    """y[n] = x[n + s] for real s (circular, Fourier phase shift)."""
    x = np.asarray(x, dtype=float)
    n = len(x)
    X = np.fft.rfft(x)
    k = np.arange(len(X))
    return np.fft.irfft(X * np.exp(2j * np.pi * k * s / n), n)


def peak_norm(x, i0, i1):
    m = float(np.max(np.abs(x[i0:i1]))) if i1 > i0 else 0.0
    return x / m if m > 0 else x


def pick_onset(stack, fs, c0):
    """First sample in ONSET_SEARCH where the stack envelope reaches ONSET_FRAC of its peak in SIG."""
    env = envelope(stack, fs, 0.25)
    w0, w1 = c0 + int(round(SIG[0] * fs)), c0 + int(round(SIG[1] * fs))
    peak = float(np.max(env[w0:w1]))
    if not peak > 0:
        return None
    s0, s1 = c0 + int(round(ONSET_SEARCH[0] * fs)), c0 + int(round(ONSET_SEARCH[1] * fs))
    for i in range(s0, s1):
        if env[i] >= ONSET_FRAC * peak:
            return i
    return None


def mkar_pair(xa, xb, fs, c0):
    """Align beam xb on beam xa (both with ak135 PKIKP at sample c0) and measure the coda
    difference. Returns a dict or None."""
    if xa is None or xb is None:
        return None
    xa = np.asarray(xa, dtype=float)
    xb = np.asarray(xb, dtype=float)
    i0, i1 = c0 + int(round(SIG[0] * fs)), c0 + int(round(SIG[1] * fs))
    L0, c_env, _, _ = lag_scan(envelope(xa, fs), envelope(xb, fs), i0, i1, int(round(COARSE_LAG * fs)))
    if L0 is None:
        return None
    y0 = fshift(xb, L0)
    onset = pick_onset(peak_norm(xa, i0, i1) + peak_norm(y0, i0, i1), fs, c0)
    if onset is None:
        return None
    j0, j1 = onset + int(round(DIRECT[0] * fs)), onset + int(round(DIRECT[1] * fs))
    L1, c1, cm, cp = lag_scan(xa, y0, j0, j1, int(round(FINE_LAG * fs)))
    if L1 is None:
        return None
    tot = L0 + L1 + parabolic(cm, c1, cp)
    y = fshift(xb, tot)
    out = dict(lag_s=tot / fs, cc_env=c_env, cc_direct=ncc(xa[j0:j1], y[j0:j1]),
               onset_minus_pred_s=(onset - c0) / fs)
    for name, (w0, w1) in (("D", CODA), ("D_c310", CODA_LATE)):
        k0, k1 = onset + int(round(w0 * fs)), onset + int(round(w1 * fs))
        out[name] = 1.0 - ncc(xa[k0:k1], y[k0:k1])
    k0, k1 = onset + int(round(CODA[0] * fs)), onset + int(round(CODA[1] * fs))
    Ld, cd, dm_, dp_ = lag_scan(xa, y, k0, k1, int(round(DT_LAG * fs)))
    out["dt_coda_s"] = (Ld + parabolic(dm_, cd, dp_)) / fs if Ld is not None else NAN
    return out


def gate_pair(snr_a, snr_b, m):
    return bool(m is not None and isf(snr_a) and isf(snr_b) and snr_a >= SNR_MIN and snr_b >= SNR_MIN
                and isf(m.get("cc_direct")) and m["cc_direct"] >= CC_DIRECT_MIN and isf(m.get("D")))


# ----------------------------------------------------------------------------- secondary statistics
def family_boot(rows, stat, nboot=NBOOT, seed=20):
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


def verified(rows, scores, s_max=1.0 / 3.0, min_st=2, key="D"):
    return [r for r in rows if isf(r.get(key)) and r["key"] in scores and scores[r["key"]]["n_st"] >= min_st
            and scores[r["key"]]["S"] <= s_max + 1e-12]


def positive_control(V):
    B = [r for r in V if r["kind"] == "before"]
    rho = spearman([r["dM"] for r in B], [r["D"] for r in B])
    bo = family_boot(B, lambda s: spearman([r["dM"] for r in s], [r["D"] for r in s]))
    p05 = pct(bo, 0.05)
    return dict(n_pairs=len(B), n_families=len({r["mult"] for r in B}), rho=rho, p05=p05, p95=pct(bo, 0.95),
                passed=bool(len(B) >= PC_MIN and isf(rho) and rho >= PC_RHO and isf(p05) and p05 > 0))


def near_far(V):
    S = [r for r in V if r["kind"] == "straddle"]
    nr = [r["D"] for r in S if r["dB"] <= NEAR + EPS]
    fr = [r["D"] for r in S if r["dB"] >= FAR - EPS]
    return dict(n_straddle=len(S), n_near=len(nr), n_far=len(fr), median_D_near=median(nr),
                median_D_far=median(fr))


def dt_stats(V):
    W = [r for r in V if isf(r.get("dt_coda_s"))]
    a = [abs(r["dt_coda_s"]) for r in W]
    return dict(n=len(W), rho_B=spearman([r["dB"] for r in W], a), rho_M=spearman([r["dM"] for r in W], a),
                median_abs_dt_s=median(a))


def stability(ev_rows):
    ok = [r for r in ev_rows if isf(r.get("beam_snr"))]
    allcc = [r["el_beam_cc"] for r in ok if isf(r.get("el_beam_cc"))]
    m = median(allcc)
    rsd = 1.4826 * median([abs(x - m) for x in allcc]) if allcc else NAN
    years = {}
    for r in ok:
        years.setdefault(int(math.floor(r["t"])), []).append(r)
    table, flags, prev = [], [], None
    for y in sorted(years):
        rs = years[y]
        med_cc = median(r.get("el_beam_cc") for r in rs)
        table.append(dict(year=y, n=len(rs), median_snr=median(r["beam_snr"] for r in rs), median_el_beam_cc=med_cc))
        if prev is not None and isf(prev) and isf(med_cc) and isf(rsd) and rsd > 0 and abs(med_cc - prev) > 3 * rsd:
            flags.append(y)
        prev = med_cc
    return dict(robust_sd_el_beam_cc=rsd, years=table, step_flags=flags)


# ----------------------------------------------------------------------------- I/O
def fetch_bytes(url, tries=6, wait=20, timeout=90):
    last = None
    for k in range(tries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "inner-core-repeaters (MV2, Munim, Inc.)"})
            with urllib.request.urlopen(req, timeout=timeout) as r:
                return r.read()
        except Exception as ex:  # noqa: BLE001
            last = ex
            time.sleep(wait * (k + 1))
    raise RuntimeError("fetch failed %s: %r" % (url, last))


def ensure_deps(here):
    info = {}
    for name, sha, must in DEPS:
        p = os.path.join(here, name)
        if not os.path.exists(p):
            body = fetch_bytes(RAW + sha + "/" + name)
            if must.encode() not in body:
                raise RuntimeError("unexpected content in %s at %s" % (name, sha))
            with open(p, "wb") as f:
                f.write(body)
        with open(p, "rb") as f:
            info[name] = dict(commit=sha, sha256=hashlib.sha256(f.read()).hexdigest())
    return info


def plan_public(wait_min):
    t0 = time.time()
    while True:
        try:
            body = fetch_bytes(RAW + "main/" + PLAN, tries=1, wait=0, timeout=60)
            if PLAN_TOKEN in body.lower():
                return True, round((time.time() - t0) / 60.0, 1)
        except Exception:  # noqa: BLE001
            pass
        if time.time() - t0 > wait_min * 60:
            return False, round((time.time() - t0) / 60.0, 1)
        time.sleep(120)


def jsonable(o):
    if isinstance(o, dict):
        return {str(k): jsonable(v) for k, v in o.items()}
    if isinstance(o, (list, tuple, set)):
        return [jsonable(v) for v in o]
    if isinstance(o, (np.floating, np.integer)):
        o = o.item()
    if isinstance(o, float) and not math.isfinite(o):
        return None
    return o


def write_csv(path, rows):
    keys = []
    for r in rows:
        for k in r:
            if k not in keys:
                keys.append(k)
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        for r in rows:
            w.writerow({k: (round(v, 5) if isinstance(v, float) and math.isfinite(v) else
                            ("" if isinstance(v, float) else (";".join(v) if isinstance(v, tuple) else v)))
                        for k, v in r.items()})


# ----------------------------------------------------------------------------- MKAR beams
def mkar_beam(cp, client, model, ev, invs, data_dir):
    """MKAR 1-2 Hz beam with chilepilot's functions. Returns (beam or None, row)."""
    from obspy import UTCDateTime

    row = dict(id=ev["id"], time=ev["time"], t=ev["t"], mag=ev["mag"], depth=ev.get("depth"),
               lat=ev["lat"], lon=ev["lon"])
    els = []
    for net, sta in MKAR:
        inv, chans = invs.get((net, sta), (None, []))
        if inv is None:
            continue
        c = cp.choose_chan(chans, ev["time"])
        if c is None:
            continue
        els.append(dict(net=net, sta=sta, cha=c[0], loc=c[1], lat=c[2], lon=c[3], inv=inv))
    row["n_el_epoch"] = len(els)
    if len(els) < cp.MIN_EL:
        row["status"] = "few-epochs"
        return None, row
    z = max(ev["depth"], 0.0) if isf(ev.get("depth")) else 10.0
    for e in els:
        dist = cp.gc_deg(ev["lat"], ev["lon"], e["lat"], e["lon"])
        arr = model.get_travel_times(source_depth_in_km=z, distance_in_degree=dist, phase_list=["PKIKP"])
        e["tdf"] = min((a.time for a in arr if a.name == "PKIKP"), default=NAN)
    els = [e for e in els if isf(e["tdf"])]
    if len(els) < cp.MIN_EL:
        row["status"] = "no-phase"
        return None, row
    o = UTCDateTime(ev["time"])
    t1 = o + min(e["tdf"] for e in els) - cp.PRE - cp.MARGIN
    t2 = o + max(e["tdf"] for e in els) + cp.POST + cp.MARGIN
    st, how = cp.fetch_array(client, "MKAR", [(e["net"], e["sta"], e["loc"], e["cha"], t1, t2) for e in els],
                             ev, data_dir)
    row["fetch"] = how
    if st is None:
        row["status"] = how
        return None, row
    xs, whys = [], {}
    for e in els:
        sub = st.select(network=e["net"], station=e["sta"], channel=e["cha"])
        sub = type(sub)([tr for tr in sub if (tr.stats.location or "") == e["loc"]])
        if len(sub) == 0:
            whys["missing"] = whys.get("missing", 0) + 1
            continue
        res, why = cp.process(sub, e["inv"], ev, e["tdf"], BAND)
        if res is None:
            whys[why[:20]] = whys.get(why[:20], 0) + 1
            continue
        xs.append(res[0])
    row["n_el_ok"] = len(xs)
    row["drop_reasons"] = ";".join("%s:%d" % kv for kv in sorted(whys.items()))
    b = cp.beam(xs)
    if b is None:
        row["status"] = "few-elements"
        return None, row
    row["beam_snr"] = cp.window_snrs(b)[0]
    i0, i1 = int(round(cp.idx(SIG[0]))), int(round(cp.idx(SIG[1])))
    row["el_beam_cc"] = median(ncc(x[i0:i1], b[i0:i1]) for x in xs)
    row["status"] = "ok"
    return b, row


# ----------------------------------------------------------------------------- plot
def plot(path, V):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(12, 5), sharey=True)
    colours = dict(straddle="C0", before="C1", after="C2")
    for ax, key, lab in ((axes[0], "dB", "backtracking-clock distance dB (yr)"),
                         (axes[1], "dM", "calendar distance dM (yr)")):
        for k, c in colours.items():
            pts = [(r[key], r["D"]) for r in V if r["kind"] == k]
            if pts:
                ax.scatter([p[0] for p in pts], [p[1] for p in pts], s=14, c=c, label=k, alpha=0.8)
        ax.set_xlabel(lab)
    axes[0].set_ylabel("MKAR coda difference D (1-6 s, 1-2 Hz)")
    axes[1].legend(fontsize=8)
    fig.suptitle("Pre-registration 20: source-verified Chile -> MKAR pairs")
    fig.tight_layout()
    fig.savefig(path, dpi=110)
    plt.close(fig)


# ----------------------------------------------------------------------------- main
def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", default="/data")
    ap.add_argument("--out-dir", default="/out")
    ap.add_argument("--wait-min", type=float, default=0.0)
    ap.add_argument("--time-budget-min", type=float, default=300.0)
    a = ap.parse_args(argv)
    t0 = time.time()
    deadline = t0 + 60.0 * a.time_budget_min
    os.makedirs(a.out_dir, exist_ok=True)
    os.makedirs(os.path.join(a.out_dir, "plots"), exist_ok=True)
    os.makedirs(a.data_dir, exist_ok=True)
    here = os.path.dirname(os.path.abspath(__file__))
    if here not in sys.path:
        sys.path.insert(0, here)
    s = dict(version=VERSION, started=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()), errors={},
             constants=dict(turn=TURN, ratio=RATIO, pair=PAIR, band=BAND, sig=SIG, onset_search=ONSET_SEARCH,
                            onset_frac=ONSET_FRAC, coarse_lag=COARSE_LAG, direct=DIRECT, fine_lag=FINE_LAG,
                            coda=CODA, coda_late=CODA_LATE, dt_lag=DT_LAG, snr_min=SNR_MIN,
                            cc_direct_min=CC_DIRECT_MIN, near=NEAR, far=FAR, pc_min=PC_MIN, pc_rho=PC_RHO,
                            mantle_stations=["%s.%s" % (n, t) for n, t, _, _ in MANTLE_STATIONS]))
    s["deps"] = ensure_deps(here)
    s["plan_public_at_start"], s["plan_wait_min"] = plan_public(a.wait_min)
    log("plan public:", s["plan_public_at_start"], s["plan_wait_min"])

    import chilepilot as cp
    import chilepower as cpow
    import mantlep3 as m3

    cat_path = os.path.join(a.out_dir, "chile_catalogue.csv")
    with open(cat_path, "wb") as f:
        f.write(fetch_bytes(CATALOGUE))
    rows = cpow.read_catalogue(cat_path)
    ev = [r for r in cpow.eligible(rows) if r["mag"] >= PAIR["min_mag"] - EPS]
    s["n_catalogue"], s["n_eligible_m5"] = len(rows), len(ev)
    pairs = candidate_pairs(ev)
    s["candidate_pairs"] = {k: sum(1 for p in pairs if p["kind"] == k) for k in ("straddle", "before", "after")}
    times = {r["id"]: r["t"] for r in ev}
    fam, frows = families(pairs, times)
    evs = {r["id"]: r for r in ev if r["id"] in fam}
    sizes = {}
    for e, f in fam.items():
        sizes[f] = sizes.get(f, 0) + 1
    s["n_events_E"], s["n_families"] = len(evs), len(sizes)
    s["family_sizes_top10"] = sorted(sizes.values(), reverse=True)[:10]
    s["family_rows"] = {k: sum(1 for p in frows if p["kind"] == k) for k in ("straddle", "before", "after")}
    log("eligible M5", len(ev), "pairs", s["candidate_pairs"], "E", len(evs), "families", len(sizes),
        "rows", s["family_rows"])

    from obspy.clients.fdsn import Client
    from obspy.taup import TauPyModel

    client = Client("https://service.earthscope.org", timeout=180, _discover_services=False)
    model = TauPyModel("ak135")
    invs = {}
    for net, sta in MKAR:
        inv = cp.get_inventory(client, net, sta, a.data_dir)
        invs[(net, sta)] = (inv, cp.chans_of(inv) if inv is not None else [])
    s["mkar_channel_epochs"] = {"%s.%s" % k: len(v[1]) for k, v in invs.items()}

    beams, ev_rows = {}, []
    mkar_deadline = t0 + 0.5 * 60.0 * a.time_budget_min
    for k, e in enumerate(sorted(evs.values(), key=lambda r: r["t"])):
        if time.time() > mkar_deadline:
            ev_rows.append(dict(id=e["id"], time=e["time"], t=e["t"], status="skipped-time-budget"))
            continue
        try:
            b, row = mkar_beam(cp, client, model, e, invs, a.data_dir)
        except Exception as ex:  # noqa: BLE001
            b, row = None, dict(id=e["id"], time=e["time"], t=e["t"], status="error:" + repr(ex)[:100])
            s["errors"]["mkar " + e["id"]] = traceback.format_exc()[-400:]
        row["family"] = fam.get(e["id"], "")
        ev_rows.append(row)
        if b is not None:
            beams[e["id"]] = b
        if k % 25 == 0:
            log("mkar", k + 1, "/", len(evs), e["id"], row.get("status"), row.get("beam_snr"))
            write_csv(os.path.join(a.out_dir, "mkar20_events.csv"), ev_rows)
    write_csv(os.path.join(a.out_dir, "mkar20_events.csv"), ev_rows)
    snr = {r["id"]: r.get("beam_snr", NAN) for r in ev_rows}
    st_counts = {}
    for r in ev_rows:
        st_counts[str(r.get("status"))[:30]] = st_counts.get(str(r.get("status"))[:30], 0) + 1
    s["mkar_status_counts"] = st_counts
    n_snr = sum(1 for r in ev_rows if isf(r.get("beam_snr")) and r["beam_snr"] >= SNR_MIN)
    s["mkar_events_snr_ge_min"] = n_snr

    c0 = int(round(cp.PRE * cp.FS))
    for r in frows:
        m = mkar_pair(beams.get(r["ev1"]), beams.get(r["ev2"]), cp.FS, c0)
        if m:
            r.update(m)
        r["snr1"], r["snr2"] = snr.get(r["ev1"], NAN), snr.get(r["ev2"], NAN)
        r["mkar_ok"] = gate_pair(r["snr1"], r["snr2"], m)
        if not r["mkar_ok"]:
            for k in ("D", "D_c310", "dt_coda_s"):
                r["raw_" + k] = r.get(k, NAN)
                r[k] = NAN
    ok = [r for r in frows if r["mkar_ok"]]
    s["mkar_pairs_ok"] = {k: sum(1 for p in ok if p["kind"] == k) for k in ("straddle", "before", "after")}
    log("mkar pairs ok", s["mkar_pairs_ok"])

    import requests
    session = requests.Session()
    session.headers["User-Agent"] = "inner-core-repeaters (MV2, Munim, Inc.)"
    mfams = {}
    for r in ok:
        f = mfams.setdefault(r["mult"], dict(members=[], need=set()))
        for e in (r["ev1"], r["ev2"]):
            if e not in f["members"]:
                f["members"].append(e)
        f["need"].add(r["key"])
    mevents = {}
    for f in mfams.values():
        for e in f["members"]:
            E = evs[e]
            mevents[e] = dict(id=e, time=E["time"], lat=E["lat"], lon=E["lon"],
                              depth=E["depth"] if isf(E.get("depth")) else 33.0)
    cache = os.path.join(a.data_dir, "mantlep")
    os.makedirs(cache, exist_ok=True)
    st_rows, mev_rows, mfam_rows, mpair_rows = [], [], [], []
    for net, sta, lat, lon in MANTLE_STATIONS:
        if time.time() > deadline:
            st_rows.append(dict(station="%s.%s" % (net, sta), skipped_time_budget=True))
            continue
        try:
            st, er, fr, pr = m3.run_station(session, model, net, sta, lat, lon, mevents, mfams, cache, deadline)
        except Exception as ex:  # noqa: BLE001
            s["errors"]["mantle %s.%s" % (net, sta)] = traceback.format_exc()[-400:]
            st, er, fr, pr = dict(station="%s.%s" % (net, sta), error=repr(ex)[:300]), [], [], []
        st_rows.append(st)
        mev_rows += er
        mfam_rows += fr
        mpair_rows += pr
        log("mantle", st)
    s["mantle_stations"] = st_rows
    write_csv(os.path.join(a.out_dir, "mkar20_mantle_pairs.csv"), mpair_rows)
    write_csv(os.path.join(a.out_dir, "mkar20_mantle_events.csv"), mev_rows)
    scores, used = m3.source_scores(mpair_rows)
    s["stations_scored"] = used
    write_csv(os.path.join(a.out_dir, "mkar20_scores.csv"),
              [dict(ev1=k[0], ev2=k[1], **v) for k, v in sorted(scores.items())])
    for r in frows:
        sc = scores.get(r["key"])
        r["S"] = sc["S"] if sc else NAN
        r["n_st"] = sc["n_st"] if sc else 0

    g = dict(primary=m3.clock_test(frows, scores, "D"))
    g["decision"] = g["primary"]["decision"]
    g["secondary_D_c310"] = m3.clock_test(frows, scores, "D_c310")
    g["sensitivity_S_le_0.25"] = m3.clock_test(frows, scores, "D", s_max=0.25)
    g["sensitivity_S_le_0.5"] = m3.clock_test(frows, scores, "D", s_max=0.5)
    V = verified(frows, scores)
    g["positive_control"] = positive_control(V)
    g["near_far"] = near_far(V)
    g["dt_coda"] = dt_stats(V)
    s["grade"] = g
    s["stability"] = stability(ev_rows)
    s["forecasts"] = dict(
        H1_frac_E_snr_ge_3=(n_snr / float(len(evs))) if evs else NAN,
        H2_straddle_pairs_mkar_ok=s["mkar_pairs_ok"]["straddle"],
        H3_testable=g["decision"] not in ("untestable", "untestable-calibration"),
        H4_decision=g["decision"],
        H5_positive_control=g["positive_control"]["passed"])
    write_csv(os.path.join(a.out_dir, "mkar20_pairs.csv"), frows)
    s["finished"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    with open(os.path.join(a.out_dir, "mkar20_summary.json"), "w", encoding="utf-8") as f:
        json.dump(jsonable(s), f, indent=1)
    try:
        plot(os.path.join(a.out_dir, "plots", "mkar20_clock.png"), V)
    except Exception as ex:  # noqa: BLE001
        log("plot failed", repr(ex)[:160])
    log("decision", g["decision"], json.dumps(jsonable(s["forecasts"])))
    return 0


if __name__ == "__main__":
    sys.exit(main())

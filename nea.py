#!/usr/bin/env python3
"""nea.py: does the 2024 inner-core "backtracking" clock show up on a second, independent
inner-core path? South Sandwich Islands (SSI) repeating earthquakes recorded near their
antipode, at long-running stations in north-east Asia (pre-registration 16).

Part of inner-core-repeaters by MV2, Munim, Inc. (made with AI agents). Licence: Apache-2.0.
Pre-registered in findings/2026-11-19-prereg16-antipodal-clock.md, written and published
before any waveform from these stations was fetched by this project.

Idea. If the inner core turns rigidly and reversed at 2008.5, retracing its path 2.5x more
slowly (Wang, Vidale et al. 2024, Nature 631:340), then two events at the same "equivalent
date" see the inner core in the same orientation. Their PKIKP (DF) waves should then look
alike on EVERY path through a laterally varying inner core, not only on the Alaska/Canada paths
where the clock was found. A local change of the inner-core surface (Wen 2006; Vidale et al.
2025) predicts no such shared clock on a path that crosses the boundary elsewhere. SSI -> NE
Asia at 155-179.5 deg is such a path: DF dives deep into the inner core and crosses the
boundary far from the Alaska paths' entry and exit points.

Source control on the same trace: PKPab travels through the mantle and outer core only. Pairs
are aligned on PKPab and must be similar there (cc_ab >= 0.8); the same clock test run on the
PKPab dissimilarity flags a source or instrument artefact.

Families: ILAR multiplets Z000-Z019 from results/m0324-v1 (pinned commit), i.e. SSI events
already known to repeat. No family membership is decided at the new stations.
"""
import argparse
import csv
import io
import json
import math
import os
import sys
import time
import traceback
import urllib.request

import numpy as np

VERSION = "0.1.0"
NAN = float("nan")

# ----------------------------------------------------------------------------- constants
REPO = "https://raw.githubusercontent.com/munimtechnologies/inner-core-repeaters/"
INPUT_COMMIT = "d0186eece9f9bae982a55a3302a0773cf8627d98"   # v0.14.6
PREREG_URL = REPO + "main/findings/2026-11-19-prereg16-antipodal-clock.md"
CENTROID = (-59.0365, -25.6939)          # results/m0324-v1 centroid
ANTIPODE = (-CENTROID[0], CENTROID[1] + 180.0 if CENTROID[1] < 0 else CENTROID[1] - 180.0)
MIN_DIST, MAX_DIST = 155.0, 179.5        # deg from the SSI centroid
NETWORKS = "IU,II,IC,G,GE,GT,CU"
MAX_STATIONS = 10
MIN_COVER = 20                           # multiplet events inside a BHZ epoch
LOC_PREF = ("00", "", "10", "20", "60")

FS = 20.0
PRE = 60.0                               # s kept before ak135 PKIKP
POST = 90.0                              # s kept after ak135 PKIKP
N = int(round((PRE + POST) * FS))
BAND = (0.6, 1.6)                        # primary band (Hz)
BAND_ALT = (1.0, 2.0)                    # secondary band
PREFILT = (0.1, 0.2, 7.0, 9.0)
NOISE = (-45.0, -10.0)                   # s, relative to ak135 PKIKP
DF_WIN = (-0.3, 5.0)                     # s, relative to the picked DF onset (primary D)
DF_CODA = (1.0, 5.0)                     # secondary D_c15
AB_WIN = (-0.5, 6.0)                     # s, relative to the picked AB onset (reference)
MIN_SEP = 8.0                            # AB - DF needed (s) so the windows do not overlap

# pre-registered gates and decision constants (pre-registration 16)
MIN_CC_AB = 0.80
MIN_SNR_DF = 3.0
MIN_SNR_AB = 4.0
MIN_ALIGN_CC = 0.30
MAX_LAG_COARSE_S = 3.0
MAX_LAG_REF_S = 0.5
MAX_LAG_DF_S = 0.4
MIN_BLOCK = 4                            # ok pairs per (station, multiplet) block
MIN_PAIRS, MIN_CROSS, MIN_CLUSTERS = 30, 5, 4
NEAR = 3.0                               # calibration: |dt| <= 3 and |d_equiv| <= 3
FAR_M, FAR_B = 10.0, 5.0                 # calibration: |dt| >= 10 and |d_equiv| >= 5
MIN_CAL = 10
CAL_RATIO = 1.25
FRAC = 0.9
NBOOT = 2000
TURN, BACK_RATIO = 2008.5, 2.5
NEW_FROM = 2019.0                        # secondary new-event best-match test


def log(*a):
    print(time.strftime("[%H:%M:%S]"), *a, flush=True)


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


def median(xs):
    v = sorted(float(x) for x in xs if isf(x))
    n = len(v)
    if not n:
        return NAN
    return v[n // 2] if n % 2 else 0.5 * (v[n // 2 - 1] + v[n // 2])


def equiv(t, T=TURN, ratio=BACK_RATIO):
    t = float(t)
    return t if t <= T else T - (t - T) / ratio


def d_b(t1, t2, T=TURN, ratio=BACK_RATIO):
    return abs(equiv(t1, T, ratio) - equiv(t2, T, ratio))


def d_m(t1, t2):
    return abs(float(t2) - float(t1))


def gc_deg(lat1, lon1, lat2, lon2):
    """Great-circle distance in degrees on a sphere (haversine)."""
    p1, p2 = math.radians(lat1), math.radians(lat2)
    dp, dl = p2 - p1, math.radians(lon2 - lon1)
    a = math.sin(dp / 2) ** 2 + math.cos(p1) * math.cos(p2) * math.sin(dl / 2) ** 2
    return math.degrees(2.0 * math.asin(min(1.0, math.sqrt(a))))


def decimal_year(iso):
    from datetime import datetime, timezone

    s = str(iso).strip().replace("Z", "+00:00")
    dt = datetime.fromisoformat(s)
    if dt.tzinfo is None:
        dt = dt.replace(tzinfo=timezone.utc)
    y0 = datetime(dt.year, 1, 1, tzinfo=timezone.utc)
    y1 = datetime(dt.year + 1, 1, 1, tzinfo=timezone.utc)
    return dt.year + (dt - y0).total_seconds() / (y1 - y0).total_seconds()


def ranks(v):
    idx = sorted(range(len(v)), key=lambda i: v[i])
    r = [0.0] * len(v)
    i = 0
    while i < len(idx):
        j = i
        while j + 1 < len(idx) and v[idx[j + 1]] == v[idx[i]]:
            j += 1
        for k in range(i, j + 1):
            r[idx[k]] = 0.5 * (i + j)
        i = j + 1
    return r


def spearman(a, b):
    pr = [(x, y) for x, y in zip(a, b) if isf(x) and isf(y)]
    if len(pr) < 3:
        return NAN
    ra, rb = ranks([p[0] for p in pr]), ranks([p[1] for p in pr])
    ma, mb = sum(ra) / len(ra), sum(rb) / len(rb)
    sab = sum((x - ma) * (y - mb) for x, y in zip(ra, rb))
    saa = sum((x - ma) ** 2 for x in ra)
    sbb = sum((y - mb) ** 2 for y in rb)
    return sab / math.sqrt(saa * sbb) if saa > 0 and sbb > 0 else NAN


def parse_members(s):
    """'id(date);id(date)' -> [id, ...]"""
    out = []
    for part in str(s).split(";"):
        part = part.strip()
        if part:
            out.append(part.split("(")[0])
    return out


def choose_loc(chans, iso, pref=LOC_PREF):
    """chans: [(loc, start_iso, end_iso)] of one station's BHZ epochs. The preferred location
    code whose epoch covers `iso` (ISO strings compare correctly), or None."""
    s = str(iso)[:19]
    cover = {loc for loc, a, b in chans if str(a)[:19] <= s < (str(b)[:19] if b else "9999")}
    for loc in pref:
        if loc in cover:
            return loc
    return sorted(cover)[0] if cover else None


def block_normalise(rows, key, min_block=MIN_BLOCK):
    """Rows with finite `key`, grouped by (station, mult). Within each block of >= min_block rows,
    u = rank / (n - 1) in [0, 1] (0 = most similar pair of the block). Returns new rows."""
    blocks = {}
    for r in rows:
        if isf(r.get(key)):
            blocks.setdefault((r["station"], r["mult"]), []).append(r)
    out = []
    for b, rs in sorted(blocks.items()):
        if len(rs) < min_block:
            continue
        rk = ranks([float(r[key]) for r in rs])
        n = len(rs)
        for r, k in zip(rs, rk):
            out.append(dict(r, u=k / (n - 1.0)))
    return out


def _rhos(rows):
    u = [r["u"] for r in rows]
    gb = [d_b(r["t1"], r["t2"]) for r in rows]
    gm = [d_m(r["t1"], r["t2"]) for r in rows]
    return spearman(gb, u), spearman(gm, u)


def _cluster_boot(rows, fn, nboot=NBOOT, seed=0):
    """Resample multiplets (clusters across stations) with replacement; fn(rows) per resample."""
    by = {}
    for r in rows:
        by.setdefault(r["mult"], []).append(r)
    keys = sorted(by)
    rng = np.random.default_rng(seed)
    out = []
    for _ in range(nboot):
        pick = rng.integers(0, len(keys), len(keys))
        rs = [r for k in pick for r in by[keys[k]]]
        out.append(fn(rs))
    return out


def clock_test(rows, key, nboot=NBOOT, seed=0):
    """P16-1: within-block ranks of `key` against the backtracking distance (rho_B) and the plain
    time distance (rho_M), pooled; multiplet bootstrap of rho_B - rho_M."""
    nr = block_normalise(rows, key)
    n_cross = sum(1 for r in nr if min(r["t1"], r["t2"]) < TURN < max(r["t1"], r["t2"]))
    clusters = sorted({r["mult"] for r in nr})
    res = dict(key=key, n_pairs=len(nr), n_cross_turn=n_cross, n_blocks=len({(r["station"], r["mult"]) for r in nr}),
               n_clusters=len(clusters), rho_B=None, rho_M=None, frac_B_better=None, label="untestable")
    if len(nr) < MIN_PAIRS or n_cross < MIN_CROSS or len(clusters) < MIN_CLUSTERS:
        return res
    rb, rm = _rhos(nr)

    def diff(rs):
        a, b = _rhos(rs)
        return (a - b) if isf(a) and isf(b) else NAN

    ds = [d for d in _cluster_boot(nr, diff, nboot, seed) if isf(d)]
    frac = float(np.mean([d > 0 for d in ds])) if ds else NAN
    res.update(rho_B=rb, rho_M=rm, frac_B_better=frac,
               diff_p05=float(np.percentile(ds, 5)) if ds else None,
               diff_p95=float(np.percentile(ds, 95)) if ds else None)
    if isf(frac) and frac >= FRAC and isf(rb) and rb > 0:
        res["label"] = "B"
    elif isf(frac) and frac <= 1.0 - FRAC and isf(rm) and rm > 0:
        res["label"] = "M"
    else:
        res["label"] = "inconclusive"
    return res


def calibration(rows, key, nboot=NBOOT, seed=1):
    """P16-0: is there a time signal at all? Uses only pairs on which both clocks agree:
    near = |dt| <= 3 and |d_equiv| <= 3; far = |dt| >= 10 and |d_equiv| >= 5. Holds if the median
    within-block rank of far pairs is >= 1.25x that of near pairs and far > near in >= 90% of
    multiplet bootstraps."""
    nr = block_normalise(rows, key)

    def split(rs):
        near = [r["u"] for r in rs if d_m(r["t1"], r["t2"]) <= NEAR and d_b(r["t1"], r["t2"]) <= NEAR]
        far = [r["u"] for r in rs if d_m(r["t1"], r["t2"]) >= FAR_M and d_b(r["t1"], r["t2"]) >= FAR_B]
        return near, far

    near, far = split(nr)
    res = dict(key=key, n_near=len(near), n_far=len(far), median_u_near=median(near), median_u_far=median(far),
               ratio=None, frac_far_gt_near=None, testable=False, holds=None)
    if len(near) < MIN_CAL or len(far) < MIN_CAL:
        return res
    mn, mf = median(near), median(far)
    if mn > 0:
        ratio = mf / mn
    elif mf > 0:
        ratio = float("inf")
    else:
        ratio = NAN

    def gt(rs):
        a, b = split(rs)
        if not a or not b:
            return NAN
        return float(median(b) > median(a))

    bs = [v for v in _cluster_boot(nr, gt, nboot, seed) if isf(v)]
    frac = float(np.mean(bs)) if bs else NAN
    ratio_ok = ratio == float("inf") or (isf(ratio) and ratio >= CAL_RATIO)
    res.update(ratio=ratio if isf(ratio) else (None if ratio != float("inf") else "inf"), frac_far_gt_near=frac,
               testable=True, holds=bool(ratio_ok and isf(frac) and frac >= FRAC))
    return res


def headline(cal, core, src):
    """Pre-registered headline of pre-registration 16."""
    if not cal.get("testable"):
        return "untestable"
    if not cal.get("holds"):
        return "no-time-signal"
    lab = core.get("label")
    if lab in ("B", "M"):
        if src.get("label") == lab:
            return "source-confounded:" + lab
        return lab + "-supported"
    return lab or "untestable"


def new_event_matches(rows, key, new_from=NEW_FROM):
    """S16-b: for each (station, event >= new_from) with >= 3 ok earlier partners, one of which is
    before TURN and one after it: does the best match (min key) equal the backtracking pick
    (min |d_equiv|) or the time-order pick (min |dt|)?"""
    by = {}
    for r in rows:
        if not isf(r.get(key)):
            continue
        if r["t2"] >= new_from and r["t1"] < r["t2"]:
            by.setdefault((r["station"], r["ev2"]), []).append(r)
    cases = []
    for (sta, ev), rs in sorted(by.items()):
        if len(rs) < 3 or not any(r["t1"] < TURN for r in rs) or not any(r["t1"] > TURN for r in rs):
            continue
        best = min(rs, key=lambda r: float(r[key]))
        bp = min(rs, key=lambda r: d_b(r["t1"], r["t2"]))
        mp = min(rs, key=lambda r: d_m(r["t1"], r["t2"]))
        cases.append(dict(station=sta, ev2=ev, t2=rs[0]["t2"], n=len(rs), best=best["ev1"], b_pick=bp["ev1"],
                          m_pick=mp["ev1"], discriminating=bp["ev1"] != mp["ev1"],
                          matches_B=best["ev1"] == bp["ev1"], matches_M=best["ev1"] == mp["ev1"]))
    disc = [c for c in cases if c["discriminating"]]
    return dict(cases=cases, n_cases=len(cases), n_discriminating=len(disc),
                n_match_B=sum(c["matches_B"] for c in disc), n_match_M=sum(c["matches_M"] for c in disc))


def clock_scan(rows, key):
    """Exploratory: best (turn, ratio) by rho on within-block ranks."""
    nr = block_normalise(rows, key)
    if len(nr) < MIN_PAIRS:
        return dict(n=len(nr))
    u = [r["u"] for r in nr]
    best = None
    for T in np.arange(2000.0, 2016.01, 0.5):
        for ratio in (1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0):
            g = [d_b(r["t1"], r["t2"], T, ratio) for r in nr]
            rho = spearman(g, u)
            if isf(rho) and (best is None or rho > best[0]):
                best = (rho, float(T), ratio)
    rm = spearman([d_m(r["t1"], r["t2"]) for r in nr], u)
    return dict(n=len(nr), best_rho=best[0] if best else None, best_T=best[1] if best else None,
                best_ratio=best[2] if best else None, rho_M=rm)


# ----------------------------------------------------------------------------- signal helpers
def _z(x):
    x = np.asarray(x, dtype=float)
    x = x - x.mean()
    n = np.linalg.norm(x)
    return x / n if n > 0 else x * 0.0


def cc(a, b):
    a, b = np.asarray(a, dtype=float), np.asarray(b, dtype=float)
    if len(a) != len(b) or len(a) < 2:
        return NAN
    return float(np.dot(_z(a), _z(b)))


def lag_scan(x, y, i0, i1, max_lag):
    """Correlate x[i0:i1] with y[i0+k:i1+k], k in [-max_lag, max_lag] (y[n+k] ~ x[n])."""
    x, y = np.asarray(x, dtype=float), np.asarray(y, dtype=float)
    i0, i1, max_lag = int(i0), int(i1), int(max_lag)
    lags = np.arange(-max_lag, max_lag + 1)
    ccs = np.full(lags.shape, np.nan)
    L = i1 - i0
    if L < 3 or i0 < 0 or i1 > len(x):
        return lags, ccs
    t = _z(x[i0:i1])
    lo, hi = max(-max_lag, -i0), min(max_lag, len(y) - i1)
    if hi < lo:
        return lags, ccs
    win = np.lib.stride_tricks.sliding_window_view(y[i0 + lo:i1 + hi], L)
    num = win @ t
    s1 = win.sum(axis=1)
    s2 = np.einsum("ij,ij->i", win, win)
    var = s2 - s1 * s1 / L
    with np.errstate(invalid="ignore", divide="ignore"):
        c = num / np.sqrt(np.where(var > 0, var, np.nan))
    ccs[lo + max_lag:hi + max_lag + 1] = c
    return lags, ccs


def best_lag(lags, ccs):
    ccs = np.asarray(ccs, dtype=float)
    if not np.any(np.isfinite(ccs)):
        return NAN, NAN
    k = int(np.nanargmax(ccs))
    c0, frac, cpk = float(ccs[k]), 0.0, float(ccs[k])
    if 0 < k < len(ccs) - 1 and np.isfinite(ccs[k - 1]) and np.isfinite(ccs[k + 1]):
        ym, yp = float(ccs[k - 1]), float(ccs[k + 1])
        den = ym - 2.0 * c0 + yp
        if den < 0:
            frac = max(-0.5, min(0.5, 0.5 * (ym - yp) / den))
            cpk = c0 - 0.25 * (ym - yp) * frac
    return float(lags[k]) + frac, float(min(cpk, 1.0))


def fshift(y, s):
    """Delay y by s samples (positive = later)."""
    y = np.asarray(y, dtype=float)
    n = len(y)
    if s == 0 or n == 0:
        return y.copy()
    Y = np.fft.rfft(y, 2 * n)
    f = np.fft.rfftfreq(2 * n)
    return np.fft.irfft(Y * np.exp(-2j * np.pi * f * s), 2 * n)[:n]


def sta_lta(x, nsta, nlta):
    e = np.asarray(x, dtype=float) ** 2
    n = len(e)
    out = np.zeros(n)
    if n < nsta + nlta + 1 or nsta < 1 or nlta < 1:
        return out
    c = np.concatenate([[0.0], np.cumsum(e)])
    idx = np.arange(nsta + nlta - 1, n)
    s = (c[idx + 1] - c[idx + 1 - nsta]) / nsta
    l = (c[idx + 1 - nsta] - c[idx + 1 - nsta - nlta]) / nlta
    with np.errstate(invalid="ignore", divide="ignore"):
        out[idx] = np.where(l > 0, s / l, 0.0)
    return out


def pick_onset(x, fs, center, search_s, thr, sta_s=0.5, lta_s=8.0):
    """First STA/LTA trigger within center +/- search_s, walked back to where the 0.1-s power
    falls to 3x the pre-trigger noise; None if nothing triggers."""
    x = np.asarray(x, dtype=float)
    nsta, nlta = int(round(sta_s * fs)), int(round(lta_s * fs))
    r = sta_lta(x, nsta, nlta)
    a, b = max(0, int(center - search_s * fs)), min(len(x), int(center + search_s * fs))
    hit = np.nonzero(r[a:b] >= thr)[0]
    if len(hit) == 0:
        return None
    trig = int(a + hit[0])
    e = x ** 2
    n0, n1 = max(0, trig - nsta - nlta), max(1, trig - nsta)
    noise = float(np.mean(e[n0:n1])) if n1 > n0 else 0.0
    w = max(1, int(round(0.1 * fs)))
    c = np.concatenate([[0.0], np.cumsum(e)])

    def env(i):
        lo = max(0, i - w + 1)
        return (c[i + 1] - c[lo]) / (i + 1 - lo)

    i, lim = trig, max(0, trig - nsta - int(round(0.5 * fs)))
    while i - 1 >= lim and env(i - 1) > 3.0 * noise:
        i -= 1
    return int(i)


def _rms(x, a, b):
    a, b = int(a), int(b)
    if a < 0 or b > len(x) or b - a < 2:
        return NAN
    return float(np.sqrt(np.mean(np.asarray(x[a:b], dtype=float) ** 2)))


def measure_family(traces, ab_off, fs=FS, pre=PRE):
    """One ILAR multiplet at one station.

    traces : {event id: noise-normalised trace, ak135 PKIKP at sample pre*fs}
    ab_off : {event id: ak135 PKPab - PKIKP (s)}
    Returns (info, {id: dict(lag, snr_df, snr_ab)}, aligned {id: trace}, onsets (df_i, ab_i)).
    """
    info = dict(n_in=len(traces), status="")
    ids = sorted(traces)
    if len(ids) < 2:
        info["status"] = "too-few"
        return info, {}, {}, None
    sep = median([ab_off[i] for i in ids])
    info["ab_minus_df_pred_s"] = sep
    if not isf(sep) or sep < MIN_SEP:
        info["status"] = "ab-too-close"
        return info, {}, {}, None
    p0 = int(round(pre * fs))
    snr0 = {i: _rms(traces[i], p0 - 2 * fs, p0 + (sep + 8) * fs) for i in ids}
    root = max(ids, key=lambda i: snr0[i] if isf(snr0[i]) else -1)
    i0, i1 = p0 - int(3 * fs), p0 + int((sep + 8) * fs)
    lag, keep = {}, []
    for i in ids:
        if i == root:
            lag[i], keep = 0.0, keep + [i]
            continue
        lg, c = best_lag(*lag_scan(traces[root], traces[i], i0, i1, int(MAX_LAG_COARSE_S * fs)))
        if isf(lg) and isf(c) and c >= MIN_ALIGN_CC:
            lag[i] = lg
            keep.append(i)
    info.update(root=root, n_aligned=len(keep))
    if len(keep) < 2:
        info["status"] = "alignment-failed"
        return info, {}, {}, None
    al = {i: fshift(traces[i], -lag[i]) for i in keep}
    stack = np.mean([al[i] for i in keep], axis=0)
    df_i = pick_onset(stack, fs, center=p0, search_s=6.0, thr=4.0)
    if df_i is None:
        info["status"] = "no-DF-pick"
        return info, {}, {}, None
    ab_pred = p0 + int(round(sep * fs))
    ab_i = pick_onset(stack, fs, center=ab_pred, search_s=4.0, thr=3.0)
    info["ab_picked"] = ab_i is not None
    if ab_i is None:
        ab_i = ab_pred + (df_i - p0)
    if ab_i - df_i < MIN_SEP * fs:
        info["status"] = "ab-too-close"
        return info, {}, {}, None
    info.update(df_pick_minus_pred_s=(df_i - p0) / fs, ab_pick_minus_pred_s=(ab_i - ab_pred) / fs,
                status="ok")
    n0, n1 = p0 + int(NOISE[0] * fs), p0 + int(NOISE[1] * fs)
    per = {}
    for i in keep:
        nz = _rms(al[i], n0, n1)
        sdf = _rms(al[i], df_i + DF_WIN[0] * fs, df_i + DF_WIN[1] * fs) / nz if nz and isf(nz) else NAN
        sab = _rms(al[i], ab_i + AB_WIN[0] * fs, ab_i + AB_WIN[1] * fs) / nz if nz and isf(nz) else NAN
        per[i] = dict(lag_s=lag[i] / fs, snr_df=sdf, snr_ab=sab)
    return info, per, al, (df_i, ab_i)


def measure_pair(xa, xb, df_i, ab_i, fs=FS, dfw=DF_WIN):
    """Align xb on xa in the PKPab window (+/- 0.5 s), then in the DF window: zero-lag and best
    (+/- 0.4 s) correlation and the residual DF delay. Returns dict or None."""
    r0, r1 = ab_i + int(round(AB_WIN[0] * fs)), ab_i + int(round(AB_WIN[1] * fs))
    l_ref, cc_ref = best_lag(*lag_scan(xa, xb, r0, r1, int(round(MAX_LAG_REF_S * fs))))
    if not isf(l_ref):
        return None
    yb = fshift(xb, -l_ref)
    a, b = df_i + int(round(dfw[0] * fs)), df_i + int(round(dfw[1] * fs))
    if a < 0 or b > len(xa) or b - a < 3:
        return None
    cc0 = cc(xa[a:b], yb[a:b])
    l_df, cc_b = best_lag(*lag_scan(xa, yb, a, b, int(round(MAX_LAG_DF_S * fs))))
    return dict(lag_ref_s=l_ref / fs, cc_ab=cc_ref, cc_df0=cc0, cc_df_best=cc_b,
                dt_df_s=l_df / fs if isf(l_df) else NAN)


def family_pairs(station, mult, evs, al, per, onsets, fs=FS):
    """All pairs of one family at one station (earlier event first)."""
    df_i, ab_i = onsets
    order = sorted(al, key=lambda i: evs[i]["t"])
    rows = []
    for k, a in enumerate(order):
        for b in order[k + 1:]:
            r = measure_pair(al[a], al[b], df_i, ab_i, fs)
            if r is None:
                continue
            rc = measure_pair(al[a], al[b], df_i, ab_i, fs, dfw=DF_CODA)
            ea, eb = evs[a], evs[b]
            D = 1.0 - r["cc_df_best"] if isf(r["cc_df_best"]) else NAN
            D_ab = 1.0 - r["cc_ab"] if isf(r["cc_ab"]) else NAN
            snrs_df = (per[a]["snr_df"], per[b]["snr_df"])
            snrs_ab = (per[a]["snr_ab"], per[b]["snr_ab"])
            ok = bool(isf(D) and isf(r["cc_ab"]) and r["cc_ab"] >= MIN_CC_AB
                      and all(isf(v) and v >= MIN_SNR_DF for v in snrs_df)
                      and all(isf(v) and v >= MIN_SNR_AB for v in snrs_ab))
            rows.append(dict(station=station, mult=mult, ev1=a, ev2=b, time1=ea["time"], time2=eb["time"],
                             t1=ea["t"], t2=eb["t"], mag1=ea["mag"], mag2=eb["mag"],
                             cc_ab=r["cc_ab"], cc_df0=r["cc_df0"], cc_df_best=r["cc_df_best"], dt_df_s=r["dt_df_s"],
                             D=D, D0=1.0 - r["cc_df0"] if isf(r["cc_df0"]) else NAN,
                             D_c15=(1.0 - rc["cc_df_best"]) if rc is not None and isf(rc["cc_df_best"]) else NAN,
                             D_ab=D_ab, Dx=D - D_ab if isf(D) and isf(D_ab) else NAN,
                             snr_df1=per[a]["snr_df"], snr_df2=per[b]["snr_df"], snr_ab1=per[a]["snr_ab"],
                             snr_ab2=per[b]["snr_ab"], dB=d_b(ea["t"], eb["t"]), dM=d_m(ea["t"], eb["t"]),
                             cross_turn=ea["t"] < TURN < eb["t"], ok=ok))
    return rows


# ----------------------------------------------------------------------------- IO
def fetch_text(url, tries=5, wait=15, timeout=90):
    for k in range(tries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "inner-core-repeaters (MV2, Munim, Inc.)"})
            with urllib.request.urlopen(req, timeout=timeout) as r:
                if r.status == 200:
                    return r.read().decode("utf-8")
        except Exception as ex:  # noqa: BLE001
            log("fetch error", url, repr(ex)[:160])
        time.sleep(wait * (k + 1))
    return None


def load_inputs(data_dir):
    """ILAR multiplets and catalogue from results/m0324-v1 at the pinned commit (cached)."""
    out = {}
    for name in ("m0324_multiplets.csv", "m0324_catalog.csv"):
        p = os.path.join(data_dir, "inputs", name)
        if not os.path.exists(p):
            body = fetch_text(REPO + INPUT_COMMIT + "/results/m0324-v1/" + name)
            if not body:
                raise RuntimeError("could not fetch " + name)
            os.makedirs(os.path.dirname(p), exist_ok=True)
            with open(p, "w", newline="") as f:
                f.write(body)
        with open(p, newline="") as f:
            out[name] = list(csv.DictReader(f))
    cat = {}
    for r in out["m0324_catalog.csv"]:
        try:
            cat[r["id"]] = dict(id=r["id"], time=r["time"], lat=float(r["lat"]), lon=float(r["lon"]),
                                depth=float(r["depth"]) if r.get("depth") not in (None, "") else 10.0,
                                mag=fnum(r.get("mag")), t=decimal_year(r["time"]))
        except (KeyError, ValueError):
            continue
    mults = {}
    for r in out["m0324_multiplets.csv"]:
        if r.get("status") == "ok":
            ids = [i for i in parse_members(r["members"]) if i in cat]
            if len(ids) >= 2:
                mults[r["mult"]] = ids
    return cat, mults


def discover(client, cat, mults):
    """BHZ stations of NETWORKS at MIN_DIST..MAX_DIST from the SSI centroid, ranked by how many
    multiplet events fall inside a BHZ epoch; at most MAX_STATIONS with >= MIN_COVER."""
    inv = client.get_stations(network=NETWORKS, channel="BHZ", latitude=ANTIPODE[0], longitude=ANTIPODE[1],
                              maxradius=180.0 - MIN_DIST, level="channel")
    events = sorted({i for ids in mults.values() for i in ids})
    found = []
    for net in inv.networks:
        for sta in net.stations:
            chans = []
            lat, lon = sta.latitude, sta.longitude
            for ch in sta.channels:
                chans.append((ch.location_code or "", str(ch.start_date)[:19],
                              str(ch.end_date)[:19] if ch.end_date else "9999"))
                lat, lon = ch.latitude or lat, ch.longitude or lon
            dist = gc_deg(CENTROID[0], CENTROID[1], lat, lon)
            if not (MIN_DIST <= dist <= MAX_DIST):
                continue
            cover = sum(1 for e in events if choose_loc(chans, cat[e]["time"].rstrip("Z")) is not None)
            found.append(dict(net=net.code, sta=sta.code, lat=lat, lon=lon, dist=dist, cover=cover, chans=chans))
    found.sort(key=lambda s: (-s["cover"], s["net"], s["sta"]))
    seen, chosen = set(), []
    for s in found:
        if s["cover"] < MIN_COVER or s["sta"] in seen:
            continue
        seen.add(s["sta"])
        chosen.append(s)
        if len(chosen) >= MAX_STATIONS:
            break
    return found, chosen


def get_response(client, s, data_dir):
    from obspy import read_inventory

    p = os.path.join(data_dir, "resp", "%s.%s.xml" % (s["net"], s["sta"]))
    if os.path.exists(p):
        try:
            return read_inventory(p)
        except Exception:  # noqa: BLE001
            os.remove(p)
    for k in range(3):
        try:
            inv = client.get_stations(network=s["net"], station=s["sta"], channel="BHZ", level="response")
            os.makedirs(os.path.dirname(p), exist_ok=True)
            inv.write(p, format="STATIONXML")
            return inv
        except Exception as ex:  # noqa: BLE001
            log("response fetch failed", s["sta"], repr(ex)[:160])
            time.sleep(10 * (k + 1))
    return None


def fetch_wave(client, s, loc, ev, tdf, data_dir):
    from obspy import UTCDateTime, read
    from obspy.clients.fdsn.header import FDSNNoDataException

    d = os.path.join(data_dir, "raw", "%s.%s" % (s["net"], s["sta"]))
    os.makedirs(d, exist_ok=True)
    p = os.path.join(d, ev["id"] + ".mseed")
    if os.path.exists(p):
        try:
            return read(p), "cached"
        except Exception:  # noqa: BLE001
            os.remove(p)
    if os.path.exists(p + ".nodata"):
        return None, "nodata"
    o = UTCDateTime(ev["time"])
    t1, t2 = o + tdf - PRE - 60, o + tdf + POST + 60
    err = ""
    for k in range(3):
        try:
            st = client.get_waveforms(s["net"], s["sta"], loc or "--", "BHZ", t1, t2)
            if len(st) == 0:
                open(p + ".nodata", "w").close()
                return None, "nodata"
            st.write(p, format="MSEED")
            return st, "fetched"
        except FDSNNoDataException:
            open(p + ".nodata", "w").close()
            return None, "nodata"
        except Exception as ex:  # noqa: BLE001
            err = repr(ex)[:160]
            time.sleep(5 * (k + 1))
    return None, "error:" + err


def process(st, inv, ev, tdf, band):
    """Velocity (response removed when possible), causal band-pass, FS Hz, ak135 PKIKP at PRE,
    normalised by the rms of the NOISE window. Returns (x, resp_ok) or (None, reason)."""
    from obspy import UTCDateTime

    o = UTCDateTime(ev["time"])
    st = st.copy()
    try:
        st.merge(method=0)
    except Exception:  # noqa: BLE001
        pass
    if len(st) == 0:
        return None, "empty"
    tr = max(st, key=lambda t: t.stats.npts)
    if np.ma.isMaskedArray(tr.data):
        if np.ma.count_masked(tr.data) > 0:
            return None, "gappy"
        tr.data = np.asarray(tr.data)
    tr.data = np.asarray(tr.data, dtype=np.float64)
    if tr.stats.npts < 120 * tr.stats.sampling_rate or np.ptp(tr.data) == 0:
        return None, "short-or-flat"
    resp_ok = False
    try:
        tr.detrend("demean")
        tr.detrend("linear")
        tr.taper(0.05)
        if inv is not None:
            try:
                nyq = 0.5 * tr.stats.sampling_rate
                pf = (PREFILT[0], PREFILT[1], min(PREFILT[2], 0.7 * nyq), min(PREFILT[3], 0.9 * nyq))
                tr.remove_response(inventory=inv, output="VEL", pre_filt=pf, water_level=None)
                resp_ok = True
            except Exception:  # noqa: BLE001
                resp_ok = False
        tr.detrend("demean")
        tr.taper(0.05)
        tr.filter("bandpass", freqmin=band[0], freqmax=band[1], corners=4, zerophase=False)
        tr.interpolate(sampling_rate=FS, method="lanczos", a=12, starttime=o + tdf - PRE, npts=N)
    except Exception as ex:  # noqa: BLE001
        return None, "process-error:" + repr(ex)[:80]
    x = np.asarray(tr.data, dtype=np.float64)
    if len(x) != N or not np.all(np.isfinite(x)):
        return None, "bad-length"
    nz = _rms(x, (PRE + NOISE[0]) * FS, (PRE + NOISE[1]) * FS)
    if not (isf(nz) and nz > 0):
        return None, "dead"
    return (x / nz, resp_ok), "ok"


# ----------------------------------------------------------------------------- plots
def plot_family(path, title, al, evs, onsets, per):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    df_i, ab_i = onsets
    order = sorted(al, key=lambda i: evs[i]["t"])
    fig, ax = plt.subplots(figsize=(10, 0.55 * len(order) + 1.6))
    t = np.arange(N) / FS - PRE
    a, b = max(0, df_i - int(8 * FS)), min(N - 1, ab_i + int(10 * FS))
    for k, i in enumerate(order):
        x = al[i]
        sc = np.max(np.abs(x[a:b])) or 1.0
        ax.plot(t[a:b], x[a:b] / sc * 0.45 - k, lw=0.6, color="k")
        ax.text(t[b] + 0.3, -k, "%s M%.1f snr %.1f/%.1f" % (evs[i]["time"][:10], evs[i]["mag"], per[i]["snr_df"],
                                                           per[i]["snr_ab"]), fontsize=6, va="center")
    for w0, w1, i0, col in ((DF_WIN[0], DF_WIN[1], df_i, "C3"), (AB_WIN[0], AB_WIN[1], ab_i, "C0")):
        ax.axvspan(t[i0] + w0, t[i0] + w1, color=col, alpha=0.15)
    ax.set_yticks([])
    ax.set_xlabel("time relative to ak135 PKIKP (s); red = DF window, blue = PKPab reference window")
    ax.set_title(title, fontsize=9)
    fig.tight_layout()
    fig.savefig(path, dpi=100)
    plt.close(fig)


def plot_clock(path, rows, key):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    nr = block_normalise(rows, key)
    if not nr:
        return
    fig, axs = plt.subplots(1, 2, figsize=(10, 4), sharey=True)
    for ax, (lab, fn) in zip(axs, (("backtracking distance |d_equiv| (yr)", d_b), ("time apart |dt| (yr)", d_m))):
        x = [fn(r["t1"], r["t2"]) for r in nr]
        ct = [r["cross_turn"] for r in nr]
        ax.scatter([a for a, c in zip(x, ct) if not c], [r["u"] for r, c in zip(nr, ct) if not c], s=6, c="C0",
                   label="same side of 2008.5")
        ax.scatter([a for a, c in zip(x, ct) if c], [r["u"] for r, c in zip(nr, ct) if c], s=6, c="C3",
                   label="across 2008.5")
        ax.set_xlabel(lab)
    axs[0].set_ylabel("within-block rank of %s (0 = most alike)" % key)
    axs[0].legend(fontsize=7)
    fig.suptitle("SSI -> stations near the SSI antipode (155-179.5 deg), ok pairs, %s" % key, fontsize=9)
    fig.tight_layout()
    fig.savefig(path, dpi=100)
    plt.close(fig)


# ----------------------------------------------------------------------------- output
def jsonable(o):
    if isinstance(o, dict):
        return {str(k): jsonable(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [jsonable(v) for v in o]
    if isinstance(o, (bool, np.bool_)):
        return bool(o)
    if isinstance(o, (float, np.floating)):
        return float(o) if math.isfinite(float(o)) else None
    if isinstance(o, np.integer):
        return int(o)
    return o


def write_csv(path, rows):
    keys = []
    for r in rows:
        for k in r:
            if k not in keys:
                keys.append(k)
    with open(path, "w", newline="") as f:
        if keys:
            w = csv.DictWriter(f, fieldnames=keys)
            w.writeheader()
            for r in rows:
                w.writerow({k: ("" if (isinstance(v, float) and not math.isfinite(v)) else v) for k, v in r.items()})


def grade(rows):
    ok = [r for r in rows if r.get("ok")]
    g = dict(n_rows=len(rows), n_ok=len(ok))
    g["P16_0_calibration"] = calibration(ok, "D")
    g["P16_1_core"] = clock_test(ok, "D")
    g["P16_2_source"] = clock_test(ok, "D_ab")
    g["headline"] = headline(g["P16_0_calibration"], g["P16_1_core"], g["P16_2_source"])
    sec = {}
    for key in ("Dx", "D0", "D_c15"):
        sec[key] = dict(calibration=calibration(ok, key), clock=clock_test(ok, key))
    sec["D_ab_calibration"] = calibration(ok, "D_ab")
    sec["new_event_matches"] = new_event_matches(ok, "D")
    sec["clock_scan"] = clock_scan(ok, "D")
    per = {}
    for s in sorted({r["station"] for r in ok}):
        rs = [r for r in ok if r["station"] == s]
        per[s] = dict(n_ok=len(rs), calibration=calibration(rs, "D", nboot=300), clock=clock_test(rs, "D", nboot=300))
    sec["per_station"] = per
    g["secondary"] = sec
    return g


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", default="/data")
    ap.add_argument("--out-dir", default="/out")
    ap.add_argument("--time-budget-min", type=float, default=150.0)
    args = ap.parse_args()
    t0 = time.time()
    deadline = t0 + 60.0 * args.time_budget_min * 0.85
    os.makedirs(args.data_dir, exist_ok=True)
    os.makedirs(args.out_dir, exist_ok=True)
    plots = os.path.join(args.out_dir, "plots")
    os.makedirs(plots, exist_ok=True)
    s = dict(version=VERSION, started=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()), input_commit=INPUT_COMMIT,
             constants=dict(min_dist=MIN_DIST, max_dist=MAX_DIST, band=BAND, band_alt=BAND_ALT, df_win=DF_WIN,
                            ab_win=AB_WIN, min_cc_ab=MIN_CC_AB, min_snr_df=MIN_SNR_DF, min_snr_ab=MIN_SNR_AB,
                            near=NEAR, far_m=FAR_M, far_b=FAR_B, cal_ratio=CAL_RATIO, frac=FRAC, turn=TURN,
                            back_ratio=BACK_RATIO))
    s["prereg_public_at_start"] = bool(fetch_text(PREREG_URL, tries=2, wait=5))
    cat, mults = load_inputs(args.data_dir)
    s["n_multiplets"] = len(mults)
    s["n_events"] = len({i for v in mults.values() for i in v})

    from obspy.clients.fdsn import Client
    from obspy.taup import TauPyModel

    client = None
    for base in ("https://service.earthscope.org", "IRIS"):
        try:
            client = Client(base, timeout=120)
            break
        except Exception as ex:  # noqa: BLE001
            log("client failed", base, repr(ex)[:160])
    if client is None:
        s["error"] = "no FDSN client"
        with open(os.path.join(args.out_dir, "nea_summary.json"), "w") as f:
            json.dump(jsonable(s), f, indent=1)
        return 1
    found, chosen = [], []
    for k in range(3):
        try:
            found, chosen = discover(client, cat, mults)
            break
        except Exception as ex:  # noqa: BLE001
            log("discover failed", repr(ex)[:200])
            time.sleep(30)
    write_csv(os.path.join(args.out_dir, "nea_stations.csv"),
              [dict({k: v for k, v in f_.items() if k != "chans"}, chosen=any(c is f_ for c in chosen),
                    epochs=";".join("%s:%s..%s" % c for c in f_["chans"])) for f_ in found])
    s["stations_chosen"] = ["%s.%s" % (c["net"], c["sta"]) for c in chosen]
    log("stations:", s["stations_chosen"])
    model = TauPyModel("ak135")
    rows, alt_rows, ev_status, fam_info = [], [], [], []
    for st_ in chosen:
        name = "%s.%s" % (st_["net"], st_["sta"])
        if time.time() > deadline:
            fam_info.append(dict(station=name, status="skipped-time-budget"))
            continue
        inv = get_response(client, st_, args.data_dir)
        for mult, ids in sorted(mults.items()):
            tr, tr_alt, ab_off, evs = {}, {}, {}, {}
            for eid in ids:
                ev = cat[eid]
                row = dict(station=name, mult=mult, id=eid, time=ev["time"], status="")
                ev_status.append(row)
                loc = choose_loc(st_["chans"], ev["time"].rstrip("Z"))
                if loc is None:
                    row["status"] = "no-epoch"
                    continue
                if time.time() > deadline:
                    row["status"] = "skipped-time-budget"
                    continue
                dist = gc_deg(ev["lat"], ev["lon"], st_["lat"], st_["lon"])
                try:
                    arr = model.get_travel_times(source_depth_in_km=max(ev["depth"], 0.0), distance_in_degree=dist,
                                                 phase_list=["PKIKP", "PKPab"])
                except Exception as ex:  # noqa: BLE001
                    row["status"] = "taup-error:" + repr(ex)[:60]
                    continue
                tdf = next((a.time for a in arr if a.name == "PKIKP"), NAN)
                tab = next((a.time for a in arr if a.name == "PKPab"), NAN)
                row.update(loc=loc, dist=dist, ab_minus_df_pred=tab - tdf if isf(tab) and isf(tdf) else NAN)
                if not (isf(tdf) and isf(tab)):
                    row["status"] = "no-phase"
                    continue
                stw, why = fetch_wave(client, st_, loc, ev, tdf, args.data_dir)
                if stw is None:
                    row["status"] = why
                    continue
                pr, why = process(stw, inv, ev, tdf, BAND)
                if pr is None:
                    row["status"] = why
                    continue
                pa, _ = process(stw, inv, ev, tdf, BAND_ALT)
                tr[eid], ab_off[eid], evs[eid] = pr[0], tab - tdf, ev
                if pa is not None:
                    tr_alt[eid] = pa[0]
                row.update(status="ok", resp_ok=pr[1])
            for band_name, traces, out_rows in (("primary", tr, rows), ("alt", tr_alt, alt_rows)):
                if len(traces) < 2:
                    if band_name == "primary":
                        fam_info.append(dict(station=name, mult=mult, band=band_name, n_in=len(traces),
                                             status="too-few"))
                    continue
                try:
                    info, per, al, onsets = measure_family(traces, {i: ab_off[i] for i in traces})
                except Exception as ex:  # noqa: BLE001
                    traceback.print_exc()
                    info, per, al, onsets = dict(status="exception:" + repr(ex)[:100]), {}, {}, None
                info.update(station=name, mult=mult, band=band_name)
                fam_info.append(info)
                if onsets is None:
                    continue
                out_rows += family_pairs(name, mult, evs, al, per, onsets)
                if band_name == "primary" and len(al) >= 3:
                    try:
                        plot_family(os.path.join(plots, "nea_%s_%s_%s.png" % (st_["net"], st_["sta"], mult)),
                                    "%s %s (ILAR multiplet), %.1f-%.1f Hz" % (name, mult, BAND[0], BAND[1]),
                                    al, evs, onsets, per)
                    except Exception as ex:  # noqa: BLE001
                        log("plot failed", repr(ex)[:120])
        log(name, "pairs so far", len(rows), "ok", sum(1 for r in rows if r["ok"]))
    s["families"] = fam_info
    s["n_pairs"], s["n_pairs_ok"] = len(rows), sum(1 for r in rows if r["ok"])
    s["grade"] = grade(rows)
    s["grade_alt_band"] = grade(alt_rows)
    try:
        plot_clock(os.path.join(plots, "nea_clock_D.png"), [r for r in rows if r["ok"]], "D")
        plot_clock(os.path.join(plots, "nea_clock_D_ab.png"), [r for r in rows if r["ok"]], "D_ab")
    except Exception as ex:  # noqa: BLE001
        log("clock plot failed", repr(ex)[:120])
    s["finished"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    s["minutes"] = round((time.time() - t0) / 60.0, 1)
    write_csv(os.path.join(args.out_dir, "nea_pairs.csv"), rows)
    write_csv(os.path.join(args.out_dir, "nea_pairs_alt.csv"), alt_rows)
    write_csv(os.path.join(args.out_dir, "nea_event_status.csv"), ev_status)
    with open(os.path.join(args.out_dir, "nea_summary.json"), "w") as f:
        json.dump(jsonable(s), f, indent=1)
    log(json.dumps(jsonable(dict(headline=s["grade"]["headline"], cal=s["grade"]["P16_0_calibration"],
                                 core=s["grade"]["P16_1_core"], src=s["grade"]["P16_2_source"])))[:3000])
    return 0


if __name__ == "__main__":
    sys.exit(main())

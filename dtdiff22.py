#!/usr/bin/env python3
# Copyright 2026 MV2, Munim, Inc.
# Licensed under the Apache License, Version 2.0 (Apache-2.0). See LICENSE.
# Made with AI agents working in MV2's name.
"""dtdiff22.py - pre-registration 22: Chile -> Makanchi (MKAR) inner-core differential travel times.

Plan (published before this job ran): findings/2026-11-20-prereg22-mkar-ddt-plan.md.

For every source-verified Chilean earthquake pair from pre-registration 21 (S <= 1/3 at >= 2 mantle-P
stations), rebuild both MKAR beams with mkar20.mkar_beam (pinned commit; waveforms come from the
shared cache), align the later event on the earlier one using the outer-core PKP arrival (BC, or AB
where BC does not exist), and measure the residual shift of the inner-core arrival PKIKP (DF):

    ddt = [t_DF - t_PKP](later event) - [t_DF - t_PKP](earlier event)     (seconds)

Then compare three published clocks for how t_DF - t_PKP evolves with calendar time t:
    L  steady change            g_L(t) = t - 2008.5
    B  backtracking (Wang, Vidale et al. 2024): change before 2008.5, reverses 2.5x slower after
                                g_B(t) = t - 2008.5 (t < 2008.5), -(t - 2008.5)/2.5 (after)
    P  pause (Yang & Song 2023): change before 2009.0, none after
                                g_P(t) = min(t, 2009.0) - 2009.0
Each model predicts ddt = a * (g(t2) - g(t1)) with one free amplitude a (s per model-year).
The decision rule (grade()) is fixed in the plan.

Top level needs only the standard library and numpy; obspy and the pinned helpers are imported
inside main().

Usage: python dtdiff22.py --data-dir /data --out-dir /out --wait-min 120 --time-budget-min 240
"""

import argparse
import csv
import hashlib
import io
import json
import math
import os
import random
import sys
import time
import urllib.request

import numpy as np

VERSION = "0.1.0"
NAN = float("nan")
EPS = 1e-9

RAW = "https://raw.githubusercontent.com/munimtechnologies/inner-core-repeaters/"
M20 = ("mkar20.py", "a7f1ceb06c422297020938e08805f2d5d81ce929", 'VERSION = "0.1.0"')
PLAN = "findings/2026-11-20-prereg22-mkar-ddt-plan.md"
PLAN_TOKEN = b"pre-registration 22"
PAIRS = "results/mkar21-v1/mkar20_pairs.csv"          # pre-registration 21's pair table
FALLBACK_PAIRS = "results/mkar20-v1/mkar20_pairs.csv"  # only if P21's table never appears

TURN, RATIO, PAUSE = 2008.5, 2.5, 2009.0
BANDS = (("primary", (1.0, 2.0)), ("secondary", (0.6, 1.6)))
S_MAX, MIN_ST = 1.0 / 3.0, 2

# measurement (seconds)
SIG = (-3.0, 7.0)            # coarse-alignment / onset-peak window, relative to ak135 PKIKP
ONSET_SEARCH = (-3.0, 4.0)
ONSET_FRAC = 0.3
COARSE_LAG = 3.0
DF_ONSET_MAX = 1.0           # onset within +1.0 s of predicted DF -> DF is the onset
DL_RANGE = (1.0, 6.0)        # accepted predicted t_PKP - t_DF
LATE_WIN = (-0.4, 1.6)       # relative to the PKP anchor
LATE_LAG = 0.4
DF_PRE, DF_GAP = -0.5, 0.4   # DF window: anchor - 0.5 s ... anchor + dl - 0.4 s
DF_LAG = 0.25
SNR_DF_MIN = 3.0
CC_LATE_MIN = 0.7
CC_DF_MIN = 0.6

# decision
SUPPLY = dict(n=20, families=6, straddle=8, lever_n=5, lever_gb=2.0)
NBOOT = 1000
SEED = 22
R2_MIN = 0.10
WIN_FRAC = 0.80
NONE_RATE = 0.005            # s per model-year
SHORT_SPAN = 0.5             # yr, noise-floor pairs


def log(*a):
    print(time.strftime("%H:%M:%S"), "[dtdiff22]", *a, flush=True)


def isf(x):
    try:
        return math.isfinite(float(x)) and not isinstance(x, bool)
    except (TypeError, ValueError):
        return False


def fnum(x):
    return float(x) if isf(x) else NAN


def median(xs):
    v = sorted(float(x) for x in xs if isf(x))
    if not v:
        return NAN
    n = len(v)
    return v[n // 2] if n % 2 else 0.5 * (v[n // 2 - 1] + v[n // 2])


def pct(xs, q):
    v = sorted(float(x) for x in xs if isf(x))
    if not v:
        return NAN
    if len(v) == 1:
        return v[0]
    p = q * (len(v) - 1)
    k = int(math.floor(p))
    k1 = min(k + 1, len(v) - 1)
    return v[k] + (p - k) * (v[k1] - v[k])


def robust_sd(xs):
    m = median(xs)
    if not isf(m):
        return NAN
    return 1.4826 * median(abs(float(x) - m) for x in xs if isf(x))


# ----------------------------------------------------------------------------- clocks
def g_L(t):
    return t - TURN


def g_B(t):
    return (t - TURN) if t < TURN else -(t - TURN) / RATIO


def g_P(t):
    return min(t, PAUSE) - PAUSE


MODELS = (("L", g_L), ("B", g_B), ("P", g_P))


def dg(g, r):
    return g(float(r["t2"])) - g(float(r["t1"]))


def fit_origin(xs, ys):
    """Least squares y = a x through the origin. Returns (a, R2 against y = 0)."""
    sxx = sum(x * x for x in xs)
    syy = sum(y * y for y in ys)
    if sxx <= 0 or not xs:
        return NAN, NAN
    a = sum(x * y for x, y in zip(xs, ys)) / sxx
    rss = sum((y - a * x) ** 2 for x, y in zip(xs, ys))
    return a, (1.0 - rss / syy) if syy > 0 else NAN


def fit_all(rows):
    ys = [float(r["ddt_s"]) for r in rows]
    return {name: fit_origin([dg(g, r) for r in rows], ys) for name, g in MODELS}


def resample(rows, rng):
    by = {}
    for r in rows:
        by.setdefault(r["mult"], []).append(r)
    fams = sorted(by)
    return [r for f in (rng.choice(fams) for _ in fams) for r in by[f]]


def supply(rows):
    return dict(n=len(rows), families=len({r["mult"] for r in rows}),
                straddle=sum(1 for r in rows if r.get("kind") == "straddle"),
                lever=sum(1 for r in rows if abs(dg(g_B, r)) >= SUPPLY["lever_gb"] - EPS))


def grade(rows, nboot=NBOOT, seed=SEED):
    """Pre-registered decision on rows with finite ddt_s (keys t1 < t2, ddt_s, mult, kind)."""
    rows = [r for r in rows if isf(r.get("ddt_s"))]
    sup = supply(rows)
    out = dict(supply=sup, models={})
    if not (sup["n"] >= SUPPLY["n"] and sup["families"] >= SUPPLY["families"]
            and sup["straddle"] >= SUPPLY["straddle"] and sup["lever"] >= SUPPLY["lever_n"]):
        out["decision"] = "untestable"
        return out
    full = fit_all(rows)
    rng = random.Random(seed)
    boots = [fit_all(resample(rows, rng)) for _ in range(nboot)]
    for name, _ in MODELS:
        a_s = [b[name][0] for b in boots]
        p05, p95 = pct(a_s, 0.05), pct(a_s, 0.95)
        out["models"][name] = dict(a=full[name][0], r2=full[name][1], a_p05=p05, a_p95=p95,
                                   significant=bool(isf(p05) and isf(p95) and (p05 > 0 or p95 < 0)))
    r2s = {n: (m["r2"] if isf(m["r2"]) else -1e9) for n, m in out["models"].items()}
    best = max(r2s, key=lambda n: r2s[n])
    out["best"] = best
    bm = out["models"][best]
    if bm["r2"] >= R2_MIN - EPS and bm["significant"]:
        wins = {}
        for other, _ in MODELS:
            if other == best:
                continue
            k = sum(1 for b in boots if isf(b[best][1]) and isf(b[other][1]) and b[best][1] > b[other][1])
            wins[other] = k / float(len(boots)) if boots else NAN
        out["win_fraction"] = wins
        out["decision"] = best if all(w >= WIN_FRAC - EPS for w in wins.values()) else "change-ambiguous"
    elif all(isf(m["a_p05"]) and isf(m["a_p95"]) and m["a_p05"] > -NONE_RATE and m["a_p95"] < NONE_RATE
             for m in out["models"].values()):
        out["decision"] = "none"
    else:
        out["decision"] = "inconclusive"
    return out


# ----------------------------------------------------------------------------- pair table
def select_pairs(csv_rows):
    """P21 pair rows (strings) -> measurable rows (mkar_ok) with verified flag, t1 < t2."""
    out = []
    for r in csv_rows:
        if str(r.get("mkar_ok", "")).strip() != "True":
            continue
        t1, t2 = fnum(r.get("t1")), fnum(r.get("t2"))
        if not (isf(t1) and isf(t2)):
            continue
        e1, e2 = str(r["ev1"]), str(r["ev2"])
        if t2 < t1:
            t1, t2, e1, e2 = t2, t1, e2, e1
        S, n_st = fnum(r.get("S")), fnum(r.get("n_st"))
        out.append(dict(ev1=e1, ev2=e2, t1=t1, t2=t2, kind=r.get("kind", ""), mult=r.get("mult", ""),
                        S=S, n_st=int(n_st) if isf(n_st) else 0, dB=fnum(r.get("dB")), dM=fnum(r.get("dM")),
                        verified=bool(isf(S) and S <= S_MAX + 1e-12 and isf(n_st) and n_st >= MIN_ST)))
    return out


# ----------------------------------------------------------------------------- signal helpers
def envelope(x, fs, smooth_s=0.5):
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


def lag_scan(a, b, i0, i1, maxlag):
    """Best integer lag L in [-maxlag, maxlag] maximising ncc(a[i0:i1], b[i0+L:i1+L]) (b later by L).
    Returns (L, cc, cc at L-1, cc at L+1) or (None, nan, nan, nan)."""
    ccs = {}
    for L in range(-maxlag - 1, maxlag + 2):
        j0, j1 = i0 + L, i1 + L
        if i0 < 0 or i1 > len(a) or j0 < 0 or j1 > len(b) or i1 - i0 < 2:
            continue
        ccs[L] = ncc(a[i0:i1], b[j0:j1])
    inside = [L for L in range(-maxlag, maxlag + 1) if isf(ccs.get(L))]
    if not inside:
        return None, NAN, NAN, NAN
    best = max(inside, key=lambda L: ccs[L])
    return best, ccs[best], ccs.get(best - 1, NAN), ccs.get(best + 1, NAN)


def parabolic(cm, c0, cp):
    if not (isf(cm) and isf(c0) and isf(cp)):
        return 0.0
    den = cm - 2.0 * c0 + cp
    if den >= 0:
        return 0.0
    return max(-0.5, min(0.5, 0.5 * (cm - cp) / den))


def fshift(x, s):
    """y[n] = x[n + s] for real s (circular Fourier shift)."""
    x = np.asarray(x, dtype=float)
    n = len(x)
    X = np.fft.rfft(x)
    k = np.arange(len(X))
    return np.fft.irfft(X * np.exp(2j * np.pi * k * s / n), n)


def peak_norm(x, i0, i1):
    m = float(np.max(np.abs(x[i0:i1]))) if i1 > i0 else 0.0
    return x / m if m > 0 else x


def rms(x):
    x = np.asarray(x, dtype=float)
    return float(np.sqrt(np.mean(x * x))) if len(x) else NAN


def pick_onset(stack, fs, c0, sig_end):
    env = envelope(stack, fs, 0.25)
    w0, w1 = c0 + int(round(SIG[0] * fs)), c0 + int(round(sig_end * fs))
    peak = float(np.max(env[w0:w1]))
    if not peak > 0:
        return None
    s0, s1 = c0 + int(round(ONSET_SEARCH[0] * fs)), c0 + int(round(ONSET_SEARCH[1] * fs))
    for i in range(s0, s1):
        if env[i] >= ONSET_FRAC * peak:
            return i
    return None


def measure_pair(xa, xb, fs, c0, dl, late_onset_ok=False):
    """xa: earlier event's beam, xb: later event's beam, both with ak135 PKIKP at sample c0 and
    noise rms 1. dl: predicted t_PKP - t_DF (s). Returns a dict with 'status' ('ok' if all gates
    pass) and, when measured, ddt_s = DF shift of xb after aligning it on the PKP arrival."""
    out = dict(dl_s=dl, status="")
    if xa is None or xb is None:
        out["status"] = "no-beam"
        return out
    if not (isf(dl) and DL_RANGE[0] - EPS <= dl <= DL_RANGE[1] + EPS):
        out["status"] = "dl-out-of-range"
        return out
    xa = np.asarray(xa, dtype=float)
    xb = np.asarray(xb, dtype=float)
    r = lambda s: int(round(s * fs))  # noqa: E731
    sig_end = max(SIG[1], dl + 2.5)
    i0, i1 = c0 + r(SIG[0]), c0 + r(sig_end)
    L0, c_env, _, _ = lag_scan(envelope(xa, fs), envelope(xb, fs), i0, i1, r(COARSE_LAG))
    if L0 is None:
        out["status"] = "no-coarse-lag"
        return out
    y0 = fshift(xb, L0)
    onset = pick_onset(peak_norm(xa, i0, i1) + peak_norm(y0, i0, i1), fs, c0, sig_end)
    if onset is None:
        out["status"] = "no-onset"
        return out
    out["onset_rel_s"] = (onset - c0) / fs
    if out["onset_rel_s"] <= DF_ONSET_MAX + EPS:
        out["onset_class"], kdf = "df-onset", onset
    elif late_onset_ok:
        out["onset_class"], kdf = "late-onset", onset - r(dl)
    else:
        out["onset_class"] = "late-onset"
        out["status"] = "df-not-onset"
        return out
    klate = kdf + r(dl)
    j0, j1 = klate + r(LATE_WIN[0]), klate + r(LATE_WIN[1])
    L1, c1, cm, cp_ = lag_scan(xa, y0, j0, j1, r(LATE_LAG))
    if L1 is None:
        out["status"] = "no-late-lag"
        return out
    tot = L0 + L1 + parabolic(cm, c1, cp_)
    y = fshift(xb, tot)
    out["lag_late_s"] = tot / fs
    out["cc_late"] = ncc(xa[j0:j1], y[j0:j1])
    d0, d1 = kdf + r(DF_PRE), kdf + r(dl - DF_GAP)
    Ld, cd, dm_, dp_ = lag_scan(xa, y, d0, d1, r(DF_LAG))
    if Ld is None:
        out["status"] = "no-df-lag"
        return out
    s = Ld + parabolic(dm_, cd, dp_)
    out["ddt_s"] = s / fs
    out["cc_df"] = ncc(xa[d0:d1], fshift(y, s)[d0:d1])
    out["snr_df_1"] = rms(xa[d0:d1])
    out["snr_df_2"] = rms(y[d0:d1])
    if abs(L1) >= r(LATE_LAG):
        out["status"] = "late-edge"
    elif abs(Ld) >= r(DF_LAG):
        out["status"] = "df-edge"
    elif not (out["snr_df_1"] >= SNR_DF_MIN and out["snr_df_2"] >= SNR_DF_MIN):
        out["status"] = "low-snr-df"
    elif not (isf(out["cc_late"]) and out["cc_late"] >= CC_LATE_MIN):
        out["status"] = "low-cc-late"
    elif not (isf(out["cc_df"]) and out["cc_df"] >= CC_DF_MIN):
        out["status"] = "low-cc-df"
    else:
        out["status"] = "ok"
    return out


def predicted_dl(model, dist, depth):
    """ak135 t(first PKP, i.e. BC or AB) - t(PKIKP) and the PKP ray parameter (s/deg)."""
    z = max(float(depth), 0.0) if isf(depth) else 10.0
    arr = model.get_travel_times(source_depth_in_km=z, distance_in_degree=dist, phase_list=["PKIKP", "PKP"])
    tdf = min((a.time for a in arr if a.name == "PKIKP"), default=NAN)
    late = [a for a in arr if a.name != "PKIKP"]
    if not late or not isf(tdf):
        return NAN, NAN
    a = min(late, key=lambda a: a.time)
    return a.time - tdf, float(a.ray_param_sec_degree)


def forecasts(primary, secondary, frac_df_onset, noise_sd):
    return dict(H1_frac_verified_df_onset_ge_half=bool(isf(frac_df_onset) and frac_df_onset >= 0.5),
                H2_primary_decision=primary.get("decision"),
                H3_noise_floor_le_0p02s=bool(isf(noise_sd) and noise_sd <= 0.02 + EPS),
                H4_secondary_testable=secondary.get("decision") != "untestable")


# ----------------------------------------------------------------------------- I/O
def fetch_bytes(url, tries=6, wait=20, timeout=90):
    last = None
    for k in range(tries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "inner-core-repeaters (MV2, Munim, Inc.)"})
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                return resp.read()
        except Exception as ex:  # noqa: BLE001
            last = ex
            time.sleep(wait * (k + 1))
    raise RuntimeError("fetch failed %s: %r" % (url, last))


def wait_for(path, must, wait_min):
    """Poll main/<path> until its bytes contain `must` (lower-cased). Returns (bytes or None, minutes)."""
    t0 = time.time()
    while True:
        try:
            body = fetch_bytes(RAW + "main/" + path, tries=1, wait=0, timeout=60)
            if must in body.lower():
                return body, round((time.time() - t0) / 60.0, 1)
        except Exception:  # noqa: BLE001
            pass
        if time.time() - t0 > wait_min * 60:
            return None, round((time.time() - t0) / 60.0, 1)
        time.sleep(120)


def ensure_m20(here):
    name, sha, must = M20
    p = os.path.join(here, name)
    if not os.path.exists(p):
        body = fetch_bytes(RAW + sha + "/" + name)
        if must.encode() not in body:
            raise RuntimeError("unexpected content in %s at %s" % (name, sha))
        with open(p, "wb") as f:
            f.write(body)
    with open(p, "rb") as f:
        return dict(file=name, commit=sha, sha256=hashlib.sha256(f.read()).hexdigest())


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
    for row in rows:
        for k in row:
            if k not in keys:
                keys.append(k)
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        for row in rows:
            w.writerow({k: (round(v, 5) if isinstance(v, float) and math.isfinite(v) else
                            ("" if isinstance(v, float) else v)) for k, v in row.items()})


def plot(path, rows, g):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.6), sharey=True)
    colours = dict(straddle="C0", before="C1", after="C2")
    for ax, (name, fn) in zip(axes, MODELS):
        for k, c in colours.items():
            pts = [(dg(fn, r), r["ddt_s"]) for r in rows if r.get("kind") == k]
            if pts:
                ax.scatter([p[0] for p in pts], [p[1] for p in pts], s=14, c=c, alpha=0.8, label=k)
        m = g.get("models", {}).get(name)
        if m and isf(m.get("a")):
            xs = [dg(fn, r) for r in rows]
            if xs:
                lo, hi = min(xs + [0.0]), max(xs + [0.0])
                ax.plot([lo, hi], [m["a"] * lo, m["a"] * hi], "k-", lw=1)
            ax.set_title("%s: a = %.4f s/yr, R2 = %.2f" % (name, m["a"], m["r2"] if isf(m["r2"]) else NAN), fontsize=9)
        ax.axhline(0, color="0.6", lw=0.6)
        ax.set_xlabel("clock difference g_%s(t2) - g_%s(t1) (yr)" % (name, name))
    axes[0].set_ylabel("ddt = shift of DF relative to PKP (s)")
    axes[0].legend(fontsize=8)
    fig.suptitle("Pre-registration 22: Chile -> MKAR, verified pairs, 1-2 Hz (decision: %s)" % g.get("decision"))
    fig.tight_layout()
    fig.savefig(path, dpi=110)
    plt.close(fig)


# ----------------------------------------------------------------------------- main
def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", default="/data")
    ap.add_argument("--out-dir", default="/out")
    ap.add_argument("--wait-min", type=float, default=0.0)
    ap.add_argument("--time-budget-min", type=float, default=240.0)
    ap.add_argument("--pairs", default=PAIRS)
    a = ap.parse_args(argv)
    t0 = time.time()
    deadline = t0 + 60.0 * a.time_budget_min
    os.makedirs(os.path.join(a.out_dir, "plots"), exist_ok=True)
    os.makedirs(a.data_dir, exist_ok=True)
    here = os.path.dirname(os.path.abspath(__file__))
    if here not in sys.path:
        sys.path.insert(0, here)
    s = dict(version=VERSION, started=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
             constants=dict(turn=TURN, ratio=RATIO, pause=PAUSE, bands=BANDS, s_max=S_MAX, min_st=MIN_ST,
                            df_onset_max=DF_ONSET_MAX, dl_range=DL_RANGE, late_win=LATE_WIN, late_lag=LATE_LAG,
                            df_pre=DF_PRE, df_gap=DF_GAP, df_lag=DF_LAG, snr_df_min=SNR_DF_MIN,
                            cc_late_min=CC_LATE_MIN, cc_df_min=CC_DF_MIN, supply=SUPPLY, nboot=NBOOT,
                            seed=SEED, r2_min=R2_MIN, win_frac=WIN_FRAC, none_rate=NONE_RATE))
    s["mkar20"] = ensure_m20(here)
    import mkar20 as m20

    s["deps"] = m20.ensure_deps(here)
    plan, s["plan_wait_min"] = wait_for(PLAN, PLAN_TOKEN, a.wait_min)
    s["plan_public_at_start"] = plan is not None
    body, s["pairs_wait_min"] = wait_for(a.pairs, b"ev1,", a.wait_min)
    s["pairs_path"] = a.pairs
    if body is None:
        s["status"] = "no-input"
        with open(os.path.join(a.out_dir, "dtdiff22_summary.json"), "w", encoding="utf-8") as f:
            json.dump(jsonable(s), f, indent=1)
        log("pair table not found on main:", a.pairs)
        return 2
    s["pairs_sha256"] = hashlib.sha256(body).hexdigest()
    with open(os.path.join(a.out_dir, "input_pairs.csv"), "wb") as f:
        f.write(body)
    pairs = select_pairs(list(csv.DictReader(io.StringIO(body.decode("utf-8")))))
    s["n_pairs_mkar_ok"] = len(pairs)
    s["n_pairs_verified"] = sum(1 for p in pairs if p["verified"])
    log("pairs mkar_ok", len(pairs), "verified", s["n_pairs_verified"])

    import chilepilot as cp
    import chilepower as cpow
    from obspy.clients.fdsn import Client
    from obspy.taup import TauPyModel

    cat_path = os.path.join(a.out_dir, "chile_catalogue.csv")
    with open(cat_path, "wb") as f:
        f.write(fetch_bytes(m20.CATALOGUE))
    evmap = {r["id"]: r for r in cpow.eligible(cpow.read_catalogue(cat_path))}
    os.remove(cat_path)
    need = sorted({e for p in pairs for e in (p["ev1"], p["ev2"]) if e in evmap}, key=lambda e: evmap[e]["t"])
    s["n_events"] = len(need)
    s["n_pair_events_missing_from_catalogue"] = len({e for p in pairs for e in (p["ev1"], p["ev2"])} - set(evmap))

    model = TauPyModel("ak135")
    dls = {}
    for e in need:
        E = evmap[e]
        try:
            dls[e] = predicted_dl(model, E["dist"], E.get("depth"))
        except Exception:  # noqa: BLE001
            dls[e] = (NAN, NAN)

    client = Client("https://service.earthscope.org", timeout=180, _discover_services=False)
    invs = {}
    for net, sta in m20.MKAR:
        inv = cp.get_inventory(client, net, sta, a.data_dir)
        invs[(net, sta)] = (inv, cp.chans_of(inv) if inv is not None else [])

    beams, ev_rows = {}, []
    for bname, band in BANDS:
        m20.BAND = band
        beams[bname] = {}
        for k, e in enumerate(need):
            if time.time() > deadline:
                ev_rows.append(dict(band=bname, id=e, status="skipped-time-budget"))
                continue
            try:
                b, row = m20.mkar_beam(cp, client, model, evmap[e], invs, a.data_dir)
            except Exception as ex:  # noqa: BLE001
                b, row = None, dict(id=e, status="error:" + repr(ex)[:100])
            row = dict(band=bname, **row)
            row["dl_s"], row["pkp_ray_param_s_deg"] = dls[e]
            row["dist_deg"] = evmap[e]["dist"]
            ev_rows.append(row)
            if b is not None:
                beams[bname][e] = b
            if k % 50 == 0:
                log(bname, k + 1, "/", len(need), e, row.get("status"))
    write_csv(os.path.join(a.out_dir, "dtdiff22_events.csv"), ev_rows)

    c0 = int(round(cp.PRE * cp.FS))
    prow = []
    for p in pairs:
        dl = dls.get(p["ev1"], (NAN, NAN))[0]
        for bname, _ in BANDS:
            bm = beams.get(bname, {})
            for strict in (True, False):
                m = measure_pair(bm.get(p["ev1"]), bm.get(p["ev2"]), cp.FS, c0, dl, late_onset_ok=not strict)
                prow.append(dict(p, band=bname, anchor="strict" if strict else "all-onset", **m))
    write_csv(os.path.join(a.out_dir, "dtdiff22_pairs.csv"), prow)

    def sel(band, anchor, verified=True):
        return [r for r in prow if r["band"] == band and r["anchor"] == anchor and r["status"] == "ok"
                and (r["verified"] or not verified)]

    g = dict(primary=grade(sel("primary", "strict")))
    g["decision"] = g["primary"]["decision"]
    g["secondary_band_0p6_1p6"] = grade(sel("secondary", "strict"))
    g["exploratory_all_onset"] = grade(sel("primary", "all-onset"))
    g["diagnostic_unverified"] = grade(sel("primary", "strict", verified=False))
    short = [r["ddt_s"] for r in sel("primary", "strict", verified=False) if r["t2"] - r["t1"] < SHORT_SPAN]
    g["noise_floor"] = dict(n=len(short), robust_sd_s=robust_sd(short), median_s=median(short))
    pre = [r for r in sel("primary", "strict") if r["kind"] == "before" and r["t2"] - r["t1"] >= SHORT_SPAN - EPS]
    a_pre, _ = fit_origin([dg(g_L, r) for r in pre], [r["ddt_s"] for r in pre])
    g["pre_turn_long_span"] = dict(n=len(pre), families=len({r["mult"] for r in pre}), a_L=a_pre)
    V = [r for r in prow if r["band"] == "primary" and r["anchor"] == "strict" and r["verified"]]
    st = {}
    for r in V:
        st[r["status"]] = st.get(r["status"], 0) + 1
    g["verified_status_counts_primary"] = st
    frac = (sum(1 for r in V if r.get("onset_class") == "df-onset") / float(len(V))) if V else NAN
    g["frac_verified_df_onset"] = frac
    s["grade"] = g
    s["forecasts"] = forecasts(g["primary"], g["secondary_band_0p6_1p6"], frac, g["noise_floor"]["robust_sd_s"])
    s["finished"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    with open(os.path.join(a.out_dir, "dtdiff22_summary.json"), "w", encoding="utf-8") as f:
        json.dump(jsonable(s), f, indent=1)
    try:
        plot(os.path.join(a.out_dir, "plots", "dtdiff22_clocks.png"), sel("primary", "strict"), g["primary"])
    except Exception as ex:  # noqa: BLE001
        log("plot failed", repr(ex)[:160])
    log("decision", g["decision"], json.dumps(jsonable(s["forecasts"])))
    return 0


if __name__ == "__main__":
    sys.exit(main())

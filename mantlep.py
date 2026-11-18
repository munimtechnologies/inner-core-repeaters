#!/usr/bin/env python3
"""mantlep.py: is the late-coda split of multiplet M03-24 (Z004) a source effect or a core effect?

Part of inner-core-repeaters by MV2, Munim, Inc. (made with AI agents). Licence: Apache-2.0.
Pre-registered in findings/2026-11-18-mantlep-prereg.md (written before this ran).

Background (exploratory, found after reading icr-m0324chk-v2): in the late PKIKP coda (3-10 s
after onset at YKA; -0.3..6 s at ILAR) the eight good Z004 events split into two groups:
  X = {2003-03-21, 2020-07-11, 2022-03-07, 2024-03-21, 2025-12-13}
  Y = {2008-02-17, 2011-02-20, 2014-12-09}
on BOTH arrays. That is exactly the split Wang et al.'s (2024) backtracking clock predicts
(turn 2008.5, back 2.5x slower: effective dates X 2001.5-2003.7, Y 2006.1-2008.1). It cuts
across the YKA hardware change, so hardware cannot make it. But a source difference (depth,
so the timing of pP/sP, location or mechanism) between two sub-clusters would make the same
split at both arrays.

The control: mantle P at nine long-running GSN stations 31-84 deg from the source, whose rays
never enter the core. If the split is a source effect it must appear there too. If it is a
core effect it must not.

Network: EarthScope FDSN dataselect + station (open data). Writes to --out-dir:
  mantlep_summary.json, mantlep_events.csv, mantlep_pairs.csv, plots/.
"""
import argparse
import csv
import io
import itertools
import json
import math
import os
import sys
import time
import traceback

import numpy as np

VERSION = "0.10.0"

# USGS ComCat (id, origin UTC, lat, lon, depth km, mb). us6000ittv (2022-10-15) is left out:
# beam SNR 6 at YKA, excluded by every earlier gate.
EVENTS = [
    ("usp000bt5x", "2003-03-21T10:47:47.620", -59.071, -25.47, 33.0, 4.8),
    ("usp000fz4d", "2008-02-17T19:32:09.100", -58.986, -25.697, 66.5, 4.9),
    ("usp000hute", "2011-02-20T17:37:37.320", -58.979, -25.857, 10.0, 5.0),
    ("usc000t4wj", "2014-12-09T19:23:36.400", -58.9931, -25.6003, 50.81, 5.1),
    ("us7000am0e", "2020-07-11T18:58:06.501", -59.0899, -25.7164, 57.13, 5.3),
    ("us6000h2wz", "2022-03-07T23:25:40.208", -59.0811, -25.7858, 35.0, 5.1),
    ("us6000mk65", "2024-03-21T02:47:21.923", -59.0115, -25.6366, 43.89, 5.2),
    ("us6000runq", "2025-12-13T14:04:33.132", -59.0011, -25.6684, 72.555, 5.0),
]
IDS = [e[0] for e in EVENTS]
GROUP_Y = ("usp000fz4d", "usp000hute", "usc000t4wj")   # the hypothesised 3-event group
TIGHT = ("us6000h2wz", "us6000mk65", "us6000runq")     # near-identical at YKA (diagnostic only)

# (net, sta, lat, lon): BHZ open 2003-2025 on EarthScope, 31-84 deg from the source (mantle P).
STATIONS = [
    ("II", "SUR", -32.3797, 20.8117),
    ("IU", "TSUM", -19.20227, 17.58376),
    ("II", "ASCN", -7.9327, -14.3601),
    ("II", "SHEL", -15.9594, -5.7455),
    ("IU", "RCBR", -5.8274, -35.9014),
    ("II", "NNA", -11.9875, -76.8422),
    ("IU", "CASY", -66.2792, 110.5354),
    ("IU", "SJG", 18.1091, -66.15),
    ("IU", "QSPA", -89.9289, 144.4382),
]
LOCS = ("00", "10", "")
DATASELECT = "https://service.earthscope.org/fdsnws/dataselect/1/query"
STATION = "https://service.earthscope.org/fdsnws/station/1/query"

FS = 20.0
BAND = (1.0, 2.0)            # same band as the YKA and ILAR measurements
PRE, POST = 30.0, 40.0       # grid around each event's own ak135 P time
FETCH_PRE, FETCH_POST = 120.0, 120.0
DIRECT = (-0.3, 0.9)         # pair alignment window (as at YKA)
WINDOWS = {"D_early": (1.0, 3.0), "D_late": (3.0, 10.0)}   # as YKA D_c13 / D_c310
NOISE = (-28.0, -6.0)
SNR_DIRECT_MIN = 5.0         # pre-registered gates, every event at a station
SNR_LATE_MIN = 2.0
HIT_RANK = 3                 # a station "hits" if the hypothesised split ranks <= 3 of 56
COARSE_S, FINE_S, MAX_PAIR_LAG_S = 3.0, 0.3, 0.25


# ----------------------------------------------------------------------------- pure helpers
def ncc(a, b):
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    if len(a) == 0 or len(a) != len(b):
        return float("nan")
    a = a - a.mean()
    b = b - b.mean()
    d = math.sqrt(float(np.dot(a, a)) * float(np.dot(b, b)))
    return float(np.dot(a, b) / d) if d > 0 else float("nan")


def shift(x, lag):
    """y[i] = x[i + lag], zero-padded."""
    x = np.asarray(x, dtype=float)
    y = np.zeros_like(x)
    n = len(x)
    if lag >= 0:
        if lag < n:
            y[:n - lag] = x[lag:]
    else:
        if -lag < n:
            y[-lag:] = x[:n + lag]
    return y


def best_lag(ref, x, i0, i1, maxlag, center=0):
    """Lag L in center +/- maxlag maximising ncc(ref[i0:i1], x[i0+L:i1+L])."""
    best = (None, -2.0)
    for lag in range(center - maxlag, center + maxlag + 1):
        a, b = i0 + lag, i1 + lag
        if a < 0 or b > len(x):
            continue
        c = ncc(ref[i0:i1], x[a:b])
        if np.isfinite(c) and c > best[1]:
            best = (lag, c)
    return best


def envelope(x, fs, smooth_s=0.25):
    w = max(1, int(round(smooth_s * fs)))
    return np.sqrt(np.convolve(np.asarray(x, dtype=float) ** 2, np.ones(w) / w, mode="same"))


def align(ref, x, fs, center, win=(-1.0, 3.0), coarse_s=COARSE_S, fine_s=FINE_S):
    """Coarse lag on the envelope (+/- coarse_s), then waveform refinement (+/- fine_s)."""
    i0 = center + int(round(win[0] * fs))
    i1 = center + int(round(win[1] * fs))
    lc, _ = best_lag(envelope(ref, fs), envelope(x, fs), i0, i1, int(round(coarse_s * fs)))
    if lc is None:
        return None, float("nan")
    lf, c = best_lag(ref, x, i0, i1, int(round(fine_s * fs)), center=lc)
    if lf is None:
        return None, float("nan")
    return lf, c


def envelope_onset(x, fs, center, half_s=3.0, frac=0.25, smooth_s=0.1):
    e = envelope(x, fs, smooth_s)
    a = max(0, int(center - half_s * fs))
    b = min(len(x), int(center + half_s * fs))
    seg = e[a:b]
    if len(seg) == 0 or not seg.max() > 0:
        return None
    return int(a + np.nonzero(seg >= frac * seg.max())[0][0])


def snr(x, fs, center, sig=(-3.0, 6.0), noise=NOISE):
    s = x[center + int(sig[0] * fs): center + int(sig[1] * fs)]
    n = x[center + int(noise[0] * fs): center + int(noise[1] * fs)]
    rn = float(np.sqrt(np.mean(n ** 2))) if len(n) else 0.0
    return float(np.max(np.abs(s)) / rn) if rn > 0 and len(s) else float("nan")


def late_snr(x, fs, onset, win=WINDOWS["D_late"], noise=NOISE):
    s = x[onset + int(win[0] * fs): onset + int(win[1] * fs)]
    n = x[onset + int(noise[0] * fs): onset + int(noise[1] * fs)]
    rn = float(np.sqrt(np.mean(n ** 2))) if len(n) else 0.0
    return float(np.sqrt(np.mean(s ** 2)) / rn) if rn > 0 and len(s) else float("nan")


def pair_D(xa, xb, onset, fs, windows=WINDOWS, direct=DIRECT, max_lag_s=MAX_PAIR_LAG_S):
    i0 = onset + int(round(direct[0] * fs))
    i1 = onset + int(round(direct[1] * fs))
    lag, c = best_lag(xa, xb, i0, i1, int(round(max_lag_s * fs)))
    if lag is None:
        return None
    out = {"cc_direct": c, "lag_s": lag / fs}
    for k, (w0, w1) in windows.items():
        j0 = onset + int(round(w0 * fs))
        j1 = onset + int(round(w1 * fs))
        out[k] = 1.0 - ncc(xa[j0:j1], xb[j0 + lag:j1 + lag])
    return out


def sep(D, ygroup):
    """mean D between the groups minus mean D within them (pooled)."""
    ys = set(ygroup)
    bw, wi = [], []
    n = D.shape[0]
    for i in range(n):
        for j in range(i + 1, n):
            v = D[i, j]
            if not np.isfinite(v):
                continue
            (wi if ((i in ys) == (j in ys)) else bw).append(v)
    if not bw or not wi:
        return float("nan")
    return float(np.mean(bw) - np.mean(wi))


def partition_rank(D, ygroup):
    """Rank (1 = best separated) of the split ygroup | rest among all splits of the same sizes."""
    k = len(ygroup)
    target = sep(D, tuple(sorted(ygroup)))
    combos = list(itertools.combinations(range(D.shape[0]), k))
    vals = [sep(D, c) for c in combos]
    rank = 1 + sum(1 for v in vals if np.isfinite(v) and v > target + 1e-12)
    return rank, len(vals), target, vals, combos


def decide(ranks, hit_rank=HIT_RANK):
    """Pre-registered network decision from the usable stations' ranks."""
    n = len(ranks)
    hits = sum(1 for r in ranks if r <= hit_rank)
    if n < 3:
        return "untestable", hits
    if hits >= int(math.ceil(n / 2.0)):
        return "S", hits
    if hits <= max(1, int(math.floor(0.2 * n))):
        return "N", hits
    return "inconclusive", hits


def matrix_from_pairs(rows, key, ids=IDS):
    idx = {e: i for i, e in enumerate(ids)}
    D = np.full((len(ids), len(ids)), np.nan)
    for r in rows:
        a, b = r.get("ev1"), r.get("ev2")
        if a in idx and b in idx:
            try:
                v = float(r[key])
            except (TypeError, ValueError, KeyError):
                continue
            D[idx[a], idx[b]] = D[idx[b], idx[a]] = v
    return D


def ygroup_idx(ids=IDS, group=GROUP_Y):
    return tuple(sorted(ids.index(g) for g in group))


def core_controls(path):
    """Ranks of the hypothesised split on the published YKA and ILAR D values."""
    with open(path, newline="") as f:
        rows = list(csv.DictReader(f))
    out = {}
    for key in ("yka_c310", "yka_c13", "ilar_long", "ilar_early"):
        D = matrix_from_pairs(rows, key)
        r, n, t, _, _ = partition_rank(D, ygroup_idx())
        out[key] = dict(rank=r, of=n, sep=t)
    return out


def network_combined(per_station_vals, hyp_index):
    """Mean over stations of the per-station z-scored separation; rank of the hypothesis."""
    Z = []
    for vals in per_station_vals:
        v = np.asarray(vals, dtype=float)
        sd = float(np.nanstd(v))
        if not sd > 0:
            continue
        Z.append((v - np.nanmean(v)) / sd)
    if not Z:
        return None
    m = np.nanmean(np.vstack(Z), axis=0)
    rank = 1 + int(np.sum(m > m[hyp_index] + 1e-12))
    return dict(rank=rank, of=len(m), z=float(m[hyp_index]), n_stations=len(Z))


# ----------------------------------------------------------------------------- data access
def get(session, url, params, tries=3):
    for k in range(tries):
        try:
            r = session.get(url, params=params, timeout=90)
            if r.status_code == 200 and r.content:
                return r.content
            if r.status_code in (204, 404):
                return None
        except Exception as ex:
            print("fetch error", url, repr(ex)[:200], flush=True)
        time.sleep(3 * (k + 1))
    return None


def iso(t):
    return t.strftime("%Y-%m-%dT%H:%M:%S")


def p_time(model, depth, dist):
    arr = model.get_travel_times(source_depth_in_km=max(0.0, depth), distance_in_degree=dist,
                                 phase_list=["P", "p", "Pdiff"])
    return min(a.time for a in arr) if arr else None


def load_trace(session, net, sta, t1, t2, cache):
    """Return (trace, inventory_or_None, loc) for the first location with data."""
    from obspy import read, read_inventory
    for loc in LOCS:
        lc = loc or "--"
        key = "%s.%s.%s.BHZ.%s" % (net, sta, lc, t1.strftime("%Y%m%dT%H%M%S"))
        mpath = os.path.join(cache, key + ".mseed")
        xpath = os.path.join(cache, key + ".xml")
        if not os.path.exists(mpath):
            body = get(session, DATASELECT, dict(net=net, sta=sta, loc=lc, cha="BHZ",
                                                 start=iso(t1), end=iso(t2)))
            if not body:
                continue
            with open(mpath, "wb") as f:
                f.write(body)
        try:
            st = read(mpath)
        except Exception:
            continue
        st.merge()
        if len(st) != 1 or np.ma.isMaskedArray(st[0].data):
            continue
        inv = None
        if not os.path.exists(xpath):
            body = get(session, STATION, dict(net=net, sta=sta, loc=lc, cha="BHZ", starttime=iso(t1),
                                               endtime=iso(t2), level="response"))
            if body:
                with open(xpath, "wb") as f:
                    f.write(body)
        if os.path.exists(xpath):
            try:
                inv = read_inventory(xpath)
            except Exception:
                inv = None
        return st[0], inv, loc
    return None, None, None


def process(tr, inv, t_start):
    tr = tr.copy()
    tr.detrend("demean")
    tr.detrend("linear")
    tr.taper(0.05)
    resp = False
    if inv is not None:
        try:
            tr.remove_response(inventory=inv, output="VEL", pre_filt=(0.2, 0.4, 7.0, 9.0))
            resp = True
        except Exception as ex:
            print("response removal failed", tr.id, repr(ex)[:150], flush=True)
    tr.filter("bandpass", freqmin=BAND[0], freqmax=BAND[1], corners=4, zerophase=False)
    n = int(round((PRE + POST) * FS))
    tr.interpolate(sampling_rate=FS, method="lanczos", a=12, starttime=t_start, npts=n)
    return np.asarray(tr.data, dtype=float), resp


# ----------------------------------------------------------------------------- main
def run_station(session, model, net, sta, slat, slon, cache, plots, ev_rows, pair_rows):
    from obspy import UTCDateTime
    from obspy.geodetics import locations2degrees
    center = int(round(PRE * FS))
    X, info = {}, {}
    for eid, ot, lat, lon, dep, mb in EVENTS:
        o = UTCDateTime(ot)
        dist = locations2degrees(lat, lon, slat, slon)
        tp = p_time(model, dep, dist)
        row = dict(station="%s.%s" % (net, sta), event=eid, time=ot, dist_deg=round(dist, 2),
                   status="", loc="", response_removed="", snr="", late_snr="", align_lag_s="", align_cc="")
        if tp is None:
            row["status"] = "no_P"
            ev_rows.append(row)
            continue
        tr, inv, loc = load_trace(session, net, sta, o + tp - FETCH_PRE, o + tp + FETCH_POST, cache)
        if tr is None:
            row["status"] = "no_data"
            ev_rows.append(row)
            continue
        try:
            x, resp = process(tr, inv, o + tp - PRE)
        except Exception as ex:
            row["status"] = "process_error:" + repr(ex)[:80]
            ev_rows.append(row)
            continue
        if not np.all(np.isfinite(x)) or not np.any(x):
            row["status"] = "bad_samples"
            ev_rows.append(row)
            continue
        X[eid] = x
        info[eid] = row
        row.update(status="ok", loc=loc or "--", response_removed=resp, snr=round(snr(x, FS, center), 2))
        ev_rows.append(row)
    res = dict(station="%s.%s" % (net, sta), n_fetched=len(X))
    if len(X) < 2:
        res["status"] = "too_few"
        return res
    ref = max(X, key=lambda e: info[e]["snr"] if np.isfinite(info[e]["snr"]) else -1)
    A = {}
    for e, x in X.items():
        if e == ref:
            A[e] = x
            info[e].update(align_lag_s=0.0, align_cc=1.0)
            continue
        lag, c = align(X[ref], x, FS, center)
        if lag is None:
            info[e]["status"] = "align_failed"
            continue
        A[e] = shift(x, lag)
        info[e].update(align_lag_s=round(lag / FS, 3), align_cc=round(c, 3))
    stack = np.mean([a / (np.max(np.abs(a[center - int(3 * FS):center + int(6 * FS)])) or 1.0)
                     for a in A.values()], axis=0)
    onset = envelope_onset(stack, FS, center)
    if onset is None:
        res["status"] = "no_onset"
        return res
    res.update(ref=ref, onset_minus_pred_s=round((onset - center) / FS, 2))
    for e, a in A.items():
        info[e]["late_snr"] = round(late_snr(a, FS, onset), 2)
    gates = {e: (e in A and np.isfinite(info[e]["snr"]) and info[e]["snr"] >= SNR_DIRECT_MIN
                 and np.isfinite(info[e]["late_snr"] or float("nan"))
                 and info[e]["late_snr"] >= SNR_LATE_MIN) for e in IDS}
    res["gates"] = gates
    n = len(IDS)
    D = {k: np.full((n, n), np.nan) for k in WINDOWS}
    for i, j in itertools.combinations(range(n), 2):
        a, b = IDS[i], IDS[j]
        if a not in A or b not in A:
            continue
        r = pair_D(A[a], A[b], onset, FS)
        if r is None:
            continue
        for k in WINDOWS:
            D[k][i, j] = D[k][j, i] = r[k]
        pair_rows.append(dict(station=res["station"], ev1=a, ev2=b, cc_direct=round(r["cc_direct"], 4),
                              lag_s=r["lag_s"], D_early=round(r["D_early"], 5), D_late=round(r["D_late"], 5),
                              group1="Y" if a in GROUP_Y else "X", group2="Y" if b in GROUP_Y else "X",
                              gated_ok=bool(gates[a] and gates[b])))
    res["usable"] = bool(all(gates.values()))
    yg = ygroup_idx()
    for k in WINDOWS:
        if np.all(np.isfinite(D[k][np.triu_indices(n, 1)])):
            r, m, t, vals, combos = partition_rank(D[k], yg)
            res[k] = dict(rank=r, of=m, sep=t, median_D=float(np.nanmedian(D[k][np.triu_indices(n, 1)])))
            if k == "D_late":
                res["_vals_late"] = vals
                res["_hyp_index"] = combos.index(yg)
    # diagnostics only
    ti = [IDS.index(e) for e in TIGHT]
    tight = [D["D_late"][a, b] for a, b in itertools.combinations(ti, 2)]
    other = [D["D_late"][a, b] for a, b in itertools.combinations(range(n), 2)
             if not (a in ti and b in ti)]
    res["tight_median_D_late"] = float(np.nanmedian(tight)) if tight else None
    res["other_median_D_late"] = float(np.nanmedian(other)) if other else None
    xs = [A[e] for e in IDS if e in A and e not in GROUP_Y]
    ys = [A[e] for e in IDS if e in A and e in GROUP_Y]
    if xs and ys:
        sx, sy = np.mean(xs, axis=0), np.mean(ys, axis=0)
        w0 = onset + int(WINDOWS["D_late"][0] * FS)
        w1 = onset + int(WINDOWS["D_late"][1] * FS)
        lag, c = best_lag(sx, sy, w0, w1, int(2.0 * FS))
        res["xy_stack_late"] = dict(cc0=ncc(sx[w0:w1], sy[w0:w1]), best_lag_s=None if lag is None else lag / FS,
                                    best_cc=c)
    try:
        plot_station(res["station"], A, onset, plots)
    except Exception as ex:
        print("plot failed", repr(ex)[:120], flush=True)
    res["status"] = "ok"
    return res


def plot_station(name, A, onset, plots):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(10, 8))
    t = (np.arange(len(next(iter(A.values())))) - onset) / FS
    for k, e in enumerate(IDS):
        if e not in A:
            continue
        a = A[e]
        m = np.max(np.abs(a[onset - int(2 * FS):onset + int(10 * FS)])) or 1.0
        col = "tab:red" if e in GROUP_Y else "tab:blue"
        ax.plot(t, a / m - 2.2 * k, color=col, lw=0.8)
        ax.text(10.5, -2.2 * k, dict((x[0], x[1][:10]) for x in EVENTS)[e], fontsize=8)
    ax.axvspan(*DIRECT, color="0.85")
    ax.axvspan(*WINDOWS["D_late"], color="orange", alpha=0.15)
    ax.set_xlim(-5, 13)
    ax.set_yticks([])
    ax.set_xlabel("s relative to stack P onset (blue X group, red Y group)")
    ax.set_title("%s: mantle P of Z004, 1-2 Hz velocity" % name)
    fig.tight_layout()
    fig.savefig(os.path.join(plots, "mantlep_%s.png" % name.replace(".", "_")), dpi=110)
    plt.close(fig)


def write_csv(path, rows):
    if not rows:
        open(path, "w").close()
        return
    keys = []
    for r in rows:
        for k in r:
            if k not in keys:
                keys.append(k)
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        w.writerows(rows)


def jsonable(o):
    if isinstance(o, dict):
        return {str(k): jsonable(v) for k, v in o.items() if not str(k).startswith("_")}
    if isinstance(o, (list, tuple)):
        return [jsonable(v) for v in o]
    if isinstance(o, (np.bool_,)):
        return bool(o)
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating, float)):
        v = float(o)
        return v if math.isfinite(v) else None
    return o


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", default="/data")
    ap.add_argument("--out-dir", default="/out")
    args = ap.parse_args()
    t0 = time.time()
    import requests
    from obspy.taup import TauPyModel
    cache = os.path.join(args.data_dir, "mantlep")
    plots = os.path.join(args.out_dir, "plots")
    os.makedirs(cache, exist_ok=True)
    os.makedirs(plots, exist_ok=True)
    here = os.path.dirname(os.path.abspath(__file__))
    s = dict(version=VERSION, started=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
             group_Y=GROUP_Y, band_hz=BAND, windows=WINDOWS, gates=dict(snr_direct=SNR_DIRECT_MIN,
             snr_late=SNR_LATE_MIN), hit_rank=HIT_RANK)
    try:
        s["core_controls"] = core_controls(os.path.join(here, "z004_core_D.csv"))
    except Exception as ex:
        s["core_controls"] = dict(error=repr(ex)[:200])
    session = requests.Session()
    session.headers["User-Agent"] = "inner-core-repeaters (MV2, Munim, Inc.)"
    model = TauPyModel("ak135")
    ev_rows, pair_rows, stations = [], [], []
    for net, sta, slat, slon in STATIONS:
        try:
            res = run_station(session, model, net, sta, slat, slon, cache, plots, ev_rows, pair_rows)
        except Exception as ex:
            traceback.print_exc()
            res = dict(station="%s.%s" % (net, sta), status="exception", error=repr(ex)[:300])
        stations.append(res)
        print(json.dumps(jsonable(res))[:1500], flush=True)
    usable = [r for r in stations if r.get("usable") and "D_late" in r]
    ranks = [r["D_late"]["rank"] for r in usable]
    decision, hits = decide(ranks)
    s.update(stations=stations, n_usable=len(usable), usable=[r["station"] for r in usable],
             ranks_D_late=ranks, hits=hits, decision=decision)
    s["ranks_D_early"] = [r["D_early"]["rank"] for r in usable if "D_early" in r]
    s["network_combined_D_late"] = network_combined([r["_vals_late"] for r in usable],
                                                    usable[0]["_hyp_index"]) if usable else None
    allr = [r for r in stations if "D_late" in r]
    s["all_stations_ungated"] = dict(ranks=[r["D_late"]["rank"] for r in allr],
                                     stations=[r["station"] for r in allr],
                                     label="includes stations failing the gates; not the decision")
    s["finished"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    s["minutes"] = round((time.time() - t0) / 60.0, 1)
    write_csv(os.path.join(args.out_dir, "mantlep_events.csv"), ev_rows)
    write_csv(os.path.join(args.out_dir, "mantlep_pairs.csv"), pair_rows)
    with open(os.path.join(args.out_dir, "mantlep_summary.json"), "w") as f:
        json.dump(jsonable(s), f, indent=1)
    print(json.dumps(jsonable({k: s.get(k) for k in ("core_controls", "n_usable", "ranks_D_late", "hits",
                                                     "decision", "network_combined_D_late")}))[:3000],
          flush=True)


if __name__ == "__main__":
    sys.exit(main())

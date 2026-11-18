#!/usr/bin/env python3
"""mantlep2.py: corrected mantle-P control of the M03-24 (Z004) late-coda split.

Part of inner-core-repeaters by MV2, Munim, Inc. (made with AI agents). Licence: Apache-2.0.
Pre-registered in findings/2026-11-18-mantlep2-prereg.md (pre-registration 11), written after
icr-mantlep-v1 came back "untestable" and before this script fetched any new station.

Why v2 (see findings/2026-11-18-mantlep-v1-results.md):
  * v1 anchored every window on the ak135 P time. Where the observed P is 2-6 s late (QSPA,
    CASY) the alignment window held little or no P, and events were misaligned by whole cycles.
  * v1's late-SNR noise window started 28 s before the stack onset, but the grid began only
    30 s before ak135 P, so when the onset sat at the -3 s search edge the noise slice was
    empty and the gate silently failed (SUR, SHEL, RCBR, CASY).
  * v1 required all eight events at a station; most long-running stations miss one or two.

v2 anchors on the observed onset of the best event, aligns with a wide envelope search then
a waveform refinement on the whole P train, measures SNR around the aligned onset, allows a
station with >= 6 gated events (>= 2 per group), and adds every IU/II/GT BHZ station 30-85 deg
from the source. The statistic and the S/N rule are those of pre-registration 10.

Writes to --out-dir: mantlep2_summary.json, mantlep2_events.csv, mantlep2_pairs.csv, plots/.
"""
import argparse
import csv
import itertools
import json
import math
import os
import sys
import time
import traceback

import numpy as np

VERSION = "0.13.0"

# USGS ComCat (id, origin UTC, lat, lon, depth km, mb). us6000ittv (2022-10-15) left out as in v1.
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
DATES = {e[0]: e[1][:10] for e in EVENTS}
GROUP_Y = ("usp000fz4d", "usp000hute", "usc000t4wj")   # hypothesised 3-event group
SRC_LAT, SRC_LON = -59.03, -25.68

# Stations whose v1 waveform plots were looked at before this pre-registration (QSPA, TSUM, CASY),
# and SPA, which sits at the same site as QSPA. Graded separately, never in the primary decision.
SEEN = ("IU.QSPA", "IU.TSUM", "IU.CASY", "IU.SPA")

NETS = "IU,II,GT"
DIST_RANGE = (30.0, 85.0)
MAX_STATIONS = 60
# Used only if the station service cannot be reached: the nine v1 stations.
FALLBACK_STATIONS = [
    ("II", "SUR", -32.3797, 20.8117), ("IU", "TSUM", -19.20227, 17.58376),
    ("II", "ASCN", -7.9327, -14.3601), ("II", "SHEL", -15.9594, -5.7455),
    ("IU", "RCBR", -5.8274, -35.9014), ("II", "NNA", -11.9875, -76.8422),
    ("IU", "CASY", -66.2792, 110.5354), ("IU", "SJG", 18.1091, -66.15),
    ("IU", "QSPA", -89.9289, 144.4382),
]
LOCS = ("00", "10", "")
DATASELECT = "https://service.earthscope.org/fdsnws/dataselect/1/query"
STATION = "https://service.earthscope.org/fdsnws/station/1/query"

FS = 20.0
BAND = (1.0, 2.0)
PRE, POST = 40.0, 50.0            # grid around each event's own ak135 P time
FETCH_PRE, FETCH_POST = 120.0, 120.0
ONSET_SEARCH = (-5.0, 15.0)       # where the reference onset is looked for (s from ak135 P)
RAW_SIG = (-5.0, 15.0)            # reference choice only
RAW_NOISE = (-38.0, -8.0)
ALIGN_WIN = (-0.5, 4.0)           # around the reference onset: the whole P train
COARSE_S, FINE_S = 6.0, 0.6
STACK_HALF_S = 3.0
SIG = (-1.0, 4.0)                 # gate SNR, around the aligned stack onset
NOISE = (-30.0, -8.0)
DIRECT = (-0.3, 0.9)
WINDOWS = {"D_early": (1.0, 3.0), "D_late": (3.0, 10.0)}
MAX_PAIR_LAG_S = 0.25
SNR_DIRECT_MIN = 5.0
SNR_LATE_MIN = 2.0
MIN_EVENTS, MIN_PER_GROUP = 6, 2
MIN_SPLITS = 15                   # so that a rank-1 hit has chance <= 1/15
HIT_P = 0.07                      # hit if rank / number of splits <= 0.07 (3/56 = 0.054 in v1)


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
    if i0 < 0 or i1 > len(ref) or i1 <= i0:
        return best
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


def onset_in(x, fs, a_s, b_s, base, frac=0.25, smooth_s=0.1):
    """First sample in [base + a_s, base + b_s) whose envelope reaches frac * the window maximum."""
    e = envelope(x, fs, smooth_s)
    a = max(0, base + int(round(a_s * fs)))
    b = min(len(x), base + int(round(b_s * fs)))
    seg = e[a:b]
    if len(seg) == 0 or not seg.max() > 0:
        return None
    return int(a + np.nonzero(seg >= frac * seg.max())[0][0])


def window_ratio(x, fs, base, sig, noise, kind="peak"):
    """Peak (or RMS) in base+sig over noise RMS in base+noise; nan if a window falls off the grid."""
    s0, s1 = base + int(round(sig[0] * fs)), base + int(round(sig[1] * fs))
    n0, n1 = base + int(round(noise[0] * fs)), base + int(round(noise[1] * fs))
    if s0 < 0 or n0 < 0 or s1 > len(x) or n1 > len(x) or s1 <= s0 or n1 <= n0:
        return float("nan")
    s, n = x[s0:s1], x[n0:n1]
    n = n[n != 0.0]                      # zero padding from an alignment shift is not noise
    if len(n) < int(2 * fs):
        return float("nan")
    rn = float(np.sqrt(np.mean(n ** 2)))
    if not rn > 0:
        return float("nan")
    v = float(np.max(np.abs(s))) if kind == "peak" else float(np.sqrt(np.mean(s ** 2)))
    return v / rn


def align2(ref, x, fs, anchor, win=ALIGN_WIN, coarse_s=COARSE_S, fine_s=FINE_S):
    """Envelope lag (+/- coarse_s) on the window around the reference onset, then waveform
    refinement (+/- fine_s) on the same window. Returns (lag in samples, cc) with x shifted by
    lag matching ref."""
    i0 = anchor + int(round(win[0] * fs))
    i1 = anchor + int(round(win[1] * fs))
    lc, _ = best_lag(envelope(ref, fs), envelope(x, fs), i0, i1, int(round(coarse_s * fs)))
    if lc is None:
        return None, float("nan")
    lf, c = best_lag(ref, x, i0, i1, int(round(fine_s * fs)), center=lc)
    if lf is None:
        return None, float("nan")
    return lf, c


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
        if j0 + lag < 0 or j1 + lag > len(xb) or j1 > len(xa):
            out[k] = float("nan")
            continue
        out[k] = 1.0 - ncc(xa[j0:j1], xb[j0 + lag:j1 + lag])
    return out


def sep(D, ygroup):
    """Mean D between the groups minus mean D within them (pooled)."""
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


def usable_subset(gates, ids=IDS, group=GROUP_Y, min_events=MIN_EVENTS, min_per_group=MIN_PER_GROUP):
    """Gated events in catalogue order, and whether the station is usable under prereg 11."""
    sub = [e for e in ids if gates.get(e)]
    ny = sum(1 for e in sub if e in group)
    nx = len(sub) - ny
    ok = (len(sub) >= min_events and ny >= min_per_group and nx >= min_per_group
          and n_unique_splits(len(sub), ny) >= MIN_SPLITS)
    return sub, ok


def n_unique_splits(n, k):
    """Distinct ways to split n events into groups of k and n - k."""
    c = math.comb(n, k)
    return c // 2 if 2 * k == n else c


def rank_subset(D, sub, ids=IDS, group=GROUP_Y, hit_p=HIT_P):
    """Rank of the hypothesised split restricted to the events in sub, among distinct splits."""
    idx = [ids.index(e) for e in sub]
    S = D[np.ix_(idx, idx)]
    yg = tuple(k for k, e in enumerate(sub) if e in group)
    r, m, t, vals, combos = partition_rank(S, yg)
    if 2 * len(yg) == len(sub):          # every split appears twice (a group and its complement)
        r, m = 1 + (r - 1) // 2, m // 2
    p = r / float(m)
    return dict(rank=r, of=m, p=p, sep=t, hit=bool(p <= hit_p + 1e-12),
                n_events=len(sub), n_y=len(yg))


def decide(hits):
    """Pre-registered network decision (rule of prereg 10) from per-station hit flags."""
    n = len(hits)
    h = sum(1 for x in hits if x)
    if n < 3:
        return "untestable", h
    if h >= int(math.ceil(n / 2.0)):
        return "S", h
    if h <= max(1, int(math.floor(0.2 * n))):
        return "N", h
    return "inconclusive", h


def rankdata(v):
    v = np.asarray(v, dtype=float)
    order = np.argsort(v, kind="mergesort")
    r = np.empty(len(v))
    r[order] = np.arange(len(v), dtype=float)
    # average ties
    sv = v[order]
    i = 0
    while i < len(sv):
        j = i
        while j + 1 < len(sv) and sv[j + 1] == sv[i]:
            j += 1
        if j > i:
            r[order[i:j + 1]] = (i + j) / 2.0
        i = j + 1
    return r


def spearman(a, b):
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    m = np.isfinite(a) & np.isfinite(b)
    if m.sum() < 5:
        return float("nan")
    return ncc(rankdata(a[m]), rankdata(b[m]))


def parse_station_text(text):
    """FDSN station text (level=station) -> sorted [(net, sta, lat, lon)], one per net.sta."""
    out = {}
    for line in (text or "").splitlines():
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        f = [s.strip() for s in line.split("|")]
        if len(f) < 4:
            continue
        try:
            lat, lon = float(f[2]), float(f[3])
        except ValueError:
            continue
        out.setdefault((f[0], f[1]), (lat, lon))
    return [(k[0], k[1], v[0], v[1]) for k, v in sorted(out.items())]


def gc_deg(lat1, lon1, lat2, lon2):
    p1, p2 = math.radians(lat1), math.radians(lat2)
    dl = math.radians(lon2 - lon1)
    c = math.sin(p1) * math.sin(p2) + math.cos(p1) * math.cos(p2) * math.cos(dl)
    return math.degrees(math.acos(max(-1.0, min(1.0, c))))


def load_core(path):
    with open(path, newline="") as f:
        rows = list(csv.DictReader(f))
    out = {}
    for r in rows:
        k = tuple(sorted((r["ev1"], r["ev2"])))
        out[k] = {c: float(r[c]) for c in ("yka_c13", "yka_c310", "ilar_early", "ilar_long")}
    return out


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


def discover(session):
    body = get(session, STATION, dict(net=NETS, cha="BHZ", lat=SRC_LAT, lon=SRC_LON,
                                      minradius=DIST_RANGE[0], maxradius=DIST_RANGE[1],
                                      starttime="2003-03-01", endtime="2025-12-31",
                                      level="station", format="text"))
    if not body:
        return None
    try:
        return parse_station_text(body.decode("utf-8", "replace"))
    except Exception:
        return None


def p_time(model, depth, dist):
    arr = model.get_travel_times(source_depth_in_km=max(0.0, depth), distance_in_degree=dist,
                                 phase_list=["P", "p", "Pdiff"])
    return min(a.time for a in arr) if arr else None


def load_trace(session, net, sta, t1, t2, cache):
    """(trace, inventory_or_None, loc) for the first location with one gap-free BHZ trace.
    Uses the same cache file names as mantlep.py v1."""
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


# ----------------------------------------------------------------------------- per station
def measure_station(X, raw_snr, fs=FS, center=None):
    """Everything after the waveforms are on the grid. X: {event: trace} (same grid, ak135 P at
    center). Returns (result dict, aligned traces, per-event info, D matrices)."""
    center = int(round(PRE * fs)) if center is None else center
    info = {e: dict(align_lag_s=None, align_cc=None, snr=None, late_snr=None) for e in X}
    res = dict(n_fetched=len(X))
    cands = [e for e in X if np.isfinite(raw_snr.get(e, float("nan")))]
    if len(cands) < 2:
        res["status"] = "too_few"
        return res, {}, info, None
    ref = max(cands, key=lambda e: raw_snr[e])
    anchor = onset_in(X[ref], fs, ONSET_SEARCH[0], ONSET_SEARCH[1], center)
    if anchor is None:
        res["status"] = "no_ref_onset"
        return res, {}, info, None
    A = {}
    for e, x in X.items():
        if e == ref:
            A[e] = x
            info[e].update(align_lag_s=0.0, align_cc=1.0)
            continue
        lag, c = align2(X[ref], x, fs, anchor)
        if lag is None:
            continue
        A[e] = shift(x, lag)
        info[e].update(align_lag_s=round(lag / fs, 3), align_cc=round(c, 3))
    norm = []
    for a in A.values():
        seg = a[anchor + int(round(ALIGN_WIN[0] * fs)):anchor + int(round(ALIGN_WIN[1] * fs))]
        m = float(np.max(np.abs(seg))) if len(seg) else 0.0
        norm.append(a / m if m > 0 else a)
    stack = np.mean(norm, axis=0)
    onset = onset_in(stack, fs, -STACK_HALF_S, STACK_HALF_S, anchor)
    if onset is None:
        res["status"] = "no_onset"
        return res, A, info, None
    res.update(ref=ref, ref_onset_minus_pred_s=round((anchor - center) / fs, 2),
               onset_minus_pred_s=round((onset - center) / fs, 2))
    gates = {}
    for e in IDS:
        if e not in A:
            gates[e] = False
            continue
        s = window_ratio(A[e], fs, onset, SIG, NOISE, "peak")
        l = window_ratio(A[e], fs, onset, WINDOWS["D_late"], NOISE, "rms")
        info[e].update(snr=round(s, 2) if np.isfinite(s) else None,
                       late_snr=round(l, 2) if np.isfinite(l) else None)
        gates[e] = bool(np.isfinite(s) and s >= SNR_DIRECT_MIN and np.isfinite(l) and l >= SNR_LATE_MIN)
    res["gates"] = gates
    n = len(IDS)
    D = {k: np.full((n, n), np.nan) for k in WINDOWS}
    D["cc_direct"] = np.full((n, n), np.nan)
    for i, j in itertools.combinations(range(n), 2):
        a, b = IDS[i], IDS[j]
        if a not in A or b not in A:
            continue
        r = pair_D(A[a], A[b], onset, fs)
        if r is None:
            continue
        for k in WINDOWS:
            D[k][i, j] = D[k][j, i] = r[k]
        D["cc_direct"][i, j] = D["cc_direct"][j, i] = r["cc_direct"]
    sub, ok = usable_subset(gates)
    res["gated_events"] = sub
    res["usable"] = bool(ok)
    if len(sub) >= 4 and 1 <= sum(1 for e in sub if e in GROUP_Y) < len(sub):
        for k in WINDOWS:
            res[k] = rank_subset(D[k], sub)
    res["onset"] = onset
    res["status"] = "ok"
    return res, A, info, D


def pair_contrast(D, sub):
    """Median D_late of gated between-group pairs over gated within-group pairs (diagnostic)."""
    bw, wi = [], []
    for a, b in itertools.combinations(sub, 2):
        v = D[IDS.index(a), IDS.index(b)]
        if not np.isfinite(v):
            continue
        if (a in GROUP_Y) == (b in GROUP_Y):
            wi.append(v)
        else:
            bw.append(v)
    if not bw or not wi:
        return None
    return dict(median_between=float(np.median(bw)), median_within=float(np.median(wi)),
                ratio=float(np.median(bw) / np.median(wi)) if np.median(wi) > 0 else None)


def core_corr(D, sub, core):
    """Spearman of mantle D_late against core D on the same gated pairs (diagnostic)."""
    out = {}
    for c in ("yka_c310", "ilar_long"):
        m, k = [], []
        for a, b in itertools.combinations(sub, 2):
            key = tuple(sorted((a, b)))
            v = D[IDS.index(a), IDS.index(b)]
            if key in core and np.isfinite(v):
                m.append(v)
                k.append(core[key][c])
        out[c] = dict(rho=spearman(m, k), n=len(m))
    return out


def run_station(session, model, net, sta, slat, slon, cache, plots, ev_rows, pair_rows, core):
    from obspy import UTCDateTime
    from obspy.geodetics import locations2degrees
    name = "%s.%s" % (net, sta)
    center = int(round(PRE * FS))
    X, raw, rows = {}, {}, {}
    for eid, ot, lat, lon, dep, mb in EVENTS:
        o = UTCDateTime(ot)
        dist = locations2degrees(lat, lon, slat, slon)
        row = dict(station=name, event=eid, date=DATES[eid], dist_deg=round(dist, 2), status="",
                   loc="", response_removed="", raw_snr="", snr="", late_snr="", align_lag_s="",
                   align_cc="", gated="")
        rows[eid] = row
        tp = p_time(model, dep, dist)
        if tp is None:
            row["status"] = "no_P"
            continue
        tr, inv, loc = load_trace(session, net, sta, o + tp - FETCH_PRE, o + tp + FETCH_POST, cache)
        if tr is None:
            row["status"] = "no_data"
            continue
        try:
            x, resp = process(tr, inv, o + tp - PRE)
        except Exception as ex:
            row["status"] = "process_error:" + repr(ex)[:80]
            continue
        if not np.all(np.isfinite(x)) or not np.any(x):
            row["status"] = "bad_samples"
            continue
        X[eid] = x
        raw[eid] = window_ratio(x, FS, center, RAW_SIG, RAW_NOISE, "peak")
        row.update(status="ok", loc=loc or "--", response_removed=resp,
                   raw_snr=round(raw[eid], 2) if np.isfinite(raw[eid]) else "")
    res, A, info, D = measure_station(X, raw)
    res["station"] = name
    res["seen_before_prereg"] = name in SEEN
    for e, row in rows.items():
        if e in info:
            for k in ("align_lag_s", "align_cc", "snr", "late_snr"):
                row[k] = "" if info[e][k] is None else info[e][k]
        row["gated"] = bool(res.get("gates", {}).get(e, False))
        ev_rows.append(row)
    if D is None:
        return res
    for i, j in itertools.combinations(range(len(IDS)), 2):
        a, b = IDS[i], IDS[j]
        if not np.isfinite(D["D_late"][i, j]):
            continue
        pair_rows.append(dict(station=name, ev1=a, ev2=b, cc_direct=round(D["cc_direct"][i, j], 4),
                              D_early=round(D["D_early"][i, j], 5), D_late=round(D["D_late"][i, j], 5),
                              group1="Y" if a in GROUP_Y else "X", group2="Y" if b in GROUP_Y else "X",
                              gated_ok=bool(res["gates"][a] and res["gates"][b])))
    sub = res.get("gated_events", [])
    if len(sub) >= 4:
        res["pair_contrast_late"] = pair_contrast(D["D_late"], sub)
        if core:
            res["core_corr"] = core_corr(D["D_late"], sub, core)
    try:
        plot_station(name, A, res["onset"], res.get("gates", {}), plots)
    except Exception as ex:
        print("plot failed", repr(ex)[:120], flush=True)
    return res


def plot_station(name, A, onset, gates, plots):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(10, 8))
    t = (np.arange(len(next(iter(A.values())))) - onset) / FS
    for k, e in enumerate(IDS):
        if e not in A:
            continue
        a = A[e]
        m = np.max(np.abs(a[max(0, onset - int(2 * FS)):onset + int(10 * FS)])) or 1.0
        col = "tab:red" if e in GROUP_Y else "tab:blue"
        ax.plot(t, a / m - 2.2 * k, color=col, lw=0.8, alpha=1.0 if gates.get(e) else 0.35)
        ax.text(10.5, -2.2 * k, DATES[e] + ("" if gates.get(e) else " (gated out)"), fontsize=8)
    ax.axvspan(*DIRECT, color="0.85")
    ax.axvspan(*WINDOWS["D_late"], color="orange", alpha=0.15)
    ax.set_xlim(-5, 13)
    ax.set_yticks([])
    ax.set_xlabel("s relative to stack P onset (blue X group, red Y group; faded = gated out)")
    ax.set_title("%s: mantle P of Z004, 1-2 Hz velocity (mantlep2)" % name)
    fig.tight_layout()
    fig.savefig(os.path.join(plots, "mantlep2_%s.png" % name.replace(".", "_")), dpi=100)
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


def grade(stations):
    usable = [r for r in stations if r.get("usable") and "D_late" in r]
    prim = [r for r in usable if not r.get("seen_before_prereg")]
    d1, h1 = decide([r["D_late"]["hit"] for r in prim])
    d2, h2 = decide([r["D_late"]["hit"] for r in usable])
    out = dict(primary=dict(decision=d1, n_usable=len(prim), hits=h1,
                            stations=[r["station"] for r in prim],
                            p=[r["D_late"]["p"] for r in prim]),
               all_stations=dict(decision=d2, n_usable=len(usable), hits=h2,
                                 stations=[r["station"] for r in usable],
                                 p=[r["D_late"]["p"] for r in usable],
                                 label="includes stations whose v1 plots were seen; secondary"))
    ratios = [r["pair_contrast_late"]["ratio"] for r in prim
              if r.get("pair_contrast_late") and r["pair_contrast_late"].get("ratio")]
    v1 = set("%s.%s" % (a, b) for a, b, _, _ in FALLBACK_STATIONS)
    new = [r for r in prim if r["station"] not in v1]
    d3, h3 = decide([r["D_late"]["hit"] for r in new])
    out["new_stations_only"] = dict(decision=d3, n_usable=len(new), hits=h3,
                                    stations=[r["station"] for r in new],
                                    label="excludes all nine v1 stations; sensitivity")
    out["primary_median_contrast_ratio"] = float(np.median(ratios)) if ratios else None
    rhos = [r["core_corr"]["yka_c310"]["rho"] for r in prim
            if r.get("core_corr") and np.isfinite(r["core_corr"]["yka_c310"]["rho"])]
    out["primary_median_rho_vs_yka_c310"] = float(np.median(rhos)) if rhos else None
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", default="/data")
    ap.add_argument("--out-dir", default="/out")
    ap.add_argument("--time-budget-min", type=float, default=300.0)
    args = ap.parse_args()
    t0 = time.time()
    import requests
    from obspy.taup import TauPyModel
    cache = os.path.join(args.data_dir, "mantlep")       # shared with v1
    plots = os.path.join(args.out_dir, "plots")
    os.makedirs(cache, exist_ok=True)
    os.makedirs(plots, exist_ok=True)
    here = os.path.dirname(os.path.abspath(__file__))
    try:
        core = load_core(os.path.join(here, "z004_core_D.csv"))
    except Exception as ex:
        print("core table not loaded", repr(ex)[:120], flush=True)
        core = {}
    session = requests.Session()
    session.headers["User-Agent"] = "inner-core-repeaters (MV2, Munim, Inc.)"
    found = discover(session)
    source = "station-service"
    if not found:
        found, source = FALLBACK_STATIONS, "fallback"
    cand = []
    for net, sta, lat, lon in found:
        d = gc_deg(SRC_LAT, SRC_LON, lat, lon)
        if DIST_RANGE[0] <= d <= DIST_RANGE[1]:
            cand.append((d, net, sta, lat, lon))
    cand.sort()
    cand = cand[:MAX_STATIONS]
    s = dict(version=VERSION, started=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
             station_source=source, n_candidates=len(cand), group_Y=GROUP_Y, seen=SEEN,
             band_hz=BAND, windows=WINDOWS, gates=dict(snr_direct=SNR_DIRECT_MIN, snr_late=SNR_LATE_MIN,
             min_events=MIN_EVENTS, min_per_group=MIN_PER_GROUP), hit_p=HIT_P)
    model = TauPyModel("ak135")
    model_ev, pair_rows, stations = [], [], []
    skipped = []
    for d, net, sta, lat, lon in cand:
        if (time.time() - t0) / 60.0 > args.time_budget_min:
            skipped.append("%s.%s" % (net, sta))
            continue
        try:
            res = run_station(session, model, net, sta, lat, lon, cache, plots, model_ev, pair_rows, core)
        except Exception as ex:
            traceback.print_exc()
            res = dict(station="%s.%s" % (net, sta), status="exception", error=repr(ex)[:300])
        res["dist_deg"] = round(d, 2)
        stations.append(res)
        print(json.dumps(jsonable(res))[:1500], flush=True)
    s.update(stations=stations, skipped_time_budget=skipped)
    s["grade"] = grade(stations)
    s["finished"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    s["minutes"] = round((time.time() - t0) / 60.0, 1)
    write_csv(os.path.join(args.out_dir, "mantlep2_events.csv"), model_ev)
    write_csv(os.path.join(args.out_dir, "mantlep2_pairs.csv"), pair_rows)
    with open(os.path.join(args.out_dir, "mantlep2_summary.json"), "w") as f:
        json.dump(jsonable(s), f, indent=1)
    print(json.dumps(jsonable(s["grade"]))[:3000], flush=True)


if __name__ == "__main__":
    sys.exit(main())

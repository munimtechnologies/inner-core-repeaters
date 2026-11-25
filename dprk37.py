#!/usr/bin/env python3
# Copyright 2026 MV2, Munim, Inc.
# Licensed under the Apache License, Version 2.0 (Apache-2.0). See LICENSE.
# Made with AI agents working in MV2's name.
"""dprk37.py - pre-registration 37: do inner-core waves from the North Korean explosions change?

Implements findings/2026-11-25-prereg37-dprk-doublets-plan.md exactly. For each pair of tests and
each gated station 150-179.5 deg from Punggye-ri, it aligns PKPab (AB, outer core only) between
the two tests, measures the PKIKP (DF, inner core) lag relative to AB (dt) and DF's decorrelation
beyond noise relative to AB's (X), then labels G1 (2016-01 / 2016-09) and G2 (2013-02 / 2016-01)
with the fixed rules, and runs P37b on stations whose DF crosses the inner-core boundary near the
points sampled by Zhang & Wen's multiplet F at YKA and ILAR.

Predicted phase times come from one fixed site (SRC) for both tests, so catalogue location noise
(kilometres) never enters dt; the real relative positions (about 100 m) give a cosine-of-azimuth
term that the T fit absorbs.

Top-level imports are stdlib only (the tests run without numpy/obspy).

Changes
- 0.1.1: no change to the method. 0.1.0 never reached the repo because its unit test of the
  P37b permutation test expected p < 0.05 from 4 vs 4 tied values, where the exact answer is
  17/70 (a median-difference permutation test cannot reach p <= 0.05 with four tied values per
  group). The test now uses 10 vs 10 separated values and checks the 4 vs 4 case against 17/70.
"""

import argparse
import csv
import json
import math
import os
import random
import statistics
import sys
import time
import traceback
import urllib.request

VERSION = "0.1.1"
NAN = float("nan")
REPO = "https://raw.githubusercontent.com/munimtechnologies/inner-core-repeaters/"
PROBE_COMMIT = "0fdf6cbdc7bd681cd5ca967a0cb0876725333fad"   # dprk37probe.py 0.1.1
PROBE_VERSION = "0.1.1"
PREREG_URL = REPO + "main/findings/2026-11-25-prereg37-dprk-doublets-plan.md"
UA = "inner-core-repeaters (MV2, Munim, Inc.)"
FDSN = "https://service.earthscope.org"

SRC_LAT, SRC_LON = 41.30, 129.05
PAIRS = [
    ("G1", "us10004bnm", "us10006n8a"),
    ("G2", "usc000f5t0", "us10004bnm"),
]
EVENT_IDS = ("usp000eurb", "usp000gxgc", "usc000f5t0", "us10004bnm", "us10006n8a", "us2000aert")
PRIMARY = (0.7, 2.0)
SECONDARY = (1.0, 2.0)
FS = 100.0
SNR_MIN = 5.0
CC_AB_MIN = 0.80
AB_WIN = (-0.5, 5.0)
AB_LAG = 1.0
DF_WIN = (-0.5, 4.0)
DF_LAG = 0.4
N_MIN = 8
N_BOOT = 2000
N_PERM = 10000
SEED = 37
W_CHANGED, W_FLOOR = 0.03, 0.01
T_SHIFT = 0.020
PATCH_RADIUS = 10.0
PATCH_MIN = 4
POS_BOUND_DEG = 0.0027        # 300 m
ICB_DEPTH = 5153.5            # ak135
ICB_TOL = 5.0
F_CENTROID = (-58.7887, -25.3637)
REF_ARRAYS = {"YKA": (62.493, -114.605), "ILAR": (64.7714, -146.8866)}


# ----------------------------------------------------------------------------- pure helpers
def isf(x):
    try:
        return math.isfinite(float(x))
    except (TypeError, ValueError):
        return False


def _unit(lat, lon):
    p, l = math.radians(lat), math.radians(lon)
    return (math.cos(p) * math.cos(l), math.cos(p) * math.sin(l), math.sin(p))


def gc_deg(lat1, lon1, lat2, lon2):
    """Great-circle distance in degrees, atan2(|a x b|, a.b): accurate at 0 and 180 deg."""
    a, b = _unit(lat1, lon1), _unit(lat2, lon2)
    cx = a[1] * b[2] - a[2] * b[1]
    cy = a[2] * b[0] - a[0] * b[2]
    cz = a[0] * b[1] - a[1] * b[0]
    dot = a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
    return math.degrees(math.atan2(math.sqrt(cx * cx + cy * cy + cz * cz), dot))


def azimuth(lat1, lon1, lat2, lon2):
    """Forward azimuth (deg, 0-360) from point 1 to point 2."""
    p1, p2 = math.radians(lat1), math.radians(lat2)
    dl = math.radians(lon2 - lon1)
    y = math.sin(dl) * math.cos(p2)
    x = math.cos(p1) * math.sin(p2) - math.sin(p1) * math.cos(p2) * math.cos(dl)
    return math.degrees(math.atan2(y, x)) % 360.0


def gc_point(lat, lon, az, dist):
    """Point reached from (lat, lon) along azimuth az (deg) after dist (deg)."""
    p1, l1 = math.radians(lat), math.radians(lon)
    th, d = math.radians(az), math.radians(dist)
    s = math.sin(p1) * math.cos(d) + math.cos(p1) * math.sin(d) * math.cos(th)
    p2 = math.asin(max(-1.0, min(1.0, s)))
    l2 = l1 + math.atan2(math.sin(th) * math.sin(d) * math.cos(p1),
                         math.cos(d) - math.sin(p1) * math.sin(p2))
    lon2 = (math.degrees(l2) + 180.0) % 360.0 - 180.0
    return math.degrees(p2), lon2


def parabolic_peak(ym1, y0, yp1):
    """Sub-sample offset (samples) and value of the parabola through three points."""
    den = ym1 - 2.0 * y0 + yp1
    if not den < 0:
        return 0.0, y0
    delta = 0.5 * (ym1 - yp1) / den
    return delta, y0 - 0.25 * (ym1 - yp1) * delta


def signal_fraction(r):
    """Share of window power that is signal for an rms SNR r (signal window includes noise)."""
    if not isf(r) or float(r) <= 0:
        return 0.0
    return max(0.0, 1.0 - 1.0 / (float(r) ** 2))


def expected_cc(ra, rb):
    """Expected best correlation of two noisy copies of one signal with SNRs ra and rb."""
    return math.sqrt(signal_fraction(ra) * signal_fraction(rb))


def excess_x(cc_df, cc_ab, snr_df_a, snr_df_b, snr_ab_a, snr_ab_b):
    """X = (D_DF - E_DF) - (D_AB - E_AB), with D = 1 - cc and E = 1 - expected_cc."""
    x_df = (1.0 - cc_df) - (1.0 - expected_cc(snr_df_a, snr_df_b))
    x_ab = (1.0 - cc_ab) - (1.0 - expected_cc(snr_ab_a, snr_ab_b))
    return x_df - x_ab


def quantile(xs, q):
    """Linear-interpolation quantile of a list (sorted inside)."""
    v = sorted(xs)
    if not v:
        return NAN
    pos = q * (len(v) - 1)
    i = int(math.floor(pos))
    if i >= len(v) - 1:
        return v[-1]
    f = pos - i
    return v[i] + f * (v[i + 1] - v[i])


def solve3(m, y):
    """Solve a 3x3 linear system with partial pivoting; None if singular."""
    a = [list(m[i]) + [y[i]] for i in range(3)]
    scale = max(abs(a[i][j]) for i in range(3) for j in range(3)) or 1.0
    for c in range(3):
        piv = max(range(c, 3), key=lambda r: abs(a[r][c]))
        if abs(a[piv][c]) < 1e-12 * scale:
            return None
        a[c], a[piv] = a[piv], a[c]
        for r in range(3):
            if r != c:
                f = a[r][c] / a[c][c]
                for k in range(c, 4):
                    a[r][k] -= f * a[c][k]
    return [a[i][3] / a[i][i] for i in range(3)]


def fit_cos(rows):
    """rows: [(baz_deg, dt)]. Least-squares dt = c + a cos(baz) + b sin(baz) -> (c, a, b) or None."""
    if len(rows) < 3:
        return None
    m = [[0.0] * 3 for _ in range(3)]
    y = [0.0] * 3
    for baz, v in rows:
        g = (1.0, math.cos(math.radians(baz)), math.sin(math.radians(baz)))
        for i in range(3):
            y[i] += g[i] * v
            for j in range(3):
                m[i][j] += g[i] * g[j]
    return solve3(m, y)


def bootstrap(items, stat, n=N_BOOT, seed=SEED):
    """Statistic over n resamples of items (with replacement); None results are skipped."""
    rng = random.Random(seed)
    out = []
    k = len(items)
    if k == 0:
        return out
    for _ in range(n):
        s = stat([items[rng.randrange(k)] for _ in range(k)])
        if s is not None and isf(s):
            out.append(float(s))
    return out


def stat_median_x(sample):
    return statistics.median([r["X"] for r in sample])


def stat_c(sample):
    f = fit_cos([(r["baz"], r["dt"]) for r in sample])
    return None if f is None else f[0]


def label_w(m, lo, hi):
    if m >= W_CHANGED and lo > W_FLOOR:
        return "changed"
    if hi < W_CHANGED:
        return "same"
    return "inconclusive"


def label_t(c, lo, hi):
    if abs(c) >= T_SHIFT and (lo > 0 or hi < 0):
        return "shift"
    if lo >= -T_SHIFT and hi <= T_SHIFT:
        return "no shift"
    return "inconclusive"


def headline(w, t):
    if w == "changed" or t == "shift":
        return "DF changed"
    if w == "same" and t == "no shift":
        return "no change at this level"
    return "inconclusive"


def label_patch(d, p):
    if d >= W_CHANGED and p <= 0.05:
        return "patch excess"
    if d <= W_FLOOR:
        return "no patch excess"
    return "inconclusive"


def perm_p(patch, other, n=N_PERM, seed=SEED):
    """Observed d = median(patch) - median(other) and the share of label permutations >= d."""
    d = statistics.median(patch) - statistics.median(other)
    pool = list(patch) + list(other)
    k = len(patch)
    rng = random.Random(seed)
    hits = 0
    for _ in range(n):
        rng.shuffle(pool)
        if statistics.median(pool[:k]) - statistics.median(pool[k:]) >= d:
            hits += 1
    return d, hits / float(n)


def icb_crossings(points, icb=ICB_DEPTH, tol=ICB_TOL):
    """points: [(dist_deg, depth_km)] along a PKIKP ray -> (first, last) ICB distance or None."""
    hit = [d for d, z in points if abs(z - icb) <= tol]
    if not hit:
        return None
    return min(hit), max(hit)


def is_patch(points, refs, radius=PATCH_RADIUS):
    """True if any (lat, lon) in points lies within radius deg of any (lat, lon) in refs."""
    for p in points:
        for r in refs:
            if gc_deg(p[0], p[1], r[0], r[1]) <= radius:
                return True
    return False


def gate(ra, rb, band_label, clear):
    """Station gates 1-4 of P37 for one pair. Returns '' if passed, else the first failed gate."""
    for r in (ra, rb):
        if r is None or r.get("status") != "ok":
            return "record not ok"
        if r.get("covers") is not True and str(r.get("covers")) != "True":
            return "record does not cover window"
    if ra.get("channel") != rb.get("channel") or ra.get("epoch_start") != rb.get("epoch_start"):
        return "channel or epoch differs"
    if not clear:
        return "windows not clear"
    for r in (ra, rb):
        for ph in ("df", "ab"):
            v = r.get("snr_%s_%s" % (ph, band_label))
            if not (isf(v) and float(v) >= SNR_MIN):
                return "snr"
    return ""


def all_pairs():
    """G1, G2, then every other pair of the six tests as report-only (R:...)."""
    out = list(PAIRS)
    have = {(a, b) for _, a, b in PAIRS}
    for i in range(len(EVENT_IDS)):
        for j in range(i + 1, len(EVENT_IDS)):
            a, b = EVENT_IDS[i], EVENT_IDS[j]
            if (a, b) not in have:
                out.append(("R:%s/%s" % (a, b), a, b))
    return out


def pick_branch_p(arrivals, t):
    """Ray parameter (s/deg) of the arrival whose time equals t (None if t is None)."""
    if t is None:
        return None
    best = None
    for name, tt, p in arrivals:
        if abs(tt - t) < 1e-6:
            best = p
            break
    return best


def summarise_pair(meas, patch_flags):
    """meas: list of dicts with X, dt, baz, station. Applies sections 6 and 7 of P37."""
    n = len(meas)
    out = dict(n=n)
    if n < N_MIN:
        out.update(headline="untestable", W="untestable", T="untestable",
                   P37b="untestable")
        return out
    mx = statistics.median([r["X"] for r in meas])
    bx = bootstrap(meas, stat_median_x)
    lo, hi = quantile(bx, 0.025), quantile(bx, 0.975)
    w = label_w(mx, lo, hi)
    f = fit_cos([(r["baz"], r["dt"]) for r in meas])
    if f is None:
        t, c, clo, chi = "untestable", NAN, NAN, NAN
        a = b = NAN
    else:
        c, a, b = f
        bc = bootstrap(meas, stat_c)
        clo, chi = quantile(bc, 0.025), quantile(bc, 0.975)
        t = label_t(c, clo, chi)
    out.update(median_X=mx, X_ci=[lo, hi], W=w, c=c, a_cos=a, b_sin=b, c_ci=[clo, chi], T=t,
               median_dt=statistics.median([r["dt"] for r in meas]), headline=headline(w, t))
    pat = [r["X"] for r in meas if patch_flags.get(r["station"])]
    oth = [r["X"] for r in meas if not patch_flags.get(r["station"])]
    out["n_patch"], out["n_nonpatch"] = len(pat), len(oth)
    if len(pat) < PATCH_MIN or len(oth) < PATCH_MIN:
        out["P37b"] = "untestable"
    else:
        d, p = perm_p(pat, oth)
        out.update(P37b=label_patch(d, p), patch_d=d, patch_p=p)
    return out


def write_csv(path, rows):
    keys = []
    for r in rows:
        for k in r:
            if k not in keys:
                keys.append(k)
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=keys or ["empty"])
        w.writeheader()
        for r in rows:
            w.writerow(r)


# ----------------------------------------------------------------------------- I/O helpers
def fetch(url, tries=4, wait=10):
    for k in range(tries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": UA})
            with urllib.request.urlopen(req, timeout=60) as r:
                if r.status == 200:
                    return r.read().decode("utf-8", "replace")
        except Exception as ex:  # noqa: BLE001
            print("fetch error", url, repr(ex)[:200], flush=True)
        time.sleep(wait * (k + 1))
    return None


def load_pinned(dest, name, commit, version):
    os.makedirs(dest, exist_ok=True)
    p = os.path.join(dest, name + ".py")
    if not os.path.exists(p):
        body = fetch(REPO + commit + "/" + name + ".py")
        if not body or ('VERSION = "%s"' % version) not in body:
            return None
        compile(body, name + ".py", "exec")
        with open(p, "w", encoding="utf-8") as f:
            f.write(body)
    if dest not in sys.path:
        sys.path.insert(0, dest)
    return __import__(name)


# ----------------------------------------------------------------------------- numerics
def prep(tr, band):
    """Processed copy: demean, detrend, 5% taper, causal 4-pole band-pass, cubic to 100 Hz."""
    x = tr.copy()
    x.data = x.data.astype("float64")
    x.detrend("demean")
    x.detrend("linear")
    x.taper(0.05)
    x.filter("bandpass", freqmin=band[0], freqmax=band[1], corners=4, zerophase=False)
    if abs(x.stats.sampling_rate - FS) > 1e-6:
        x.interpolate(sampling_rate=FS, method="cubic")
    return x


def xcorr(xa, xb, ta0, tb0, win, k0, kmax):
    """Normalised cc of xa[ta0+win] with xb[tb0+win+lag], lag = (k0 + k)/FS for |k| <= kmax.
    Returns (lag_s, cc_peak, at_edge) or None when a window falls outside a record."""
    import numpy as np
    n = int(round((win[1] - win[0]) * FS))
    ia = int(round(((ta0 + win[0]) - xa.stats.starttime) * FS))
    ib0 = int(round(((tb0 + win[0]) - xb.stats.starttime) * FS))
    if ia < 0 or ia + n > len(xa.data):
        return None
    a = xa.data[ia:ia + n]
    ea = float(np.dot(a, a))
    ccs = []
    for k in range(k0 - kmax, k0 + kmax + 1):
        ib = ib0 + k
        if ib < 0 or ib + n > len(xb.data):
            return None
        b = xb.data[ib:ib + n]
        den = math.sqrt(ea * float(np.dot(b, b)))
        ccs.append(float(np.dot(a, b)) / den if den > 0 else 0.0)
    j = max(range(len(ccs)), key=lambda i: ccs[i])
    edge = j == 0 or j == len(ccs) - 1
    delta, peak = (0.0, ccs[j]) if edge else parabolic_peak(ccs[j - 1], ccs[j], ccs[j + 1])
    return (k0 - kmax + j + delta) / FS, min(1.0, peak), edge


def measure(sta, ra, rb, ea, eb, tt, band, band_label, cache):
    """One station, one pair: dt, cc, X or a drop reason."""
    from obspy import read
    out = dict(station=sta, channel=ra["channel"])
    xs = []
    for r, ev in ((ra, ea), (rb, eb)):
        path = os.path.join(cache, ev["id"], "%s.%s.mseed" % (sta, r["channel"]))
        if not os.path.exists(path):
            out["drop"] = "cached record missing"
            return out
        st = read(path)
        st.merge(method=0)
        if len(st) != 1:
            out["drop"] = "merge left %d traces" % len(st)
            return out
        tr = st[0]
        import numpy as np
        if np.ma.isMaskedArray(tr.data):
            if np.ma.getmaskarray(tr.data).any():
                out["drop"] = "gap"
                return out
            tr.data = np.asarray(tr.data)
        xs.append(prep(tr, band))
    ta_ab, tb_ab = ea["origin"] + tt["ab"], eb["origin"] + tt["ab"]
    ab = xcorr(xs[0], xs[1], ta_ab, tb_ab, AB_WIN, 0, int(round(AB_LAG * FS)))
    if ab is None:
        out["drop"] = "AB window outside record"
        return out
    out["tau_ab"], out["cc_ab"] = round(ab[0], 4), round(ab[1], 4)
    if ab[2]:
        out["drop"] = "AB lag at edge"
        return out
    if ab[1] < CC_AB_MIN:
        out["drop"] = "cc_AB < 0.80"
        return out
    ta_df, tb_df = ea["origin"] + tt["df"], eb["origin"] + tt["df"]
    df = xcorr(xs[0], xs[1], ta_df, tb_df, DF_WIN, int(round(ab[0] * FS)),
               int(round(DF_LAG * FS)))
    if df is None:
        out["drop"] = "DF window outside record"
        return out
    out["tau_df"], out["cc_df"] = round(df[0], 4), round(df[1], 4)
    if df[2]:
        out["drop"] = "DF lag at edge"
        return out
    out["dt"] = df[0] - ab[0]
    out["X"] = excess_x(df[1], ab[1], ra["snr_df_" + band_label], rb["snr_df_" + band_label],
                        ra["snr_ab_" + band_label], rb["snr_ab_" + band_label])
    out["drop"] = ""
    return out


def pierce_points(model, dist, src, az):
    """ICB entry/exit (lat, lon) of PKIKP for a surface source at src along azimuth az."""
    arr = model.get_pierce_points(source_depth_in_km=0.0, distance_in_degree=dist,
                                  phase_list=["PKIKP"])
    arr = [a for a in arr if a.name in ("PKIKP", "PKPdf")]
    if not arr:
        return None
    a = min(arr, key=lambda q: q.time)
    pts = [(math.degrees(float(q["dist"])), float(q["depth"])) for q in a.pierce]
    cr = icb_crossings(pts)
    if cr is None:
        return None
    return [gc_point(src[0], src[1], az, cr[0]), gc_point(src[0], src[1], az, cr[1])]


def run(args):
    t0 = time.time()
    os.makedirs(args.out_dir, exist_ok=True)
    s = dict(version=VERSION, started=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
             prereg=PREREG_URL, probe_commit=PROBE_COMMIT)
    pre = fetch(PREREG_URL, tries=2, wait=5)
    s["prereg_public_at_start"] = bool(pre and "Pre-registration 37" in pre)
    pr = load_pinned(os.path.join(args.out_dir, "_pinned"), "dprk37probe", PROBE_COMMIT,
                     PROBE_VERSION)
    if pr is None:
        s["error"] = "could not load dprk37probe.py at the pinned commit"
        json.dump(s, open(os.path.join(args.out_dir, "dprk37_summary.json"), "w"), indent=1)
        return 1
    from obspy import UTCDateTime
    from obspy.clients.fdsn import Client
    from obspy.taup import TauPyModel
    client = Client(base_url=FDSN, _discover_services=False, user_agent=UA, timeout=180)
    model = TauPyModel("ak135")
    cache = os.path.join(args.data_dir, "dprk37")
    os.makedirs(cache, exist_ok=True)
    evs = {e[0]: dict(id=e[0], origin=UTCDateTime(e[1]), date=e[1][:10]) for e in pr.EVENTS}
    rows, coords = {}, {}
    for e in pr.EVENTS:
        try:
            for r in pr.probe_event(client, model, e, cache, {}, s):
                rows[(e[0], r["station"])] = r
        except Exception as ex:  # noqa: BLE001
            traceback.print_exc()
            s.setdefault("event_errors", {})[e[0]] = repr(ex)[:300]
        try:
            o = UTCDateTime(e[1])
            inv = client.get_stations(latitude=SRC_LAT, longitude=SRC_LON,
                                      minradius=pr.DIST_RANGE[0], maxradius=pr.DIST_RANGE[1],
                                      channel="BHZ,HHZ", starttime=o - 60, endtime=o + 3600,
                                      level="station")
            for net in inv:
                for sta in net:
                    coords.setdefault((e[0], "%s.%s" % (net.code, sta.code)),
                                      (sta.latitude, sta.longitude))
        except Exception as ex:  # noqa: BLE001
            s.setdefault("station_errors", {})[e[0]] = repr(ex)[:200]
    # reference ICB points (Zhang & Wen multiplet F to YKA and ILAR)
    refs = []
    for name, (la, lo) in REF_ARRAYS.items():
        d = gc_deg(F_CENTROID[0], F_CENTROID[1], la, lo)
        az = azimuth(F_CENTROID[0], F_CENTROID[1], la, lo)
        pp = pierce_points(model, d, F_CENTROID, az)
        s.setdefault("reference_icb_points", {})[name] = dict(dist=d, points=pp)
        if pp:
            refs.extend(pp)
    # per-station geometry from the fixed site
    geo, ttc = {}, {}
    for (eid, sta), (la, lo) in coords.items():
        key = (sta, round(la, 4), round(lo, 4))
        if key in geo:
            continue
        d = gc_deg(SRC_LAT, SRC_LON, la, lo)
        az = azimuth(SRC_LAT, SRC_LON, la, lo)
        dk = round(d, 2)
        if dk not in ttc:
            arr = model.get_travel_times(source_depth_in_km=0.0, distance_in_degree=dk,
                                         phase_list=["PKIKP", "PKP"])
            al = [(a.name, a.time, a.ray_param_sec_degree) for a in arr]
            b = pr.branches([(n, t) for n, t, _ in al], dk)
            b["p_df"] = pick_branch_p(al, b["df"])
            b["p_ab"] = pick_branch_p(al, b["ab"])
            b["clear"] = pr.windows_clear(b)
            ttc[dk] = b
        pp = None
        try:
            pp = pierce_points(model, d, (SRC_LAT, SRC_LON), az)
        except Exception as ex:  # noqa: BLE001
            s.setdefault("pierce_errors", []).append("%s %s" % (sta, repr(ex)[:80]))
        geo[key] = dict(station=sta, lat=la, lon=lo, dist=d, az=az,
                        baz=azimuth(la, lo, SRC_LAT, SRC_LON), tt=ttc[dk], icb=pp,
                        patch=bool(pp and refs and is_patch(pp, refs)))
    geo_rows = [dict(station=g["station"], lat=g["lat"], lon=g["lon"], dist=round(g["dist"], 3),
                     baz=round(g["baz"], 2), icb=json.dumps(g["icb"]), patch=g["patch"],
                     clear=g["tt"]["clear"]) for g in geo.values()]
    write_csv(os.path.join(args.out_dir, "dprk37_geometry.csv"), geo_rows)
    results, meas_rows = {}, []
    for band, tag in ((PRIMARY, "primary"), (SECONDARY, "secondary")):
        bl = "%.1f-%.1f" % band
        for label, a, b in all_pairs():
            stations = sorted({st for (e, st) in rows if e == a} & {st for (e, st) in rows if e == b})
            meas, drops = [], {}
            flags = {}
            for st in stations:
                ca, cb = coords.get((a, st)), coords.get((b, st))
                if ca is None or cb is None:
                    drops["no coordinates"] = drops.get("no coordinates", 0) + 1
                    continue
                g = geo.get((st, round(ca[0], 4), round(ca[1], 4)))
                if g is None or g["tt"]["df"] is None or g["tt"]["ab"] is None:
                    drops["no DF/AB prediction"] = drops.get("no DF/AB prediction", 0) + 1
                    continue
                why = gate(rows[(a, st)], rows[(b, st)], bl, g["tt"]["clear"])
                if why:
                    drops[why] = drops.get(why, 0) + 1
                    continue
                try:
                    m = measure(st, rows[(a, st)], rows[(b, st)], evs[a], evs[b], g["tt"], band,
                                bl, cache)
                except Exception as ex:  # noqa: BLE001
                    m = dict(station=st, drop="error " + repr(ex)[:80])
                m.update(pair=label, band=tag, dist=round(g["dist"], 3), baz=g["baz"],
                         patch=g["patch"])
                if g["tt"]["p_df"] is not None and g["tt"]["p_ab"] is not None:
                    m["pos_budget_s"] = round(abs(g["tt"]["p_ab"] - g["tt"]["p_df"])
                                              * POS_BOUND_DEG, 4)
                meas_rows.append(m)
                if m.get("drop"):
                    drops[m["drop"]] = drops.get(m["drop"], 0) + 1
                    continue
                meas.append(m)
                flags[st] = g["patch"]
            res = summarise_pair(meas, flags)
            res.update(a=a, b=b, a_date=evs[a]["date"], b_date=evs[b]["date"], drops=drops,
                       graded=(tag == "primary" and not label.startswith("R:")))
            results["%s|%s" % (tag, label)] = res
            print("PAIR", tag, label, res.get("n"), res.get("headline"), flush=True)
    write_csv(os.path.join(args.out_dir, "dprk37_measurements.csv"), meas_rows)
    s["pairs"] = results
    g1 = results.get("primary|G1", {})
    s["P37_headline"] = g1.get("headline", "untestable")
    s["G2_headline"] = results.get("primary|G2", {}).get("headline", "untestable")
    s["P37b_G1"] = g1.get("P37b", "untestable")
    s["minutes"] = round((time.time() - t0) / 60.0, 1)
    s["finished"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    with open(os.path.join(args.out_dir, "dprk37_summary.json"), "w", encoding="utf-8") as f:
        json.dump(s, f, indent=1, default=str)
    print("P37", s["P37_headline"], "G2", s["G2_headline"], "P37b", s["P37b_G1"], flush=True)
    return 0


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--data-dir", default="/data")
    ap.add_argument("--out-dir", default="/out")
    return run(ap.parse_args(argv))


if __name__ == "__main__":
    sys.exit(main())

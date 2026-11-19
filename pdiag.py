#!/usr/bin/env python3
"""pdiag.py: why do many global stations record no mantle P from the M03-24 events?

Part of inner-core-repeaters by MV2, Munim, Inc. (made with AI agents). Licence: Apache-2.0.
Plan: findings/2026-11-19-pdiag-plan.md (published before this script ran). Exploratory
diagnostic; it uses only the eight M03-24 events already used in mantlep2-v1 and does not
touch any event or station measurement of pre-registration 13.

In mantlep2-v1, 22 of the 37 candidate stations with data (30-85 deg) gated at most one of
the eight M03-24 South Sandwich Islands events, many with raw SNR near 3 (noise only). Example: II.SUR at 41 deg
(raw SNR 2.4-4.7) while IU.LCO at 43 deg had 35-52. Before any wider station set is used for
mantle-P source scores, this asks why. Explanations and their tests:

  R  radiation: those stations sit near a P nodal plane. Test: predicted P amplitude
     A = gamma.M.gamma / M0 from the GCMT moment tensor vs observed SNR.
  F  frequency: 1-2 Hz is the wrong band there. Test: SNR at 0.5-1 and 2-4 Hz.
  L  location code / sensor: the first location with data is a poor sensor. Test: SNR at
     every location code (00, 10, empty).
  X  processing: response removal fails or misbehaves. Test: SNR without response removal.
  D  depth phase: where direct P is weak, an envelope picker may lock on pP or sP (this
     could explain the 2.5-7 s "late P" seen at some stations). Test: onset time vs ak135
     pP-P and sP-P.

Writes to --out-dir: pdiag_rows.csv, pdiag_stations.csv, pdiag_events.csv,
pdiag_summary.json, plots/pdiag_<NET>_<STA>.png, plots/pdiag_overview.png.
"""
import argparse
import calendar
import csv
import json
import math
import os
import re
import sys
import time
import traceback

import numpy as np

VERSION = "0.1.0"

# Identical to mantlep2.EVENTS (v0.13.0): ComCat id, origin UTC, lat, lon, depth km, mb.
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
SRC_LAT, SRC_LON = -59.03, -25.68
NETS = "IU,II,GT"
DIST_RANGE = (30.0, 85.0)
LOCS = ("00", "10", "")
DATASELECT = "https://service.earthscope.org/fdsnws/dataselect/1/query"
STATION = "https://service.earthscope.org/fdsnws/station/1/query"

FS = 20.0
FETCH_PRE, FETCH_POST = 120.0, 120.0     # as mantlep2, so its cached files are reused
GRID_PRE, GRID_POST = 100.0, 100.0       # analysis grid around each event's ak135 P
BANDS = (("0.5-1", 0.5, 1.0), ("1-2", 1.0, 2.0), ("2-4", 2.0, 4.0))
MAIN_BAND = "1-2"
NOISE_W = (-90.0, -15.0)
SIG_W = (-5.0, 15.0)
RAW_SIG, RAW_NOISE = (-5.0, 15.0), (-38.0, -8.0)   # mantlep2's raw SNR windows, for comparison
PICK_W = (-15.0, 30.0)
PICK_K = 4.0                             # onset = first envelope sample >= PICK_K x noise RMS
GOOD_SNR = 5.0
CLEAR_SNR = 10.0
LATE_S = 2.0

GCMT_BASE = "https://www.ldeo.columbia.edu/~gcmt/projects/CMT/catalog/"
GCMT_FILES = ("jan76_dec20.ndk", "NEW_MONTHLY/2022/mar22.ndk", "NEW_MONTHLY/2024/mar24.ndk",
              "NEW_MONTHLY/2025/dec25.ndk", "NEW_QUICK/qcmt.ndk")
BOX = dict(lat=(-61.5, -54.0), lon=(-31.0, -22.0))
NAN = float("nan")


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
    xs = sorted(float(x) for x in xs if isf(x))
    n = len(xs)
    if n == 0:
        return NAN
    m = n // 2
    return xs[m] if n % 2 else 0.5 * (xs[m - 1] + xs[m])


def rankdata(v):
    """Ranks 0..n-1, ties averaged."""
    v = np.asarray(v, dtype=float)
    order = np.argsort(v, kind="mergesort")
    sv = v[order]
    r = np.empty(len(v), dtype=float)
    i = 0
    while i < len(v):
        j = i
        while j + 1 < len(v) and sv[j + 1] == sv[i]:
            j += 1
        r[order[i:j + 1]] = 0.5 * (i + j)
        i = j + 1
    return r


def pearson(a, b):
    a = np.asarray(a, dtype=float) - np.mean(a)
    b = np.asarray(b, dtype=float) - np.mean(b)
    d = math.sqrt(float(np.dot(a, a)) * float(np.dot(b, b)))
    return float(np.dot(a, b) / d) if d > 0 else NAN


def spearman(a, b):
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    m = np.isfinite(a) & np.isfinite(b)
    if m.sum() < 5:
        return NAN
    return pearson(rankdata(a[m]), rankdata(b[m]))


def envelope(x, fs, smooth_s=0.25):
    w = max(1, int(round(smooth_s * fs)))
    return np.sqrt(np.convolve(np.asarray(x, dtype=float) ** 2, np.ones(w) / w, mode="same"))


def measure(x, fs=FS, c=None):
    """SNR, pick and envelope-maximum time of one trace on the grid (ak135 P at sample c)."""
    x = np.asarray(x, dtype=float)
    c = int(round(GRID_PRE * fs)) if c is None else c

    def idx(t):
        return c + int(round(t * fs))

    n0, n1 = idx(NOISE_W[0]), idx(NOISE_W[1])
    s0, s1 = idx(SIG_W[0]), idx(SIG_W[1])
    p0, p1 = idx(PICK_W[0]), idx(PICK_W[1])
    r0, r1 = idx(RAW_NOISE[0]), idx(RAW_NOISE[1])
    if min(n0, s0, p0, r0) < 0 or max(n1, s1, p1, r1) > len(x):
        return None
    rn = float(np.sqrt(np.mean(x[n0:n1] ** 2)))
    if not rn > 0:
        return None
    peak = float(np.max(np.abs(x[s0:s1])))
    env = envelope(x, fs)
    seg = env[p0:p1]
    hit = np.nonzero(seg >= PICK_K * rn)[0]
    rr = float(np.sqrt(np.mean(x[r0:r1] ** 2)))
    q0, q1 = idx(RAW_SIG[0]), idx(RAW_SIG[1])
    raw = float(np.max(np.abs(x[q0:q1]))) / rr if rr > 0 else NAN
    return dict(noise_rms=rn, peak=peak, snr=peak / rn,
                t_pick=(p0 + int(hit[0]) - c) / fs if hit.size else NAN,
                t_envmax=(p0 + int(np.argmax(seg)) - c) / fs, raw_snr_m2win=raw)


# ----------------------------------------------------------------------------- moment tensors
def sdr_to_ned(strike, dip, rake, m0=1.0):
    """Aki & Richards Box 4.4: double couple -> M in x=north, y=east, z=down."""
    s, d, r = math.radians(strike), math.radians(dip), math.radians(rake)
    mxx = -m0 * (math.sin(d) * math.cos(r) * math.sin(2 * s) + math.sin(2 * d) * math.sin(r) * math.sin(s) ** 2)
    mxy = m0 * (math.sin(d) * math.cos(r) * math.cos(2 * s) + 0.5 * math.sin(2 * d) * math.sin(r) * math.sin(2 * s))
    mxz = -m0 * (math.cos(d) * math.cos(r) * math.cos(s) + math.cos(2 * d) * math.sin(r) * math.sin(s))
    myy = m0 * (math.sin(d) * math.cos(r) * math.sin(2 * s) - math.sin(2 * d) * math.sin(r) * math.cos(s) ** 2)
    myz = -m0 * (math.cos(d) * math.cos(r) * math.sin(s) - math.cos(2 * d) * math.sin(r) * math.cos(s))
    mzz = m0 * math.sin(2 * d) * math.sin(r)
    return np.array([[mxx, mxy, mxz], [mxy, myy, myz], [mxz, myz, mzz]])


def mt_ned(mrr, mtt, mpp, mrt, mrp, mtp):
    """GCMT (r=up, t=south, p=east) -> x=north, y=east, z=down."""
    return np.array([[mtt, -mtp, mrt], [-mtp, mpp, -mrp], [mrt, -mrp, mrr]], dtype=float)


def ned_to_rtp(M):
    return dict(mrr=M[2, 2], mtt=M[0, 0], mpp=M[1, 1], mrt=M[0, 2], mrp=-M[1, 2], mtp=-M[0, 1])


def scalar_moment(M):
    return float(math.sqrt(0.5 * float(np.sum(np.asarray(M) ** 2))))


def p_amp(M, takeoff_deg, az_deg):
    """Far-field P radiation gamma.M.gamma / M0 for a ray leaving at takeoff_deg from the
    downward vertical toward azimuth az_deg (clockwise from north). |A| <= 1 for a pure DC."""
    m0 = scalar_moment(M)
    if not m0 > 0:
        return NAN
    i, a = math.radians(takeoff_deg), math.radians(az_deg)
    g = np.array([math.sin(i) * math.cos(a), math.sin(i) * math.sin(a), math.cos(i)])
    return float(g @ np.asarray(M, dtype=float) @ g) / m0


_NUM = re.compile(r"[-+]?\d+\.\d*")


def _epoch(date, hms):
    y, mo, d = (int(v) for v in date.split("/"))
    h, mi, s = hms.split(":")
    return calendar.timegm((y, mo, d, int(h), int(mi), 0, 0, 0, 0)) + float(s)


def iso_epoch(s):
    s = s.strip().rstrip("Z")
    date, hms = s.split("T")
    return _epoch(date.replace("-", "/"), hms)


def parse_ndk(text, box=BOX, quick=False):
    """GCMT NDK (5 lines per event) -> moment tensors of events inside box."""
    lines = [ln for ln in text.splitlines() if ln.strip()]
    out = []
    for i in range(0, len(lines) - 4, 5):
        l1, l2, _, l4, _ = lines[i:i + 5]
        try:
            t = l1[4:].split()
            ep = _epoch(t[0], t[1])
            lat, lon, dep = float(t[2]), float(t[3]), float(t[4])
        except (IndexError, ValueError):
            continue
        if not (box["lat"][0] <= lat <= box["lat"][1] and box["lon"][0] <= lon <= box["lon"][1]):
            continue
        try:
            ex = int(l4[:2])
            v = [float(s) for s in _NUM.findall(l4[2:])]
        except ValueError:
            continue
        if len(v) < 12:
            continue
        k = 10.0 ** ex
        out.append(dict(name=l2[:16].strip(), epoch=ep, lat=lat, lon=lon, depth=dep, quick=quick,
                        mrr=v[0] * k, mtt=v[2] * k, mpp=v[4] * k, mrt=v[6] * k, mrp=v[8] * k, mtp=v[10] * k))
    return out


def match_gcmt(epoch, lat, lon, cmts, tol_s=20.0, tol_deg=1.0):
    """Closest GCMT entry in time within tolerances; monthly/final entries before quick ones."""
    best = None
    for c in cmts:
        dt = abs(c["epoch"] - epoch)
        if dt <= tol_s and abs(c["lat"] - lat) <= tol_deg and abs(c["lon"] - lon) <= tol_deg:
            key = (bool(c.get("quick")), dt)
            if best is None or key < best[0]:
                best = (key, c)
    return best[1] if best else None


def mw_from_dyncm(m0):
    return 2.0 / 3.0 * (math.log10(m0) - 16.1) if m0 > 0 else NAN


# ----------------------------------------------------------------------------- analysis
def _main_rows(rows):
    return [r for r in rows if r.get("band") == MAIN_BAND and r.get("resp_requested") is True
            and r.get("first_loc") is True and isf(r.get("snr"))]


def station_table(rows):
    by = {}
    for r in rows:
        if isf(r.get("snr")):
            by.setdefault(r["station"], []).append(r)
    out = []
    for name, v in sorted(by.items()):
        def med(filt, key="snr"):
            return median([r[key] for r in v if filt(r)])
        main = [r for r in v if r["band"] == MAIN_BAND and r["resp_requested"] is True and r["first_loc"] is True]
        d = dict(station=name, dist_deg=v[0].get("dist_deg"), az_deg=v[0].get("az_deg"),
                 n_events=len(main), n_good=sum(1 for r in main if r["snr"] >= GOOD_SNR),
                 median_absA=median([abs(r["A_pred"]) for r in main if isf(r.get("A_pred"))]))
        for b, _, _ in BANDS:
            d["snr_" + b] = med(lambda r, b=b: r["band"] == b and r["resp_requested"] is True and r["first_loc"] is True)
        d["snr_noresp"] = med(lambda r: r["band"] == MAIN_BAND and r["resp_requested"] is False and r["first_loc"] is True)
        for loc in LOCS:
            lc = loc or "--"
            d["snr_loc_" + lc] = med(lambda r, lc=lc: r["band"] == MAIN_BAND and r["resp_requested"] is True
                                     and r["loc"] == lc)
        clear = [r for r in main if r["snr"] >= CLEAR_SNR and isf(r.get("t_pick"))]
        d["median_t_pick_clear"] = median([r["t_pick"] for r in clear])
        d["rho_tpick_pP"] = spearman([r["t_pick"] for r in clear], [fnum(r.get("pP_minus_P_s")) for r in clear])
        d["median_noise_ms"] = median([r["noise_rms"] for r in main if r.get("resp_removed") is True])
        d["median_peak_ms"] = median([r["peak"] for r in main if r.get("resp_removed") is True])
        out.append(d)
    return out


def analyse(rows):
    out = {}
    main = _main_rows(rows)
    ra = [r for r in main if isf(r.get("A_pred")) and r["snr"] > 0]
    good = [abs(r["A_pred"]) for r in ra if r["snr"] >= GOOD_SNR]
    bad = [abs(r["A_pred"]) for r in ra if r["snr"] < GOOD_SNR]
    st = {}
    for r in ra:
        st.setdefault(r["station"], []).append(r)
    sa = [(median([abs(r["A_pred"]) for r in v]), median([r["snr"] for r in v])) for v in st.values()]
    out["R_radiation"] = dict(
        n_rows=len(ra), rho_rows_absA_logsnr=spearman([abs(r["A_pred"]) for r in ra],
                                                      [math.log10(r["snr"]) for r in ra]),
        n_stations=len(sa), rho_stations=spearman([a for a, _ in sa], [s for _, s in sa]),
        n_good=len(good), n_bad=len(bad), median_absA_good=median(good), median_absA_bad=median(bad))
    tab = station_table(rows)
    fails = [d for d in tab if isf(d["snr_" + MAIN_BAND]) and d["snr_" + MAIN_BAND] < GOOD_SNR]
    out["n_stations"] = len(tab)
    out["failing_stations"] = [d["station"] for d in fails]
    other = [b for b, _, _ in BANDS if b != MAIN_BAND]
    out["F_band_rescues"] = [d["station"] for d in fails if any(isf(d["snr_" + b]) and d["snr_" + b] >= GOOD_SNR
                                                                for b in other)]
    out["L_location_rescues"] = [d["station"] for d in fails
                                 if any(isf(d["snr_loc_" + (l or "--")]) and d["snr_loc_" + (l or "--")] >= GOOD_SNR
                                        for l in LOCS)]
    out["X_response_changes"] = [d["station"] for d in tab
                                 if isf(d["snr_noresp"]) and isf(d["snr_" + MAIN_BAND]) and d["snr_" + MAIN_BAND] > 0
                                 and abs(math.log10(d["snr_noresp"] / d["snr_" + MAIN_BAND])) > math.log10(1.5)]
    clear = [r for r in main if r["snr"] >= CLEAR_SNR and isf(r.get("t_pick"))]
    late = [r for r in clear if r["t_pick"] > LATE_S]
    out["D_depth_phase"] = dict(
        n_clear=len(clear), median_t_pick=median([r["t_pick"] for r in clear]),
        n_late=len(late), frac_late=(len(late) / float(len(clear))) if clear else NAN,
        late_median_abs_minus_pP=median([abs(r["t_pick"] - fnum(r.get("pP_minus_P_s"))) for r in late]),
        late_median_abs_minus_sP=median([abs(r["t_pick"] - fnum(r.get("sP_minus_P_s"))) for r in late]),
        late_median_abs_minus_P=median([abs(r["t_pick"]) for r in late]),
        rho_tpick_pP_clear=spearman([r["t_pick"] for r in clear], [fnum(r.get("pP_minus_P_s")) for r in clear]))
    return out, tab


# ----------------------------------------------------------------------------- data access
def get(session, url, params, tries=3, timeout=90, verify=True):
    for k in range(tries):
        try:
            r = session.get(url, params=params, timeout=timeout, verify=verify)
            if r.status_code == 200 and r.content:
                return r.content
            if r.status_code in (204, 404):
                return None
        except Exception as ex:  # noqa: BLE001
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
    out = {}
    for line in (body or b"").decode("utf-8", "replace").splitlines():
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        f = [s.strip() for s in line.split("|")]
        try:
            out.setdefault((f[0], f[1]), (float(f[2]), float(f[3])))
        except (IndexError, ValueError):
            continue
    return [(k[0], k[1], v[0], v[1]) for k, v in sorted(out.items())]


def fetch_gcmt(session, cache):
    os.makedirs(cache, exist_ok=True)
    cmts, status = [], {}
    for f in GCMT_FILES:
        path = os.path.join(cache, "gcmt_" + f.replace("/", "_"))
        body, how = None, ""
        if os.path.exists(path) and os.path.getsize(path) > 0:
            with open(path, "rb") as fh:
                body = fh.read()
            how = "cache"
        else:
            for base, verify in ((GCMT_BASE, True), (GCMT_BASE.replace("https://", "http://"), True),
                                 (GCMT_BASE, False)):
                body = get(session, base + f, {}, tries=2, timeout=300, verify=verify)
                if body:
                    how = ("https" if base.startswith("https") else "http") + ("" if verify else "-unverified-tls")
                    with open(path, "wb") as fh:
                        fh.write(body)
                    break
        if not body:
            status[f] = "unavailable"
            continue
        got = parse_ndk(body.decode("latin-1"), quick=f.startswith("NEW_QUICK"))
        status[f] = "%s, %d bytes, %d in box" % (how, len(body), len(got))
        cmts += got
    return cmts, status


def load_loc(session, net, sta, loc, t1, t2, cache):
    """(trace, inventory or None, status) for one location code. mantlep2's cache names."""
    from obspy import read, read_inventory
    lc = loc or "--"
    key = "%s.%s.%s.BHZ.%s" % (net, sta, lc, t1.strftime("%Y%m%dT%H%M%S"))
    mpath = os.path.join(cache, key + ".mseed")
    xpath = os.path.join(cache, key + ".xml")
    if not os.path.exists(mpath):
        body = get(session, DATASELECT, dict(net=net, sta=sta, loc=lc, cha="BHZ", start=iso(t1), end=iso(t2)))
        if not body:
            return None, None, "no_data"
        with open(mpath, "wb") as f:
            f.write(body)
    try:
        st = read(mpath)
    except Exception:  # noqa: BLE001
        return None, None, "unreadable"
    st.merge()
    if len(st) != 1 or np.ma.isMaskedArray(st[0].data):
        return None, None, "gaps"
    if not os.path.exists(xpath):
        body = get(session, STATION, dict(net=net, sta=sta, loc=lc, cha="BHZ", starttime=iso(t1),
                                           endtime=iso(t2), level="response"))
        if body:
            with open(xpath, "wb") as f:
                f.write(body)
    inv = None
    if os.path.exists(xpath):
        try:
            inv = read_inventory(xpath)
        except Exception:  # noqa: BLE001
            inv = None
    return st[0], inv, "ok"


def process(tr, inv, t_start, band, remove=True):
    tr = tr.copy()
    tr.detrend("demean")
    tr.detrend("linear")
    tr.taper(0.05)
    resp = False
    if remove and inv is not None:
        try:
            tr.remove_response(inventory=inv, output="VEL", pre_filt=(0.2, 0.4, 7.0, 9.0))
            resp = True
        except Exception:  # noqa: BLE001
            resp = False
    tr.filter("bandpass", freqmin=band[0], freqmax=band[1], corners=4, zerophase=True)
    tr.interpolate(sampling_rate=FS, method="lanczos", a=12, starttime=t_start,
                   npts=int(round((GRID_PRE + GRID_POST) * FS)))
    return np.asarray(tr.data, dtype=float), resp


def p_time(model, depth, dist):
    """Exactly mantlep2.p_time (so the cache keys match)."""
    arr = model.get_travel_times(source_depth_in_km=max(0.0, depth), distance_in_degree=dist,
                                 phase_list=["P", "p", "Pdiff"])
    return min(a.time for a in arr) if arr else None


def phase_info(model, depth, dist):
    out = dict(takeoff_deg=NAN, pP_minus_P_s=NAN, sP_minus_P_s=NAN)
    arr = model.get_travel_times(source_depth_in_km=max(0.0, depth), distance_in_degree=dist,
                                 phase_list=["P", "p", "Pdiff", "pP", "sP"])
    first = {}
    for a in arr:
        first.setdefault(a.name, a)
    direct = [first[n] for n in ("P", "p", "Pdiff") if n in first]
    if not direct:
        return out
    p = min(direct, key=lambda a: a.time)
    out["takeoff_deg"] = float(p.takeoff_angle)
    for n in ("pP", "sP"):
        if n in first:
            out[n + "_minus_P_s"] = float(first[n].time - p.time)
    return out


def run_station(session, model, net, sta, slat, slon, cache, cmt_by_event, deadline):
    from obspy import UTCDateTime
    from obspy.geodetics import gps2dist_azimuth, locations2degrees
    name = "%s.%s" % (net, sta)
    rows, traces, first_loc = [], {}, {}
    for eid, ot, lat, lon, dep, mb in EVENTS:
        if time.time() > deadline:
            break
        o = UTCDateTime(ot)
        dist = locations2degrees(lat, lon, slat, slon)
        az = gps2dist_azimuth(lat, lon, slat, slon)[1]
        tp = p_time(model, dep, dist)
        if tp is None:
            continue
        info = phase_info(model, dep, dist)
        M = cmt_by_event.get(eid)
        A = p_amp(M, info["takeoff_deg"], az) if (M is not None and isf(info["takeoff_deg"])) else NAN
        for loc in LOCS:
            base = dict(station=name, event=eid, date=ot[:10], mb=mb, depth_km=dep, dist_deg=round(dist, 2),
                        az_deg=round(az, 1), A_pred=round(A, 4) if isf(A) else "", loc=loc or "--",
                        first_loc=False, band="", resp_requested="", resp_removed="", status="",
                        **{k: round(v, 2) if isf(v) else "" for k, v in info.items()})
            try:
                tr, inv, status = load_loc(session, net, sta, loc, o + tp - FETCH_PRE, o + tp + FETCH_POST, cache)
            except Exception as ex:  # noqa: BLE001
                tr, inv, status = None, None, "error:" + repr(ex)[:60]
            if tr is None:
                rows.append(dict(base, status=status))
                continue
            if eid not in first_loc:
                first_loc[eid] = loc
            base["first_loc"] = first_loc[eid] == loc
            base["sps"] = float(tr.stats.sampling_rate)
            for bname, f1, f2 in BANDS:
                for remove in ((True, False) if bname == MAIN_BAND else (True,)):
                    row = dict(base, band=bname, resp_requested=remove)
                    try:
                        x, resp = process(tr, inv, o + tp - GRID_PRE, (f1, f2), remove)
                        m = measure(x)
                    except Exception as ex:  # noqa: BLE001
                        rows.append(dict(row, status="process_error:" + repr(ex)[:60]))
                        continue
                    if m is None or not np.all(np.isfinite(x)):
                        rows.append(dict(row, status="bad_grid"))
                        continue
                    row.update(status="ok", resp_removed=resp, **m)
                    rows.append(row)
                    if bname == MAIN_BAND and remove and base["first_loc"]:
                        traces[eid] = x / m["noise_rms"]
    return rows, traces


# ----------------------------------------------------------------------------- output
def round_row(r):
    out = {}
    for k, v in r.items():
        if isinstance(v, float):
            out[k] = "" if not math.isfinite(v) else (float("%.4g" % v) if abs(v) < 1e-3 else round(v, 4))
        else:
            out[k] = v
    return out


def write_csv(path, rows):
    keys = []
    for r in rows:
        for k in r:
            if k not in keys:
                keys.append(k)
    with open(path, "w", newline="") as f:
        if not keys:
            return
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        for r in rows:
            w.writerow(round_row(r))


def jsonable(o):
    if isinstance(o, dict):
        return {str(k): jsonable(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [jsonable(v) for v in o]
    if isinstance(o, np.bool_):
        return bool(o)
    if isinstance(o, np.integer):
        return int(o)
    if isinstance(o, (np.floating, float)):
        v = float(o)
        return v if math.isfinite(v) else None
    return o


def plot_station(out_dir, name, traces, rows):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    if not traces:
        return
    c = int(round(GRID_PRE * FS))
    t = (np.arange(len(next(iter(traces.values())))) - c) / FS
    keep = (t >= -20) & (t <= 40)
    fig, ax = plt.subplots(figsize=(9, 7))
    info = {}
    for r in rows:
        info.setdefault(r["event"], r)
    for k, (eid, x) in enumerate(sorted(traces.items(), key=lambda kv: kv[0])):
        y = x[keep] / 12.0 - 2.0 * k
        ax.plot(t[keep], y, lw=0.6, color="k")
        r = info.get(eid, {})
        lab = "%s %s A=%s" % (r.get("date", ""), eid, r.get("A_pred", ""))
        ax.text(-19.5, -2.0 * k + 0.5, lab, fontsize=7)
        for key, col in (("pP_minus_P_s", "tab:red"), ("sP_minus_P_s", "tab:blue")):
            v = r.get(key)
            if isf(v):
                ax.plot([v, v], [-2.0 * k - 0.8, -2.0 * k + 0.8], color=col, lw=1)
    ax.axvline(0, color="tab:green", lw=0.8)
    r0 = rows[0] if rows else {}
    ax.set_title("%s  %.1f deg  az %.0f  (1-2 Hz, /noise RMS /12; green P, red pP, blue sP)"
                 % (name, float(r0.get("dist_deg", NAN)), float(r0.get("az_deg", NAN))), fontsize=9)
    ax.set_xlabel("s after ak135 P")
    ax.set_yticks([])
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, "plots", "pdiag_%s.png" % name.replace(".", "_")), dpi=80)
    plt.close(fig)


def plot_overview(out_dir, tab, rows):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 2, figsize=(13, 4.8))
    az = [d["az_deg"] for d in tab]
    sn = [d["snr_" + MAIN_BAND] for d in tab]
    aa = [d["median_absA"] if isf(d["median_absA"]) else -0.1 for d in tab]
    sc = ax[0].scatter(az, sn, c=aa, cmap="viridis", vmin=-0.1, vmax=1.0, s=30)
    for d in tab:
        ax[0].annotate(d["station"].split(".")[1], (d["az_deg"], d["snr_" + MAIN_BAND]), fontsize=6)
    ax[0].set_yscale("log")
    ax[0].axhline(GOOD_SNR, color="k", lw=0.7)
    ax[0].set_xlabel("azimuth from source (deg)")
    ax[0].set_ylabel("median 1-2 Hz P SNR")
    fig.colorbar(sc, ax=ax[0], label="median |predicted P amplitude| (-0.1 = no GCMT)")
    main = [r for r in _main_rows(rows) if isf(r.get("A_pred"))]
    ax[1].scatter([abs(r["A_pred"]) for r in main], [r["snr"] for r in main], s=8, alpha=0.6)
    ax[1].set_yscale("log")
    ax[1].axhline(GOOD_SNR, color="k", lw=0.7)
    ax[1].set_xlabel("|predicted P amplitude| (GCMT, gamma.M.gamma/M0)")
    ax[1].set_ylabel("1-2 Hz P SNR (event-station)")
    fig.suptitle("pdiag %s: mantle P of the eight M03-24 events at 30-85 deg" % VERSION)
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, "plots", "pdiag_overview.png"), dpi=90)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", default="/data")
    ap.add_argument("--out-dir", default="/out")
    ap.add_argument("--time-budget-min", type=float, default=100.0)
    args = ap.parse_args()
    t0 = time.time()
    deadline = t0 + 60.0 * args.time_budget_min
    import requests
    from obspy.taup import TauPyModel
    cache = os.path.join(args.data_dir, "mantlep")
    os.makedirs(cache, exist_ok=True)
    os.makedirs(os.path.join(args.out_dir, "plots"), exist_ok=True)
    session = requests.Session()
    session.headers["User-Agent"] = "inner-core-repeaters (MV2, Munim, Inc.)"
    s = dict(version=VERSION, started=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
             bands=BANDS, main_band=MAIN_BAND, windows=dict(noise=NOISE_W, sig=SIG_W, pick=PICK_W, pick_k=PICK_K),
             filter="4-pole Butterworth, zero phase (mantlep2 used causal)")
    cmts, s["gcmt_files"] = fetch_gcmt(session, os.path.join(args.data_dir, "gcmt"))
    cmt_by_event, ev_rows = {}, []
    for eid, ot, lat, lon, dep, mb in EVENTS:
        c = match_gcmt(iso_epoch(ot), lat, lon, cmts)
        row = dict(event=eid, origin=ot, mb=mb, depth_km=dep, gcmt="", quick="", mw="", gcmt_depth="")
        if c is not None:
            M = mt_ned(c["mrr"], c["mtt"], c["mpp"], c["mrt"], c["mrp"], c["mtp"])
            cmt_by_event[eid] = M
            row.update(gcmt=c["name"], quick=c["quick"], mw=round(mw_from_dyncm(scalar_moment(M)), 2),
                       gcmt_depth=c["depth"], **{k: "%.4g" % c[k] for k in ("mrr", "mtt", "mpp", "mrt", "mrp", "mtp")})
        ev_rows.append(row)
    s["n_events_with_gcmt"] = len(cmt_by_event)
    print("gcmt", s["gcmt_files"], "matched", len(cmt_by_event), flush=True)
    stations = discover(session)
    cand = []
    for net, sta, lat, lon in stations:
        from obspy.geodetics import locations2degrees
        d = locations2degrees(SRC_LAT, SRC_LON, lat, lon)
        if DIST_RANGE[0] <= d <= DIST_RANGE[1]:
            cand.append((d, net, sta, lat, lon))
    cand.sort()
    s["n_candidates"] = len(cand)
    model = TauPyModel("ak135")
    rows, skipped = [], []
    for d, net, sta, lat, lon in cand:
        name = "%s.%s" % (net, sta)
        if time.time() > deadline:
            skipped.append(name)
            continue
        try:
            r, traces = run_station(session, model, net, sta, lat, lon, cache, cmt_by_event, deadline)
        except Exception:  # noqa: BLE001
            traceback.print_exc()
            continue
        rows += r
        try:
            plot_station(args.out_dir, name, traces, [x for x in r if x.get("first_loc")])
        except Exception as ex:  # noqa: BLE001
            print("plot failed", name, repr(ex)[:120], flush=True)
        print(name, "%.1f deg" % d, "rows", len(r), "%.1f min" % ((time.time() - t0) / 60.0), flush=True)
    s["skipped_time_budget"] = skipped
    res, tab = analyse(rows)
    s["analysis"] = res
    s["finished"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    s["minutes"] = round((time.time() - t0) / 60.0, 1)
    write_csv(os.path.join(args.out_dir, "pdiag_rows.csv"), rows)
    write_csv(os.path.join(args.out_dir, "pdiag_stations.csv"), tab)
    write_csv(os.path.join(args.out_dir, "pdiag_events.csv"), ev_rows)
    with open(os.path.join(args.out_dir, "pdiag_summary.json"), "w") as f:
        json.dump(jsonable(s), f, indent=1)
    try:
        plot_overview(args.out_dir, tab, rows)
    except Exception as ex:  # noqa: BLE001
        print("overview plot failed", repr(ex)[:200], flush=True)
    print(json.dumps(jsonable(res))[:4000], flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())

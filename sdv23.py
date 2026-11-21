#!/usr/bin/env python3
# Copyright 2026 MV2, Munim, Inc.
# Licensed under the Apache License, Version 2.0 (Apache-2.0). See LICENSE.
# Made with AI agents working in MV2's name.
"""sdv23.py - pre-registration 23: is the reported decadal ULVZ change at IU.SDV
robust to the station's 2009-09-18 sensor re-orientation, and does the source
cluster hold a third member that dates the change?

Plan: findings/2026-11-21-prereg23-sdv-ulvz-plan.md (written and published before
this code downloads any waveform). Amendment 1 (0.1.1, also before any waveform):
findings/2026-11-21-prereg23-amendment-1.md - body-wave lags are measured after a
coarse alignment of the two events on the transverse S window (+/- 5 s), because a
catalogue origin-time difference larger than the 0.6 s fine search would otherwise
push every lag to the edge. Decision rules are unchanged.

Claim under test: Zhang & Wen 2025, Nat. Commun. 16:8754 (doi 10.1038/s41467-025-63814-3):
for the Aleutian doublet 2000-07-07 / 2009-12-17 recorded at IU.SDV (Venezuela,
94.9 deg), SKS on the radial component arrives 160 ms earlier in 2009 and radial
S-Scd-ScS shifts 80 -> 30 -> 0 ms, while the transverse component is unchanged.
IU.SDV station metadata change the horizontal azimuths of location 00 on 2009-09-18
(BHN 0 -> 357, BHE 90 -> 87; location 10: 0 -> 343, 90 -> 73), between the two events.
The instrument responses at 0.25-4 Hz are identical at both event times.

Parts (all decision rules in PARAMS and grade_* functions):
 A. relative sensor rotation between the doublet events from long-period waveforms
    (delta_wf) and from Rayleigh-wave polarisation (Rayleigh-Z method);
 B. SKS (radial) minus S (transverse) differential time change, 0.4-2 Hz, under three
    rotations: station metadata ("meta"), nominal N/E ("naive"), physically
    consistent with the measured delta ("cons"); location-corrected with a relative
    relocation from mantle P at IU/II stations;
 C. membership (mantle P) and the same measurement for every other catalogue event
    within 12 km (history);
 S2. sensor azimuth time series from large global earthquakes before and after
    2009-09-18 at locations 00 and 10 (was the metadata step a physical rotation?).

Top level imports are stdlib + numpy only; ObsPy is imported inside functions.
"""

import argparse
import csv
import hashlib
import json
import math
import os
import sys
import time
import urllib.request

import numpy as np

VERSION = "0.1.1"
REPO_RAW = "https://raw.githubusercontent.com/munimtechnologies/inner-core-repeaters/main/"
PLAN_PATH = "findings/2026-11-21-prereg23-sdv-ulvz-plan.md"
AMEND_PATH = "findings/2026-11-21-prereg23-amendment-1.md"
BASE = "https://service.earthscope.org"
USGS = "https://earthquake.usgs.gov"
R_EARTH = 6371.0
SDV = dict(net="IU", sta="SDV", lat=8.8839, lon=-70.634)
STEP_TIME = "2009-09-18T00:00:00"
META_STEP = {"00": -3.0, "10": -17.0}   # metadata azimuth change of the N channel

# USGS ComCat, M >= 4.5 within 12 km of (51.405 N, 179.97 E), from 1994-08-19 (IU.SDV
# start) to 2026-11-21, queried 2026-11-21. Columns: id, time, lat, lon, depth, mag, type.
EVENTS = [
    ("E1996a", "1996-01-20T07:16:37.84", 51.418, -179.990, 33.0, 5.0, "mb"),
    ("E1998a", "1998-11-08T23:48:47.86", 51.334, -179.925, 33.0, 5.1, "mwc"),
    ("E1998b", "1998-11-09T03:00:41.10", 51.417, -179.884, 33.0, 4.5, "mb"),
    ("E2000a", "2000-07-07T15:46:44.56", 51.411, 179.978, 31.0, 5.9, "mwc"),
    ("E2000b", "2000-07-11T02:10:33.69", 51.482, -179.974, 61.7, 4.6, "mb"),
    ("E2000c", "2000-12-25T20:43:11.81", 51.301, 179.967, 55.1, 5.0, "mb"),
    ("E2007a", "2007-08-02T03:21:42.82", 51.307, -179.971, 21.0, 6.7, "mwc"),
    ("E2009a", "2009-12-17T20:01:21.61", 51.398, 179.963, 35.0, 5.8, "mwc"),
]
MASTER = "E2000a"
PARTNER = "E2009a"

PARAMS = dict(
    sr_up=200.0,              # Hz, common rate after FFT upsampling (body waves)
    band_body=(0.4, 2.0),     # Hz, as Zhang & Wen for S phases
    band_p=(0.8, 2.0),        # Hz, mantle P membership (Zhang & Wen P band)
    band_lp=[(0.01, 0.05), (0.02, 0.06)],  # Hz, delta_wf bands
    band_rz=(0.02, 0.05),     # Hz, Rayleigh-Z orientation
    pre_filt_body=(0.1, 0.2, 7.0, 9.0),
    pre_filt_lp=(0.003, 0.005, 0.5, 1.0),
    maxlag_body_s=0.6,
    coarse_maxlag_s=5.0,      # amendment 1: coarse S_T alignment before the fine lags
    maxlag_p_s=3.0,
    maxlag_lp_s=10.0,
    peak_search_s=5.0,
    p_member_cc=0.90, p_member_n=3, p_strict_cc=0.95, p_strict_n=4,
    p_dist=(30.0, 90.0), p_max_stations=30, reloc_min=4,
    st_cc_min=0.90, sks_cc_min=0.80, sks_snr_min=5.0,
    dwf_cc_min=0.90, dwf_band_spread_max=1.0, dwf_rz_agree_max=2.0,
    rz_c_min=0.80, rz_half=45.0, pair_half=30.0, grid_step=0.1,
    reproduce_s=0.10, explain_s=0.04,
    hist_stable_s=0.03, hist_change_s=0.06,
    s2_min_before=5, s2_min_after=3, s2_nboot=2000, s2_seed=23,
    s2_dist=(25.0, 150.0), s2_maxdepth=150.0,
)


def log(*a):
    print(time.strftime("%H:%M:%S"), "[sdv23]", *a, flush=True)


# ----------------------------------------------------------------- pure numerics

def wrap180(x):
    """Wrap angle(s) in degrees to [-180, 180)."""
    return (np.asarray(x, float) + 180.0) % 360.0 - 180.0


def upsample(x, f):
    """Band-limited (FFT) upsampling of a real series by integer factor f."""
    x = np.asarray(x, float)
    n = len(x)
    f = int(f)
    if f <= 1 or n < 2:
        return x.copy()
    X = np.fft.rfft(x)
    m = n * f
    Y = np.zeros(m // 2 + 1, complex)
    Y[:len(X)] = X
    if n % 2 == 0:
        Y[n // 2] *= 0.5
    return np.fft.irfft(Y, m) * f


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


def shift_left(b, s):
    """c[i] = b[i + s] (zeros where i + s is outside b). s may be negative."""
    b = np.asarray(b, float)
    s = int(s)
    c = np.zeros_like(b)
    n = len(b)
    if s >= 0:
        if s < n:
            c[:n - s] = b[s:]
    else:
        if -s < n:
            c[-s:] = b[:n + s]
    return c


def rotate_to_ne(c1, c2, az1, az2):
    """Horizontal channels c1, c2 pointing at azimuths az1, az2 (deg) -> (north, east)."""
    a1, a2 = math.radians(az1), math.radians(az2)
    det = math.sin(a2 - a1)
    if abs(det) < 1e-6:
        raise ValueError("channels are parallel")
    c1 = np.asarray(c1, float)
    c2 = np.asarray(c2, float)
    n = (math.sin(a2) * c1 - math.sin(a1) * c2) / det
    e = (-math.cos(a2) * c1 + math.cos(a1) * c2) / det
    return n, e


def ne_to_rt(n, e, baz):
    """ObsPy convention: R positive away from the source, baz = station -> source azimuth."""
    b = math.radians(baz)
    r = -np.asarray(e) * math.sin(b) - np.asarray(n) * math.cos(b)
    t = -np.asarray(e) * math.cos(b) + np.asarray(n) * math.sin(b)
    return r, t


def to_rt(c1, c2, az1, az2, baz):
    n, e = rotate_to_ne(c1, c2, az1, az2)
    return ne_to_rt(n, e, baz)


def window_lag(a, b, i1, i2, maxlag):
    """Normalised cross-correlation of a[i1:i2] with b[i1+k:i2+k], |k| <= maxlag samples.
    Returns (lag, cc, edge). lag > 0 means b is delayed relative to a. Parabolic refine."""
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


def snr(x, n1, n2, s1, s2):
    x = np.asarray(x, float)
    noise = x[max(0, int(n1)):max(0, int(n2))]
    sig = x[max(0, int(s1)):max(0, int(s2))]
    if len(noise) < 3 or len(sig) < 1:
        return float("nan")
    rms = math.sqrt(float(np.mean(noise ** 2)))
    return float(np.max(np.abs(sig)) / rms) if rms > 0 else float("inf")


def orient_rz(c1, c2, z, daz, baz, center, half=45.0, step=0.1):
    """Rayleigh-Z sensor orientation: azimuth a of channel 1 (channel 2 at a + daz) that
    maximises |corr(Hilbert(R), Z)|, searched within center +/- half. Returns (a, |c|)."""
    h1, h2 = hilbert(c1), hilbert(c2)
    z = np.asarray(z, float)
    nz = math.sqrt(float(np.dot(z, z)))
    best = (float("nan"), -1.0)
    for a in np.arange(center - half, center + half + step / 2, step):
        rh, _ = to_rt(h1, h2, a, a + daz, baz)
        nr = math.sqrt(float(np.dot(rh, rh)))
        if nr <= 0 or nz <= 0:
            continue
        c = abs(float(np.dot(rh, z)) / (nr * nz))
        if c > best[1]:
            best = (float(a % 360.0), c)
    return best


def orient_pair(a1, a2, b1, b2, daz_a, daz_b, half=30.0, step=0.1):
    """Relative physical rotation delta (deg) of the sensor between two recordings of the
    same wavefield: b's channel 1 points at delta in a's frame. Returns (delta, cc)."""
    an, ae = rotate_to_ne(a1, a2, 0.0, daz_a)
    na = float(np.dot(an, an) + np.dot(ae, ae))
    best = (float("nan"), -2.0)
    for d in np.arange(-half, half + step / 2, step):
        bn, be = rotate_to_ne(b1, b2, d, d + daz_b)
        nb = float(np.dot(bn, bn) + np.dot(be, be))
        if na <= 0 or nb <= 0:
            continue
        c = float(np.dot(an, bn) + np.dot(ae, be)) / math.sqrt(na * nb)
        if c > best[1]:
            best = (float(d), c)
    return best


def relocate(lags, p_skm, az_deg):
    """Relative relocation from P lags (b - a, s): lag_i = dt0 - p_i (dxN cos az_i +
    dxE sin az_i). Returns dict(dt0, dxN, dxE, rms, n) or None if fewer than 4."""
    lags = np.asarray(lags, float)
    p = np.asarray(p_skm, float)
    az = np.radians(np.asarray(az_deg, float))
    if len(lags) < 4:
        return None
    G = np.column_stack([np.ones(len(lags)), -p * np.cos(az), -p * np.sin(az)])
    m, *_ = np.linalg.lstsq(G, lags, rcond=None)
    res = lags - G @ m
    return dict(dt0=float(m[0]), dxN=float(m[1]), dxE=float(m[2]),
                rms=float(math.sqrt(float(np.mean(res ** 2)))), n=int(len(lags)))


def loc_pred(dxN, dxE, az_deg, dp_skm):
    """Predicted change of a differential time (phase1 - phase2) for a source shift
    (dxN, dxE) km, station azimuth az (from source), dp = p1 - p2 (s/km)."""
    a = math.radians(az_deg)
    return -dp_skm * (dxN * math.cos(a) + dxE * math.sin(a))


def grade_primary(gates, x_meta, x_naive, x_cons, P=PARAMS):
    """Pre-registered primary grade for E2000a/E2009a (x = SKS_R - S_T change, s,
    location-corrected; negative = SKS earlier in 2009)."""
    if not gates or not all(bool(v) for v in gates.values()):
        return "untestable"
    if not all(np.isfinite([x_meta, x_naive, x_cons])):
        return "untestable"
    thr = P["reproduce_s"]
    if not (x_meta <= -thr or x_naive <= -thr):
        return "not-reproduced"
    if abs(x_cons) <= P["explain_s"]:
        return "orientation-explains"
    if x_cons <= -thr:
        return "change-survives"
    return "partly-orientation"


def grade_history(dt, ok, P=PARAMS):
    if not ok or dt is None or not np.isfinite(dt):
        return "untestable"
    if abs(dt) <= P["hist_stable_s"]:
        return "stable"
    if abs(dt) >= P["hist_change_s"]:
        return "changed"
    return "unclear"


def grade_rotation(before, after, step, P=PARAMS):
    """S2: was the metadata azimuth step a physical rotation? Bootstrap the difference of
    medians (after - before) of wrapped channel-1 azimuths."""
    b = wrap180(np.asarray(before, float))
    a = wrap180(np.asarray(after, float))
    if len(b) < P["s2_min_before"] or len(a) < P["s2_min_after"]:
        return dict(grade="untestable", diff=float("nan"), lo=float("nan"), hi=float("nan"),
                    n_before=int(len(b)), n_after=int(len(a)))
    rng = np.random.default_rng(P["s2_seed"])
    diffs = np.empty(P["s2_nboot"])
    for k in range(P["s2_nboot"]):
        diffs[k] = (np.median(rng.choice(a, len(a))) - np.median(rng.choice(b, len(b))))
    diff = float(np.median(a) - np.median(b))
    lo, hi = (float(v) for v in np.percentile(diffs, [5, 95]))
    has_step = lo <= step <= hi
    has_zero = lo <= 0.0 <= hi
    if has_step and not has_zero:
        g = "physical-rotation"
    elif has_zero and not has_step:
        g = "no-physical-rotation"
    else:
        g = "unclear"
    return dict(grade=g, diff=diff, lo=lo, hi=hi, n_before=int(len(b)), n_after=int(len(a)))


def jsonable(o):
    """json.dump default: numpy scalars/arrays -> Python types, anything else -> str."""
    if isinstance(o, np.bool_):
        return bool(o)
    if isinstance(o, np.integer):
        return int(o)
    if isinstance(o, np.floating):
        return float(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    return str(o)


def sdv_location(t_iso):
    return "--" if t_iso < "1999-11-06" else "00"


# ----------------------------------------------------------------- data access

def sha256(data):
    return hashlib.sha256(data).hexdigest()


def fetch_url(url, tries=3, timeout=60):
    err = ""
    for k in range(tries):
        try:
            with urllib.request.urlopen(url, timeout=timeout) as r:
                return r.read()
        except Exception as ex:  # noqa: BLE001
            err = repr(ex)[:160]
            time.sleep(10 * (k + 1))
    raise RuntimeError("fetch failed %s: %s" % (url, err))


def check_plan(wait_min, plan_path=PLAN_PATH):
    deadline = time.time() + wait_min * 60
    while True:
        try:
            data = fetch_url(REPO_RAW + plan_path, tries=1)
            return True, sha256(data)
        except Exception:  # noqa: BLE001
            if time.time() > deadline:
                return False, ""
            time.sleep(300)


class Data:
    """FDSN access with an on-disk cache in data_dir (MiniSEED / StationXML)."""

    def __init__(self, data_dir):
        from obspy.clients.fdsn import Client
        self.dir = os.path.join(data_dir, "sdv23")
        os.makedirs(self.dir, exist_ok=True)
        self.wf = Client(BASE, timeout=180, _discover_services=False)
        self.ev = Client(USGS, timeout=180, _discover_services=False)

    def _retry(self, fn, tries=3):
        from obspy.clients.fdsn.header import FDSNNoDataException
        err = None
        for k in range(tries):
            try:
                return fn()
            except FDSNNoDataException:
                return None
            except Exception as ex:  # noqa: BLE001
                err = ex
                time.sleep(20 * (k + 1))
        raise err

    def waveforms(self, net, sta, loc, cha, t1, t2):
        from obspy import read
        key = "%s.%s.%s.%s.%d.%d" % (net, sta, loc.replace("-", "_").replace("*", "x"),
                                     cha.replace("?", "x"), int(t1.timestamp), int(t2.timestamp))
        f = os.path.join(self.dir, key + ".mseed")
        if os.path.exists(f):
            return read(f)
        if os.path.exists(f + ".nodata"):
            return None
        st = self._retry(lambda: self.wf.get_waveforms(net, sta, loc, cha, t1, t2))
        if st is None or len(st) == 0:
            open(f + ".nodata", "w").close()
            return None
        for tr in st:
            if tr.data.dtype.kind not in "if" or isinstance(tr.data, np.ma.MaskedArray):
                tr.data = np.asarray(tr.data, dtype=np.float64)
        st.write(f, format="MSEED")
        return st

    def stations(self, key, **kw):
        from obspy import read_inventory
        f = os.path.join(self.dir, "inv_" + key + ".xml")
        if os.path.exists(f):
            return read_inventory(f)
        inv = self._retry(lambda: self.wf.get_stations(**kw))
        if inv is not None:
            inv.write(f, format="STATIONXML")
        return inv

    def events(self, **kw):
        return self._retry(lambda: self.ev.get_events(**kw))


def phase_info(model, depth, dist_deg):
    """First P-type, S-type and SKS arrival: name -> (time s, horizontal slowness s/km)."""
    names = {"P": "P", "Pdiff": "P", "S": "S", "Sdiff": "S", "SKS": "SKS"}
    arr = model.get_travel_times(source_depth_in_km=depth, distance_in_degree=dist_deg,
                                 phase_list=["P", "Pdiff", "S", "Sdiff", "SKS"])
    out = {}
    for a in arr:
        k = names.get(a.name)
        if k and k not in out:
            out[k] = (float(a.time), float(a.ray_param) / (R_EARTH - depth))
    return out


def process(st, inv, t0, tb, te, output, pre_filt, band, sr_up=None):
    """Remove response, band-pass (zero phase), trim to [t0+tb, t0+te]; optional FFT
    upsampling. Returns comp -> dict(x, sr, off, az, dip, id)."""
    st = st.copy()
    for tr in st:
        if tr.stats.sampling_rate not in (20.0, 40.0, 50.0, 100.0, 200.0):
            tr.resample(20.0)
    st.merge(method=1, fill_value=0)
    st.detrend("demean")
    st.detrend("linear")
    st.taper(0.05)
    if inv is not None:
        st.remove_response(inventory=inv, output=output, pre_filt=pre_filt)
    else:
        st.remove_response(output=output, pre_filt=pre_filt)
    st.filter("bandpass", freqmin=band[0], freqmax=band[1], corners=2, zerophase=True)
    out = {}
    for tr in st:
        tr.trim(t0 + tb, t0 + te, pad=True, fill_value=0.0, nearest_sample=True)
        x = np.asarray(tr.data, float)
        sr = float(tr.stats.sampling_rate)
        if sr_up:
            f = int(round(sr_up / sr))
            x = upsample(x, f)
            sr = sr * f
        c = tr.stats.channel[-1]
        c = {"N": "1", "E": "2"}.get(c, c)
        az, dip = float("nan"), float("nan")
        if inv is not None:
            try:
                o = inv.get_orientation(tr.id, datetime=t0)
                az, dip = float(o["azimuth"]), float(o["dip"])
            except Exception:  # noqa: BLE001
                pass
        out[c] = dict(x=x, sr=sr, off=float(tr.stats.starttime - (t0 + tb)), az=az, dip=dip,
                      id=tr.id)
    return out


# ----------------------------------------------------------------- analysis steps

def ev_dict(row):
    from obspy import UTCDateTime
    return dict(id=row[0], iso=row[1], t=UTCDateTime(row[1]), lat=row[2], lon=row[3],
                depth=row[4], mag=row[5], magtype=row[6])


def membership(D, model, M, C, P, sta_cache):
    """Mantle-P cc of candidate C against master M at IU/II stations 30-90 deg."""
    from obspy.geodetics import gps2dist_azimuth, locations2degrees

    def sta_set(E):
        if E["id"] in sta_cache:
            return sta_cache[E["id"]]
        inv = D.stations("p_" + E["id"], network="IU,II", channel="BHZ", starttime=E["t"],
                         endtime=E["t"] + 600, level="channel")
        out = {}
        if inv is not None:
            for net in inv:
                for s in net:
                    out[(net.code, s.code)] = (s.latitude, s.longitude)
        sta_cache[E["id"]] = out
        return out

    sm, sc = sta_set(M), sta_set(C)
    rows = []
    for key in sorted(set(sm) & set(sc)):
        lat, lon = sm[key]
        dist = locations2degrees(M["lat"], M["lon"], lat, lon)
        if P["p_dist"][0] <= dist <= P["p_dist"][1]:
            az = gps2dist_azimuth(M["lat"], M["lon"], lat, lon)[1]
            rows.append(dict(net=key[0], sta=key[1], dist=dist, az=az))
    rows.sort(key=lambda r: r["dist"])
    rows = rows[:P["p_max_stations"]]
    out = []
    for r in rows:
        ph = phase_info(model, M["depth"], r["dist"])
        if "P" not in ph:
            continue
        tP, pP = ph["P"]
        rec = dict(cand=C["id"], net=r["net"], sta=r["sta"], dist=round(r["dist"], 3),
                   az=round(r["az"], 3), lag=float("nan"), cc=float("nan"), status="")
        try:
            arrs = []
            for E in (M, C):
                st = D.waveforms(r["net"], r["sta"], "*", "BHZ", E["t"] + tP - 120,
                                 E["t"] + tP + 120)
                if st is None:
                    raise RuntimeError("nodata " + E["id"])
                locs = sorted({tr.stats.location for tr in st},
                              key=lambda l: {"00": 0, "": 1, "10": 2}.get(l, 3))
                st = st.select(location=locs[0])
                inv = D.stations("p_%s_%s_%s" % (r["net"], r["sta"], E["id"]), network=r["net"],
                                 station=r["sta"], location=locs[0] or "--", channel="BHZ",
                                 starttime=E["t"] + tP - 120, endtime=E["t"] + tP + 120,
                                 level="response")
                pr = process(st, inv, E["t"], tP - 60, tP + 60, "VEL", P["pre_filt_body"],
                             P["band_p"], P["sr_up"])
                arrs.append(pr["Z"])
            a, b = arrs
            sr = a["sr"]
            i0 = lambda t, o: int(round((t - (tP - 60) - o) * sr))  # noqa: E731
            s1, s2 = i0(tP - 3, a["off"]), i0(tP + 5, a["off"])
            pk = s1 + int(np.argmax(np.abs(a["x"][s1:s2])))
            k, cc, edge = window_lag(a["x"], b["x"], pk - int(1 * sr), pk + int(5 * sr),
                                     int(P["maxlag_p_s"] * sr))
            rec["lag"] = k / sr + b["off"] - a["off"]
            rec["cc"] = cc
            rec["status"] = "edge" if edge else "ok"
            rec["p_skm"] = pP
        except Exception as ex:  # noqa: BLE001
            rec["status"] = "error: " + repr(ex)[:120]
        out.append(rec)
    good = [r for r in out if r["status"] == "ok" and np.isfinite(r["cc"])]
    n_mem = sum(r["cc"] >= P["p_member_cc"] for r in good)
    n_strict = sum(r["cc"] >= P["p_strict_cc"] for r in good)
    passing = [r for r in good if r["cc"] >= P["p_member_cc"]]
    reloc = relocate([r["lag"] for r in passing], [r["p_skm"] for r in passing],
                     [r["az"] for r in passing]) if len(passing) >= P["reloc_min"] else None
    summ = dict(cand=C["id"], n_stations=len(out), n_measured=len(good), n_cc090=int(n_mem),
                n_cc095=int(n_strict), member=bool(n_mem >= P["p_member_n"]),
                strict=bool(n_strict >= P["p_strict_n"]), reloc=reloc)
    return summ, out


def sdv_event(D, E, tb, te):
    loc = sdv_location(E["iso"])
    st = D.waveforms("IU", "SDV", loc, "BH?", E["t"] + tb - 300, E["t"] + te + 300)
    if st is None:
        return None, None, loc
    inv = D.stations("sdv_%s" % E["id"], network="IU", station="SDV", location=loc,
                     channel="BH?", starttime=E["t"] + tb - 300, endtime=E["t"] + te + 300,
                     level="response")
    return st, inv, loc


def measure_sdv_pair(A, B, rots, baz, wins, tb, P):
    """A, B: processed body-wave dicts (comps '1','2','Z'). rots: name -> (azA1, azA2,
    azB1, azB2). wins: name -> (comp, t1, t2), must include 'S_T'. Returns list of row
    dicts. Amendment 1: per rotation, B is first aligned to A on the S_T window
    (+/- coarse_maxlag_s, integer samples); every fine lag (+/- maxlag_body_s) is then
    measured on the shifted B and reported as the total lag (coarse + fine)."""
    rows = []
    sr = A["Z"]["sr"]
    offA = float(np.mean([A[c]["off"] for c in ("1", "2", "Z")]))
    offB = float(np.mean([B[c]["off"] for c in ("1", "2", "Z")]))
    ml = int(P["maxlag_body_s"] * sr)
    idx = lambda t: int(round((t - tb - offA) * sr))  # noqa: E731
    for rname, (a1, a2, b1, b2) in rots.items():
        RA, TA = to_rt(A["1"]["x"], A["2"]["x"], a1, a2, baz)
        RB, TB = to_rt(B["1"]["x"], B["2"]["x"], b1, b2, baz)
        _, tT1, tT2 = wins["S_T"]
        kc, _, edc = window_lag(TA, TB, idx(tT1), idx(tT2), int(P["coarse_maxlag_s"] * sr))
        sh = int(round(kc)) if (np.isfinite(kc) and not edc) else 0
        comps = {"R": (RA, shift_left(RB, sh)), "T": (TA, shift_left(TB, sh)),
                 "Z": (A["Z"]["x"], shift_left(B["Z"]["x"], sh))}
        for wname, (comp, t1, t2) in wins.items():
            a, b = comps[comp]
            k, cc, edge = window_lag(a, b, idx(t1), idx(t2), ml)
            rows.append(dict(rotation=rname, window=wname, comp=comp, t1=round(t1, 3),
                             t2=round(t2, 3), lag_s=(k + sh) / sr + offB - offA, cc=cc,
                             edge=edge, coarse_s=sh / sr, coarse_ok=bool(np.isfinite(kc)
                                                                         and not edc)))
    return rows


def run(args):
    from obspy import UTCDateTime
    from obspy.geodetics import gps2dist_azimuth, locations2degrees
    from obspy.taup import TauPyModel

    P = PARAMS
    t_start = time.time()
    budget = args.time_budget_min * 60
    os.makedirs(args.out_dir, exist_ok=True)
    summary = dict(version=VERSION, started=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
                   params=P, errors=[])
    with open(os.path.abspath(__file__), "rb") as f:
        summary["code_sha256"] = sha256(f.read())
    ok, psha = check_plan(args.wait_min)
    summary["plan_public_at_start"] = ok
    summary["plan_sha256"] = psha
    log("plan public:", ok, psha[:16])
    ok2, asha = check_plan(args.wait_min, AMEND_PATH)
    summary["amendment_public_at_start"] = ok2
    summary["amendment_sha256"] = asha
    log("amendment 1 public:", ok2, asha[:16])

    D = Data(args.data_dir)
    model = TauPyModel("iasp91")
    EV = {r[0]: ev_dict(r) for r in EVENTS}
    M = EV[MASTER]
    dist_deg = locations2degrees(M["lat"], M["lon"], SDV["lat"], SDV["lon"])
    dist_km, az_src, baz = gps2dist_azimuth(M["lat"], M["lon"], SDV["lat"], SDV["lon"])
    dist_km /= 1000.0
    ph = phase_info(model, M["depth"], dist_deg)
    tP, tS, tSKS = ph["P"][0], ph["S"][0], ph["SKS"][0]
    pP, pS, pSKS = ph["P"][1], ph["S"][1], ph["SKS"][1]
    summary["geometry"] = dict(dist_deg=dist_deg, dist_km=dist_km, az_src_to_sdv=az_src,
                               baz=baz, tP=tP, tS=tS, tSKS=tSKS, pP=pP, pS=pS, pSKS=pSKS)
    log("SDV dist %.2f deg baz %.2f tP %.1f tSKS %.1f tS %.1f" % (dist_deg, baz, tP, tSKS, tS))

    tb, te = tP - 300.0, dist_km / 2.6 + 150.0
    tbb, teb = tP - 150.0, tS + 120.0        # body-wave trim window
    body, lp, sdvmeta, ev_rows = {}, {}, {}, []
    for eid, E in EV.items():
        row = dict(id=eid, time=E["iso"], lat=E["lat"], lon=E["lon"], depth=E["depth"],
                   mag=E["mag"], magtype=E["magtype"], sdv_loc=sdv_location(E["iso"]),
                   sdv_status="", az1_meta=float("nan"), az2_meta=float("nan"))
        try:
            st, inv, loc = sdv_event(D, E, tb, te)
            if st is None:
                row["sdv_status"] = "nodata"
            else:
                chans = sorted({tr.stats.channel for tr in st})
                b = process(st, inv, E["t"], tbb, teb, "DISP", P["pre_filt_body"],
                            P["band_body"], P["sr_up"])
                if not all(c in b for c in ("1", "2", "Z")):
                    row["sdv_status"] = "missing components " + ",".join(chans)
                else:
                    body[eid] = b
                    lp[eid] = {bb: process(st, inv, E["t"], tb, te, "VEL", P["pre_filt_lp"], bb)
                               for bb in P["band_lp"] + [P["band_rz"]]}
                    row["az1_meta"], row["az2_meta"] = b["1"]["az"], b["2"]["az"]
                    sdvmeta[eid] = (b["1"]["az"], b["2"]["az"])
                    row["sdv_status"] = "ok"
        except Exception as ex:  # noqa: BLE001
            row["sdv_status"] = "error: " + repr(ex)[:160]
        ev_rows.append(row)
        log("SDV", eid, row["sdv_status"])

    # ---- windows from the master (metadata rotation)
    wins = None
    if MASTER in body:
        A = body[MASTER]
        sr = A["Z"]["sr"]
        off = A["Z"]["off"]
        RA, TA = to_rt(A["1"]["x"], A["2"]["x"], sdvmeta[MASTER][0], sdvmeta[MASTER][1], baz)
        idx = lambda t: int(round((t - tbb - off) * sr))  # noqa: E731
        tt = lambda i: tbb + off + i / sr  # noqa: E731

        def peak(x, tc, lo=-P["peak_search_s"], hi=P["peak_search_s"]):
            i1, i2 = idx(tc + lo), idx(tc + hi)
            return tt(i1 + int(np.argmax(np.abs(x[i1:i2]))))
        kR = peak(RA, tSKS)
        kZ = peak(A["Z"]["x"], tSKS)
        kT = peak(TA, tS, hi=8.0)
        kRs = peak(RA, tS, hi=8.0)
        kP = peak(A["Z"]["x"], tP)
        wins = {"SKS_R": ("R", kR - 2.0, kR + 3.0), "postSKS_R": ("R", kR + 5.0, kR + 11.0),
                "SKS_Z": ("Z", kZ - 2.0, kZ + 3.0), "S_T": ("T", kT - 3.0, kT + 10.0),
                "S_R_early": ("R", kRs - 2.0, kRs + 1.5), "S_R_mid": ("R", kRs + 1.5, kRs + 5.0),
                "P_Z": ("Z", kP - 1.0, kP + 4.0)}
        summary["windows"] = {k: [v[0], round(v[1], 3), round(v[2], 3)] for k, v in wins.items()}
        noise = (idx(tP - 120.0), idx(tP - 20.0))
    slow = {"SKS_R": pSKS, "postSKS_R": pSKS, "SKS_Z": pSKS, "S_T": pS, "S_R_early": pS,
            "S_R_mid": pS, "P_Z": pP}

    # ---- membership + relocation (mantle P)
    sta_cache, mem, mem_rows = {}, {}, []
    for eid, E in EV.items():
        if eid == MASTER:
            continue
        if time.time() - t_start > budget * 0.5:
            summary["errors"].append("membership stopped by time budget at " + eid)
            break
        try:
            s, rows = membership(D, model, M, E, P, sta_cache)
            mem[eid] = s
            mem_rows += rows
            log("member", eid, s["n_cc090"], s["n_cc095"], s["member"], s["reloc"])
        except Exception as ex:  # noqa: BLE001
            mem[eid] = dict(cand=eid, member=False, strict=False, reloc=None,
                            error=repr(ex)[:160])
    summary["membership"] = mem

    # ---- per pair: orientation + SDV differential times
    pair_rows, pairs = [], {}
    for eid in EV:
        if eid == MASTER or MASTER not in body or eid not in body or wins is None:
            continue
        A, B = body[MASTER], body[eid]
        res = dict(pair="%s-%s" % (MASTER, eid))
        # A. delta from long-period waveforms (two bands) and Rayleigh-Z
        dw = []
        for bb in P["band_lp"]:
            la, lb = lp[MASTER][bb], lp[eid][bb]
            srl = la["Z"]["sr"]
            i1 = int(round((tS - 60.0 - tb - la["Z"]["off"]) * srl))
            i2 = int(round((dist_km / 2.8 - tb - la["Z"]["off"]) * srl))
            ml = int(P["maxlag_lp_s"] * srl)
            k, ccz, _ = window_lag(la["Z"]["x"], lb["Z"]["x"], i1, i2, ml)
            ki = int(round(k)) if np.isfinite(k) else 0
            daz_a = la["2"]["az"] - la["1"]["az"]
            daz_b = lb["2"]["az"] - lb["1"]["az"]
            d, cc = orient_pair(la["1"]["x"][i1:i2], la["2"]["x"][i1:i2],
                                lb["1"]["x"][i1 + ki:i2 + ki], lb["2"]["x"][i1 + ki:i2 + ki],
                                daz_a, daz_b, P["pair_half"], P["grid_step"])
            dw.append(dict(band=list(bb), delta=d, cc=cc, lag_z_s=k / srl, cc_z=ccz))
        res["delta_wf"] = dw
        rz = {}
        for e2 in (MASTER, eid):
            L = lp[e2][P["band_rz"]]
            srl = L["Z"]["sr"]
            i1 = int(round((dist_km / 4.2 - tb - L["Z"]["off"]) * srl))
            i2 = int(round((dist_km / 3.0 - tb - L["Z"]["off"]) * srl))
            daz = L["2"]["az"] - L["1"]["az"]
            a_hat, c = orient_rz(L["1"]["x"][i1:i2], L["2"]["x"][i1:i2], L["Z"]["x"][i1:i2],
                                 daz, baz, L["1"]["az"], P["rz_half"], P["grid_step"])
            rz[e2] = dict(a_hat=a_hat, c=c)
        res["rz"] = rz
        dvals = [x["delta"] for x in dw if np.isfinite(x["delta"])]
        delta = float(np.mean(dvals)) if dvals else float("nan")
        spread = float(max(dvals) - min(dvals)) if len(dvals) == 2 else float("nan")
        rz_ok = all(rz[e]["c"] >= P["rz_c_min"] for e in rz)
        d_rz = float(wrap180(rz[eid]["a_hat"] - rz[MASTER]["a_hat"])) if rz_ok else float("nan")
        dwf_ok = (len(dvals) == 2 and all(x["cc"] >= P["dwf_cc_min"] for x in dw)
                  and spread <= P["dwf_band_spread_max"])
        agree_ok = (not rz_ok) or abs(delta - d_rz) <= P["dwf_rz_agree_max"]
        res.update(delta=delta, delta_spread=spread, delta_rz=d_rz, dwf_ok=bool(dwf_ok),
                   rz_ok=bool(rz_ok), agree_ok=bool(agree_ok))
        aA1, aA2 = sdvmeta[MASTER]
        aB1, aB2 = sdvmeta[eid]
        rots = {"meta": (aA1, aA2, aB1, aB2), "naive": (0.0, 90.0, 0.0, 90.0)}
        if dwf_ok:
            rots["cons"] = (aA1, aA2, aA1 + delta, aA2 + delta)
        rows = measure_sdv_pair(A, B, rots, baz, wins, tbb, P)
        sr = A["Z"]["sr"]
        R_A = {r: to_rt(A["1"]["x"], A["2"]["x"], v[0], v[1], baz)[0] for r, v in rots.items()}
        R_B = {r: to_rt(B["1"]["x"], B["2"]["x"], v[2], v[3], baz)[0] for r, v in rots.items()}
        iS1 = int(round((wins["SKS_R"][1] - tbb - A["Z"]["off"]) * sr))
        iS2 = int(round((wins["SKS_R"][2] - tbb - A["Z"]["off"]) * sr))
        res["snr_sks_r"] = {r: [snr(R_A[r], noise[0], noise[1], iS1, iS2),
                                snr(R_B[r], noise[0], noise[1], iS1, iS2)] for r in rots}
        reloc = mem.get(eid, {}).get("reloc")
        lagST = {r["rotation"]: r["lag_s"] for r in rows if r["window"] == "S_T"}
        lagPZ = {r["rotation"]: r["lag_s"] for r in rows if r["window"] == "P_Z"}
        for r in rows:
            r["pair"] = res["pair"]
            r["dt_vs_ST"] = r["lag_s"] - lagST[r["rotation"]]
            r["dt_vs_PZ"] = r["lag_s"] - lagPZ[r["rotation"]]
            if reloc:
                r["locpred_vs_ST"] = loc_pred(reloc["dxN"], reloc["dxE"], az_src,
                                              slow[r["window"]] - pS)
                r["locpred_vs_PZ"] = loc_pred(reloc["dxN"], reloc["dxE"], az_src,
                                              slow[r["window"]] - pP)
            else:
                r["locpred_vs_ST"] = r["locpred_vs_PZ"] = 0.0
            r["dtc_vs_ST"] = r["dt_vs_ST"] - r["locpred_vs_ST"]
            r["dtc_vs_PZ"] = r["dt_vs_PZ"] - r["locpred_vs_PZ"]
        res["relocated"] = bool(reloc)
        get = lambda rot, w, k="dtc_vs_ST": next((x[k] for x in rows  # noqa: E731
                                                   if x["rotation"] == rot and x["window"] == w),
                                                  float("nan"))
        res["x"] = {rot: get(rot, "SKS_R") for rot in rots}
        res["x_uncorrected"] = {rot: get(rot, "SKS_R", "dt_vs_ST") for rot in rots}
        res["cc_ST_meta"] = get("meta", "S_T", "cc")
        res["cc_SKS_R"] = {rot: get(rot, "SKS_R", "cc") for rot in rots}
        res["coarse_s"] = {rot: get(rot, "S_T", "coarse_s") for rot in rots}
        res["x_Z"] = get("meta", "SKS_Z")
        res["x_postSKS"] = {rot: get(rot, "postSKS_R") for rot in rots}
        res["x_S_R_early"] = {rot: get(rot, "S_R_early") for rot in rots}
        res["x_S_R_mid"] = {rot: get(rot, "S_R_mid") for rot in rots}
        res["x_vs_PZ"] = {rot: get(rot, "SKS_R", "dtc_vs_PZ") for rot in rots}
        pairs[eid] = res
        pair_rows += rows
        log("pair", res["pair"], "delta %.2f spread %.2f rz %.2f" % (delta, spread, d_rz),
            "x", {k: round(v, 4) for k, v in res["x"].items()})
    summary["pairs"] = pairs

    # ---- primary grade
    prim = dict(grade="untestable")
    if PARTNER in pairs:
        r = pairs[PARTNER]
        mres = mem.get(PARTNER, {})
        gates = dict(
            G1_partner_member=bool(mres.get("member")),
            G2_sdv_data=True,
            G3_ST_cc=bool(np.isfinite(r["cc_ST_meta"]) and r["cc_ST_meta"] >= P["st_cc_min"]),
            G4_sks_snr=all(np.isfinite(v) and v >= P["sks_snr_min"]
                           for vv in r["snr_sks_r"].values() for v in vv),
            G4_sks_cc=("cons" in r["cc_SKS_R"]) and all(
                np.isfinite(v) and v >= P["sks_cc_min"] for v in r["cc_SKS_R"].values()),
            G5_delta=bool(r["dwf_ok"] and r["agree_ok"]),
        )
        xm, xn = r["x"].get("meta", float("nan")), r["x"].get("naive", float("nan"))
        xc = r["x"].get("cons", float("nan"))
        prim = dict(grade=grade_primary(gates, xm, xn, xc), gates=gates, x_meta=xm, x_naive=xn,
                    x_cons=xc, relocated=r["relocated"], delta=r["delta"],
                    delta_rz=r["delta_rz"], meta_step=META_STEP["00"])
    else:
        prim["gates"] = dict(G2_sdv_data=False)
    summary["primary"] = prim
    log("PRIMARY", prim.get("grade"), prim)

    # ---- history (secondary)
    hist = {}
    for eid, r in pairs.items():
        if eid == PARTNER:
            continue
        rot = "cons" if "cons" in r["x"] else "meta"
        okh = (bool(mem.get(eid, {}).get("member")) and np.isfinite(r["cc_ST_meta"])
               and r["cc_ST_meta"] >= P["st_cc_min"]
               and all(np.isfinite(v) and v >= P["sks_snr_min"] for v in r["snr_sks_r"][rot])
               and r["cc_SKS_R"][rot] >= P["sks_cc_min"])
        hist[eid] = dict(rotation=rot, x=r["x"][rot], member=bool(mem.get(eid, {}).get("member")),
                         grade=grade_history(r["x"][rot], okh))
    summary["history"] = hist

    # ---- S2: orientation time series from global earthquakes
    orient_rows = []
    s2 = {}
    windows = [("era2000", "1999-11-06", "2001-06-30", 6.8, ["00"]),
               ("before", "2007-09-18", "2009-09-17T23:59:59", 6.5, ["00", "10"]),
               ("after", STEP_TIME, "2010-01-16T23:59:59", 6.0, ["00", "10"])]
    for wname, t1, t2, mmin, locs in windows:
        try:
            cat = D.events(starttime=UTCDateTime(t1), endtime=UTCDateTime(t2),
                           minmagnitude=mmin, latitude=SDV["lat"], longitude=SDV["lon"],
                           minradius=P["s2_dist"][0], maxradius=P["s2_dist"][1],
                           maxdepth=P["s2_maxdepth"], orderby="time-asc")
        except Exception as ex:  # noqa: BLE001
            summary["errors"].append("events %s: %s" % (wname, repr(ex)[:120]))
            continue
        for ev in (cat or []):
            if time.time() - t_start > budget:
                summary["errors"].append("S2 stopped by time budget in " + wname)
                break
            o = ev.preferred_origin() or ev.origins[0]
            dkm, _, bz = gps2dist_azimuth(o.latitude, o.longitude, SDV["lat"], SDV["lon"])
            dkm /= 1000.0
            for loc in locs:
                row = dict(window=wname, time=str(o.time), lat=o.latitude, lon=o.longitude,
                           loc=loc, dist_km=round(dkm, 1), baz=round(bz, 2), a_hat=float("nan"),
                           c=float("nan"), az1_meta=float("nan"), status="")
                try:
                    w1, w2 = dkm / 4.2 - 150.0, dkm / 3.0 + 150.0
                    st = D.waveforms("IU", "SDV", loc, "BH?", o.time + w1 - 300, o.time + w2 + 300)
                    if st is None:
                        row["status"] = "nodata"
                    else:
                        inv = D.stations("s2_%s_%d" % (loc, int(o.time.timestamp)), network="IU",
                                         station="SDV", location=loc, channel="BH?",
                                         starttime=o.time + w1 - 300, endtime=o.time + w2 + 300,
                                         level="response")
                        L = process(st, inv, o.time, w1, w2, "VEL", P["pre_filt_lp"],
                                    P["band_rz"])
                        srl = L["Z"]["sr"]
                        i1 = int(round((dkm / 4.2 - w1 - L["Z"]["off"]) * srl))
                        i2 = int(round((dkm / 3.0 - w1 - L["Z"]["off"]) * srl))
                        daz = L["2"]["az"] - L["1"]["az"]
                        a_hat, c = orient_rz(L["1"]["x"][i1:i2], L["2"]["x"][i1:i2],
                                             L["Z"]["x"][i1:i2], daz, bz, L["1"]["az"],
                                             P["rz_half"], 0.5)
                        row.update(a_hat=a_hat, c=c, az1_meta=L["1"]["az"],
                                   status="ok" if c >= P["rz_c_min"] else "low-c")
                except Exception as ex:  # noqa: BLE001
                    row["status"] = "error: " + repr(ex)[:120]
                orient_rows.append(row)
        log("S2 window", wname, "rows", len(orient_rows))
    for loc in ("00", "10"):
        bef = [r["a_hat"] for r in orient_rows if r["loc"] == loc and r["window"] == "before"
               and r["status"] == "ok"]
        aft = [r["a_hat"] for r in orient_rows if r["loc"] == loc and r["window"] == "after"
               and r["status"] == "ok"]
        s2[loc] = grade_rotation(bef, aft, META_STEP[loc])
        s2[loc]["meta_step"] = META_STEP[loc]
        s2[loc]["median_before"] = float(np.median(wrap180(bef))) if bef else float("nan")
        s2[loc]["median_after"] = float(np.median(wrap180(aft))) if aft else float("nan")
    e2k = [r["a_hat"] for r in orient_rows if r["loc"] == "00" and r["window"] == "era2000"
           and r["status"] == "ok"]
    s2["00"]["median_era2000"] = float(np.median(wrap180(e2k))) if e2k else float("nan")
    s2["00"]["n_era2000"] = len(e2k)
    summary["s2"] = s2
    log("S2", s2)

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
    wcsv("sdv23_events.csv", ev_rows)
    wcsv("sdv23_members.csv", mem_rows)
    wcsv("sdv23_pairs.csv", pair_rows)
    wcsv("sdv23_orient.csv", orient_rows)
    summary["finished"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    summary["runtime_min"] = round((time.time() - t_start) / 60.0, 1)
    with open(os.path.join(args.out_dir, "sdv23_summary.json"), "w") as f:
        json.dump(summary, f, indent=1, default=jsonable)
    try:
        plot(args.out_dir, body, sdvmeta, baz, wins, tbb, pairs, orient_rows)
    except Exception as ex:  # noqa: BLE001
        log("plot failed", repr(ex)[:200])
    return summary


def plot(out_dir, body, sdvmeta, baz, wins, tbb, pairs, orient_rows):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    os.makedirs(os.path.join(out_dir, "plots"), exist_ok=True)
    fig, ax = plt.subplots(3, 1, figsize=(9, 10))
    if MASTER in body and PARTNER in body and wins and PARTNER in pairs:
        A, B = body[MASTER], body[PARTNER]
        sr = A["Z"]["sr"]
        t1, t2 = wins["SKS_R"][1] - 3.0, wins["postSKS_R"][2]
        for k, (rot, lab) in enumerate([("meta", "station metadata"), ("cons", "consistent")]):
            r = pairs[PARTNER]
            if rot == "cons" and "cons" not in r["x"]:
                continue
            aA1, aA2 = sdvmeta[MASTER]
            if rot == "meta":
                b1, b2 = sdvmeta[PARTNER]
            else:
                b1, b2 = aA1 + r["delta"], aA2 + r["delta"]
            RA = to_rt(A["1"]["x"], A["2"]["x"], aA1, aA2, baz)[0]
            RB = to_rt(B["1"]["x"], B["2"]["x"], b1, b2, baz)[0]
            for x, off, lab2 in ((RA, A["Z"]["off"], "2000-07-07"), (RB, B["Z"]["off"],
                                                                       "2009-12-17")):
                i1 = int(round((t1 - tbb - off) * sr))
                i2 = int(round((t2 - tbb - off) * sr))
                seg = x[i1:i2]
                tt = t1 + np.arange(len(seg)) / sr
                ax[k].plot(tt, seg / (np.max(np.abs(seg)) or 1), label=lab2)
            ax[k].set_title("IU.SDV radial SKS, 0.4-2 Hz, rotation: " + lab)
            ax[k].legend()
            ax[k].set_xlabel("time after origin (s)")
    for loc, col in (("00", "C0"), ("10", "C1")):
        rr = [r for r in orient_rows if r["loc"] == loc and r["status"] == "ok"]
        if rr:
            from obspy import UTCDateTime
            ax[2].plot([UTCDateTime(r["time"]).decimal_year for r in rr],
                       wrap180([r["a_hat"] for r in rr]), "o", color=col, ms=3,
                       label="loc %s Rayleigh-Z azimuth of N channel" % loc)
    ax[2].axvline(2009.71, color="k", lw=0.8)
    ax[2].set_ylabel("deg")
    ax[2].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, "plots", "sdv23_overview.png"), dpi=110)


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", default="/data")
    ap.add_argument("--out-dir", default="/out")
    ap.add_argument("--wait-min", type=float, default=90.0)
    ap.add_argument("--time-budget-min", type=float, default=200.0)
    a = ap.parse_args(argv)
    run(a)
    return 0


if __name__ == "__main__":
    sys.exit(main())

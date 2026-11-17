#!/usr/bin/env python3
"""YKA: a second, independent array for inner-core-repeaters.

The Yellowknife array (network CN; legacy YKB0-9/YKR1-9 at 20 sps from 1989 to 2014,
current YKAB0-9/YKAR1-9 at 40 sps from 2013-06-26) sits about 133-142 deg from the South
Sandwich Islands (SSI). There, PKIKP (DF) has no PKPbc/PKPab companion, so repeaters
cannot be found from outer-core phases at YKA itself. This script therefore

  stage A (always): downloads YKA records around PKIKP for every catalogue event,
      makes 1-2 Hz DF beams with per-element statics, and writes diagnostics:
      archive coverage by year, a census of per-element static offsets (Wang et al.
      2024 report YKA clock errors of 0.125/0.25/0.375 s in 2013-2020) and quick-look
      plots of large events with ak135 phases marked;
  stage B (with --ilar-multiplets): takes the repeater groups found at ILAR from
      PKPbc+PKPab (outer core only, so independent of the inner core), and at YKA
      aligns each pair on the direct PKIKP pulse and measures how similar the
      following PKIKP coda is. The pre-registered predictions P1-P5 of model.py are
      then graded on YKA alone. See findings/2026-11-17-yka-prereg.md.

By MV2, Munim, Inc. (made with AI agents). Licence: Apache-2.0.
"""
import argparse
import csv
import io
import json
import math
import os
import time
import traceback
from itertools import combinations

import numpy as np

import model
from phases import DEPTHS, Phases, gc_deg
from xcorr import best_lag, cc, fshift, lag_scan, measure_pair, pick_onset

VERSION = "0.2.0"
FS = 40.0
PRE = 40.0
POST = 40.0
N = int(round((PRE + POST) * FS))
BAND = (1.0, 2.0)
LONG_PRE, LONG_POST = 70.0, 170.0  # download window around predicted PKIKP (s)
NET = "CN"
STATIONS = "YKB*,YKR*,YKAB*,YKAR*"
YKA_CENTER = (62.4931, -114.6061)
USGS = "https://earthquake.usgs.gov/fdsnws/event/1/query"
BOX = dict(minlatitude=-61.5, maxlatitude=-54.0, minlongitude=-31.0, maxlongitude=-22.0)
PROC_TAG = "yka_proc_v1"
# EarthScope station service, 2026-11-17. Legacy and current codes share positions, so both
# map to one canonical element name (YKB3 and YKAB3 -> B3, YKR7 and YKAR7 -> R7).
COORDS = {
    "B0": (62.6059, -114.6060), "B1": (62.4023, -114.6063), "B2": (62.4247, -114.6064),
    "B3": (62.4485, -114.6061), "B4": (62.4709, -114.6057), "B6": (62.5164, -114.6057),
    "B7": (62.5389, -114.6061), "B8": (62.5614, -114.6054), "B9": (62.5829, -114.6047),
    "R1": (62.4928, -114.9445), "R2": (62.4928, -114.8959), "R3": (62.4929, -114.8477),
    "R4": (62.4927, -114.7998), "R5": (62.4931, -114.7504), "R6": (62.4932, -114.7019),
    "R7": (62.4932, -114.6544), "R8": (62.4931, -114.6061), "R9": (62.4930, -114.5564),
}
ELEMENTS = sorted(COORDS)

# Pre-registered YKA windows, in seconds relative to the picked PKIKP onset.
DIRECT = (-0.3, 0.9)      # direct PKIKP pulse: used only to align the pair
CODA = (1.0, 6.0)         # PKIKP coda: the measured window (D)
CODA_SENS = {"D_c13": (1.0, 3.0), "D_c310": (3.0, 10.0), "D_full": (-0.3, 6.0)}
T0 = time.time()


def log(*a):
    print("[%7.1f min]" % ((time.time() - T0) / 60.0), *a, flush=True)


def elapsed_min():
    return (time.time() - T0) / 60.0


def canonical(sta):
    """YKAB3 / YKB3 -> 'B3'; YKAR7 / YKR7 -> 'R7'; anything else -> None."""
    s = str(sta).upper()
    if s.startswith("YKA") and len(s) == 5:
        c = s[3:]
    elif s.startswith("YK") and len(s) == 4:
        c = s[2:]
    else:
        return None
    return c if c in COORDS else None


def generation(sta):
    return "new" if str(sta).upper().startswith("YKA") else "legacy"


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


def align_member(root, beam, fs, pre_s, coarse_lag_s=5.0, fine_lag_s=0.4, allow_flip=False):
    """Offset (samples) of `beam` relative to `root` around the predicted PKIKP at pre_s.

    Coarse: envelope correlation over pre_s-4 .. pre_s+14 s (no cycle skipping, polarity-blind).
    Fine: waveform correlation over pre_s-2 .. pre_s+4 s within +/- fine_lag_s. If allow_flip,
    the negated trace is also tried (for pairs across the 2013-2014 hardware change) and kept
    only if it beats the original by more than 0.1.
    """
    ea, eb = envelope(root, fs), envelope(beam, fs)
    L, C = lag_scan(ea, eb, int((pre_s - 4) * fs), int((pre_s + 14) * fs), int(coarse_lag_s * fs))
    l0, c0 = best_lag(L, C)
    out = dict(off=float("nan"), cc_env=c0, cc_fine=float("nan"), flip=False)
    if not math.isfinite(l0):
        return out
    y = fshift(beam, -l0)
    i0, i1, ml = int((pre_s - 2) * fs), int((pre_s + 4) * fs), int(fine_lag_s * fs)
    l1, c1 = best_lag(*lag_scan(root, y, i0, i1, ml))
    flip = False
    if allow_flip:
        l2, c2 = best_lag(*lag_scan(root, -y, i0, i1, ml))
        if math.isfinite(c2) and (not math.isfinite(c1) or c2 > c1 + 0.1):
            l1, c1, flip = l2, c2, True
    if not math.isfinite(l1):
        out.update(off=l0)
        return out
    out.update(off=l0 + l1, cc_fine=c1, flip=flip)
    return out


def yka_pair(xa, xb, onset, off_a, off_b, fs, direct=DIRECT, window=CODA, max_lag_ref_s=0.3,
             max_lag_df_s=0.25):
    """Align event b on event a using the direct PKIKP pulse only, then compare `window`.

    onset is a sample index in the multiplet (root) frame; event m's arrivals sit at
    index + off_m in its own trace. Returns measure_pair's dict (cc_ref = direct-pulse
    similarity, cc_df0 = window similarity after alignment) or None.
    """
    sa = int(round(off_a))
    ref = (onset + int(round(direct[0] * fs)) + sa, onset + int(round(direct[1] * fs)) + sa)
    dfw = (onset + int(round(window[0] * fs)) + sa, onset + int(round(window[1] * fs)) + sa)
    if ref[0] < 0 or dfw[0] < 0 or dfw[1] > len(xa) or ref[1] > len(xa):
        return None
    return measure_pair(xa, xb, ref, dfw, fs, max_lag_ref_s=max_lag_ref_s, max_lag_df_s=max_lag_df_s,
                        lag0=off_b - off_a)


def parse_members(s):
    """'us7000abcd(2005-06-01);usp000xyz(2001-02-03)' -> ['us7000abcd', 'usp000xyz']."""
    out = []
    for part in str(s or "").split(";"):
        part = part.strip()
        if not part:
            continue
        out.append(part.split("(")[0].strip())
    return [p for p in out if p]


# ----------------------------------------------------------------------------- I/O helpers
def write_csv(path, rows, keys=None):
    if not rows:
        with open(path, "w") as f:
            f.write("")
        return
    keys = keys or list(rows[0].keys())
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)


def jsonable(o):
    if isinstance(o, dict):
        return {str(k): jsonable(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [jsonable(v) for v in o]
    if isinstance(o, (np.floating, float)):
        return None if not math.isfinite(float(o)) else float(o)
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.bool_,)):
        return bool(o)
    return o


def fetch_catalog(data_dir, end_iso, min_mag_old, min_mag_new, new_from):
    """USGS ComCat SSI events (same query as run.py); cached copy used if the service fails."""
    import requests

    path = os.path.join(data_dir, "catalog_usgs.csv")
    rows, complete = {}, True
    for s, e, m in [("1990-01-01", new_from, min_mag_old), (new_from, end_iso, min_mag_new)]:
        params = dict(format="csv", starttime=s, endtime=e, minmagnitude=m, orderby="time-asc", **BOX)
        text = None
        for attempt in range(4):
            try:
                r = requests.get(USGS, params=params, timeout=180)
                if r.status_code == 200:
                    text = r.text
                    break
                log("USGS status", r.status_code)
            except Exception as ex:
                log("USGS error", repr(ex)[:200])
            time.sleep(10 * (attempt + 1))
        if text is None:
            complete = False
            continue
        for row in csv.DictReader(io.StringIO(text)):
            if row.get("type", "earthquake") == "earthquake":
                rows[row["id"]] = row
    if not complete and os.path.exists(path):
        with open(path, newline="") as f:
            for row in csv.DictReader(f):
                rows.setdefault(row["id"], row)
    if rows:
        keys = list(next(iter(rows.values())).keys())
        with open(path, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=keys, extrasaction="ignore")
            w.writeheader()
            for r in rows.values():
                w.writerow(r)
    evs = []
    for r in rows.values():
        try:
            evs.append(dict(id=r["id"], time=r["time"], lat=float(r["latitude"]), lon=float(r["longitude"]),
                            depth=float(r["depth"]) if r.get("depth") not in (None, "") else 10.0,
                            mag=float(r["mag"]) if r.get("mag") not in (None, "") else float("nan"),
                            t=model.decimal_year(r["time"])))
        except Exception:
            continue
    evs.sort(key=lambda e: e["t"])
    return evs


# ----------------------------------------------------------------------------- travel times
def build_grid(dists=None, depths=DEPTHS):
    from obspy.taup import TauPyModel

    if dists is None:
        dists = np.round(np.arange(128.0, 146.0 + 1e-9, 0.1), 2)
    m = TauPyModel(model="ak135")
    g = {k: np.full((len(depths), len(dists)), np.nan) for k in ("df", "ikp", "p_df")}
    for a, z in enumerate(depths):
        for b, d in enumerate(dists):
            try:
                arr = m.get_travel_times(source_depth_in_km=float(z), distance_in_degree=float(d),
                                         phase_list=["PKIKP", "PKiKP"])
            except Exception:
                continue
            df = [x for x in arr if x.name == "PKIKP"]
            ik = [x for x in arr if x.name == "PKiKP"]
            if df:
                g["df"][a, b] = df[0].time
                g["p_df"][a, b] = df[0].ray_param_sec_degree
            if ik:
                g["ikp"][a, b] = ik[0].time
    return dict(dists=np.asarray(dists, dtype=float), depths=np.asarray(depths, dtype=float), **g)


def load_or_build(path):
    if os.path.exists(path):
        z = np.load(path)
        return {k: z[k] for k in z.files}
    g = build_grid()
    np.savez(path, **g)
    return g


# ----------------------------------------------------------------------------- waveforms
# Checked 2026-11-17: EarthScope serves the current YKAR/YKAB channels (data for 2025) but
# returned nothing for legacy YKR1 in 1993 or 2005; NRCan (Canadian National Seismograph
# Network) served legacy YKR1 for both 1993 and 2005. So legacy-era events go to NRCan
# first and current-era events to EarthScope first, each falling back to the other.
EARTHSCOPE = "https://service.earthscope.org/fdsnws/dataselect/1/query"
NRCAN = "https://earthquakescanada.nrcan.gc.ca/fdsnws/dataselect/1/query"
NEW_ERA = 2013.45  # YKAR/YKAB from 2013-06-26 (overlap with legacy until 2014-04-02)


def providers_for(t):
    return (EARTHSCOPE, NRCAN) if t >= NEW_ERA else (NRCAN, EARTHSCOPE)


def dataselect(url, t1, t2, session=None):
    """One FDSN dataselect request. Returns (Stream | None, 'ok' | 'nodata' | 'error:...')."""
    import requests
    from obspy import read

    params = dict(net=NET, sta=STATIONS, loc="*", cha="SHZ", start=t1.strftime("%Y-%m-%dT%H:%M:%S"),
                  end=t2.strftime("%Y-%m-%dT%H:%M:%S"))
    s = session or requests
    err = ""
    for attempt in range(3):
        try:
            r = s.get(url, params=params, timeout=120)
            if r.status_code in (204, 404):
                return None, "nodata"
            if r.status_code == 200 and r.content:
                st = read(io.BytesIO(r.content))
                return (st, "ok") if len(st) else (None, "nodata")
            err = "status %d" % r.status_code
        except Exception as ex:
            err = repr(ex)[:160]
        time.sleep(5 * (attempt + 1))
    return None, "error:" + err


def fetch_event(session, ev, ph, data_dir):
    from obspy import UTCDateTime, read

    rawdir = os.path.join(data_dir, "yka_raw")
    os.makedirs(rawdir, exist_ok=True)
    p = os.path.join(rawdir, ev["id"] + ".mseed")
    pn = p + ".nodata"
    if os.path.exists(p):
        try:
            return read(p), "cached"
        except Exception:
            os.remove(p)
    if os.path.exists(pn):
        return None, "nodata"
    dist = float(gc_deg(ev["lat"], ev["lon"], *YKA_CENTER))
    tdf = ph.t("df", dist, ev["depth"])
    if not math.isfinite(tdf):
        return None, "no-df-prediction(dist=%.2f)" % dist
    o = UTCDateTime(ev["time"])
    t1, t2 = o + tdf - LONG_PRE, o + tdf + LONG_POST
    results, best = [], None
    for url in providers_for(ev["t"]):
        st, s = dataselect(url, t1, t2, session)
        results.append(s)
        if st is not None:
            n_el = len({canonical(tr.stats.station) for tr in st} - {None})
            if best is None or n_el > best[1]:
                best = (st, n_el, url)
            if n_el >= 5:
                break  # enough elements; otherwise also try the other archive
    if best is not None:
        st, _, url = best
        try:
            st.write(p, format="MSEED")
        except Exception as ex:
            log("could not cache", ev["id"], repr(ex)[:120])
        return st, "fetched:" + ("earthscope" if url == EARTHSCOPE else "nrcan")
    if all(s == "nodata" for s in results):
        open(pn, "w").close()
        return None, "nodata"
    return None, ";".join(results)


def process_event(st, ev, ph):
    """1-2 Hz causal filter, per-element cut on its own ak135 PKIKP time, noise normalisation,
    per-element statics (+/-0.5 s, enough for the reported 0.375-s clock errors) and a beam."""
    from obspy import UTCDateTime

    o = UTCDateTime(ev["time"])
    depth = ev["depth"]
    st = st.copy()
    try:
        st.merge(method=0)
    except Exception:
        pass
    by_gen = {"legacy": {}, "new": {}}
    for tr in st:
        c = canonical(tr.stats.station)
        if c is None or tr.stats.channel != "SHZ":
            continue
        by_gen[generation(tr.stats.station)].setdefault(c, tr)
    gen = "new" if len(by_gen["new"]) >= len(by_gen["legacy"]) else "legacy"
    names, traces = [], []
    for c in ELEMENTS:
        tr = by_gen[gen].get(c)
        if tr is None:
            continue
        data = tr.data
        if np.ma.isMaskedArray(data):
            if np.ma.count_masked(data) > 0:
                continue
            data = np.asarray(data)
        sr = float(tr.stats.sampling_rate)
        if tr.stats.npts < 100 * sr:
            continue
        tr = tr.copy()
        tr.data = np.asarray(data, dtype=np.float64)
        try:
            tr.detrend("demean")
            tr.detrend("linear")
            tr.taper(0.05)
            tr.filter("bandpass", freqmin=BAND[0], freqmax=BAND[1], corners=4, zerophase=False)
            lat, lon = COORDS[c]
            tdf = ph.t("df", float(gc_deg(ev["lat"], ev["lon"], lat, lon)), depth)
            if not math.isfinite(tdf):
                continue
            tr.interpolate(sampling_rate=FS, method="lanczos", a=12, starttime=o + tdf - PRE, npts=N)
        except Exception:
            continue
        x = np.asarray(tr.data, dtype=np.float64)
        if len(x) != N or not np.all(np.isfinite(x)):
            continue
        names.append(c)
        traces.append(x)
    if len(traces) < 5:
        return None, "few-stations(%d,%s)" % (len(traces), gen)
    X = np.vstack(traces)
    n0, n1 = int((PRE - 35) * FS), int((PRE - 5) * FS)
    nrms = np.sqrt(np.mean(X[:, n0:n1] ** 2, axis=1))
    med = np.median(nrms)
    keep = (nrms > 0.2 * med) & (nrms < 5 * med) & (nrms > 0)
    X, nrms = X[keep], nrms[keep]
    names = [n for n, k in zip(names, keep) if k]
    if len(X) < 5:
        return None, "few-good-stations(%d,%s)" % (len(X), gen)
    X = X / nrms[:, None]
    g0, g1 = int((PRE - 3) * FS), int((PRE + 12) * FS)
    beam0 = X.mean(axis=0)
    lags = np.full(len(X), np.nan)
    ccs = np.full(len(X), np.nan)
    for k in range(len(X)):
        l, c = best_lag(*lag_scan(beam0, X[k], g0, g1, int(0.5 * FS)))
        ccs[k] = c
        if math.isfinite(l) and c >= 0.5:
            lags[k] = l
    good = np.isfinite(lags)
    if good.sum() >= 3:
        lags[good] -= np.median(lags[good])
    for k in range(len(X)):
        if good[k] and lags[k] != 0:
            X[k] = fshift(X[k], -lags[k])
    beam = X.mean(axis=0)
    nb = np.sqrt(np.mean(beam[n0:n1] ** 2))
    if not nb > 0:
        return None, "dead-beam"
    beam = beam / nb
    snr = float(np.sqrt(np.mean(beam[int((PRE - 1) * FS):int((PRE + 6) * FS)] ** 2)))
    dist_c = float(gc_deg(ev["lat"], ev["lon"], *YKA_CENTER))
    return dict(names=np.array(names), X=X.astype(np.float16), beam=beam, dist=dist_c, snr=snr,
                gen=gen, statics_s=lags / FS, static_cc=ccs), "ok"


def save_proc(path, d):
    np.savez(path, names=d["names"], X=d["X"], beam=d["beam"], statics_s=d["statics_s"],
             static_cc=d["static_cc"], meta=np.array([d["dist"], d["snr"], 1.0 if d["gen"] == "new" else 0.0]))


def load_proc(path, with_X=False):
    z = np.load(path, allow_pickle=False)
    m = z["meta"]
    d = dict(names=z["names"], beam=z["beam"], statics_s=z["statics_s"], static_cc=z["static_cc"],
             dist=float(m[0]), snr=float(m[1]), gen="new" if m[2] > 0.5 else "legacy")
    if with_X:
        d["X"] = z["X"].astype(np.float64)
    return d


# ----------------------------------------------------------------------------- plots
def quicklook(path, ev, st, d):
    """Long single-element record with ak135 phases marked, plus the DF beam zoom."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from obspy import UTCDateTime
    from obspy.taup import TauPyModel

    m = TauPyModel(model="ak135")
    dist = float(gc_deg(ev["lat"], ev["lon"], *YKA_CENTER))
    arr = m.get_travel_times(source_depth_in_km=max(0.0, float(ev["depth"])), distance_in_degree=dist,
                             phase_list=["PKIKP", "PKiKP", "pPKIKP", "sPKIKP", "Pdiff", "PP", "SKP", "PKS",
                                         "SKIKP", "PKIKS"])
    tdf = [a.time for a in arr if a.name == "PKIKP"]
    if not tdf:
        return
    tdf = tdf[0]
    fig, axs = plt.subplots(2, 1, figsize=(11, 6))
    o = UTCDateTime(ev["time"])
    tr = None
    for want in ("R5", "B4", "R3", "B6"):
        for t in st:
            if canonical(t.stats.station) == want:
                tr = t.copy()
                break
        if tr is not None:
            break
    if tr is not None:
        try:
            tr.detrend("demean")
            tr.taper(0.05)
            tr.filter("bandpass", freqmin=BAND[0], freqmax=BAND[1], corners=4, zerophase=False)
            tt = tr.times(reftime=o + tdf)
            y = tr.data / (np.max(np.abs(tr.data)) or 1.0)
            axs[0].plot(tt, y, lw=0.5, color="k")
            for a in arr:
                rel = a.time - tdf
                if -LONG_PRE < rel < LONG_POST:
                    axs[0].axvline(rel, color="C3", lw=0.6)
                    axs[0].text(rel, 1.02, a.name, rotation=90, fontsize=6, va="bottom")
            axs[0].set_xlim(-LONG_PRE, LONG_POST)
            axs[0].set_title("YKA %s %s M%.1f %.1f km, %.2f deg (1-2 Hz, element %s); times rel. ak135 PKIKP"
                             % (ev["id"], ev["time"][:19], ev["mag"], ev["depth"], dist,
                                canonical(tr.stats.station)), fontsize=8)
        except Exception:
            pass
    t = np.arange(N) / FS - PRE
    axs[1].plot(t, d["beam"], lw=0.7)
    for a in arr:
        rel = a.time - tdf
        if -10 < rel < 40:
            axs[1].axvline(rel, color="C3", lw=0.6)
            axs[1].text(rel, axs[1].get_ylim()[1], a.name, rotation=90, fontsize=6, va="top")
    axs[1].set_xlim(-10, 40)
    axs[1].set_xlabel("s relative to predicted PKIKP (beam, %d elements, %s)" % (len(d["names"]), d["gen"]))
    fig.tight_layout()
    fig.savefig(path, dpi=100)
    plt.close(fig)


def plot_statics(path, rows):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(11, 7))
    for k, el in enumerate(ELEMENTS):
        rr = [r for r in rows if r["element"] == el and r["static_s"] is not None
              and math.isfinite(r["static_s"])]
        if rr:
            ax.scatter([r["t"] for r in rr], [k + 2.0 * r["static_s"] for r in rr], s=3)
    ax.set_yticks(range(len(ELEMENTS)))
    ax.set_yticklabels(ELEMENTS)
    ax.set_xlabel("year")
    ax.set_ylabel("element (offset = 2 x static in s; 0.125 s = a quarter row)")
    ax.set_title("YKA per-element static offsets vs PKIKP-group beam (1-2 Hz), per event")
    fig.tight_layout()
    fig.savefig(path, dpi=100)
    plt.close(fig)


def plot_multiplet(path, mid, members, evs, data, off, flips, onset):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    order = sorted(members, key=lambda m: evs[m]["t"])
    if len(order) > 25:
        order = order[:12] + order[-13:]
    t = (np.arange(N) - onset) / FS
    fig, ax = plt.subplots(figsize=(10, 1.1 * len(order) + 1.5))
    for k, m in enumerate(order):
        b = fshift(data[m]["beam"], -off[m]) * (-1.0 if flips.get(m) else 1.0)
        seg = b[max(0, onset - int(3 * FS)):onset + int(10 * FS)]
        sc = np.max(np.abs(seg)) if len(seg) else 1.0
        ax.plot(t, b / (sc if sc > 0 else 1.0) * 0.45 - k, lw=0.8)
        ax.text(10.3, -k, "%s M%.1f %s" % (evs[m]["time"][:10], evs[m]["mag"], data[m]["gen"][:3]),
                fontsize=8, va="center")
    ax.axvspan(DIRECT[0], DIRECT[1], color="grey", alpha=0.2, label="direct PKIKP (alignment)")
    ax.axvspan(CODA[0], CODA[1], color="orange", alpha=0.2, label="PKIKP coda (measured)")
    ax.set_xlim(-3, 13)
    ax.set_yticks([])
    ax.set_xlabel("s relative to picked PKIKP onset")
    ax.set_title("ILAR multiplet %s seen at YKA (1-2 Hz beams)" % mid)
    ax.legend(loc="lower left", fontsize=7)
    fig.tight_layout()
    fig.savefig(path, dpi=110)
    plt.close(fig)


# ----------------------------------------------------------------------------- stage B
def analyze_multiplet(mid, members, evs, data, proc_dir, args):
    info = dict(mult=mid, n_ilar=len(members), status="ok")
    members = [m for m in members if m in data]
    info["n_yka"] = len(members)
    if len(members) < 2:
        info["status"] = "fewer-than-2-at-YKA"
        return [], info, None
    root = max(members, key=lambda m: data[m]["snr"])
    off, flips, align = {root: 0.0}, {root: False}, {}
    for m in members:
        if m == root:
            continue
        a = align_member(data[root]["beam"], data[m]["beam"], FS, PRE,
                         allow_flip=data[m]["gen"] != data[root]["gen"])
        align[m] = a
        if math.isfinite(a["off"]) and math.isfinite(a["cc_fine"]) and a["cc_fine"] >= 0.5:
            off[m], flips[m] = a["off"], a["flip"]
    members = [m for m in members if m in off]
    if len(members) < 2:
        info["status"] = "alignment-failed"
        return [], info, None
    stack = np.zeros(N)
    for m in members:
        stack += fshift(data[m]["beam"], -off[m]) * (-1.0 if flips[m] else 1.0)
    stack /= len(members)
    onset = pick_onset(stack, FS, center=int(PRE * FS), search_s=6.0, sta_s=0.5, lta_s=8.0, thr=6.0)
    info.update(root=evs[root]["id"], n_aligned=len(members), n_flipped=sum(bool(flips[m]) for m in members))
    if onset is None:
        info["status"] = "no-DF-pick"
        return [], info, None
    info["pick_minus_pred_s"] = round((onset - PRE * FS) / FS, 3)
    n0, n1 = int((PRE - 35) * FS), int((PRE - 5) * FS)

    def snr(m, w):
        b = data[m]["beam"]
        s = int(round(off[m]))
        a, c = onset + int(round(w[0] * FS)) + s, onset + int(round(w[1] * FS)) + s
        if a < 0 or c > N:
            return float("nan")
        nb = np.sqrt(np.mean(b[n0:n1] ** 2))
        return float(np.sqrt(np.mean(b[a:c] ** 2)) / nb) if nb > 0 else float("nan")

    snr_dir = {m: snr(m, DIRECT) for m in members}
    snr_coda = {m: snr(m, CODA) for m in members}
    Xc = {}
    rows = []
    order = sorted(members, key=lambda m: evs[m]["t"])
    for a, b in combinations(order, 2):
        xa = data[a]["beam"]
        sgn = -1.0 if flips[a] != flips[b] else 1.0
        xb = data[b]["beam"] * sgn
        r = yka_pair(xa, xb, onset, off[a], off[b], FS)
        if r is None:
            continue
        sens = {}
        for key, w in CODA_SENS.items():
            rs = yka_pair(xa, xb, onset, off[a], off[b], FS, window=w)
            sens[key] = (1.0 - rs["cc_df0"]) if rs is not None else float("nan")
        st_cc = []
        if len(members) <= args.max_station_level:
            for m in (a, b):
                if m not in Xc:
                    Xc[m] = load_proc(os.path.join(proc_dir, m + ".npz"), with_X=True)
            na, nb_ = list(Xc[a]["names"]), list(Xc[b]["names"])
            for s in sorted(set(na) & set(nb_)):
                rs = yka_pair(Xc[a]["X"][na.index(s)], Xc[b]["X"][nb_.index(s)] * sgn, onset, off[a], off[b], FS)
                if rs is not None and rs["cc_ref"] >= 0.6:
                    st_cc.append(rs["cc_df0"])
        ea, eb = evs[a], evs[b]
        D = 1.0 - r["cc_df0"]
        q = (snr_dir[a], snr_dir[b], snr_coda[a], snr_coda[b])
        okq = (all(math.isfinite(v) for v in q) and math.isfinite(D) and math.isfinite(r["cc_ref"])
               and r["cc_ref"] >= args.pair_thr and min(snr_dir[a], snr_dir[b]) >= args.min_snr_direct
               and min(snr_coda[a], snr_coda[b]) >= args.min_snr_coda)
        rows.append(dict(
            array="YKA", mult=mid, ev1=ea["id"], ev2=eb["id"], time1=ea["time"], time2=eb["time"],
            t1=ea["t"], t2=eb["t"], mag1=ea["mag"], mag2=eb["mag"], gen1=data[a]["gen"], gen2=data[b]["gen"],
            polarity_flipped=sgn < 0, cc_ref=r["cc_ref"], cc_df0=r["cc_df0"], cc_df_best=r["cc_df_best"],
            dt_df_s=r["dt_df_s"], D=D, Dx=D - (1.0 - r["cc_ref"]),
            Ds=(1.0 - float(np.median(st_cc))) if len(st_cc) >= 5 else float("nan"), n_sta=len(st_cc),
            snr_dir1=snr_dir[a], snr_dir2=snr_dir[b], snr_coda1=snr_coda[a], snr_coda2=snr_coda[b],
            dM=model.delta_M(ea["t"], eb["t"]), dB=model.delta_B(ea["t"], eb["t"]),
            cross_turn=ea["t"] < model.T_TURN < eb["t"], new_pair=eb["t"] >= model.NEW_CUT, ok=bool(okq),
            **sens))
    return rows, info, (off, flips, onset, members)


def stage_b(args, evs, data, proc_dir, plots_dir, summary):
    mults = []
    with open(args.ilar_multiplets, newline="") as f:
        for r in csv.DictReader(f):
            ids = [i for i in parse_members(r.get("members")) if i in evs]
            if len(ids) >= 2:
                mults.append((r.get("mult") or "M%03d" % len(mults), ids))
    log("ILAR multiplets read:", len(mults))
    rows, infos, frames = [], [], {}
    for mid, members in mults:
        if elapsed_min() > args.time_budget_min * 0.95:
            infos.append(dict(mult=mid, n_ilar=len(members), status="skipped-time-budget"))
            continue
        try:
            r, info, fr = analyze_multiplet(mid, members, evs, data, proc_dir, args)
        except Exception as ex:
            traceback.print_exc()
            r, info, fr = [], dict(mult=mid, n_ilar=len(members), status="exception:" + repr(ex)[:150]), None
        info["span"] = "%.1f-%.1f" % (min(evs[m]["t"] for m in members), max(evs[m]["t"] for m in members))
        info["has_new"] = any(evs[m]["t"] >= model.NEW_CUT for m in members)
        rows += r
        infos.append(info)
        if fr is not None:
            frames[mid] = fr
    write_csv(os.path.join(args.out_dir, "yka_multiplets.csv"), infos,
              keys=["mult", "n_ilar", "n_yka", "n_aligned", "n_flipped", "status", "root", "span", "has_new",
                    "pick_minus_pred_s"])
    write_csv(os.path.join(args.out_dir, "yka_pairs.csv"), rows)
    log("YKA pairs:", len(rows), "ok:", sum(r["ok"] for r in rows),
        "out-of-sample ok:", sum(r["ok"] and r["new_pair"] for r in rows))
    try:
        res = model.evaluate(rows)
        res["sensitivity_windows"] = {}
        ok = [r for r in rows if r["ok"]]
        for key in CODA_SENS:
            rr = [dict(r, D=r[key]) for r in ok if math.isfinite(r[key])]
            res["sensitivity_windows"][key] = dict(
                P1=model.grade_set([r for r in rr if r["t2"] < model.NEW_CUT], nboot=300),
                P2=model.grade_set([r for r in rr if r["t2"] >= model.NEW_CUT], min_n=6, nboot=300),
                P4=model.p4(rr))
    except Exception as ex:
        traceback.print_exc()
        res = dict(error=repr(ex))
    if args.ilar_pairs and os.path.exists(args.ilar_pairs):
        try:
            il = {}
            with open(args.ilar_pairs, newline="") as f:
                for r in csv.DictReader(f):
                    if str(r.get("ok")) == "True":
                        il[(r["ev1"], r["ev2"])] = float(r["D"])
            both = [(il[(r["ev1"], r["ev2"])], r["D"]) for r in rows if r["ok"] and (r["ev1"], r["ev2"]) in il]
            rho = model.spearman([b[0] for b in both], [b[1] for b in both])
            res["cross_array_D_spearman"] = dict(n=rho[2], rho=rho[0], p=rho[1])
        except Exception as ex:
            res["cross_array_D_spearman"] = dict(error=repr(ex))
    summary["yka_predictions"] = res
    with open(os.path.join(args.out_dir, "yka_predictions.json"), "w") as f:
        json.dump(jsonable(res), f, indent=1)
    shown = 0
    for mid in sorted(frames, key=lambda m: (not any(evs[x]["t"] >= model.NEW_CUT for x in frames[m][3]),
                                            -len(frames[m][3]))):
        if shown >= args.plots:
            break
        off, flips, onset, members = frames[mid]
        try:
            plot_multiplet(os.path.join(plots_dir, "yka_multiplet_%s.png" % mid), mid, members, evs, data, off,
                           flips, onset)
            shown += 1
        except Exception:
            traceback.print_exc()
    return res, infos


# ----------------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", default="/data")
    ap.add_argument("--out-dir", default="/out")
    ap.add_argument("--end", default=time.strftime("%Y-%m-%dT00:00:00", time.gmtime()))
    ap.add_argument("--new-from", default="2024-07-01")
    ap.add_argument("--min-mag-old", type=float, default=5.0)
    ap.add_argument("--min-mag-new", type=float, default=4.7)
    ap.add_argument("--ilar-multiplets", default="", help="multiplets.csv from run.py (enables stage B)")
    ap.add_argument("--ilar-pairs", default="", help="pairs.csv from run.py (cross-array comparison)")
    ap.add_argument("--pair-thr", type=float, default=0.90)
    ap.add_argument("--min-snr-direct", type=float, default=5.0)
    ap.add_argument("--min-snr-coda", type=float, default=2.0)
    ap.add_argument("--max-station-level", type=int, default=40)
    ap.add_argument("--time-budget-min", type=float, default=220.0)
    ap.add_argument("--fetch-frac", type=float, default=0.85)
    ap.add_argument("--quicklooks", type=int, default=14)
    ap.add_argument("--plots", type=int, default=30)
    ap.add_argument("--max-events", type=int, default=0)
    args = ap.parse_args()
    for d in (args.data_dir, args.out_dir):
        os.makedirs(d, exist_ok=True)
    plots_dir = os.path.join(args.out_dir, "plots")
    os.makedirs(plots_dir, exist_ok=True)
    proc_dir = os.path.join(args.data_dir, PROC_TAG)
    os.makedirs(proc_dir, exist_ok=True)
    summary = dict(version=VERSION, args=vars(args), started=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()))

    ph = Phases(load_or_build(os.path.join(args.data_dir, "ak135_yka_grid_v1.npz")))
    tab = []
    for dd in (132.0, 134.0, 136.0, 138.0, 140.0, 142.0):
        for z in (10.0, 100.0):
            df = ph.t("df", dd, z)
            tab.append(dict(dist=dd, depth=z, df=df, pkikp_minus_pkikp=ph.t("ikp", dd, z) - df))
    log("ak135 PKIKP and PKiKP-PKIKP:", json.dumps(jsonable(tab)))
    summary["phase_table"] = tab

    evs_list = fetch_catalog(args.data_dir, args.end, args.min_mag_old, args.min_mag_new, args.new_from)
    evs = {e["id"]: e for e in evs_list}
    dists = [float(gc_deg(e["lat"], e["lon"], *YKA_CENTER)) for e in evs_list]
    log("catalog events:", len(evs_list), "YKA distance range %.2f-%.2f deg" % (min(dists), max(dists)))
    summary["dist_range"] = [min(dists), max(dists)]

    wanted = set()
    if args.ilar_multiplets and os.path.exists(args.ilar_multiplets):
        with open(args.ilar_multiplets, newline="") as f:
            for r in csv.DictReader(f):
                wanted.update(parse_members(r.get("members")))

    def prio(e):
        g = 0 if e["id"] in wanted else (1 if e["t"] >= 2024.5 else (2 if e["t"] < 2010.0 else 3))
        return (g, -(e["mag"] if math.isfinite(e["mag"]) else 0))

    order = sorted(evs_list, key=prio)
    if args.max_events:
        order = order[:args.max_events]
    client, data, status, statics_rows, streams = None, {}, [], [], {}
    stop_at = args.time_budget_min * args.fetch_frac
    for k, e in enumerate(order):
        pp = os.path.join(proc_dir, e["id"] + ".npz")
        if os.path.exists(pp):
            try:
                data[e["id"]] = load_proc(pp)
                status.append(dict(id=e["id"], t=e["t"], mag=e["mag"], status="ok-cached",
                                   n_sta=len(data[e["id"]]["names"]), gen=data[e["id"]]["gen"],
                                   snr=round(data[e["id"]]["snr"], 2)))
                continue
            except Exception:
                os.remove(pp)
        if os.path.exists(os.path.join(args.data_dir, "yka_raw", e["id"] + ".mseed.nodata")):
            status.append(dict(id=e["id"], t=e["t"], mag=e["mag"], status="nodata"))
            continue
        if elapsed_min() > stop_at:
            status.append(dict(id=e["id"], t=e["t"], mag=e["mag"], status="skipped-time-budget"))
            continue
        if client is None:
            import requests

            client = requests.Session()
            client.headers["User-Agent"] = "inner-core-repeaters/%s (MV2, Munim, Inc.; mv2@munimtech.com)" % VERSION
        try:
            st, s = fetch_event(client, e, ph, args.data_dir)
            if st is None:
                status.append(dict(id=e["id"], t=e["t"], mag=e["mag"], status=s))
                continue
            d, s2 = process_event(st, e, ph)
            if d is None:
                status.append(dict(id=e["id"], t=e["t"], mag=e["mag"], status=s2))
                continue
            save_proc(pp, d)
            data[e["id"]] = dict((kk, v) for kk, v in d.items() if kk != "X")
            status.append(dict(id=e["id"], t=e["t"], mag=e["mag"], status="ok", n_sta=len(d["names"]),
                               gen=d["gen"], snr=round(d["snr"], 2), source=s))
        except Exception as ex:
            status.append(dict(id=e["id"], t=e["t"], mag=e["mag"], status="exception:" + repr(ex)[:150]))
        if (k + 1) % 50 == 0:
            log("events tried %d/%d, processed %d" % (k + 1, len(order), len(data)))
        time.sleep(0.2)
    write_csv(os.path.join(args.out_dir, "yka_event_status.csv"), status,
              keys=["id", "t", "mag", "status", "n_sta", "gen", "snr", "source"])
    by_year = {}
    for s in status:
        b = by_year.setdefault(int(s["t"]), dict(tried=0, ok=0))
        b["tried"] += 1
        b["ok"] += int(str(s["status"]).startswith("ok"))
    summary["coverage_by_year"] = by_year
    summary["n_processed"] = len(data)
    log("YKA processed events:", len(data))

    # statics census (clock offsets)
    for i, d in data.items():
        for el, s, c in zip(list(d["names"]), list(d["statics_s"]), list(d["static_cc"])):
            statics_rows.append(dict(id=i, t=evs[i]["t"], gen=d["gen"], snr=d["snr"], element=str(el),
                                     static_s=float(s) if math.isfinite(s) else None,
                                     cc=float(c) if math.isfinite(c) else None))
    statics_rows.sort(key=lambda r: (r["t"], r["element"]))
    write_csv(os.path.join(args.out_dir, "yka_statics.csv"), statics_rows)
    try:
        hi = [r for r in statics_rows if r["snr"] >= 4.0]
        plot_statics(os.path.join(plots_dir, "yka_statics.png"), hi)
        summary["statics_n_rows_snr4"] = len(hi)
        per = {}
        for r in hi:
            if r["static_s"] is None:
                continue
            per.setdefault((int(r["t"]), r["element"]), []).append(r["static_s"])
        summary["statics_abs_median_ge_0p1s"] = sorted("%d %s %.3f n=%d" % (y, el, float(np.median(v)), len(v))
                                                       for (y, el), v in per.items()
                                                       if len(v) >= 3 and abs(np.median(v)) >= 0.1)
    except Exception:
        traceback.print_exc()

    # quick looks: the largest events per era
    try:
        eras = [(1989, 2000), (2000, 2008.5), (2008.5, 2014), (2014, 2021), (2021, 2025), (2025, 2100)]
        picks = []
        for a, b in eras:
            cand = sorted([i for i in data if a <= evs[i]["t"] < b], key=lambda i: -data[i]["snr"])
            picks += cand[:max(1, args.quicklooks // len(eras))]
        from obspy import read

        for i in picks:
            st = streams.get(i)
            if st is None:
                p = os.path.join(args.data_dir, "yka_raw", i + ".mseed")
                st = read(p) if os.path.exists(p) else None
            if st is None:
                continue
            d = dict(data[i])
            quicklook(os.path.join(plots_dir, "yka_quicklook_%s_%s.png" % (evs[i]["time"][:10], i)), evs[i], st, d)
    except Exception:
        traceback.print_exc()

    if args.ilar_multiplets and os.path.exists(args.ilar_multiplets):
        stage_b(args, evs, data, proc_dir, plots_dir, summary)
    else:
        log("no ILAR multiplets given: stage A only (download, beams, diagnostics)")

    summary["finished"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    summary["minutes"] = round(elapsed_min(), 1)
    with open(os.path.join(args.out_dir, "yka_run_summary.json"), "w") as f:
        json.dump(jsonable(summary), f, indent=1)
    cov = summary["coverage_by_year"]
    with open(os.path.join(args.out_dir, "yka_report.md"), "w") as f:
        f.write("# inner-core-repeaters %s: YKA run report\n\n" % VERSION)
        f.write("Processed %d of %d catalogue events in %.1f min.\n\n" % (len(data), len(evs_list),
                                                                       summary["minutes"]))
        f.write("Coverage by year (tried/ok): " + ", ".join("%s: %d/%d" % (y, cov[y]["tried"], cov[y]["ok"])
                                                            for y in sorted(cov)) + "\n\n")
        f.write("Element-years with |median static| >= 0.1 s (SNR >= 4):\n\n")
        for s in summary.get("statics_abs_median_ge_0p1s", []):
            f.write("* %s\n" % s)
        if "yka_predictions" in summary:
            f.write("\n## Predictions on YKA\n\n```\n%s\n```\n" % json.dumps(jsonable(summary["yka_predictions"]),
                                                                       indent=1)[:20000])
    log("done")


if __name__ == "__main__":
    main()

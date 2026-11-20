#!/usr/bin/env python3
# Copyright 2026 MV2, Munim, Inc.
# Licensed under the Apache License, Version 2.0 (Apache-2.0). See LICENSE.
# Made with AI agents working in MV2's name.
"""chilepilot.py - is the Chile -> Central Asia path usable for a source-controlled
inner-core clock test?

Design pilot (plan: findings/2026-11-20-chile-path-plan.md). The SSI -> ILAR/YKA path is close
to its data limit. Chile -> Kazakhstan/Kyrgyzstan crosses the inner core on an independent path
at the same distances (about 143-165 deg) as SSI -> ILAR/YKA, and is recorded by three open
arrays: Makanchi (KZ MK01-MK09), Karatau (KZ KK01-KK09) and the Kyrgyz network KNET (KN).
For a time-spread sample of Chilean earthquakes in four magnitude bins (half before, half after
the published 2008.5 turn) this pilot measures, per array, whether waveforms exist and the
rms signal-to-noise ratio of an array beam aligned on ak135 PKIKP, in the 1-2 Hz band used by
the ILAR/YKA work, with a noise-only control. It also counts close event pairs in the ComCat
catalogue as a rough upper-bound proxy for repeaters straddling 2008.5. No temporal-change
claim is made.

Usage: python chilepilot.py --data-dir /data --out-dir /out [--wait-min 120]
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

REPO = "https://raw.githubusercontent.com/munimtechnologies/inner-core-repeaters/"
PLAN_URL = REPO + "main/findings/2026-11-20-chile-path-plan.md"
COMCAT = "https://earthquake.usgs.gov/fdsnws/event/1/query"
BOX = dict(minlatitude=-42.0, maxlatitude=-18.0, minlongitude=-76.0, maxlongitude=-66.0)
MAX_DEPTH = 300.0
START, END = "1995-01-01", "2026-11-20"
MIN_MAG = 4.5
BINS = ((4.5, 5.0), (5.0, 5.5), (5.5, 6.0), (6.0, 9.9))
PER_SIDE = 6                              # events per bin before and after the turn
TURN = 2008.5
DIST_RANGE = (143.0, 165.0)               # array-centre distance accepted (deg)

ARRAYS = (
    ("MKAR", tuple(("KZ", "MK%02d" % i) for i in range(1, 10))),
    ("KKAR", tuple(("KZ", "KK%02d" % i) for i in range(1, 10))),
    ("KNET", tuple(("KN", s) for s in ("AAK", "AML", "CHM", "EKS2", "KBK", "KZA", "TKM2",
                                        "UCH", "ULHL", "USP"))),
)
CHA_PREF = ("BHZ", "SHZ", "HHZ")
LOC_PREF = ("", "00", "10", "20")
MIN_EL = 3                                # elements needed for a beam

FS = 20.0
PRE, POST = 60.0, 90.0                    # s kept before / after the element's ak135 PKIKP
N = int(round((PRE + POST) * FS))
MARGIN = 60.0                             # extra s downloaded on each side
BANDS = (("primary", (1.0, 2.0)), ("alt", (0.6, 1.6)))
PREFILT = (0.3, 0.5, 7.0, 9.0)
NOISE = (-45.0, -10.0)                    # normalisation window, s rel. ak135 PKIKP
SIG = (-3.0, 7.0)                         # PKP group window (DF, and BC at these distances)
CODA = (1.0, 6.0)                         # the ILAR/YKA coda window
NULL = (-58.0, -48.0)                     # noise-only control, outside NOISE
SNR_REQ = 3.0
MIN_BEAMS = 4                             # beams needed in a bin for its median to count
PAIR_KM, PAIR_DM, PAIR_DZ, PAIR_MIN_MAG = 5.0, 0.2, 10.0, 5.0


def log(*a):
    print(time.strftime("%H:%M:%S"), *a, flush=True)


def isf(x):
    return isinstance(x, (int, float)) and not isinstance(x, bool) and math.isfinite(x)


def fnum(x):
    try:
        v = float(x)
    except (TypeError, ValueError):
        return NAN
    return v if math.isfinite(v) else NAN


def median(xs):
    v = sorted(x for x in xs if isf(x))
    if not v:
        return NAN
    n = len(v)
    return v[n // 2] if n % 2 else 0.5 * (v[n // 2 - 1] + v[n // 2])


def gc_deg(lat1, lon1, lat2, lon2):
    p1, p2 = math.radians(lat1), math.radians(lat2)
    dl = math.radians(lon2 - lon1)
    c = math.sin(p1) * math.sin(p2) + math.cos(p1) * math.cos(p2) * math.cos(dl)
    return math.degrees(math.acos(max(-1.0, min(1.0, c))))


def km_between(lat1, lon1, lat2, lon2):
    return 6371.0 * math.radians(gc_deg(lat1, lon1, lat2, lon2))


def decimal_year(iso):
    from datetime import datetime
    s = str(iso).rstrip("Z")[:19]
    d = datetime.strptime(s, "%Y-%m-%dT%H:%M:%S")
    y0 = datetime(d.year, 1, 1)
    y1 = datetime(d.year + 1, 1, 1)
    return d.year + (d - y0).total_seconds() / (y1 - y0).total_seconds()


def rms(x, a, b):
    a, b = int(round(a)), int(round(b))
    if a < 0 or b > len(x) or b - a < 2:
        return NAN
    return float(np.sqrt(np.mean(np.asarray(x[a:b], dtype=float) ** 2)))


def idx(rel_s, pre=PRE, fs=FS):
    """Sample index of a time rel_s seconds after ak135 PKIKP (which sits at `pre`)."""
    return (pre + rel_s) * fs


def window_snrs(x):
    """x: trace normalised by its NOISE-window rms. Returns (snr_sig, snr_coda, snr_null)."""
    return (rms(x, idx(SIG[0]), idx(SIG[1])), rms(x, idx(CODA[0]), idx(CODA[1])),
            rms(x, idx(NULL[0]), idx(NULL[1])))


def normalise(x):
    nz = rms(x, idx(NOISE[0]), idx(NOISE[1]))
    if not (isf(nz) and nz > 0):
        return None
    return np.asarray(x, dtype=float) / nz


def beam(traces):
    """Mean of noise-normalised, already aligned traces, renormalised by its own noise rms.
    Returns None with fewer than MIN_EL traces."""
    xs = [t for t in traces if t is not None]
    if len(xs) < MIN_EL:
        return None
    b = np.mean(np.vstack(xs), axis=0)
    return normalise(b)


def pick_events(rows, bins=BINS, per_side=PER_SIDE, turn=TURN):
    """rows: dicts with 'time' (ISO), 't' (decimal year) and 'mag'. Per bin [lo, hi) and per
    side of the turn: events in time order, at most per_side, evenly spaced through the list
    (deterministic)."""
    per_side = max(2, int(per_side))
    out = []
    for lo, hi in bins:
        for side in ("before", "after"):
            sel = sorted((r for r in rows if isf(r.get("mag")) and lo <= r["mag"] < hi and isf(r.get("t"))
                          and ((r["t"] < turn) == (side == "before"))), key=lambda r: r["time"])
            if len(sel) > per_side:
                keep = sorted({int(round(i * (len(sel) - 1) / (per_side - 1))) for i in range(per_side)})
                sel = [sel[i] for i in keep]
            for r in sel:
                d = dict(r)
                d["bin"] = "%.1f-%.1f" % (lo, hi)
                d["side"] = side
                out.append(d)
    return out


def bin_counts(rows, bins=BINS, turn=TURN):
    out = []
    for lo, hi in bins:
        sel = [r for r in rows if isf(r.get("mag")) and lo <= r["mag"] < hi]
        after = sum(1 for r in sel if isf(r.get("t")) and r["t"] >= turn)
        out.append(dict(bin="%.1f-%.1f" % (lo, hi), n=len(sel), before_turn=len(sel) - after,
                        after_turn=after))
    return out


def close_pairs(rows, max_km=PAIR_KM, max_dm=PAIR_DM, max_dz=PAIR_DZ, min_mag=PAIR_MIN_MAG, turn=TURN):
    """Catalogue pairs (both M >= min_mag) within max_km epicentral distance, |dM| <= max_dm and
    |dz| <= max_dz. A loose upper-bound proxy for repeaters (ComCat locations are uncertain by
    5-20 km). Returns counts overall and straddling the turn, and for M >= 5.5."""
    ev = [r for r in rows if isf(r.get("mag")) and r["mag"] >= min_mag and isf(r.get("lat"))
          and isf(r.get("lon")) and isf(r.get("t"))]
    res = dict(n_events=len(ev), pairs=0, straddling=0, straddling_m55=0, events_in_straddling=0)
    used = set()
    deg = max_km / 111.0
    for i in range(len(ev)):
        a = ev[i]
        for j in range(i + 1, len(ev)):
            b = ev[j]
            if abs(a["lat"] - b["lat"]) > deg or abs(a["mag"] - b["mag"]) > max_dm:
                continue
            za, zb = a.get("depth"), b.get("depth")
            if isf(za) and isf(zb) and abs(za - zb) > max_dz:
                continue
            if km_between(a["lat"], a["lon"], b["lat"], b["lon"]) > max_km:
                continue
            res["pairs"] += 1
            if (a["t"] < turn) != (b["t"] < turn):
                res["straddling"] += 1
                used.add(i)
                used.add(j)
                if min(a["mag"], b["mag"]) >= 5.5:
                    res["straddling_m55"] += 1
    res["events_in_straddling"] = len(used)
    return res


def choose_chan(chans, iso, cha_pref=CHA_PREF, loc_pref=LOC_PREF):
    """chans: [(cha, loc, start_iso, end_iso, lat, lon)]. The preferred covering
    (cha, loc, lat, lon), or None."""
    s = str(iso)[:19]
    cover = [c for c in chans if str(c[2])[:19] <= s < (str(c[3])[:19] if c[3] else "9999")]
    if not cover:
        return None

    def rank(c):
        ci = cha_pref.index(c[0]) if c[0] in cha_pref else len(cha_pref)
        li = loc_pref.index(c[1]) if c[1] in loc_pref else len(loc_pref)
        return (ci, li, c[1])
    c = sorted(cover, key=rank)[0]
    return c[0], c[1], c[4], c[5]


def array_verdict(rows, key, bins=BINS, req=SNR_REQ, min_beams=MIN_BEAMS):
    """rows: per event-array rows of one array. Lowest bin with >= min_beams beams whose median
    row[key] >= req."""
    table = []
    first = None
    for lo, hi in bins:
        name = "%.1f-%.1f" % (lo, hi)
        v = [r[key] for r in rows if r.get("bin") == name and isf(r.get(key))]
        med = median(v)
        ok = len(v) >= min_beams and isf(med) and med >= req
        table.append(dict(bin=name, n=len(v), median=med,
                          frac_ge=(sum(1 for x in v if x >= req) / len(v)) if v else NAN, passes=ok))
        if ok and first is None:
            first = name
    return dict(table=table, first_bin=first,
                verdict=("testable-from-M%s" % first.split("-")[0]) if first else "none")


def availability(rows):
    """Fraction of in-range picked events with a beam, before and after the turn."""
    out = {}
    for side in ("before", "after"):
        sel = [r for r in rows if r.get("side") == side and r.get("status") != "out-of-range"]
        ok = sum(1 for r in sel if isf(r.get("beam_sig_primary")))
        out[side] = dict(n=len(sel), with_beam=ok, frac=(ok / len(sel)) if sel else NAN)
    return out


def beam_gain_frac(rows, min_el=5):
    sel = [r for r in rows if isf(r.get("beam_sig_primary")) and isf(r.get("el_sig_primary_median"))
           and (r.get("n_el_ok") or 0) >= min_el]
    if not sel:
        return NAN, 0
    return sum(1 for r in sel if r["beam_sig_primary"] > r["el_sig_primary_median"]) / len(sel), len(sel)


def decide(verdicts, avail, pairs, min_avail=0.6):
    """Pre-registered decision (plan section 'Decision rule')."""
    order = {"testable-from-M4.5": 0, "testable-from-M5.0": 1, "testable-from-M5.5": 2,
             "testable-from-M6.0": 3, "none": 4}
    good = []
    for name, v in verdicts.items():
        a = avail.get(name, {})
        fb, fa = a.get("before", {}).get("frac", NAN), a.get("after", {}).get("frac", NAN)
        if isf(fb) and isf(fa) and fb >= min_avail and fa >= min_avail:
            good.append((order.get(v["verdict"], 4), name))
    if not good:
        return dict(best_array=None, decision="close-chile-path", reason="no array with data in both epochs")
    rank, name = sorted(good)[0]
    if rank <= 1:
        return dict(best_array=name, decision="build-catalogue", reason=verdicts[name]["verdict"])
    if rank == 2:
        if pairs.get("straddling_m55", 0) >= 30:
            return dict(best_array=name, decision="build-catalogue-m5.5",
                        reason="testable-from-M5.5 and >= 30 M5.5 close pairs straddle the turn")
        return dict(best_array=name, decision="close-chile-path", reason="M5.5 floor, too few M5.5 pairs")
    return dict(best_array=name, decision="close-chile-path", reason=verdicts[name]["verdict"])


# ----------------------------------------------------------------------------- I/O
def fetch_text(url, tries=5, wait=15, timeout=120):
    err = ""
    for k in range(tries):
        try:
            with urllib.request.urlopen(url, timeout=timeout) as r:
                return r.read().decode("utf-8", "replace")
        except Exception as ex:  # noqa: BLE001
            err = repr(ex)[:160]
            time.sleep(wait * (k + 1))
    raise RuntimeError("fetch failed %s: %s" % (url, err))


def plan_public(wait_min):
    t0 = time.time()
    while True:
        try:
            with urllib.request.urlopen(PLAN_URL, timeout=60) as r:
                if r.status == 200 and b"chilepilot" in r.read():
                    return True
        except Exception:  # noqa: BLE001
            pass
        if time.time() - t0 > wait_min * 60:
            return False
        time.sleep(120)


def load_catalogue():
    q = ("?format=csv&starttime=%s&endtime=%s&minlatitude=%.1f&maxlatitude=%.1f&minlongitude=%.1f"
         "&maxlongitude=%.1f&maxdepth=%.0f&minmagnitude=%.1f&orderby=time-asc"
         % (START, END, BOX["minlatitude"], BOX["maxlatitude"], BOX["minlongitude"], BOX["maxlongitude"],
            MAX_DEPTH, MIN_MAG))
    txt = fetch_text(COMCAT + q)
    rows = []
    for r in csv.DictReader(io.StringIO(txt)):
        rows.append(dict(id=r["id"], time=r["time"], lat=fnum(r["latitude"]), lon=fnum(r["longitude"]),
                         depth=fnum(r["depth"]), mag=fnum(r["mag"]), magType=r.get("magType", ""),
                         t=decimal_year(r["time"])))
    return rows, txt


def get_inventory(client, net, sta, data_dir):
    from obspy import read_inventory
    d = os.path.join(data_dir, "resp")
    os.makedirs(d, exist_ok=True)
    p = os.path.join(d, "%s.%s.xml" % (net, sta))
    if os.path.exists(p):
        try:
            return read_inventory(p)
        except Exception:  # noqa: BLE001
            os.remove(p)
    for k in range(3):
        try:
            inv = client.get_stations(network=net, station=sta, channel=",".join(CHA_PREF), level="response")
            inv.write(p, format="STATIONXML")
            return inv
        except Exception as ex:  # noqa: BLE001
            log("inventory failed", net, sta, repr(ex)[:120])
            time.sleep(10 * (k + 1))
    return None


def chans_of(inv):
    out = []
    for n in inv.networks:
        for s in n.stations:
            for ch in s.channels:
                if ch.code not in CHA_PREF:
                    continue
                lat = ch.latitude if ch.latitude is not None else s.latitude
                lon = ch.longitude if ch.longitude is not None else s.longitude
                out.append((ch.code, ch.location_code or "", str(ch.start_date)[:19],
                            str(ch.end_date)[:19] if ch.end_date else "9999", lat, lon))
    return out


def fetch_array(client, aname, lines, ev, data_dir):
    """lines: [(net, sta, loc, cha, t1, t2)]. One bulk request per event and array, cached."""
    from obspy import read
    from obspy.clients.fdsn.header import FDSNNoDataException

    d = os.path.join(data_dir, "raw", aname)
    os.makedirs(d, exist_ok=True)
    p = os.path.join(d, ev["id"] + ".mseed")
    if os.path.exists(p):
        try:
            return read(p), "cached"
        except Exception:  # noqa: BLE001
            os.remove(p)
    if os.path.exists(p + ".nodata"):
        return None, "nodata"
    bulk = [(n, s, (l if l else "--"), c, t1, t2) for n, s, l, c, t1, t2 in lines]
    err = ""
    for k in range(3):
        try:
            st = client.get_waveforms_bulk(bulk)
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
            time.sleep(10 * (k + 1))
    return None, "error:" + err


def process(st, inv, ev, tdf, band):
    """Velocity (response removed when possible), causal 4-pole band-pass, FS Hz, the element's
    ak135 PKIKP at PRE, normalised by the NOISE-window rms."""
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
    if tr.stats.starttime > o + tdf - PRE - 5 or tr.stats.endtime < o + tdf + POST + 5:
        return None, "window-not-covered"
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
    xn = normalise(x)
    if xn is None:
        return None, "dead"
    return (xn, resp_ok), "ok"


def jsonable(o):
    if isinstance(o, dict):
        return {str(k): jsonable(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
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
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in keys})


def plot(path, rows):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(13, 5.5), sharey=True)
    for ax, (bname, band) in zip(axes, BANDS):
        for aname, _ in ARRAYS:
            pts = [(r["mag"], r["beam_sig_" + bname]) for r in rows
                   if r["array"] == aname and isf(r.get("beam_sig_" + bname)) and isf(r.get("mag"))]
            if pts:
                ax.scatter([p[0] for p in pts], [p[1] for p in pts], s=18, label=aname + " beam", alpha=0.8)
        nul = [(r["mag"], r["beam_null_" + bname]) for r in rows if isf(r.get("beam_null_" + bname))]
        if nul:
            ax.scatter([p[0] for p in nul], [p[1] for p in nul], s=6, c="0.6", marker="x",
                       label="noise-only control")
        ax.axhline(SNR_REQ, color="k", lw=1, ls="--")
        ax.set_yscale("log")
        ax.set_xlabel("ComCat magnitude")
        ax.set_title("Beam rms SNR, PKP window, %s band %.1f-%.1f Hz" % (bname, band[0], band[1]))
    axes[0].set_ylabel("rms SNR (window / noise)")
    axes[1].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=110)
    plt.close(fig)


# ----------------------------------------------------------------------------- main
def run_event_array(client, model, ev, aname, members, invs, data_dir):
    row = dict(id=ev["id"], time=ev["time"], t=ev["t"], mag=ev["mag"], magType=ev["magType"],
               depth=ev["depth"], lat=ev["lat"], lon=ev["lon"], bin=ev["bin"], side=ev["side"], array=aname)
    from obspy import UTCDateTime

    els = []
    for net, sta in members:
        inv, chans = invs.get((net, sta), (None, []))
        if inv is None:
            continue
        c = choose_chan(chans, ev["time"])
        if c is None:
            continue
        els.append(dict(net=net, sta=sta, cha=c[0], loc=c[1], lat=c[2], lon=c[3], inv=inv))
    row["n_el_epoch"] = len(els)
    if len(els) < MIN_EL:
        row["status"] = "few-epochs"
        return row
    clat = sum(e["lat"] for e in els) / len(els)
    clon = sum(e["lon"] for e in els) / len(els)
    row["dist"] = gc_deg(ev["lat"], ev["lon"], clat, clon)
    if not (DIST_RANGE[0] <= row["dist"] <= DIST_RANGE[1]):
        row["status"] = "out-of-range"
        return row
    z = max(ev["depth"], 0.0) if isf(ev["depth"]) else 10.0
    for e in els:
        dist = gc_deg(ev["lat"], ev["lon"], e["lat"], e["lon"])
        arr = model.get_travel_times(source_depth_in_km=z, distance_in_degree=dist, phase_list=["PKIKP"])
        e["tdf"] = min((a.time for a in arr if a.name == "PKIKP"), default=NAN)
    els = [e for e in els if isf(e["tdf"])]
    if len(els) < MIN_EL:
        row["status"] = "no-phase"
        return row
    arr = model.get_travel_times(source_depth_in_km=z, distance_in_degree=row["dist"],
                                 phase_list=["PKIKP", "PKiKP", "PKP"])
    tdf_c = min((a.time for a in arr if a.name == "PKIKP"), default=NAN)
    tbc = [a.time for a in arr if a.name == "PKP"]
    row["bc_minus_df_pred"] = (min(tbc) - tdf_c) if (tbc and isf(tdf_c)) else NAN
    row["ab_minus_df_pred"] = (max(tbc) - tdf_c) if (tbc and isf(tdf_c)) else NAN
    o = UTCDateTime(ev["time"])
    t1 = o + min(e["tdf"] for e in els) - PRE - MARGIN
    t2 = o + max(e["tdf"] for e in els) + POST + MARGIN
    st, how = fetch_array(client, aname, [(e["net"], e["sta"], e["loc"], e["cha"], t1, t2) for e in els],
                          ev, data_dir)
    row["fetch"] = how
    if st is None:
        row["status"] = how
        return row
    whys = {}
    for bname, band in BANDS:
        xs, el_sig = [], []
        for e in els:
            sub = st.select(network=e["net"], station=e["sta"], channel=e["cha"])
            sub = type(sub)([tr for tr in sub if (tr.stats.location or "") == e["loc"]])
            if len(sub) == 0:
                whys["missing"] = whys.get("missing", 0) + 1
                continue
            res, why = process(sub, e["inv"], ev, e["tdf"], band)
            if res is None:
                whys[why[:20]] = whys.get(why[:20], 0) + 1
                continue
            xs.append(res[0])
            el_sig.append(window_snrs(res[0])[0])
        row["n_el_ok_" + bname] = len(xs)
        row["el_sig_%s_median" % bname] = median(el_sig)
        b = beam(xs)
        if b is not None:
            s, c, nl = window_snrs(b)
            row["beam_sig_" + bname], row["beam_coda_" + bname], row["beam_null_" + bname] = s, c, nl
    row["n_el_ok"] = row.get("n_el_ok_primary", 0)
    row["drop_reasons"] = ";".join("%s:%d" % kv for kv in sorted(whys.items()))
    row["status"] = "ok" if isf(row.get("beam_sig_primary")) else "few-elements"
    return row


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", default="/data")
    ap.add_argument("--out-dir", default="/out")
    ap.add_argument("--wait-min", type=float, default=0.0)
    args = ap.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)
    os.makedirs(os.path.join(args.out_dir, "plots"), exist_ok=True)
    os.makedirs(args.data_dir, exist_ok=True)
    summary = dict(version=VERSION, started=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
                   constants=dict(box=BOX, max_depth=MAX_DEPTH, start=START, end=END, bins=BINS,
                                  per_side=PER_SIDE, turn=TURN, dist_range=DIST_RANGE,
                                  arrays={a: [n + "." + s for n, s in m] for a, m in ARRAYS},
                                  bands=dict(BANDS), sig=SIG, coda=CODA, noise=NOISE, null=NULL,
                                  snr_req=SNR_REQ, min_el=MIN_EL, min_beams=MIN_BEAMS,
                                  pairs=dict(km=PAIR_KM, dm=PAIR_DM, dz=PAIR_DZ, min_mag=PAIR_MIN_MAG)),
                   errors={})
    summary["plan_public_at_start"] = plan_public(args.wait_min)
    log("plan public:", summary["plan_public_at_start"])

    from obspy.clients.fdsn import Client
    from obspy.taup import TauPyModel

    cat, raw = load_catalogue()
    with open(os.path.join(args.out_dir, "chile_catalogue.csv"), "w", newline="") as f:
        f.write(raw)
    summary["catalogue_n"] = len(cat)
    summary["catalogue_bins"] = bin_counts(cat)
    summary["close_pairs"] = close_pairs(cat)
    log("catalogue", len(cat), "close pairs", json.dumps(summary["close_pairs"]))
    events = pick_events(cat)
    summary["n_events"] = len(events)

    client = None
    for base in ("https://service.earthscope.org", "IRIS"):
        try:
            client = Client(base, timeout=180)
            summary["client"] = base
            break
        except Exception as ex:  # noqa: BLE001
            summary["errors"]["client-" + base] = repr(ex)[:160]
    if client is None:
        raise SystemExit("no FDSN client")
    model = TauPyModel("ak135")

    invs = {}
    inv_status = {}
    for aname, members in ARRAYS:
        for net, sta in members:
            inv = get_inventory(client, net, sta, args.data_dir)
            invs[(net, sta)] = (inv, chans_of(inv) if inv is not None else [])
            inv_status[net + "." + sta] = len(invs[(net, sta)][1])
    summary["inventory_channel_epochs"] = inv_status

    rows = []
    for k, ev in enumerate(events):
        for aname, members in ARRAYS:
            try:
                rows.append(run_event_array(client, model, ev, aname, members, invs, args.data_dir))
            except Exception as ex:  # noqa: BLE001
                rows.append(dict(id=ev["id"], time=ev["time"], mag=ev["mag"], bin=ev["bin"], side=ev["side"],
                                 array=aname, status="error:" + repr(ex)[:100]))
                summary["errors"]["%s %s" % (ev["id"], aname)] = traceback.format_exc()[-400:]
        log("event", k + 1, "/", len(events), ev["id"], ev["mag"], [r.get("status") for r in rows[-len(ARRAYS):]])
        if k % 4 == 0:
            write_csv(os.path.join(args.out_dir, "chile_rows.csv"), rows)

    write_csv(os.path.join(args.out_dir, "chile_rows.csv"), rows)
    st_counts = {}
    for r in rows:
        s = "%s %s" % (r.get("array"), str(r.get("status", ""))[:40])
        st_counts[s] = st_counts.get(s, 0) + 1
    summary["status_counts"] = st_counts

    verdicts, avail, per_array = {}, {}, {}
    for aname, _ in ARRAYS:
        ar = [r for r in rows if r.get("array") == aname]
        verdicts[aname] = array_verdict(ar, "beam_sig_primary")
        avail[aname] = availability(ar)
        g, ng = beam_gain_frac(ar)
        per_array[aname] = dict(
            verdict_primary=verdicts[aname], verdict_alt=array_verdict(ar, "beam_sig_alt"),
            verdict_coda_primary=array_verdict(ar, "beam_coda_primary"),
            availability=avail[aname], beam_gain_frac=g, beam_gain_n=ng,
            dist_min=min((r["dist"] for r in ar if isf(r.get("dist"))), default=NAN),
            dist_max=max((r["dist"] for r in ar if isf(r.get("dist"))), default=NAN),
            null_median=median(r.get("beam_null_primary") for r in ar))
    summary["arrays"] = per_array
    summary["decision"] = decide(verdicts, avail, summary["close_pairs"])
    allnull = median(r.get("beam_null_primary") for r in rows)
    order = ["testable-from-M4.5", "testable-from-M5.0", "testable-from-M5.5", "testable-from-M6.0", "none"]
    best = min((order.index(v["verdict"]) for v in verdicts.values()), default=4)
    gains = [(per_array[a]["beam_gain_frac"], per_array[a]["beam_gain_n"]) for a, _ in ARRAYS]
    gn = sum(n for _, n in gains)
    summary["forecasts"] = dict(
        F1_some_array_avail_ge_0p7_both=any(
            isf(a["before"]["frac"]) and isf(a["after"]["frac"]) and a["before"]["frac"] >= 0.7
            and a["after"]["frac"] >= 0.7 for a in avail.values()),
        F2_best_verdict=order[best],
        F3_beam_gain_frac=(sum(g * n for g, n in gains if isf(g)) / gn) if gn else NAN,
        F4_null_median=allnull,
        F5_straddling_pairs=summary["close_pairs"]["straddling"])
    summary["finished"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    with open(os.path.join(args.out_dir, "chile_summary.json"), "w") as f:
        json.dump(jsonable(summary), f, indent=1)
    try:
        plot(os.path.join(args.out_dir, "plots", "chile_snr_vs_mag.png"), rows)
    except Exception as ex:  # noqa: BLE001
        log("plot failed", repr(ex)[:160])
    log("done", json.dumps(jsonable(summary["forecasts"])), json.dumps(jsonable(summary["decision"])))


if __name__ == "__main__":
    sys.exit(main())

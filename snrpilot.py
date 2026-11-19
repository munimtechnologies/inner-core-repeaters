#!/usr/bin/env python3
# Copyright 2026 MV2, Munim, Inc.
# Licensed under the Apache License, Version 2.0 (Apache-2.0). See LICENSE.
# Made with AI agents working in MV2's name.
"""snrpilot.py - how loud is PKIKP (DF) from South Sandwich Islands earthquakes at
near-antipodal stations (155-175 deg), as a function of magnitude?

Design pilot for a large-event multiplet catalogue (plan:
findings/2026-11-19-snr-pilot-plan.md). Pre-registration 16 found the existing SSI
multiplets (mostly mb 4.3-5.2) untestable at these stations: DF sat at noise. This pilot
measures, event by event and without any family alignment, the rms signal-to-noise ratio of
DF and of PKPab in fixed windows around the ak135 predictions, plus a noise-only control
window, for a time-spread sample of events in four magnitude bins. It reports the magnitude
from which DF is reliably above the pre-registration-16 gate (rms SNR 3).

Usage: python snrpilot.py --data-dir /data --out-dir /out [--wait-min 120]
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
PLAN_URL = REPO + "main/findings/2026-11-19-snr-pilot-plan.md"
COMCAT = "https://earthquake.usgs.gov/fdsnws/event/1/query"
CENTROID = (-59.0365, -25.6939)          # results/m0324-v1 centroid (same as nea.py)
RADIUS_KM = 300.0
START, END = "1990-01-01", "2026-11-19"
MIN_MAG = 4.5
BINS = ((4.5, 5.0), (5.0, 5.5), (5.5, 6.0), (6.0, 9.9))
PER_BIN = 15
TURN = 2008.5
# the ten stations pre-registration 16 chose (results/nea-v3/nea_summary.json)
STATIONS = (("IC", "HIA"), ("IC", "MDJ"), ("II", "ERM"), ("IU", "BILL"), ("IU", "MAJO"),
            ("IU", "PET"), ("IU", "TIXI"), ("IU", "YAK"), ("IU", "YSS"), ("IU", "ADK"))
LOC_PREF = ("00", "", "10", "20", "60")

FS = 20.0
PRE, POST = 60.0, 90.0                   # s kept before / after ak135 PKIKP
N = int(round((PRE + POST) * FS))
BANDS = (("primary", (0.6, 1.6)), ("alt", (1.0, 2.0)))   # the two nea.py bands
PREFILT = (0.1, 0.2, 7.0, 9.0)
NOISE = (-45.0, -10.0)                   # normalisation window, s rel. ak135 PKIKP (as nea.py)
SIG_DF = (-3.0, 7.0)                     # DF window, s rel. ak135 PKIKP
SIG_AB = (-3.0, 7.0)                     # AB window, s rel. ak135 PKPab
NULL = (-58.0, -48.0)                    # noise-only control, outside NOISE
SNR_REQ = 3.0                            # pre-registration 16 DF gate
MIN_N = 3                                # rows needed for a station-bin median
MIN_STA_PASS = 3                         # stations whose bin median must reach SNR_REQ


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


def window_snrs(x, ab_off, pre=PRE, fs=FS):
    """x: trace normalised by its NOISE-window rms. Returns (snr_df, snr_ab, snr_null)."""
    df = rms(x, idx(SIG_DF[0], pre, fs), idx(SIG_DF[1], pre, fs))
    nl = rms(x, idx(NULL[0], pre, fs), idx(NULL[1], pre, fs))
    ab = NAN
    if isf(ab_off):
        ab = rms(x, idx(ab_off + SIG_AB[0], pre, fs), idx(ab_off + SIG_AB[1], pre, fs))
    return df, ab, nl


def pick_events(rows, bins=BINS, per_bin=PER_BIN):
    """rows: dicts with 'time' (ISO) and 'mag'. Per bin [lo, hi): events in time order,
    at most per_bin of them, evenly spaced through the list (deterministic)."""
    per_bin = max(2, int(per_bin))
    out = []
    for lo, hi in bins:
        sel = sorted((r for r in rows if isf(r.get("mag")) and lo <= r["mag"] < hi),
                     key=lambda r: r["time"])
        if len(sel) > per_bin:
            keep = sorted({int(round(i * (len(sel) - 1) / (per_bin - 1))) for i in range(per_bin)})
            sel = [sel[i] for i in keep]
        for r in sel:
            d = dict(r)
            d["bin"] = "%.1f-%.1f" % (lo, hi)
            out.append(d)
    return out


def bin_counts(rows, bins=BINS, turn=TURN):
    """Catalogue events per bin, before and after the 2008.5 turn."""
    out = []
    for lo, hi in bins:
        sel = [r for r in rows if isf(r.get("mag")) and lo <= r["mag"] < hi]
        after = sum(1 for r in sel if isf(r.get("t")) and r["t"] >= turn)
        out.append(dict(bin="%.1f-%.1f" % (lo, hi), n=len(sel), before_turn=len(sel) - after,
                        after_turn=after))
    return out


def fe_slope(rows, key):
    """Slope of log10(row[key]) on row['mag'] with a fixed effect per station."""
    by = {}
    for r in rows:
        v, m = r.get(key), r.get("mag")
        if isf(v) and v > 0 and isf(m):
            by.setdefault(r["station"], []).append((m, math.log10(v)))
    num = den = 0.0
    n = 0
    for pts in by.values():
        if len(pts) < 3:
            continue
        mm = sum(p[0] for p in pts) / len(pts)
        ss = sum(p[1] for p in pts) / len(pts)
        for m, s in pts:
            num += (m - mm) * (s - ss)
            den += (m - mm) ** 2
            n += 1
    return (num / den if den > 0 else NAN), n


def station_fit(rows, key, req=SNR_REQ):
    """Per station: OLS log10(snr) = a + b*mag and the magnitude where the fit reaches req."""
    by = {}
    for r in rows:
        v, m = r.get(key), r.get("mag")
        if isf(v) and v > 0 and isf(m):
            by.setdefault(r["station"], []).append((m, math.log10(v)))
    out = {}
    for sta, pts in sorted(by.items()):
        res = dict(n=len(pts), a=NAN, b=NAN, m_req=NAN)
        if len(pts) >= 3:
            mm = sum(p[0] for p in pts) / len(pts)
            ss = sum(p[1] for p in pts) / len(pts)
            den = sum((p[0] - mm) ** 2 for p in pts)
            if den > 0:
                b = sum((p[0] - mm) * (p[1] - ss) for p in pts) / den
                a = ss - b * mm
                res.update(a=a, b=b)
                if b > 0:
                    res["m_req"] = (math.log10(req) - a) / b
        out[sta] = res
    return out


def bin_table(rows, key, bins=BINS, req=SNR_REQ, min_n=MIN_N):
    """Per bin and station: n, median and fraction >= req of row[key]."""
    table = []
    stations = sorted({r["station"] for r in rows})
    for lo, hi in bins:
        name = "%.1f-%.1f" % (lo, hi)
        for sta in stations:
            v = [r[key] for r in rows if r["station"] == sta and r.get("bin") == name and isf(r.get(key))]
            table.append(dict(bin=name, station=sta, n=len(v), median=median(v) if v else NAN,
                              frac_ge=(sum(1 for x in v if x >= req) / len(v)) if v else NAN,
                              usable=len(v) >= min_n))
    return table


def verdict(table, bins=BINS, req=SNR_REQ, min_sta=MIN_STA_PASS):
    """Lowest bin in which at least min_sta stations (with >= MIN_N rows) have median >= req."""
    passing = {}
    for lo, hi in bins:
        name = "%.1f-%.1f" % (lo, hi)
        good = sorted(t["station"] for t in table
                      if t["bin"] == name and t["usable"] and isf(t["median"]) and t["median"] >= req)
        passing[name] = good
    first = next(("%.1f-%.1f" % b for b in bins if len(passing["%.1f-%.1f" % b]) >= min_sta), None)
    return dict(first_bin=first, verdict=("testable-from-M%s" % first.split("-")[0]) if first else "none",
                stations_passing=passing)


# ----------------------------------------------------------------------------- I/O
def fetch_text(url, tries=5, wait=15, timeout=90):
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
                if r.status == 200 and b"snrpilot" in r.read():
                    return True
        except Exception:  # noqa: BLE001
            pass
        if time.time() - t0 > wait_min * 60:
            return False
        time.sleep(120)


def load_catalogue():
    q = ("?format=csv&starttime=%s&endtime=%s&latitude=%.4f&longitude=%.4f&maxradiuskm=%.0f"
         "&minmagnitude=%.1f&orderby=time-asc" % (START, END, CENTROID[0], CENTROID[1], RADIUS_KM, MIN_MAG))
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
            inv = client.get_stations(network=net, station=sta, channel="BHZ", level="response")
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
                if ch.code != "BHZ":
                    continue
                lat = ch.latitude if ch.latitude is not None else s.latitude
                lon = ch.longitude if ch.longitude is not None else s.longitude
                out.append((ch.location_code or "", str(ch.start_date)[:19],
                            str(ch.end_date)[:19] if ch.end_date else "9999", lat, lon))
    return out


def choose_loc(chans, iso, pref=LOC_PREF):
    """chans: [(loc, start_iso, end_iso, lat, lon)]. The preferred covering location code
    and its (lat, lon), or (None, None)."""
    s = str(iso)[:19]
    cover = {}
    for loc, a, b, lat, lon in chans:
        if str(a)[:19] <= s < (str(b)[:19] if b else "9999"):
            cover.setdefault(loc, (lat, lon))
    for loc in pref:
        if loc in cover:
            return loc, cover[loc]
    if cover:
        k = sorted(cover)[0]
        return k, cover[k]
    return None, None


def fetch_wave(client, net, sta, loc, ev, tdf, data_dir):
    from obspy import UTCDateTime, read
    from obspy.clients.fdsn.header import FDSNNoDataException

    d = os.path.join(data_dir, "raw", "%s.%s" % (net, sta))
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
            st = client.get_waveforms(net, sta, loc or "--", "BHZ", t1, t2)
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
    """As nea.py process(): velocity (response removed when possible), causal 4-pole
    band-pass, FS Hz, ak135 PKIKP at PRE, normalised by the NOISE-window rms."""
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
    nz = rms(x, idx(NOISE[0]), idx(NOISE[1]))
    if not (isf(nz) and nz > 0):
        return None, "dead"
    return (x / nz, resp_ok), "ok"


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
        key = "snr_df_" + bname
        for sta in sorted({r["station"] for r in rows}):
            pts = [(r["mag"], r[key]) for r in rows if r["station"] == sta and isf(r.get(key)) and isf(r.get("mag"))]
            if pts:
                ax.scatter([p[0] for p in pts], [p[1] for p in pts], s=14, label=sta, alpha=0.8)
        nul = [r["mag"] for r in rows if isf(r.get("snr_null_" + bname))]
        nv = [r["snr_null_" + bname] for r in rows if isf(r.get("snr_null_" + bname))]
        if nul:
            ax.scatter(nul, nv, s=6, c="0.6", marker="x", label="noise-only control")
        ax.axhline(SNR_REQ, color="k", lw=1, ls="--")
        ax.set_yscale("log")
        ax.set_xlabel("ComCat magnitude")
        ax.set_title("DF rms SNR, %s band %.1f-%.1f Hz" % (bname, band[0], band[1]))
    axes[0].set_ylabel("rms SNR (window / noise)")
    axes[1].legend(fontsize=7, ncol=2)
    fig.tight_layout()
    fig.savefig(path, dpi=110)
    plt.close(fig)


# ----------------------------------------------------------------------------- main
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
                   constants=dict(bins=BINS, per_bin=PER_BIN, stations=[n + "." + s for n, s in STATIONS],
                                  bands=dict(BANDS), sig_df=SIG_DF, sig_ab=SIG_AB, noise=NOISE, null=NULL,
                                  snr_req=SNR_REQ, min_n=MIN_N, min_sta_pass=MIN_STA_PASS,
                                  radius_km=RADIUS_KM, centroid=CENTROID, start=START, end=END),
                   errors={})
    summary["plan_public_at_start"] = plan_public(args.wait_min)
    log("plan public:", summary["plan_public_at_start"])

    from obspy import UTCDateTime  # noqa: F401
    from obspy.clients.fdsn import Client
    from obspy.taup import TauPyModel

    cat, raw = load_catalogue()
    with open(os.path.join(args.out_dir, "snr_catalogue.csv"), "w", newline="") as f:
        f.write(raw)
    summary["catalogue_n"] = len(cat)
    summary["catalogue_bins"] = bin_counts(cat)
    events = pick_events(cat)
    summary["n_events"] = len(events)
    log("catalogue", len(cat), "picked", len(events))

    client = None
    for base in ("https://service.earthscope.org", "IRIS"):
        try:
            client = Client(base, timeout=120)
            summary["client"] = base
            break
        except Exception as ex:  # noqa: BLE001
            summary["errors"]["client-" + base] = repr(ex)[:160]
    if client is None:
        raise SystemExit("no FDSN client")
    model = TauPyModel("ak135")

    invs = {}
    for net, sta in STATIONS:
        inv = get_inventory(client, net, sta, args.data_dir)
        invs[(net, sta)] = (inv, chans_of(inv) if inv is not None else [])

    rows = []
    for k, ev in enumerate(events):
        for net, sta in STATIONS:
            row = dict(id=ev["id"], time=ev["time"], t=ev["t"], mag=ev["mag"], magType=ev["magType"],
                       depth=ev["depth"], bin=ev["bin"], station=net + "." + sta)
            rows.append(row)
            try:
                inv, chans = invs[(net, sta)]
                if inv is None:
                    row["status"] = "no-inventory"
                    continue
                loc, ll = choose_loc(chans, ev["time"])
                if loc is None:
                    row["status"] = "no-epoch"
                    continue
                dist = gc_deg(ev["lat"], ev["lon"], ll[0], ll[1])
                row.update(loc=loc, dist=dist)
                arr = model.get_travel_times(source_depth_in_km=max(ev["depth"], 0.0) if isf(ev["depth"]) else 10.0,
                                             distance_in_degree=dist, phase_list=["PKIKP", "PKP"])
                tdf = min((a.time for a in arr if a.name == "PKIKP"), default=NAN)
                tab = max((a.time for a in arr if a.name == "PKP"), default=NAN)
                if not isf(tdf):
                    row["status"] = "no-phase"
                    continue
                ab_off = tab - tdf if isf(tab) else NAN
                if isf(ab_off) and ab_off + SIG_AB[1] > POST - 1:
                    ab_off = NAN
                row.update(tdf=tdf, ab_minus_df_pred=ab_off)
                st, how = fetch_wave(client, net, sta, loc, ev, tdf, args.data_dir)
                row["fetch"] = how
                if st is None:
                    row["status"] = how
                    continue
                status = []
                for bname, band in BANDS:
                    res, why = process(st, inv, ev, tdf, band)
                    if res is None:
                        status.append(bname + ":" + why)
                        continue
                    x, resp_ok = res
                    df, ab, nl = window_snrs(x, ab_off)
                    row["snr_df_" + bname] = df
                    row["snr_ab_" + bname] = ab
                    row["snr_null_" + bname] = nl
                    row["resp_ok"] = resp_ok
                    status.append(bname + ":ok")
                row["status"] = "ok" if all(s.endswith(":ok") for s in status) else ";".join(status)
            except Exception as ex:  # noqa: BLE001
                row["status"] = "error:" + repr(ex)[:100]
                summary["errors"]["%s %s.%s" % (ev["id"], net, sta)] = traceback.format_exc()[-400:]
        if k % 5 == 0:
            log("event", k + 1, "/", len(events), ev["id"], ev["mag"])
            write_csv(os.path.join(args.out_dir, "snr_rows.csv"), rows)

    write_csv(os.path.join(args.out_dir, "snr_rows.csv"), rows)
    st_counts = {}
    for r in rows:
        s = str(r.get("status", ""))[:48]
        st_counts[s] = st_counts.get(s, 0) + 1
    summary["status_counts"] = st_counts

    res = {}
    for bname, _ in BANDS:
        kdf, kab, knl = "snr_df_" + bname, "snr_ab_" + bname, "snr_null_" + bname
        table = bin_table(rows, kdf)
        slope_df, n_df = fe_slope(rows, kdf)
        slope_ab, n_ab = fe_slope(rows, kab)
        res[bname] = dict(
            n_rows=sum(1 for r in rows if isf(r.get(kdf))),
            null_median=median(r.get(knl) for r in rows),
            null_frac_ge_req=(lambda v: sum(1 for x in v if x >= SNR_REQ) / len(v) if v else NAN)(
                [r[knl] for r in rows if isf(r.get(knl))]),
            df_slope_fe=slope_df, df_slope_n=n_df, ab_slope_fe=slope_ab, ab_slope_n=n_ab,
            df_station_fit=station_fit(rows, kdf), ab_station_fit=station_fit(rows, kab),
            df_bin_table=table, ab_bin_table=bin_table(rows, kab),
            df_verdict=verdict(table))
    summary["results"] = res
    p = res["primary"]
    summary["forecasts"] = dict(
        F1_lowest_bin_no_station_passes=not p["df_verdict"]["stations_passing"].get("4.5-5.0"),
        F2_verdict=p["df_verdict"]["verdict"],
        F3_df_slope_fe=p["df_slope_fe"],
        F4_null_median=p["null_median"],
        F5_ab_beats_df=median(r.get("snr_ab_primary") for r in rows) > median(r.get("snr_df_primary") for r in rows))
    summary["finished"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    with open(os.path.join(args.out_dir, "snr_summary.json"), "w") as f:
        json.dump(jsonable(summary), f, indent=1)
    try:
        plot(os.path.join(args.out_dir, "plots", "snr_vs_mag.png"), rows)
    except Exception as ex:  # noqa: BLE001
        log("plot failed", repr(ex)[:160])
    log("done", json.dumps(jsonable(summary["forecasts"])))


if __name__ == "__main__":
    sys.exit(main())

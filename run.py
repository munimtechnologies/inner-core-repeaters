#!/usr/bin/env python3
"""inner-core-repeaters: an open, objective, out-of-sample test of inner-core "backtracking".

South Sandwich Islands (SSI) earthquakes recorded on the Eielson array (ILAR, network IM,
IL01-IL19 SHZ) in Alaska, ~150 deg away, where PKIKP (DF, through the inner core) arrives a
few seconds before PKPbc and PKPab (outer core only). Repeating earthquakes are found by
cross-correlating the PKPbc+PKPab window. For every pair, the PKIKP waveform is compared
after aligning on that reference window. Pre-registered tests of the backtracking model
(Wang et al. 2024) vs a monotone model are in model.py and findings/2026-11-17-prereg.md.

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
from cluster import components, max_spanning_offsets
from phases import Phases, gc_deg, load_or_build
from xcorr import best_lag, fshift, lag_scan, measure_pair, pick_onset

VERSION = "0.1.0"
FS = 40.0
PRE = 40.0
POST = 40.0
N = int(round((PRE + POST) * FS))
BAND = (1.0, 2.0)
ARRAY_CENTER = (64.7714, -146.8866)
USGS = "https://earthquake.usgs.gov/fdsnws/event/1/query"
BOX = dict(minlatitude=-61.5, maxlatitude=-54.0, minlongitude=-31.0, maxlongitude=-22.0)
PROC_TAG = "proc_v1"
# ILAR element coordinates (EarthScope station service, 2026-11-17), used if the service fails.
FALLBACK_COORDS = {
    "IL01": (64.771599, -146.886093), "IL02": (64.784698, -146.864304),
    "IL03": (64.7714, -146.851196), "IL04": (64.757004, -146.876099),
    "IL05": (64.773102, -146.922897), "IL06": (64.779198, -146.904007),
    "IL07": (64.799301, -146.839294), "IL08": (64.790298, -146.796906),
    "IL09": (64.768097, -146.783203), "IL10": (64.752899, -146.843094),
    "IL11": (64.741501, -146.8974), "IL12": (64.744698, -146.943604),
    "IL13": (64.747902, -146.986496), "IL14": (64.775002, -146.979401),
    "IL15": (64.777702, -146.942795), "IL16": (64.793297, -146.921494),
    "IL17": (64.807198, -146.889801), "IL18": (64.7575, -146.776794),
    "IL19": (64.746101, -146.797394),
}
T0 = time.time()


def log(*a):
    print("[%7.1f min]" % ((time.time() - T0) / 60.0), *a, flush=True)


def elapsed_min():
    return (time.time() - T0) / 60.0


# ----------------------------------------------------------------------------- catalog
def fetch_catalog(data_dir, end_iso, min_mag_old, min_mag_new, new_from):
    import requests

    path = os.path.join(data_dir, "catalog_usgs.csv")
    rows = {}
    complete = True
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
            if row.get("type", "earthquake") != "earthquake":
                continue
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
                            magType=r.get("magType", ""), t=model.decimal_year(r["time"])))
        except Exception:
            continue
    evs.sort(key=lambda e: e["t"])
    return evs


# ----------------------------------------------------------------------------- waveforms
def get_client():
    from obspy.clients.fdsn import Client

    last = None
    for base in ("https://service.earthscope.org", "IRIS"):
        try:
            c = Client(base, timeout=120)
            log("FDSN client:", base)
            return c
        except Exception as ex:
            last = ex
            log("FDSN client failed:", base, repr(ex)[:200])
    raise last


def get_coords(client, data_dir):
    path = os.path.join(data_dir, "ilar_coords.json")
    if os.path.exists(path):
        with open(path) as f:
            return {k: tuple(v) for k, v in json.load(f).items()}
    coords = {}
    try:
        inv = client.get_stations(network="IM", station="IL*", location="--", channel="SHZ", level="station")
        for net in inv:
            for sta in net:
                coords[sta.code] = (float(sta.latitude), float(sta.longitude))
    except Exception as ex:
        log("station service failed, using fallback coordinates:", repr(ex)[:200])
    coords = {k: v for k, v in coords.items() if k in FALLBACK_COORDS} or {}
    if len(coords) < 10:
        coords = dict(FALLBACK_COORDS)
    with open(path, "w") as f:
        json.dump(coords, f)
    return coords


def fetch_event(client, ev, ph, data_dir):
    from obspy import UTCDateTime, read
    from obspy.clients.fdsn.header import FDSNNoDataException

    rawdir = os.path.join(data_dir, "raw")
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
    dist = float(gc_deg(ev["lat"], ev["lon"], *ARRAY_CENTER))
    tdf = ph.t("df", dist, ev["depth"])
    if not math.isfinite(tdf):
        return None, "no-df-prediction(dist=%.2f)" % dist
    o = UTCDateTime(ev["time"])
    t1 = o + tdf - PRE - 30
    t2 = o + tdf + POST + 30
    err = ""
    for attempt in range(3):
        try:
            st = client.get_waveforms("IM", "IL*", "--", "SHZ", t1, t2)
            if len(st) == 0:
                open(pn, "w").close()
                return None, "nodata"
            try:
                st.write(p, format="MSEED")
            except Exception as ex:
                log("could not cache", ev["id"], repr(ex)[:120])
            return st, "fetched"
        except FDSNNoDataException:
            open(pn, "w").close()
            return None, "nodata"
        except Exception as ex:
            err = repr(ex)[:200]
            time.sleep(5 * (attempt + 1))
    return None, "error:" + err


def process_event(st, ev, ph, coords):
    """Filter 1-2 Hz (causal), cut each element on its own predicted PKIKP time, noise-normalise,
    apply per-event station statics from the whole PKP group, and stack to a beam."""
    from obspy import UTCDateTime

    o = UTCDateTime(ev["time"])
    depth = ev["depth"]
    dist_c = float(gc_deg(ev["lat"], ev["lon"], *ARRAY_CENTER))
    st = st.copy()
    try:
        st.merge(method=0)
    except Exception:
        pass
    names, traces = [], []
    for tr in st:
        sta = tr.stats.station
        if sta not in coords or sta in names:
            continue
        data = tr.data
        if np.ma.isMaskedArray(data):
            if np.ma.count_masked(data) > 0:
                continue
            data = np.asarray(data)
        sr = float(tr.stats.sampling_rate)
        if tr.stats.npts < 60 * sr:
            continue
        tr.data = np.asarray(data, dtype=np.float64)
        try:
            tr.detrend("demean")
            tr.detrend("linear")
            tr.taper(0.05)
            tr.filter("bandpass", freqmin=BAND[0], freqmax=BAND[1], corners=4, zerophase=False)
            lat, lon = coords[sta]
            d = float(gc_deg(ev["lat"], ev["lon"], lat, lon))
            tdf = ph.t("df", d, depth)
            if not math.isfinite(tdf):
                continue
            tr.interpolate(sampling_rate=FS, method="lanczos", a=12, starttime=o + tdf - PRE, npts=N)
        except Exception:
            continue
        x = np.asarray(tr.data, dtype=np.float64)
        if len(x) != N or not np.all(np.isfinite(x)):
            continue
        names.append(sta)
        traces.append(x)
    if len(traces) < 5:
        return None, "few-stations(%d)" % len(traces)
    X = np.vstack(traces)
    n0, n1 = int((PRE - 35) * FS), int((PRE - 5) * FS)
    nrms = np.sqrt(np.mean(X[:, n0:n1] ** 2, axis=1))
    med = np.median(nrms)
    keep = (nrms > 0.2 * med) & (nrms < 5 * med) & (nrms > 0)
    X = X[keep]
    names = [n for n, k in zip(names, keep) if k]
    nrms = nrms[keep]
    if len(X) < 5:
        return None, "few-good-stations(%d)" % len(X)
    X = X / nrms[:, None]
    tdf_c = ph.t("df", dist_c, depth)
    dt_bc = ph.t("bc", dist_c, depth) - tdf_c
    dt_ab = ph.t("ab", dist_c, depth) - tdf_c
    if not math.isfinite(dt_ab):
        return None, "no-PKPab(dist=%.2f)" % dist_c
    g0, g1 = int((PRE - 4) * FS), int((PRE + dt_ab + 12) * FS)
    beam0 = X.mean(axis=0)
    lags = np.zeros(len(X))
    for k in range(len(X)):
        L, C = lag_scan(beam0, X[k], g0, g1, int(0.3 * FS))
        l, c = best_lag(L, C)
        if math.isfinite(l) and c >= 0.5:
            lags[k] = l
    nz = lags != 0
    if nz.any():
        lags[nz] -= lags[nz].mean()
    for k in range(len(X)):
        if lags[k] != 0:
            X[k] = fshift(X[k], -lags[k])
    beam = X.mean(axis=0)
    nb = np.sqrt(np.mean(beam[n0:n1] ** 2))
    if not nb > 0:
        return None, "dead-beam"
    beam = beam / nb
    snr_group = float(np.sqrt(np.mean(beam[g0:g1] ** 2)))
    return dict(names=np.array(names), X=X.astype(np.float32), beam=beam, dist=dist_c, dt_bc=dt_bc,
                dt_ab=dt_ab, snr_group=snr_group, statics_s=lags / FS), "ok"


def save_proc(path, d):
    np.savez(path, names=d["names"], X=d["X"], beam=d["beam"],
             meta=np.array([d["dist"], d["dt_bc"], d["dt_ab"], d["snr_group"]]))


def load_proc(path):
    z = np.load(path, allow_pickle=False)
    m = z["meta"]
    return dict(names=z["names"], X=z["X"], beam=z["beam"], dist=float(m[0]), dt_bc=float(m[1]),
                dt_ab=float(m[2]), snr_group=float(m[3]))


# ----------------------------------------------------------------------------- analysis
def detect_pairs(ids, evs, data, max_km, thr, max_lag_s, min_snr):
    """Edges between events whose PKPbc+PKPab beam windows correlate >= thr."""
    use = [i for i in ids if data[i]["snr_group"] >= min_snr]
    lat = np.array([evs[i]["lat"] for i in use])
    lon = np.array([evs[i]["lon"] for i in use])
    edges = {}
    ncomp = 0
    for a_pos, a in enumerate(use):
        da = data[a]
        start = da["dt_bc"] + 0.5 if math.isfinite(da["dt_bc"]) else da["dt_ab"] - 0.5
        r0 = int((PRE + start) * FS)
        r1 = int((PRE + da["dt_ab"] + 10) * FS)
        if a_pos + 1 >= len(use):
            break
        dkm = gc_deg(lat[a_pos], lon[a_pos], lat[a_pos + 1:], lon[a_pos + 1:]) * 111.195
        for off in np.nonzero(dkm <= max_km)[0]:
            b = use[a_pos + 1 + int(off)]
            L, C = lag_scan(da["beam"], data[b]["beam"], r0, r1, int(max_lag_s * FS))
            l, c = best_lag(L, C)
            ncomp += 1
            if math.isfinite(c) and c >= thr:
                edges[(a, b)] = (float(c), float(l))
    return edges, ncomp, len(use)


def analyze_multiplet(mid, members, evs, data, edges, args):
    info = dict(mult=mid, n=len(members), status="ok")
    root = max(members, key=lambda m: data[m]["snr_group"])
    off = max_spanning_offsets(members, edges, root)
    members = [m for m in members if m in off]
    n0, n1 = int((PRE - 35) * FS), int((PRE - 5) * FS)
    stack = np.zeros(N)
    for m in members:
        b = data[m]["beam"]
        stack += fshift(b, -off[m]) / max(np.sqrt(np.mean(b[n0:n1] ** 2)), 1e-12)
    stack /= len(members)
    trig = pick_onset(stack, FS, center=int(PRE * FS), search_s=6.0, sta_s=0.5, lta_s=8.0, thr=6.0)
    info.update(root=evs[root]["id"])
    if trig is None:
        info["status"] = "no-DF-pick"
        return [], info, None
    o = trig
    dt_bc, dt_ab = data[root]["dt_bc"], data[root]["dt_ab"]
    df_end = min(dt_bc - 0.35, 2.5) if math.isfinite(dt_bc) else 2.5
    if df_end < 0.7:
        info["status"] = "DF-window-too-short(%.2f s)" % df_end
        return [], info, None
    dfw = (o - int(0.3 * FS), o + int(df_end * FS))
    ref_start = (dt_bc - 0.35) if math.isfinite(dt_bc) else (dt_ab - 0.5)
    refw = (o + int(ref_start * FS), o + int((dt_ab + 8.0) * FS))
    info.update(pick_minus_pred_s=round((o - PRE * FS) / FS, 3), df_window_s=round(df_end + 0.3, 2),
                dt_bc=round(dt_bc, 3) if math.isfinite(dt_bc) else None, dt_ab=round(dt_ab, 3))

    def snr(m, w):
        b = data[m]["beam"]
        s = int(round(off[m]))
        a, c = w[0] + s, w[1] + s
        if a < 0 or c > N:
            return float("nan")
        nb = np.sqrt(np.mean(b[n0:n1] ** 2))
        return float(np.sqrt(np.mean(b[a:c] ** 2)) / nb) if nb > 0 else float("nan")

    snr_df = {m: snr(m, dfw) for m in members}
    snr_ref = {m: snr(m, refw) for m in members}
    rows = []
    order = sorted(members, key=lambda m: evs[m]["t"])
    for a, b in combinations(order, 2):
        sa = int(round(off[a]))
        ref = (refw[0] + sa, refw[1] + sa)
        dwa = (dfw[0] + sa, dfw[1] + sa)
        if ref[0] < 0 or ref[1] > N or dwa[0] < 0:
            continue
        lag0 = off[b] - off[a]
        r = measure_pair(data[a]["beam"], data[b]["beam"], ref, dwa, FS, max_lag_ref_s=0.5,
                         max_lag_df_s=0.25, lag0=lag0)
        if r is None:
            continue
        # element-by-element (robust to station clock problems and stacking artefacts)
        na, nb_ = list(data[a]["names"]), list(data[b]["names"])
        st_cc, st_dt = [], []
        common = sorted(set(na) & set(nb_)) if len(members) <= args.max_station_level else []
        for s in common:
            xs = data[a]["X"][na.index(s)].astype(float)
            ys = data[b]["X"][nb_.index(s)].astype(float)
            rs = measure_pair(xs, ys, ref, dwa, FS, max_lag_ref_s=0.25, max_lag_df_s=0.25,
                              lag0=r["lag_total"])
            if rs is not None and rs["cc_ref"] >= 0.6:
                st_cc.append(rs["cc_df0"])
                st_dt.append(rs["dt_df_s"])
        ea, eb = evs[a], evs[b]
        D = 1.0 - r["cc_df0"]
        snrs = (snr_ref[a], snr_ref[b], snr_df[a], snr_df[b])
        okq = (all(math.isfinite(v) for v in snrs) and math.isfinite(D) and math.isfinite(r["cc_ref"])
               and r["cc_ref"] >= args.pair_thr and min(snr_ref[a], snr_ref[b]) >= args.min_snr_ref
               and min(snr_df[a], snr_df[b]) >= args.min_snr_df)
        rows.append(dict(
            mult=mid, ev1=ea["id"], ev2=eb["id"], time1=ea["time"], time2=eb["time"],
            t1=ea["t"], t2=eb["t"], mag1=ea["mag"], mag2=eb["mag"],
            sep_km_catalog=float(gc_deg(ea["lat"], ea["lon"], eb["lat"], eb["lon"]) * 111.195),
            cc_ref=r["cc_ref"], cc_df0=r["cc_df0"], cc_df_best=r["cc_df_best"], dt_df_s=r["dt_df_s"],
            D=D, Dx=D - (1.0 - r["cc_ref"]),
            Ds=(1.0 - float(np.median(st_cc))) if len(st_cc) >= 5 else float("nan"),
            st_median_dt_df_s=float(np.median(st_dt)) if len(st_dt) >= 5 else float("nan"),
            n_sta=len(st_cc), snr_ref1=snr_ref[a], snr_ref2=snr_ref[b], snr_df1=snr_df[a],
            snr_df2=snr_df[b], dM=model.delta_M(ea["t"], eb["t"]), dB=model.delta_B(ea["t"], eb["t"]),
            cross_turn=ea["t"] < model.T_TURN < eb["t"], new_pair=eb["t"] >= model.NEW_CUT, ok=bool(okq)))
    frame = dict(off=off, dfw=dfw, refw=refw, order=order)
    return rows, info, frame


def plot_multiplet(path, mid, frame, evs, data, rows):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    off, dfw, refw, order = frame["off"], frame["dfw"], frame["refw"], frame["order"]
    if len(order) > 25:  # keep the newest and oldest members readable
        order = order[:12] + order[-13:]
    t = (np.arange(N) - dfw[0]) / FS - 0.3
    fig, ax = plt.subplots(figsize=(10, 1.1 * len(order) + 1.5))
    for k, m in enumerate(order):
        b = fshift(data[m]["beam"], -off[m])
        seg = b[dfw[0] - int(3 * FS):refw[1]]
        sc = np.max(np.abs(seg)) if len(seg) else 1.0
        ax.plot(t, b / (sc if sc > 0 else 1.0) * 0.45 - k, lw=0.8)
        ax.text(t[refw[1]] + 0.2, -k, "%s M%.1f" % (evs[m]["time"][:10], evs[m]["mag"]), fontsize=8,
                va="center")
    ax.axvspan(t[dfw[0]], t[dfw[1]], color="orange", alpha=0.2, label="PKIKP (DF) window")
    ax.axvspan(t[refw[0]], t[refw[1]], color="grey", alpha=0.15, label="reference (PKPbc+ab)")
    ax.set_xlim(t[dfw[0]] - 3, t[refw[1]] + 4)
    ax.set_xlabel("time relative to picked PKIKP onset (s)")
    ax.set_yticks([])
    ax.set_title("SSI multiplet %s at ILAR (1-2 Hz beams aligned on reference phases)" % mid)
    ax.legend(loc="lower left", fontsize=7)
    fig.tight_layout()
    fig.savefig(path, dpi=110)
    plt.close(fig)


def plot_summary(path, rows):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    ok = [r for r in rows if r["ok"]]
    if not ok:
        return
    fig, axs = plt.subplots(1, 2, figsize=(11, 4.2), sharey=True)
    for ax, key, lab in ((axs[0], "dM", "|t2 - t1| (years): monotone model"),
                         (axs[1], "dB", "|e(t2) - e(t1)| (years): backtracking model, T=2008.5, ratio 2.5")):
        for sel, col, name in ((lambda r: not r["new_pair"], "C0", "both events < 2025"),
                               (lambda r: r["new_pair"], "C3", "later event >= 2025 (out of sample)")):
            rr = [r for r in ok if sel(r)]
            ax.scatter([r[key] for r in rr], [r["D"] for r in rr], s=12, c=col, label=name, alpha=0.7)
        ax.set_xlabel(lab, fontsize=8)
    axs[0].set_ylabel("PKIKP dissimilarity D = 1 - cc (after reference alignment)")
    axs[0].legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(path, dpi=110)
    plt.close(fig)


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


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", default="/data")
    ap.add_argument("--out-dir", default="/out")
    ap.add_argument("--end", default=time.strftime("%Y-%m-%dT00:00:00", time.gmtime()))
    ap.add_argument("--new-from", default="2024-07-01")
    ap.add_argument("--min-mag-old", type=float, default=5.0)
    ap.add_argument("--min-mag-new", type=float, default=4.7)
    ap.add_argument("--max-km", type=float, default=60.0)
    ap.add_argument("--detect-thr", type=float, default=0.90)
    ap.add_argument("--pair-thr", type=float, default=0.93)
    ap.add_argument("--min-snr-group", type=float, default=3.0)
    ap.add_argument("--min-snr-ref", type=float, default=5.0)
    ap.add_argument("--min-snr-df", type=float, default=3.0)
    ap.add_argument("--max-lag-s", type=float, default=8.0)
    ap.add_argument("--max-events", type=int, default=0)
    ap.add_argument("--time-budget-min", type=float, default=200.0)
    ap.add_argument("--fetch-frac", type=float, default=0.7)
    ap.add_argument("--plots", type=int, default=30)
    ap.add_argument("--max-station-level", type=int, default=40,
                    help="skip element-by-element measurements in multiplets larger than this")
    args = ap.parse_args()
    os.makedirs(args.data_dir, exist_ok=True)
    os.makedirs(args.out_dir, exist_ok=True)
    plots_dir = os.path.join(args.out_dir, "plots")
    os.makedirs(plots_dir, exist_ok=True)
    proc_dir = os.path.join(args.data_dir, PROC_TAG)
    os.makedirs(proc_dir, exist_ok=True)
    summary = dict(version=VERSION, args=vars(args), started=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()))

    grid = load_or_build(os.path.join(args.data_dir, "ak135_pkp_grid_v1.npz"))
    ph = Phases(grid)
    tab = []
    for d in (145.0, 147.0, 149.0, 150.0, 151.0, 152.0, 153.0, 154.0, 156.0):
        df = ph.t("df", d, 15.0)
        tab.append(dict(dist=d, df=df, bc_minus_df=ph.t("bc", d, 15.0) - df, ab_minus_df=ph.t("ab", d, 15.0) - df,
                        p_df=ph.t("p_df", d, 15.0)))
    log("ak135 at 15 km depth:", json.dumps(jsonable(tab)))
    summary["phase_table_15km"] = tab

    evs_list = fetch_catalog(args.data_dir, args.end, args.min_mag_old, args.min_mag_new, args.new_from)
    log("catalog events:", len(evs_list))
    write_csv(os.path.join(args.out_dir, "catalog.csv"), evs_list)
    evs = {e["id"]: e for e in evs_list}

    def prio(e):
        g = 0 if e["t"] >= 2024.5 else (1 if 1998.0 <= e["t"] < 2010.0 else 2)
        return (g, -(e["mag"] if math.isfinite(e["mag"]) else 0))

    order = sorted(evs_list, key=prio)
    if args.max_events:
        order = order[:args.max_events]

    client = None
    coords = None
    data = {}
    status = []
    stop_at = args.time_budget_min * args.fetch_frac
    for k, e in enumerate(order):
        pp = os.path.join(proc_dir, e["id"] + ".npz")
        if os.path.exists(pp):
            try:
                data[e["id"]] = load_proc(pp)
                status.append(dict(id=e["id"], t=e["t"], mag=e["mag"], status="ok-cached"))
                continue
            except Exception:
                os.remove(pp)
        if os.path.exists(os.path.join(args.data_dir, "raw", e["id"] + ".mseed.nodata")):
            status.append(dict(id=e["id"], t=e["t"], mag=e["mag"], status="nodata"))
            continue
        if elapsed_min() > stop_at:
            status.append(dict(id=e["id"], t=e["t"], mag=e["mag"], status="skipped-time-budget"))
            continue
        if client is None:
            try:
                client = get_client()
            except Exception as ex:
                log("no FDSN client, stopping downloads:", repr(ex)[:200])
                status.append(dict(id=e["id"], t=e["t"], mag=e["mag"], status="no-client"))
                stop_at = -1.0
                continue
            coords = get_coords(client, args.data_dir)
            log("ILAR elements:", sorted(coords))
        try:
            st, s = fetch_event(client, e, ph, args.data_dir)
            if st is None:
                status.append(dict(id=e["id"], t=e["t"], mag=e["mag"], status=s))
                continue
            d, s2 = process_event(st, e, ph, coords)
            if d is None:
                status.append(dict(id=e["id"], t=e["t"], mag=e["mag"], status=s2))
                continue
            save_proc(pp, d)
            data[e["id"]] = d
            status.append(dict(id=e["id"], t=e["t"], mag=e["mag"], status="ok", n_sta=len(d["names"]),
                               snr_group=round(d["snr_group"], 2)))
        except Exception as ex:
            status.append(dict(id=e["id"], t=e["t"], mag=e["mag"], status="exception:" + repr(ex)[:150]))
        if (k + 1) % 50 == 0:
            log("events tried %d/%d, processed %d" % (k + 1, len(order), len(data)))
        time.sleep(0.2)
    write_csv(os.path.join(args.out_dir, "event_status.csv"), status, keys=["id", "t", "mag", "status", "n_sta",
                                                                          "snr_group"])
    by_year = {}
    for s in status:
        y = int(s["t"])
        b = by_year.setdefault(y, dict(tried=0, ok=0))
        b["tried"] += 1
        b["ok"] += int(str(s["status"]).startswith("ok"))
    summary["coverage_by_year"] = by_year
    summary["n_processed"] = len(data)
    log("processed events:", len(data))

    ids = sorted(data, key=lambda i: evs[i]["t"])
    edges, ncomp, nuse = detect_pairs(ids, evs, data, args.max_km, args.detect_thr, args.max_lag_s,
                                      args.min_snr_group)
    log("detection: %d events with snr>=%.1f, %d comparisons, %d edges >= %.2f"
        % (nuse, args.min_snr_group, ncomp, len(edges), args.detect_thr))
    write_csv(os.path.join(args.out_dir, "edges.csv"),
              [dict(ev1=evs[a]["id"], ev2=evs[b]["id"], t1=evs[a]["t"], t2=evs[b]["t"], cc=c, lag_s=l / FS)
               for (a, b), (c, l) in edges.items()])
    comps = [c for c in components(ids, edges) if len(c) >= 2]
    comps.sort(key=lambda c: -len(c))
    log("multiplets:", len(comps), "sizes:", [len(c) for c in comps[:20]])
    # analyse multiplets with a 2025+ member first, then smaller ones (a giant chained
    # component, if any, goes last so it cannot eat the time budget)
    work = sorted(enumerate(comps), key=lambda kc: (not any(evs[m]["t"] >= model.NEW_CUT for m in kc[1]),
                                                     len(kc[1])))

    rows, infos, frames = [], [], {}
    for k, members in work:
        mid = "M%03d" % k
        if elapsed_min() > args.time_budget_min * 0.95:
            infos.append(dict(mult=mid, n=len(members), status="skipped-time-budget",
                              members=";".join(evs[m]["id"] for m in members), has_new=None, span=""))
            continue
        try:
            r, info, fr = analyze_multiplet(mid, members, evs, data, edges, args)
        except Exception as ex:
            r, info, fr = [], dict(mult=mid, n=len(members), status="exception:" + repr(ex)[:150]), None
            traceback.print_exc()
        info["members"] = ";".join("%s(%s)" % (evs[m]["id"], evs[m]["time"][:10]) for m in members)
        info["has_new"] = any(evs[m]["t"] >= model.NEW_CUT for m in members)
        info["span"] = "%.1f-%.1f" % (min(evs[m]["t"] for m in members), max(evs[m]["t"] for m in members))
        rows += r
        infos.append(info)
        if fr is not None:
            frames[mid] = (fr, members)
    write_csv(os.path.join(args.out_dir, "multiplets.csv"), infos,
              keys=["mult", "n", "status", "root", "span", "has_new", "pick_minus_pred_s", "df_window_s", "dt_bc",
                    "dt_ab", "members"])
    write_csv(os.path.join(args.out_dir, "pairs.csv"), rows)
    log("pairs measured:", len(rows), "ok:", sum(r["ok"] for r in rows),
        "out-of-sample ok:", sum(r["ok"] and r["new_pair"] for r in rows),
        "cross-turn ok:", sum(r["ok"] and r["cross_turn"] for r in rows))

    try:
        res = model.evaluate(rows)
        sens = {}
        for thr in (0.90, 0.95):
            rr = [dict(r, ok=bool(all(math.isfinite(r[k]) for k in ("snr_ref1", "snr_ref2", "snr_df1",
                                                                     "snr_df2", "D", "cc_ref"))
                                  and r["cc_ref"] >= thr and min(r["snr_ref1"], r["snr_ref2"]) >= args.min_snr_ref
                                  and min(r["snr_df1"], r["snr_df2"]) >= args.min_snr_df))
                  for r in rows]
            e2 = model.evaluate(rr, nboot=300)
            sens["pair_thr_%.2f" % thr] = dict(P1=e2["P1_in_sample"], P2=e2["P2_out_of_sample"],
                                               P3=dict((k, v) for k, v in e2["P3_new_event_best_match"].items()
                                                       if k != "cases"))
        res["sensitivity"] = sens
    except Exception as ex:
        traceback.print_exc()
        res = dict(error=repr(ex))
    summary["predictions"] = res
    with open(os.path.join(args.out_dir, "predictions.json"), "w") as f:
        json.dump(jsonable(res), f, indent=1)

    try:
        plot_summary(os.path.join(plots_dir, "D_vs_time_separation.png"), rows)
        shown = 0
        pri = sorted(frames, key=lambda m: (not any(evs[x]["t"] >= model.NEW_CUT for x in frames[m][1]),
                                            -len(frames[m][1])))
        for mid in pri:
            if shown >= args.plots:
                break
            fr, members = frames[mid]
            plot_multiplet(os.path.join(plots_dir, "multiplet_%s.png" % mid), mid, fr, evs, data,
                           [r for r in rows if r["mult"] == mid])
            shown += 1
    except Exception:
        traceback.print_exc()

    summary["finished"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    summary["minutes"] = round(elapsed_min(), 1)
    summary["n_edges"] = len(edges)
    summary["n_multiplets"] = len(comps)
    summary["n_pairs"] = len(rows)
    with open(os.path.join(args.out_dir, "run_summary.json"), "w") as f:
        json.dump(jsonable(summary), f, indent=1)
    write_report(os.path.join(args.out_dir, "report.md"), summary, res, infos)
    log("done")


def write_report(path, summary, res, infos):
    L = ["# inner-core-repeaters %s run report" % VERSION, "",
         "Started %s, finished %s (%.1f min)." % (summary.get("started"), summary.get("finished"),
                                                  summary.get("minutes", 0)), "",
         "* events processed: %s" % summary.get("n_processed"),
         "* repeater edges: %s, multiplets: %s, pairs measured: %s" % (summary.get("n_edges"),
                                                                      summary.get("n_multiplets"),
                                                                      summary.get("n_pairs")), ""]
    if isinstance(res, dict) and "P1_in_sample" in res:
        for k in ("P1_in_sample", "P2_out_of_sample", "P4_calibration"):
            v = res[k]
            L.append("## %s" % k)
            L.append("```")
            L.append(json.dumps(jsonable(v), indent=1))
            L.append("```")
        p3 = res["P3_new_event_best_match"]
        L.append("## P3_new_event_best_match")
        L.append("discriminating new events: %d, best match = backtracking pick: %d, = monotone pick: %d, holds: %s"
                 % (p3["n_discriminating"], p3["n_match_B"], p3["n_match_M"], p3["holds"]))
        for c in p3["cases"]:
            L.append("* %s" % json.dumps(jsonable(c)))
        L.append("")
        L.append("## P5 exploratory fit")
        L.append("```")
        L.append(json.dumps(jsonable(dict(in_sample=res["P5_fit_in_sample"], all=res["P5_fit_all"])), indent=1))
        L.append("```")
    L.append("## Multiplets with a member from 2025 on")
    for i in infos:
        if i.get("has_new"):
            L.append("* %s n=%s %s span %s: %s" % (i["mult"], i["n"], i["status"], i["span"], i["members"]))
    L.append("")
    L.append("## Coverage by year (tried / processed)")
    cov = summary.get("coverage_by_year", {})
    L.append(", ".join("%s: %d/%d" % (y, cov[y]["tried"], cov[y]["ok"]) for y in sorted(cov)))
    with open(path, "w") as f:
        f.write("\n".join(L) + "\n")


if __name__ == "__main__":
    main()

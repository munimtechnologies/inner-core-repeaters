#!/usr/bin/env python3
"""BC-DF times of single South Sandwich Islands events at College, Alaska (COL/COLA), 1982-2026.

Part of inner-core-repeaters by MV2, Munim, Inc. (made with AI agents). Licence: Apache-2.0.

The classic inner-core rotation observable (Song & Richards 1996, Nature 382:221) on the
station that defined it. Wu, Song & Yang 2026 (GRL 53, doi:10.1029/2025GL119295) measured
151 SSI->COLA single-event residuals from 1982 to 2024. This script rebuilds that series
from the open archive with an objective, published measurement and extends it to the
2025-2026 events that no published study has used. Pre-registration:
findings/2026-11-17-cola-prereg.md (written before any COL/COLA waveform was downloaded).

The College site has one continuous record in the EarthScope archive:
  DW.COL  (blank) SHZ 20 Hz  1982-01-08 .. 1991-03-04
  IU.COL  (blank) BHZ 20 Hz  1991-03-04 .. 1996-05-11
  IU.COLA (blank) BHZ 20 Hz  1996-06-14 .. 1999-02-05
  IU.COLA 00      BHZ 20/40  1999-02-05 .. now   (primary)
  IU.COLA 10      BHZ 40 Hz  1999-02-05 .. now   (second sensor, used only as a cross-check)
BC-DF is a differential time on one trace, so station clock errors cancel.
"""
import argparse
import csv
import io
import json
import math
import os
import time

import numpy as np

import bcdf
import model
from phases import Phases, gc_deg, load_or_build

VERSION = "0.4.0"
FS = bcdf.FS
PRE = bcdf.PRE
POST = 40.0
N = int(round((PRE + POST) * FS))
BAND = (1.0, 2.0)        # primary (same as the ILAR beams)
BAND_ALT = (0.6, 1.6)    # robustness band
USGS = "https://earthquake.usgs.gov/fdsnws/event/1/query"
BOX = dict(minlatitude=-61.5, maxlatitude=-54.0, minlongitude=-31.0, maxlongitude=-22.0)
PROC_TAG = "proc_cola_v1"
T_EARLY_END = 2008.5     # C1 calibration window: events before the Wang et al. turn
# (start, end, net, sta, loc, cha); ISO date strings compare correctly with ComCat times.
EPOCHS = [
    ("1982-01-08", "1991-03-04", "DW", "COL", "", "SHZ"),
    ("1991-03-04", "1996-05-11", "IU", "COL", "", "BHZ"),
    ("1996-06-14", "1999-02-05", "IU", "COLA", "", "BHZ"),
    ("1999-02-05", "9999-12-31", "IU", "COLA", "00", "BHZ"),
]
SECONDARY = [("1999-02-05", "9999-12-31", "IU", "COLA", "10", "BHZ")]
# EarthScope station service, 2026-11-17
COORDS = {("DW", "COL"): (64.900002, -147.793304), ("IU", "COL"): (64.8993, -147.7954),
          ("IU", "COLA"): (64.873566, -147.861772)}
T0 = time.time()


def log(*a):
    print("[%7.1f min]" % ((time.time() - T0) / 60.0), *a, flush=True)


def channel_for(iso, epochs=EPOCHS):
    """The (net, sta, loc, cha) recording at College at ISO time `iso`, or None."""
    s = str(iso)[:19]
    for a, b, net, sta, loc, cha in epochs:
        if a <= s < b:
            return net, sta, loc, cha
    return None


def chan_id(ch):
    return "%s.%s.%s.%s" % ch


# ----------------------------------------------------------------------------- catalog
def fetch_catalog(data_dir, start, end_iso, min_mag_old, min_mag_new, new_from):
    import requests

    path = os.path.join(data_dir, "catalog_cola.csv")
    rows = {}
    complete = True
    for s, e, m in [(start, new_from, min_mag_old), (new_from, end_iso, min_mag_new)]:
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
    return parse_catalog(rows.values())


def parse_catalog(rows):
    evs = []
    for r in rows:
        try:
            evs.append(dict(id=r["id"], time=r["time"], lat=float(r["latitude"]), lon=float(r["longitude"]),
                            depth=float(r["depth"]) if r.get("depth") not in (None, "") else 10.0,
                            mag=float(r["mag"]) if r.get("mag") not in (None, "") else float("nan"),
                            t=model.decimal_year(r["time"])))
        except (KeyError, ValueError, TypeError):
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


def fetch(client, ev, ch, ph, data_dir):
    from obspy import UTCDateTime, read
    from obspy.clients.fdsn.header import FDSNNoDataException

    d = os.path.join(data_dir, "raw", chan_id(ch))
    os.makedirs(d, exist_ok=True)
    p = os.path.join(d, ev["id"] + ".mseed")
    pn = p + ".nodata"
    if os.path.exists(p):
        try:
            return read(p), "cached"
        except Exception:
            os.remove(p)
    if os.path.exists(pn):
        return None, "nodata"
    lat, lon = COORDS[(ch[0], ch[1])]
    dist = float(gc_deg(ev["lat"], ev["lon"], lat, lon))
    tdf = ph.t("df", dist, ev["depth"])
    if not math.isfinite(tdf):
        return None, "no-df-prediction(dist=%.2f)" % dist
    o = UTCDateTime(ev["time"])
    t1, t2 = o + tdf - PRE - 30, o + tdf + POST + 30
    err = ""
    for attempt in range(3):
        try:
            st = client.get_waveforms(ch[0], ch[1], ch[2] or "--", ch[3], t1, t2)
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


def process_trace(st, ev, ph, coords, band=BAND):
    """One filtered (causal), noise-normalised 40 Hz trace with ak135 DF at sample PRE*FS.

    Returns (dict(x, dist, dt_bc, dt_ab, snr_group), "ok") or (None, reason).
    """
    from obspy import UTCDateTime

    o = UTCDateTime(ev["time"])
    st = st.copy()
    try:
        st.merge(method=0)
    except Exception:
        pass
    if len(st) == 0:
        return None, "empty"
    tr = max(st, key=lambda t: t.stats.npts)
    data = tr.data
    if np.ma.isMaskedArray(data):
        if np.ma.count_masked(data) > 0:
            return None, "gappy"
        data = np.asarray(data)
    sr = float(tr.stats.sampling_rate)
    if tr.stats.npts < 60 * sr:
        return None, "short"
    raw = np.asarray(data, dtype=np.float64)
    if np.ptp(raw) == 0:
        return None, "flat"
    lat, lon = coords
    dist = float(gc_deg(ev["lat"], ev["lon"], lat, lon))
    depth = ev["depth"]
    tdf = ph.t("df", dist, depth)
    if not math.isfinite(tdf):
        return None, "no-df-prediction(dist=%.2f)" % dist
    dt_bc = ph.t("bc", dist, depth) - tdf
    dt_ab = ph.t("ab", dist, depth) - tdf
    tr.data = raw
    try:
        tr.detrend("demean")
        tr.detrend("linear")
        tr.taper(0.05)
        tr.filter("bandpass", freqmin=band[0], freqmax=band[1], corners=4, zerophase=False)
        tr.interpolate(sampling_rate=FS, method="lanczos", a=12, starttime=o + tdf - PRE, npts=N)
    except Exception as ex:
        return None, "process-error:" + repr(ex)[:80]
    x = np.asarray(tr.data, dtype=np.float64)
    if len(x) != N or not np.all(np.isfinite(x)):
        return None, "bad-length"
    n0, n1 = int((PRE - 35) * FS), int((PRE - 5) * FS)
    nr = float(np.sqrt(np.mean(x[n0:n1] ** 2)))
    if not nr > 0:
        return None, "dead"
    x = x / nr
    g1 = int((PRE + (dt_ab if math.isfinite(dt_ab) else 5.0) + 6) * FS)
    snr_group = float(np.sqrt(np.mean(x[int((PRE - 2) * FS):g1] ** 2)))
    return dict(x=x, dist=dist, dt_bc=dt_bc, dt_ab=dt_ab, snr_group=snr_group), "ok"


def measure_row(ev, d, ch, x_alt=None):
    """BC-DF measurement on a processed trace (bcdf.measure_bcdf), as a row; None if not measurable."""
    r = bcdf.measure_bcdf(d["x"], d["dt_bc"], d["dt_ab"])
    if r is None:
        return None
    ok = r["cc"] >= bcdf.MIN_CC and r["snr_df"] >= bcdf.MIN_SNR_DF and abs(r["onset_s"]) <= 3.0
    row = dict(id=ev["id"], time=ev["time"], t=ev["t"], lat=ev["lat"], lon=ev["lon"], depth=ev["depth"],
               mag=ev["mag"], chan=chan_id(ch), dist=d["dist"], dt_bc_pred=d["dt_bc"], dt_ab_pred=d["dt_ab"],
               onset_s=r["onset_s"], dt_obs=r["dt_obs_s"], res=r["dt_obs_s"] - d["dt_bc"], cc=r["cc"],
               snr_df=r["snr_df"], snr_group=d["snr_group"], ok=bool(ok), new=ev["t"] >= model.NEW_CUT,
               res_alt=float("nan"), cc_alt=float("nan"), res_sec=float("nan"), cc_sec=float("nan"))
    if x_alt is not None:
        ra = bcdf.measure_bcdf(x_alt, d["dt_bc"], d["dt_ab"], onset=r["onset"])
        if ra is not None:
            row["res_alt"], row["cc_alt"] = ra["dt_obs_s"] - d["dt_bc"], ra["cc"]
    return row


# ----------------------------------------------------------------------------- statistics
def screen(rows, key="res", ok_key="ok", max_dev=1.0):
    """Pre-registered usable set: quality flags, then |res - median| <= 1 s (cycle skips)."""
    ok = [r for r in rows if r.get(ok_key) and math.isfinite(r.get(key, float("nan")))]
    if not ok:
        return []
    med = float(np.median([r[key] for r in ok]))
    return [r for r in ok if abs(r[key] - med) <= max_dev]


def compare(ins, outs, k1, k2, key="res", nboot=2000, seed=0):
    """Fit k1 and k2 on `ins`, predict `outs`; the pre-registered rule of bcdf.out_of_sample, generalised.

    decisive: >= 8 new events and |mean(pred_k2 - pred_k1)| >= 2 * sigma_in / sqrt(n_new).
    verdict: k2 if k2 has the lower out-of-sample MSE in >= 90% of bootstrap resamples, k1 if
    in <= 10%, otherwise "neither"; "not decisive" otherwise.
    """
    res = dict(models=[k1, k2], n_in=len(ins), n_out=len(outs))
    if len(ins) < 15 or len(outs) < 1:
        res.update(testable=False, decisive=False, verdict="not testable",
                   reason="need >= 15 in-sample and >= 1 new events")
        return res
    f1, f2 = bcdf.fit(ins, k1, key), bcdf.fit(ins, k2, key)
    y = np.array([r[key] for r in outs], dtype=float)
    p1, p2 = bcdf.predict(f1, outs), bcdf.predict(f2, outs)
    e1, e2 = y - p1, y - p2
    sig = math.sqrt(min(f1["rss"], f2["rss"]) / max(1, len(ins) - 5))
    gap = float(np.mean(p2 - p1))
    se = sig / math.sqrt(len(outs))
    rng = np.random.default_rng(seed)
    wins = []
    for _ in range(nboot):
        k = rng.integers(0, len(outs), len(outs))
        wins.append(float(np.mean(e2[k] ** 2) < np.mean(e1[k] ** 2)))
    frac = float(np.mean(wins))
    decisive = len(outs) >= 8 and abs(gap) >= 2 * se
    names = ["c", "dlat", "dlon", "ddepth100", "b"]
    res.update(testable=True, decisive=bool(decisive), sigma_in=sig, pred_gap=gap, se_new_mean=se,
               coef_1=dict(zip(names, map(float, f1["coef"]))), coef_2=dict(zip(names, map(float, f2["coef"]))),
               rss_in_1=f1["rss"], rss_in_2=f2["rss"], mse_out_1=float(np.mean(e1 ** 2)),
               mse_out_2=float(np.mean(e2 ** 2)), mean_err_out_1=float(np.mean(e1)),
               mean_err_out_2=float(np.mean(e2)), frac_boot_2_better=frac,
               verdict=(k2 if frac >= 0.9 else (k1 if frac <= 0.1 else "neither")) if decisive else "not decisive")
    return res


def calibration_c1(rows, key="res", t_end=T_EARLY_END, nboot=1000, seed=1):
    """C1: on events before 2008.5, the L slope b (s/yr) is positive in >= 90% of event bootstraps."""
    early = [r for r in rows if r["t"] < t_end]
    if len(early) < 20 or (max(r["t"] for r in early) - min(r["t"] for r in early)) < 8:
        return dict(n=len(early), testable=False, holds=None)
    f = bcdf.fit(early, "L", key)
    rng = np.random.default_rng(seed)
    bs = []
    for _ in range(nboot):
        k = rng.integers(0, len(early), len(early))
        bs.append(float(bcdf.fit([early[i] for i in k], "L", key)["coef"][4]))
    frac = float(np.mean(np.array(bs) > 0))
    return dict(n=len(early), span=[min(r["t"] for r in early), max(r["t"] for r in early)],
                b_s_per_yr=float(f["coef"][4]), b_p05=float(np.percentile(bs, 5)), b_p95=float(np.percentile(bs, 95)),
                frac_boot_positive=frac, testable=True, holds=frac >= 0.9)


def secondary_consistency(rows):
    d = [r["res"] - r["res_sec"] for r in rows if r.get("ok") and math.isfinite(r.get("res_sec", float("nan")))
         and r.get("cc_sec", 0) >= bcdf.MIN_CC]
    if not d:
        return dict(n=0)
    d = np.array(d)
    return dict(n=int(len(d)), median_diff=float(np.median(d)), mad_diff=float(np.median(np.abs(d - np.median(d)))),
                frac_within_0p05=float(np.mean(np.abs(d) <= 0.05)))


def run_tests(use, key="res"):
    ins = [r for r in use if not r["new"]]
    outs = [r for r in use if r["new"]]
    out = dict(n_use=len(use), n_in=len(ins), n_new=len(outs))
    out["C1_calibration"] = calibration_c1(ins, key)
    if len(ins) >= 15:
        fL, fB, fP = bcdf.fit(ins, "L", key), bcdf.fit(ins, "Bk", key), bcdf.fit(ins, "P", key)
        out["Q1_in_sample"] = dict(rss_L=fL["rss"], rss_Bk=fB["rss"], rss_P=fP["rss"], b_L=float(fL["coef"][4]),
                                   b_Bk=float(fB["coef"][4]), b_P=float(fP["coef"][4]),
                                   best=min((("L", fL["rss"]), ("Bk", fB["rss"]), ("P", fP["rss"])),
                                            key=lambda kv: kv[1])[0])
        out["Q3_turn_scan_in_sample"] = bcdf.scan_turn(ins, key)
    out["Q2_L_vs_Bk"] = compare(ins, outs, "L", "Bk", key)
    out["Q2b_P_vs_Bk"] = compare(ins, outs, "P", "Bk", key)
    out["Q2c_L_vs_P"] = compare(ins, outs, "L", "P", key)
    return out


def yearly(rows, vals):
    out = {}
    for r, v in zip(rows, vals):
        out.setdefault(int(r["t"]), []).append(float(v))
    return {y: dict(n=len(v), median=float(np.median(v))) for y, v in sorted(out.items())}


# ----------------------------------------------------------------------------- output
def plot(path, use, fits):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fL = fits["L"]
    corr = bcdf.spatially_corrected(fL, use)
    t = np.array([r["t"] for r in use])
    new = t >= model.NEW_CUT
    fig, ax = plt.subplots(figsize=(10, 4.5))
    ax.scatter(t[~new], corr[~new], s=10, c="C0", label="events < 2025 (fit)")
    ax.scatter(t[new], corr[new], s=18, c="C3", label="events >= 2025 (out of sample)")
    tt = np.linspace(min(t.min(), 1982), max(t.max(), 2027), 500)
    for kind, col in (("L", "k"), ("Bk", "C2"), ("P", "C1")):
        f = fits[kind]
        c = f["coef"]
        if kind == "L":
            g = tt - f["T"]
        elif kind == "Bk":
            g = np.array([model.equiv_date(x) for x in tt]) - f["T"]
        else:
            g = np.minimum(tt, f["T"]) - f["T"]
        ax.plot(tt, c[0] + c[4] * g, color=col, lw=1.2, label="%s fit (< 2025)" % kind)
    ax.set_xlabel("year")
    ax.set_ylabel("BC-DF residual minus spatial terms (s)")
    ax.set_title("SSI -> College, Alaska (DW.COL, IU.COL, IU.COLA): single-event BC-DF, 1-2 Hz")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=110)
    plt.close(fig)


def plot_examples(path, items):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    if not items:
        return
    fig, ax = plt.subplots(figsize=(9, 0.9 * len(items) + 1.5))
    tt = np.arange(N) / FS - PRE
    for k, (row, x) in enumerate(items):
        sc = np.max(np.abs(x[int((PRE - 3) * FS):int((PRE + 10) * FS)])) or 1.0
        ax.plot(tt, x / sc * 0.45 - k, lw=0.7)
        ax.plot([row["onset_s"]] * 2, [-k - 0.4, -k + 0.4], "r-", lw=1)
        ax.plot([row["onset_s"] + row["dt_obs"]] * 2, [-k - 0.4, -k + 0.4], "b-", lw=1)
        ax.text(10.3, -k, "%s M%.1f res %+.2f cc %.2f" % (row["time"][:10], row["mag"], row["res"], row["cc"]),
                fontsize=7, va="center")
    ax.set_xlim(-4, 16)
    ax.set_yticks([])
    ax.set_xlabel("time relative to ak135 PKIKP (s); red = DF onset, blue = DF onset + measured BC-DF")
    fig.tight_layout()
    fig.savefig(path, dpi=110)
    plt.close(fig)


def jsonable(o):
    return bcdf.jsonable(o)


def write_csv(path, rows, keys):
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)


def write_report(path, s):
    L = ["# COL/COLA single-event BC-DF run (cola.py %s)" % VERSION, "",
         "Catalogue events: %s; measured: %s; usable (primary): %s (before 2025: %s, 2025 on: %s)."
         % (s.get("n_catalog"), s.get("n_measured"), s["primary"].get("n_use"), s["primary"].get("n_in"),
            s["primary"].get("n_new")), "",
         "## Coverage by year (tried / processed / usable)", ""]
    cov = s.get("coverage_by_year", {})
    L.append(", ".join("%s: %d/%d/%d" % (y, c["tried"], c["processed"], c["usable"]) for y, c in sorted(cov.items())))
    L += ["", "## Channels used", "```", json.dumps(jsonable(s.get("by_channel", {})), indent=1), "```"]
    for name in ("primary", "alt_band"):
        L += ["", "## %s" % name, "```", json.dumps(jsonable(s.get(name, {})), indent=1), "```"]
    L += ["", "## IU.COLA 10 vs 00 (same events)", "```", json.dumps(jsonable(s.get("secondary_consistency")), indent=1),
          "```", "", "## Status counts", "```", json.dumps(jsonable(s.get("status_counts")), indent=1), "```"]
    with open(path, "w") as f:
        f.write("\n".join(L) + "\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", default="/data")
    ap.add_argument("--out-dir", default="/out")
    ap.add_argument("--start", default="1982-01-08")
    ap.add_argument("--end", default=time.strftime("%Y-%m-%dT00:00:00", time.gmtime()))
    ap.add_argument("--new-from", default="2024-07-01")
    ap.add_argument("--min-mag-old", type=float, default=5.0)
    ap.add_argument("--min-mag-new", type=float, default=4.7)
    ap.add_argument("--max-events", type=int, default=0)
    ap.add_argument("--time-budget-min", type=float, default=190.0)
    ap.add_argument("--fetch-frac", type=float, default=0.9)
    ap.add_argument("--no-secondary", action="store_true")
    args = ap.parse_args()
    os.makedirs(args.data_dir, exist_ok=True)
    os.makedirs(args.out_dir, exist_ok=True)
    proc_dir = os.path.join(args.data_dir, PROC_TAG)
    os.makedirs(proc_dir, exist_ok=True)
    s = dict(version=VERSION, args=vars(args), started=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()))

    ph = Phases(load_or_build(os.path.join(args.data_dir, "ak135_pkp_grid_v1.npz")))
    lat, lon = COORDS[("IU", "COLA")]
    tab = []
    for dd in (146.0, 148.0, 150.0, 151.0, 152.0, 153.0, 154.0):
        df = ph.t("df", dd, 15.0)
        tab.append(dict(dist=dd, bc_minus_df=ph.t("bc", dd, 15.0) - df, ab_minus_df=ph.t("ab", dd, 15.0) - df))
    s["phase_table_15km"] = tab
    log("ak135:", json.dumps(jsonable(tab)))

    evs = fetch_catalog(args.data_dir, args.start, args.end, args.min_mag_old, args.min_mag_new, args.new_from)
    s["n_catalog"] = len(evs)
    log("catalogue events:", len(evs))
    dists = [float(gc_deg(e["lat"], e["lon"], lat, lon)) for e in evs]
    if dists:
        s["dist_range_cola"] = [min(dists), max(dists)]

    def prio(e):
        return (0 if e["t"] >= 2024.5 else 1, -(e["mag"] if math.isfinite(e["mag"]) else 0))

    order = sorted(evs, key=prio)
    if args.max_events:
        order = order[:args.max_events]
    client = None
    rows, status, examples = [], [], []
    stop_at = args.time_budget_min * args.fetch_frac
    for k, e in enumerate(order):
        ch = channel_for(e["time"])
        st_row = dict(id=e["id"], t=e["t"], mag=e["mag"], chan=chan_id(ch) if ch else "", status="")
        status.append(st_row)
        if ch is None:
            st_row["status"] = "no-channel-epoch"
            continue
        try:
            client = handle_event(e, ch, st_row, client, ph, args, proc_dir, rows, examples, stop_at)
        except Exception as ex:
            st_row["status"] = "exception:" + repr(ex)[:120]
        if isinstance(client, str):  # "STOP": no FDSN client could be made
            client, stop_at = None, -1.0
        if (k + 1) % 100 == 0:
            log("events tried %d/%d, measured %d" % (k + 1, len(order), len(rows)))
    _finish(args, s, rows, status, examples)


def handle_event(e, ch, st_row, client, ph, args, proc_dir, rows, examples, stop_at):
    """Fetch (or load), process and measure one event; appends to rows/examples. Returns the client."""
    pp = os.path.join(proc_dir, e["id"] + ".npz")
    d = None
    x_alt = None
    if os.path.exists(pp):
        try:
            z = np.load(pp, allow_pickle=False)
            m = z["meta"]
            d = dict(x=z["x"], dist=float(m[0]), dt_bc=float(m[1]), dt_ab=float(m[2]), snr_group=float(m[3]))
            x_alt = z["x_alt"] if "x_alt" in z.files and len(z["x_alt"]) == N else None
        except Exception:
            os.remove(pp)
            d = None
    if d is None:
        if (time.time() - T0) / 60.0 > stop_at:
            st_row["status"] = "skipped-time-budget"
            return client
        if client is None:
            try:
                client = get_client()
            except Exception as ex:
                log("no FDSN client:", repr(ex)[:200])
                st_row["status"] = "no-client"
                return "STOP"
        st, why = fetch(client, e, ch, ph, args.data_dir)
        if st is None:
            st_row["status"] = why
            return client
        d, why = process_trace(st, e, ph, COORDS[(ch[0], ch[1])], BAND)
        if d is None:
            st_row["status"] = why
            return client
        da, _ = process_trace(st, e, ph, COORDS[(ch[0], ch[1])], BAND_ALT)
        x_alt = da["x"] if da is not None else None
        np.savez(pp, x=d["x"].astype(np.float32),
                 x_alt=(x_alt if x_alt is not None else np.zeros(1)).astype(np.float32),
                 meta=np.array([d["dist"], d["dt_bc"], d["dt_ab"], d["snr_group"]]))
        time.sleep(0.1)
    row = measure_row(e, d, ch, x_alt)
    if row is None:
        st_row["status"] = "no-measurement(dt_bc=%.2f)" % d["dt_bc"]
        return client
    st_row["status"] = "ok" if row["ok"] else "low-quality"
    # second sensor (IU.COLA 10), same measurement, only for usable primary rows
    sec = channel_for(e["time"], SECONDARY)
    if row["ok"] and sec is not None and not args.no_secondary and (time.time() - T0) / 60.0 <= stop_at:
        try:
            if client is None:
                client = get_client()
            st2, _ = fetch(client, e, sec, ph, args.data_dir)
            if st2 is not None:
                d2, _ = process_trace(st2, e, ph, COORDS[(sec[0], sec[1])], BAND)
                if d2 is not None:
                    r2 = bcdf.measure_bcdf(d2["x"], d2["dt_bc"], d2["dt_ab"])
                    if r2 is not None:
                        row["res_sec"], row["cc_sec"] = r2["dt_obs_s"] - d2["dt_bc"], r2["cc"]
        except Exception as ex:
            log("secondary failed", e["id"], repr(ex)[:120])
    rows.append(row)
    if row["ok"] and (row["new"] or len(examples) < 12) and len(examples) < 24:
        examples.append((row, np.asarray(d["x"], dtype=float)))
    return client


def _finish(args, s, rows, status, examples):
    rows.sort(key=lambda r: r["t"])
    keys = ["id", "time", "t", "lat", "lon", "depth", "mag", "chan", "dist", "dt_bc_pred", "dt_ab_pred", "onset_s",
            "dt_obs", "res", "cc", "snr_df", "snr_group", "ok", "new", "res_alt", "cc_alt", "res_sec", "cc_sec"]
    write_csv(os.path.join(args.out_dir, "cola_events.csv"), rows, keys)
    write_csv(os.path.join(args.out_dir, "cola_event_status.csv"), status, ["id", "t", "mag", "chan", "status"])
    s["n_measured"] = len(rows)
    counts = {}
    for st_row in status:
        key = str(st_row["status"]).split("(")[0].split(":")[0]
        counts[key] = counts.get(key, 0) + 1
    s["status_counts"] = counts

    use = screen(rows)
    s["primary"] = run_tests(use)
    alt_rows = [dict(r, res=r["res_alt"], ok=bool(r["ok"] and r["cc_alt"] >= bcdf.MIN_CC)) for r in rows
                if math.isfinite(r["res_alt"])]
    s["alt_band"] = run_tests(screen(alt_rows))
    s["secondary_consistency"] = secondary_consistency(rows)
    usable_ids = {r["id"] for r in use}
    cov, bych = {}, {}
    for st_row in status:
        y = int(st_row["t"])
        c = cov.setdefault(y, dict(tried=0, processed=0, usable=0))
        c["tried"] += 1
        c["processed"] += int(st_row["status"] in ("ok", "low-quality") or st_row["status"].startswith("no-measurement"))
        c["usable"] += int(st_row["id"] in usable_ids)
        b = bych.setdefault(st_row["chan"] or "none", dict(tried=0, usable=0))
        b["tried"] += 1
        b["usable"] += int(st_row["id"] in usable_ids)
    s["coverage_by_year"] = cov
    s["by_channel"] = bych
    if len(use) >= 15:
        ins = [r for r in use if not r["new"]]
        fits = {k: bcdf.fit(ins, k) for k in ("L", "Bk", "P")}
        s["yearly_corrected_median"] = yearly(use, bcdf.spatially_corrected(fits["L"], use))
        try:
            plot(os.path.join(args.out_dir, "cola_bcdf_vs_time.png"), use, fits)
        except Exception as ex:
            log("plot failed", repr(ex))
    try:
        plot_examples(os.path.join(args.out_dir, "cola_examples.png"), examples)
    except Exception as ex:
        log("example plot failed", repr(ex))
    s["finished"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    s["minutes"] = round((time.time() - T0) / 60.0, 1)
    with open(os.path.join(args.out_dir, "cola_summary.json"), "w") as f:
        json.dump(jsonable(s), f, indent=1)
    write_report(os.path.join(args.out_dir, "cola_report.md"), s)
    log(json.dumps(jsonable(dict(primary=s["primary"], secondary=s["secondary_consistency"])), indent=1)[:6000])
    log("done")


if __name__ == "__main__":
    main()

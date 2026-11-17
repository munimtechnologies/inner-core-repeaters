#!/usr/bin/env python3
"""BC-DF differential travel times of single SSI events at ILAR, 1990s to 2026.

Part of inner-core-repeaters by MV2, Munim, Inc. (made with AI agents). Licence: Apache-2.0.

A second observable that needs no repeaters. For every South Sandwich Islands (SSI) event
processed by run.py (1-2 Hz ILAR beams cached in <data-dir>/proc_v1), measure the time from
the PKIKP (DF, through the inner core) onset to PKPbc (BC, outer core only) by correlating the
DF pulse with the BC pulse. Its residual against ak135 changes when the inner core's
heterogeneity moves under the path, i.e. when the inner core rotates relative to the mantle
(Song & Richards 1996; Yang & Song 2023; Wu, Song & Yang 2026, GRL, doi:10.1029/2025GL119295,
which used single events of SSI sequences at COLA, 40 km from ILAR, to 2024).

Events sit at different places, so a linear spatial gradient (latitude, longitude, depth) is
fitted together with a time term. Two pre-registered time terms are compared (see
findings/2026-11-17-bcdf-prereg.md):
  L  (monotone):     g(t) = b * (t - 2008.5)
  Bk (backtracking): g(t) = b * (e(t) - 2008.5), e(t) = Wang et al. 2024's equivalent date
Both are fitted on events before 2025.0 and then predict the residuals of 2025-2026 events,
which no published study has used. A third term, P (pause: g = b * (min(t, 2008.5) - 2008.5),
change stops at the turn), is used by cola.py (findings/2026-11-17-cola-prereg.md).
"""
import argparse
import csv
import json
import math
import os
import time

import numpy as np

import model
from xcorr import best_lag, fshift, lag_scan, pick_onset

VERSION = "0.4.0"
FS = 40.0
PRE = 40.0
N0, N1 = int((PRE - 35) * FS), int((PRE - 5) * FS)  # noise window (beam is noise-normalised there)

# Pre-registered settings (findings/2026-11-17-bcdf-prereg.md)
MIN_DT_BC = 1.0        # s, ak135 BC-DF at the array centre; shorter separations are not measured
DF_PRE = 0.2           # s before the DF onset in the DF template
DF_LEN = (0.8, 1.6)    # s, template length = clip(dt_bc_pred - 0.4, 0.8, 1.6)
SEARCH = (-0.6, 1.5)   # s, BC searched at onset + dt_bc_pred + [-0.6, +1.5] ...
AB_GUARD = 0.8         # ... but never later than dt_ab_pred - 0.8 s
MIN_CC = 0.75          # beam DF/BC correlation for a usable measurement
MIN_SNR_DF = 3.0       # rms of the DF window over noise rms
EL_MIN_CC, EL_MIN_N = 0.6, 5
NEW_CUT = model.NEW_CUT


def log(*a):
    print(*a, flush=True)


def measure_bcdf(beam, dt_bc_pred, dt_ab_pred=float("nan"), fs=FS, pre=PRE, onset=None):
    """DF onset and BC-DF time from one noise-normalised beam (DF predicted at sample pre*fs).

    Returns dict(onset_s, dt_obs_s, cc, snr_df, L_s) or None. onset_s is relative to the
    ak135 DF prediction; dt_obs_s is the BC-DF time measured from the DF onset template.
    """
    beam = np.asarray(beam, dtype=float)
    if not (math.isfinite(dt_bc_pred) and dt_bc_pred >= MIN_DT_BC):
        return None
    c0 = int(round(pre * fs))
    if onset is None:
        onset = pick_onset(beam, fs, center=c0, search_s=4.0, sta_s=0.5, lta_s=8.0, thr=5.0)
        if onset is None:
            return None
    onset = int(onset)
    L = float(np.clip(dt_bc_pred - 0.4, DF_LEN[0], DF_LEN[1]))
    i0, i1 = onset - int(round(DF_PRE * fs)), onset + int(round(L * fs))
    if i0 < 0 or i1 >= len(beam):
        return None
    lo, hi = dt_bc_pred + SEARCH[0], dt_bc_pred + SEARCH[1]
    if math.isfinite(dt_ab_pred):
        hi = min(hi, dt_ab_pred - AB_GUARD)
    if hi - lo < 0.3:
        return None
    mid = 0.5 * (lo + hi)
    half = int(round(0.5 * (hi - lo) * fs))
    cs = int(round(mid * fs))
    y = fshift(beam, -cs)  # y[n] = beam[n + cs]
    L_, C_ = lag_scan(beam, y, i0, i1, half)
    lag, c = best_lag(L_, C_)
    if not math.isfinite(lag):
        return None
    nrms = float(np.sqrt(np.mean(beam[N0:N1] ** 2))) if len(beam) > N1 else 1.0
    snr = float(np.sqrt(np.mean(beam[i0:i1] ** 2)) / nrms) if nrms > 0 else float("nan")
    return dict(onset_s=(onset - c0) / fs, dt_obs_s=(cs + lag) / fs, cc=float(c), snr_df=snr, L_s=L,
                onset=onset)


def measure_elements(X, beam_res, dt_bc_pred, dt_ab_pred=float("nan"), fs=FS):
    """Element-by-element BC-DF with the beam's onset and a +/-0.3 s search around the beam value."""
    vals = []
    if beam_res is None:
        return float("nan"), 0
    on = beam_res["onset"]
    for x in np.asarray(X, dtype=float):
        nr = float(np.sqrt(np.mean(x[N0:N1] ** 2))) if len(x) > N1 else 0.0
        if not nr > 0:
            continue
        L = beam_res["L_s"]
        i0, i1 = on - int(round(DF_PRE * fs)), on + int(round(L * fs))
        cs = int(round(beam_res["dt_obs_s"] * fs))
        y = fshift(x, -cs)
        lag, c = best_lag(*lag_scan(x, y, i0, i1, int(round(0.3 * fs))))
        if math.isfinite(lag) and c >= EL_MIN_CC:
            vals.append((cs + lag) / fs)
    if len(vals) < EL_MIN_N:
        return float("nan"), len(vals)
    return float(np.median(vals)), len(vals)


# ----------------------------------------------------------------------------- fitting
def time_term(t, kind, T=model.T_TURN, ratio=model.RATE_RATIO):
    t = np.asarray(t, dtype=float)
    if kind == "L":
        return t - T
    if kind == "Bk":
        return np.array([model.equiv_date(x, T, ratio) for x in t]) - T
    if kind == "P":
        return np.minimum(t, T) - T
    if kind == "0":
        return np.zeros_like(t)
    raise ValueError(kind)


def design(rows, kind, ref=None, T=model.T_TURN, ratio=model.RATE_RATIO):
    """Columns: 1, dlat, dlon, ddepth/100, time term. ref = (lat0, lon0, z0) centring."""
    lat = np.array([r["lat"] for r in rows], dtype=float)
    lon = np.array([r["lon"] for r in rows], dtype=float)
    z = np.array([r["depth"] for r in rows], dtype=float)
    t = np.array([r["t"] for r in rows], dtype=float)
    if ref is None:
        ref = (float(np.mean(lat)), float(np.mean(lon)), float(np.mean(z)))
    g = time_term(t, kind, T, ratio)
    A = np.column_stack([np.ones_like(t), lat - ref[0], (lon - ref[1]) * math.cos(math.radians(ref[0])),
                         (z - ref[2]) / 100.0, g])
    return A, ref


def fit(rows, kind, key="res", ref=None, T=model.T_TURN, ratio=model.RATE_RATIO):
    A, ref = design(rows, kind, ref, T, ratio)
    y = np.array([r[key] for r in rows], dtype=float)
    if kind == "0":
        A = A[:, :4]
    coef, *_ = np.linalg.lstsq(A, y, rcond=None)
    if kind == "0":
        coef = np.append(coef, 0.0)
        A = np.column_stack([A, np.zeros(len(y))])
    pred = A @ coef
    rss = float(np.sum((y - pred) ** 2))
    return dict(coef=coef, ref=ref, rss=rss, n=len(y), kind=kind, T=T, ratio=ratio)


def predict(f, rows):
    A, _ = design(rows, f["kind"], f["ref"], f["T"], f["ratio"])
    return A @ f["coef"]


def spatially_corrected(f, rows, key="res"):
    """Residual minus the fitted spatial terms (keeps the constant and the time signal)."""
    A, _ = design(rows, "0", f["ref"])
    c = f["coef"]
    return np.array([r[key] for r in rows]) - A[:, 1:4] @ c[1:4]


def out_of_sample(ins, outs, key="res", nboot=2000, seed=0):
    """Fit L and Bk on `ins`, predict `outs`. Pre-registered decision quantities."""
    res = dict(n_in=len(ins), n_out=len(outs), n_cross_in=sum(1 for r in ins if r["t"] > model.T_TURN))
    if len(ins) < 15 or len(outs) < 1:
        res.update(testable=False, reason="need >= 15 in-sample and >= 1 new events")
        return res
    fL, fB, f0 = fit(ins, "L", key), fit(ins, "Bk", key), fit(ins, "0", key)
    y = np.array([r[key] for r in outs])
    eL, eB, e0 = y - predict(fL, outs), y - predict(fB, outs), y - predict(f0, outs)
    # in-sample scatter around the better model: per-event measurement + location noise
    sig = math.sqrt(min(fL["rss"], fB["rss"]) / max(1, len(ins) - 5))
    gap = float(np.mean(predict(fB, outs) - predict(fL, outs)))
    se_new = sig / math.sqrt(len(outs))
    rng = np.random.default_rng(seed)
    wins = []
    for _ in range(nboot):
        k = rng.integers(0, len(outs), len(outs))
        wins.append(float(np.mean(eB[k] ** 2) < np.mean(eL[k] ** 2)))
    decisive = len(outs) >= 8 and abs(gap) >= 2 * se_new
    frac = float(np.mean(wins))
    res.update(
        coef_L=dict(zip(["c", "dlat", "dlon", "ddepth100", "b"], map(float, fL["coef"]))),
        coef_Bk=dict(zip(["c", "dlat", "dlon", "ddepth100", "b"], map(float, fB["coef"]))),
        rss_in_L=fL["rss"], rss_in_Bk=fB["rss"], rss_in_0=f0["rss"], sigma_in=sig,
        mse_out_L=float(np.mean(eL ** 2)), mse_out_Bk=float(np.mean(eB ** 2)), mse_out_0=float(np.mean(e0 ** 2)),
        mean_err_out_L=float(np.mean(eL)), mean_err_out_Bk=float(np.mean(eB)),
        pred_gap_Bk_minus_L=gap, se_new_mean=se_new, decisive=bool(decisive),
        frac_boot_Bk_better=frac,
        verdict=("Bk" if frac >= 0.9 else ("L" if frac <= 0.1 else "neither")) if decisive else "not decisive",
        testable=True)
    return res


def scan_turn(rows, key="res", Ts=None, ratios=(1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0)):
    """Exploratory: (T, ratio) that minimises the in-sample RSS of the Bk-type model."""
    if Ts is None:
        Ts = np.arange(1998.0, 2022.01, 0.5)
    best = None
    for T in Ts:
        for q in ratios:
            f = fit(rows, "Bk", key, T=float(T), ratio=float(q))
            if best is None or f["rss"] < best[0]:
                best = (f["rss"], float(T), float(q), float(f["coef"][4]))
    fL = fit(rows, "L", key)
    return dict(best_T=best[1], best_ratio=best[2], best_b=best[3], rss_best=best[0], rss_L=fL["rss"],
                n=len(rows)) if best else {}


def yearly(rows, vals):
    out = {}
    for r, v in zip(rows, vals):
        out.setdefault(int(r["t"]), []).append(v)
    return {y: dict(n=len(v), median=float(np.median(v)), mad=float(np.median(np.abs(np.array(v) - np.median(v)))))
            for y, v in sorted(out.items())}


# ----------------------------------------------------------------------------- main
def read_catalog(path):
    evs = {}
    with open(path, newline="") as f:
        for r in csv.DictReader(f):
            try:
                evs[r["id"]] = dict(id=r["id"], time=r["time"], lat=float(r["latitude"]), lon=float(r["longitude"]),
                                    depth=float(r["depth"]) if r.get("depth") else 10.0,
                                    mag=float(r["mag"]) if r.get("mag") else float("nan"),
                                    t=model.decimal_year(r["time"]))
            except (KeyError, ValueError):
                continue
    return evs


def plot(path, rows, corr, fits):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    t = np.array([r["t"] for r in rows])
    new = t >= NEW_CUT
    fig, ax = plt.subplots(figsize=(10, 4.5))
    ax.scatter(t[~new], corr[~new], s=10, c="C0", label="events < 2025 (fit)")
    ax.scatter(t[new], corr[new], s=18, c="C3", label="events >= 2025 (out of sample)")
    tt = np.linspace(min(t.min(), 1995), max(t.max(), 2027), 400)
    for kind, col in (("L", "k"), ("Bk", "C2")):
        f = fits.get(kind)
        if f is None:
            continue
        c = f["coef"]
        ax.plot(tt, c[0] + c[4] * time_term(tt, kind, f["T"], f["ratio"]), color=col, lw=1.2,
                label="%s fit (< 2025)" % kind)
    ax.set_xlabel("year")
    ax.set_ylabel("BC-DF residual minus spatial terms (s)")
    ax.set_title("SSI -> ILAR single-event BC-DF residuals (1-2 Hz beams)")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=110)
    plt.close(fig)


def jsonable(o):
    if isinstance(o, dict):
        return {str(k): jsonable(v) for k, v in o.items()}
    if isinstance(o, (list, tuple, np.ndarray)):
        return [jsonable(v) for v in o]
    if isinstance(o, (np.floating, float)):
        return None if not math.isfinite(float(o)) else float(o)
    if isinstance(o, np.integer):
        return int(o)
    if isinstance(o, np.bool_):
        return bool(o)
    return o


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", default="/data")
    ap.add_argument("--out-dir", default="/out")
    ap.add_argument("--proc", default="proc_v1")
    args = ap.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)
    t0 = time.time()
    cat_path = os.path.join(args.data_dir, "catalog_usgs.csv")
    pdir = os.path.join(args.data_dir, args.proc)
    if not (os.path.exists(cat_path) and os.path.isdir(pdir)):
        log("no cached ILAR data (%s, %s): run run.py first with the same data_key" % (cat_path, pdir))
        with open(os.path.join(args.out_dir, "bcdf_summary.json"), "w") as f:
            json.dump(dict(error="no cached ILAR data"), f)
        return
    evs = read_catalog(cat_path)
    rows, status = [], []
    for fn in sorted(os.listdir(pdir)):
        if not fn.endswith(".npz"):
            continue
        eid = fn[:-4]
        e = evs.get(eid)
        if e is None:
            continue
        try:
            z = np.load(os.path.join(pdir, fn), allow_pickle=False)
            beam, X, meta = z["beam"], z["X"], z["meta"]
        except Exception as ex:
            status.append(dict(id=eid, t=e["t"], status="load-error " + repr(ex)[:80]))
            continue
        dist, dt_bc, dt_ab = float(meta[0]), float(meta[1]), float(meta[2])
        r = measure_bcdf(beam, dt_bc, dt_ab)
        if r is None:
            status.append(dict(id=eid, t=e["t"], dist=dist, status="no-measurement"))
            continue
        el, nel = measure_elements(X, r, dt_bc, dt_ab)
        ok = r["cc"] >= MIN_CC and r["snr_df"] >= MIN_SNR_DF and abs(r["onset_s"]) <= 3.0
        row = dict(id=eid, time=e["time"], t=e["t"], lat=e["lat"], lon=e["lon"], depth=e["depth"], mag=e["mag"],
                   dist=dist, dt_bc_pred=dt_bc, onset_s=r["onset_s"], dt_obs=r["dt_obs_s"], res=r["dt_obs_s"] - dt_bc,
                   cc=r["cc"], snr_df=r["snr_df"], dt_obs_el=el, res_el=el - dt_bc if math.isfinite(el) else el,
                   n_el=nel, ok=bool(ok), new=e["t"] >= NEW_CUT)
        status.append(dict(id=eid, t=e["t"], dist=dist, status="ok" if ok else "low-quality"))
        rows.append(row)
    keys = ["id", "time", "t", "lat", "lon", "depth", "mag", "dist", "dt_bc_pred", "onset_s", "dt_obs", "res", "cc",
            "snr_df", "dt_obs_el", "res_el", "n_el", "ok", "new"]
    with open(os.path.join(args.out_dir, "bcdf_events.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)
    ok = [r for r in rows if r["ok"]]
    # robust outlier screen, pre-registered: drop |res - median| > 1.0 s (cycle skips)
    if ok:
        med = float(np.median([r["res"] for r in ok]))
        ok = [r for r in ok if abs(r["res"] - med) <= 1.0]
    ins = [r for r in ok if not r["new"]]
    outs = [r for r in ok if r["new"]]
    log("measured %d, usable %d (in-sample %d, new %d)" % (len(rows), len(ok), len(ins), len(outs)))
    summary = dict(version=VERSION, n_npz=len(status), n_measured=len(rows), n_ok=len(ok), n_in=len(ins),
                   n_new=len(outs))
    summary["Q1_in_sample"] = {}
    if len(ins) >= 15:
        fL, fB = fit(ins, "L"), fit(ins, "Bk")
        summary["Q1_in_sample"] = dict(rss_L=fL["rss"], rss_Bk=fB["rss"], b_L=float(fL["coef"][4]),
                                       b_Bk=float(fB["coef"][4]), n=len(ins),
                                       n_after_turn=sum(1 for r in ins if r["t"] > model.T_TURN),
                                       better="Bk" if fB["rss"] < fL["rss"] else "L")
        summary["Q3_turn_scan_in_sample"] = scan_turn(ins)
        corr = spatially_corrected(fL, ok)
        summary["yearly_corrected_median"] = yearly(ok, corr)
        try:
            plot(os.path.join(args.out_dir, "bcdf_vs_time.png"), ok, corr, dict(L=fL, Bk=fB))
        except Exception as ex:
            log("plot failed", repr(ex))
    summary["Q2_out_of_sample"] = out_of_sample(ins, outs)
    el_rows = [dict(r, res=r["res_el"]) for r in ok if math.isfinite(r["res_el"])]
    summary["robust_element_median"] = out_of_sample([r for r in el_rows if not r["new"]],
                                                     [r for r in el_rows if r["new"]])
    summary["minutes"] = round((time.time() - t0) / 60, 2)
    with open(os.path.join(args.out_dir, "bcdf_summary.json"), "w") as f:
        json.dump(jsonable(summary), f, indent=1)
    log(json.dumps(jsonable(summary), indent=1)[:6000])


if __name__ == "__main__":
    main()

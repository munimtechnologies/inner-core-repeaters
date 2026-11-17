#!/usr/bin/env python3
"""m0324_check: is the 2025-12-13 anomaly in multiplet M03-24 real?

In icr-m0324-v1 (findings/2026-11-17-m0324-results.md) the first gated 2025 member of the
M03-24 repeating-earthquake family (USGS us6000runq, 2025-12-13, mb 5.0) has a PKIKP (DF)
waveform at ILAR that is almost identical to the 2011-02-20 member (D = 0.033) and clearly
different from the 2003-03-21 member (D = 0.167). The backtracking model of Wang et al.
(2024) predicts the opposite: a late-2025 event should resemble 2001-2002, i.e. the 2003
member most. This script runs the pre-registered checks of
findings/2026-11-17-m0324-check-prereg.md:

  R0  reproduce the published D values of multiplet Z004 from the cached ILAR beams;
  C1  element bootstrap and leave-one-element-out at ILAR;
  C2  three alternative DF windows at ILAR (early, late, long DF + coda);
  Y0/Y1  the same family at a second, independent array (YKA, Canada) with the
         pre-registered YKA coda metric of yka.py (v0.2.0, downloaded from the repo at a
         pinned commit).

It reads the ILAR beams that icr-m0324-v1 cached in /data/proc_m0324_v1 (data_key icr-m0324)
and downloads YKA records for the nine Z004 events only.

By MV2, Munim, Inc. (made with AI agents). Licence: Apache-2.0.
"""
import argparse
import csv
import json
import math
import os
import sys
import time
import traceback
from itertools import combinations
from types import SimpleNamespace

import numpy as np

import model
from cluster import max_spanning_offsets
from xcorr import fshift, measure_pair, pick_onset

VERSION = "0.8.0"
FS = 40.0
PRE = 40.0
POST = 40.0
N = int(round((PRE + POST) * FS))
N0, N1 = int((PRE - 35) * FS), int((PRE - 5) * FS)
PROC_TAG = "proc_m0324_v1"  # ILAR beams written by icr-m0324-v1 (m0324.py 0.5.0, run.py save format)

NEW = "us6000runq"    # 2025-12-13
E2003 = "usp000bt5x"  # 2003-03-21
E2008 = "usp000fz4d"  # 2008-02-17
E2011 = "usp000hute"  # 2011-02-20
E2020 = "us7000am0e"  # 2020-07-11
E2024 = "us6000mk65"  # 2024-03-21
KEY_A = (E2003, NEW)  # backtracking says: most similar of the pre-turn partners
KEY_B = (E2011, NEW)  # observed: most similar

ALT_WINDOWS = {"W_early": (-0.3, 1.2), "W_late": (1.0, 2.5), "W_long": (-0.3, 6.0)}
YKA_COMMIT = "314a42d1b609d30d35a5aa51e001cf48b11948c6"
YKA_URL = ("https://raw.githubusercontent.com/munimtechnologies/inner-core-repeaters/%s/yka.py"
           % YKA_COMMIT)
T0 = time.time()


def log(*a):
    print("[%7.1f min]" % ((time.time() - T0) / 60.0), *a, flush=True)


# ----------------------------------------------------------------------------- inputs
def read_events(path):
    evs = {}
    with open(path, newline="") as f:
        for r in csv.DictReader(f):
            evs[r["id"]] = dict(id=r["id"], time=r["time"], lat=float(r["lat"]), lon=float(r["lon"]),
                                depth=float(r["depth"]), mag=float(r["mag"]), t=float(r["t"]))
    return evs


def read_edges(path, fs=FS):
    """{(ev1, ev2): (cc, lag in samples)} from an edges csv with lag in seconds."""
    edges = {}
    with open(path, newline="") as f:
        for r in csv.DictReader(f):
            edges[(r["ev1"], r["ev2"])] = (float(r["cc"]), float(r["lag_s"]) * fs)
    return edges


def read_published(path):
    out = {}
    with open(path, newline="") as f:
        for r in csv.DictReader(f):
            out[(r["ev1"], r["ev2"])] = dict(D=float(r["D"]), ok=str(r["ok"]) == "True")
    return out


def load_ilar(path):
    """Same format as run.load_proc (icr-m0324-v1)."""
    z = np.load(path, allow_pickle=False)
    m = z["meta"]
    return dict(names=np.array([str(x) for x in z["names"]]), X=np.asarray(z["X"], dtype=np.float64),
                beam=np.asarray(z["beam"], dtype=np.float64), dist=float(m[0]), dt_bc=float(m[1]),
                dt_ab=float(m[2]), snr_group=float(m[3]))


# ----------------------------------------------------------------------------- ILAR geometry
def ilar_frame(members, data, edges):
    """Alignment, DF onset and windows exactly as run.analyze_multiplet (inner-core-repeaters v0.1)."""
    root = max(members, key=lambda m: data[m]["snr_group"])
    off = max_spanning_offsets(members, edges, root)
    members = [m for m in members if m in off]
    stack = np.zeros(N)
    for m in members:
        b = data[m]["beam"]
        stack += fshift(b, -off[m]) / max(np.sqrt(np.mean(b[N0:N1] ** 2)), 1e-12)
    stack /= len(members)
    trig = pick_onset(stack, FS, center=int(PRE * FS), search_s=6.0, sta_s=0.5, lta_s=8.0, thr=6.0)
    if trig is None:
        return None
    o = trig
    dt_bc, dt_ab = data[root]["dt_bc"], data[root]["dt_ab"]
    df_end = min(dt_bc - 0.35, 2.5) if math.isfinite(dt_bc) else 2.5
    dfw = (o - int(0.3 * FS), o + int(df_end * FS))
    ref_start = (dt_bc - 0.35) if math.isfinite(dt_bc) else (dt_ab - 0.5)
    refw = (o + int(ref_start * FS), o + int((dt_ab + 8.0) * FS))
    return dict(root=root, off=off, members=members, onset=o, dfw=dfw, refw=refw, df_end=df_end,
                pick_minus_pred_s=(o - PRE * FS) / FS)


def window(onset, w, fs=FS):
    return (onset + int(round(w[0] * fs)), onset + int(round(w[1] * fs)))


def pair_D(xa, xb, off_a, off_b, refw, dfw, n=N):
    """D = 1 - zero-lag DF correlation after reference alignment (run.analyze_multiplet)."""
    sa = int(round(off_a))
    ref = (refw[0] + sa, refw[1] + sa)
    dwa = (dfw[0] + sa, dfw[1] + sa)
    if ref[0] < 0 or ref[1] > n or dwa[0] < 0:
        return float("nan"), float("nan")
    r = measure_pair(xa, xb, ref, dwa, FS, max_lag_ref_s=0.5, max_lag_df_s=0.25, lag0=off_b - off_a)
    if r is None:
        return float("nan"), float("nan")
    return 1.0 - r["cc_df0"], r["cc_ref"]


def ordered_pairs(members, evs):
    order = sorted(members, key=lambda m: evs[m]["t"])
    return list(combinations(order, 2))


def all_D(beams, frame, evs, dfw=None):
    dfw = dfw or frame["dfw"]
    out = {}
    for a, b in ordered_pairs(frame["members"], evs):
        if beams.get(a) is None or beams.get(b) is None:
            out[(a, b)] = float("nan")
            continue
        out[(a, b)] = pair_D(beams[a], beams[b], frame["off"][a], frame["off"][b], frame["refw"], dfw)[0]
    return out


# ----------------------------------------------------------------------------- element resampling
def common_elements(data, members):
    s = None
    for m in members:
        names = set(str(x) for x in data[m]["names"])
        s = names if s is None else (s & names)
    return sorted(s or [])


def beam_from(X, idx):
    """Mean of the chosen element rows, normalised by the beam's own noise rms."""
    b = np.asarray(X, dtype=np.float64)[list(idx)].mean(axis=0)
    nb = math.sqrt(float(np.mean(b[N0:N1] ** 2)))
    return b / nb if nb > 0 else None


def element_beams(data, members, names):
    out = {}
    for m in members:
        pos = {str(n): k for k, n in enumerate(data[m]["names"])}
        idx = [pos[n] for n in names if n in pos]
        out[m] = beam_from(data[m]["X"], idx) if len(idx) >= 3 else None
    return out


def element_bootstrap(data, frame, evs, nboot=500, seed=0):
    common = common_elements(data, frame["members"])
    rng = np.random.default_rng(seed)
    pairs = ordered_pairs(frame["members"], evs)
    reps = {p: [] for p in pairs}
    for _ in range(nboot):
        pick = [common[i] for i in rng.integers(0, len(common), len(common))]
        D = all_D(element_beams(data, frame["members"], pick), frame, evs)
        for p in pairs:
            reps[p].append(D[p])
    return common, {p: np.asarray(v, dtype=float) for p, v in reps.items()}


def leave_one_out(data, frame, evs, common):
    out = {}
    for el in common:
        names = [n for n in common if n != el]
        out[el] = all_D(element_beams(data, frame["members"], names), frame, evs)
    return out


def frac_greater(x, y):
    """Fraction of paired replicates with x > y (finite pairs only)."""
    x, y = np.asarray(x, dtype=float), np.asarray(y, dtype=float)
    m = np.isfinite(x) & np.isfinite(y)
    return float((x[m] > y[m]).mean()) if m.any() else float("nan")


def grade_c1(boot, loo, a=KEY_A, b=KEY_B, frac_bar=0.95):
    """C1: D(a) > D(b) in >= 95% of element bootstraps AND in every leave-one-element-out case."""
    f = frac_greater(boot[a], boot[b])
    diffs = [d[a] - d[b] for d in loo.values() if math.isfinite(d[a]) and math.isfinite(d[b])]
    loo_all = bool(diffs) and all(x > 0 for x in diffs)
    testable = math.isfinite(f) and len(diffs) >= 5
    return dict(frac_boot=f, n_boot=int(np.isfinite(boot[a]).sum()), n_loo=len(diffs),
                loo_min_diff=float(min(diffs)) if diffs else float("nan"),
                loo_max_diff=float(max(diffs)) if diffs else float("nan"),
                testable=bool(testable), robust=bool(testable and f >= frac_bar and loo_all) if testable else None)


def summarize(v):
    v = np.asarray(v, dtype=float)
    v = v[np.isfinite(v)]
    if not len(v):
        return dict(n=0, median=float("nan"), p05=float("nan"), p95=float("nan"))
    return dict(n=int(len(v)), median=float(np.median(v)), p05=float(np.percentile(v, 5)),
                p95=float(np.percentile(v, 95)))


# ----------------------------------------------------------------------------- YKA
def load_yka_module(cache_dir):
    """yka.py v0.2.0 from the public repo at a pinned commit (needs model, phases, xcorr here)."""
    p = os.path.join(cache_dir, "yka.py")
    if not os.path.exists(p):
        import requests

        os.makedirs(cache_dir, exist_ok=True)
        for attempt in range(5):
            try:
                r = requests.get(YKA_URL, timeout=60)
                if r.status_code == 200 and "def analyze_multiplet" in r.text and 'VERSION = "0.2.0"' in r.text:
                    with open(p, "w") as f:
                        f.write(r.text)
                    break
            except Exception:
                pass
            time.sleep(10 * (attempt + 1))
    if not os.path.exists(p):
        return None
    if cache_dir not in sys.path:
        sys.path.append(cache_dir)
    import yka  # noqa: E402

    return yka


def yka_key_values(rows, a, b):
    for r in rows:
        if (r["ev1"], r["ev2"]) == (a, b):
            return r
    return None


def grade_y(rows):
    """Y1 (replication) and Y0 (in-sample pattern) on YKA rows; ungated values always reported."""
    ra, rb = yka_key_values(rows, *KEY_A), yka_key_values(rows, *KEY_B)
    y1 = dict(D_2003_new=ra["D"] if ra else None, D_2011_new=rb["D"] if rb else None,
              ok_2003_new=bool(ra and ra["ok"]), ok_2011_new=bool(rb and rb["ok"]))
    y1["testable"] = bool(ra and rb and ra["ok"] and rb["ok"])
    y1["replicated"] = (ra["D"] > rb["D"]) if y1["testable"] else None
    if ra and rb and math.isfinite(ra["D"]) and math.isfinite(rb["D"]):
        y1["ungated_same_direction"] = bool(ra["D"] > rb["D"])
    r311 = yka_key_values(rows, E2003, E2011)
    r320 = yka_key_values(rows, E2003, E2020)
    r324 = yka_key_values(rows, E2003, E2024)
    trio = [r311, r320, r324]
    y0 = dict(D_2003_2011=r311["D"] if r311 else None, D_2003_2020=r320["D"] if r320 else None,
              D_2003_2024=r324["D"] if r324 else None)
    y0["testable"] = all(r is not None and r["ok"] for r in trio)
    y0["holds"] = (r311["D"] > max(r320["D"], r324["D"])) if y0["testable"] else None
    return dict(Y1=y1, Y0=y0)


def run_yka(evs, members, data_dir, out_dir, plots_dir, s):
    yka = load_yka_module("/tmp/icr_mods")
    if yka is None:
        s["yka"] = dict(status="yka.py download failed")
        return
    import requests

    ph = yka.Phases(yka.load_or_build(os.path.join(data_dir, "ak135_yka_grid_v1.npz")))
    proc_dir = os.path.join(data_dir, yka.PROC_TAG)
    os.makedirs(proc_dir, exist_ok=True)
    session = requests.Session()
    session.headers["User-Agent"] = "inner-core-repeaters/%s (MV2, Munim, Inc.; mv2@munimtech.com)" % VERSION
    data, status = {}, []
    for m in sorted(members, key=lambda i: evs[i]["t"]):
        e = evs[m]
        pp = os.path.join(proc_dir, m + ".npz")
        try:
            if os.path.exists(pp):
                data[m] = yka.load_proc(pp)
                status.append(dict(id=m, time=e["time"], status="ok-cached", gen=data[m]["gen"],
                                   snr=round(data[m]["snr"], 2), n_sta=len(data[m]["names"])))
                continue
            st, why = yka.fetch_event(session, e, ph, data_dir)
            if st is None:
                status.append(dict(id=m, time=e["time"], status=why))
                continue
            d, why2 = yka.process_event(st, e, ph)
            if d is None:
                status.append(dict(id=m, time=e["time"], status=why2))
                continue
            yka.save_proc(pp, d)
            data[m] = dict((k, v) for k, v in d.items() if k != "X")
            status.append(dict(id=m, time=e["time"], status="ok", gen=d["gen"], snr=round(d["snr"], 2),
                               n_sta=len(d["names"]), source=why))
        except Exception as ex:
            status.append(dict(id=m, time=e["time"], status="exception:" + repr(ex)[:150]))
        time.sleep(0.2)
    write_csv(os.path.join(out_dir, "yka_z004_status.csv"), status,
              keys=["id", "time", "status", "gen", "snr", "n_sta", "source"])
    args = SimpleNamespace(pair_thr=0.90, min_snr_direct=5.0, min_snr_coda=2.0, max_station_level=40)
    rows, info, fr = yka.analyze_multiplet("Z004", list(data), evs, data, proc_dir, args)
    write_csv(os.path.join(out_dir, "yka_z004_pairs.csv"), rows)
    s["yka"] = dict(status="ran", n_processed=len(data), multiplet=info, n_pairs=len(rows),
                    n_ok=sum(1 for r in rows if r["ok"]), tests=grade_y(rows))
    if fr is not None:
        off, flips, onset, mem = fr
        try:
            yka.plot_multiplet(os.path.join(plots_dir, "yka_z004.png"), "Z004", mem, evs, data, off, flips, onset)
        except Exception:
            traceback.print_exc()


# ----------------------------------------------------------------------------- output helpers
def write_csv(path, rows, keys=None):
    if not rows:
        open(path, "w").close()
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
    if isinstance(o, np.bool_):
        return bool(o)
    if isinstance(o, np.ndarray):
        return jsonable(o.tolist())
    return o


def plot_matrix(path, members, evs, D, title):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    order = sorted(members, key=lambda m: evs[m]["t"])
    n = len(order)
    M = np.full((n, n), np.nan)
    for i, a in enumerate(order):
        M[i, i] = 0.0
        for j, b in enumerate(order):
            if (a, b) in D:
                M[i, j] = M[j, i] = D[(a, b)]
    fig, ax = plt.subplots(figsize=(7.5, 6.5))
    im = ax.imshow(M, cmap="viridis_r", vmin=0, vmax=max(0.25, float(np.nanmax(M))))
    labels = [evs[m]["time"][:10] for m in order]
    ax.set_xticks(range(n))
    ax.set_xticklabels(labels, rotation=60, ha="right", fontsize=8)
    ax.set_yticks(range(n))
    ax.set_yticklabels(labels, fontsize=8)
    for i in range(n):
        for j in range(n):
            if i != j and math.isfinite(M[i, j]):
                ax.text(j, i, "%.2f" % M[i, j], ha="center", va="center", fontsize=6.5,
                        color="white" if M[i, j] > 0.12 else "black")
    fig.colorbar(im, ax=ax, label="D = 1 - cc (PKIKP window, after PKPbc+ab alignment)")
    ax.set_title(title, fontsize=9)
    fig.tight_layout()
    fig.savefig(path, dpi=110)
    plt.close(fig)


def plot_boot(path, boot):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(7.5, 4))
    for key, lab in ((KEY_A, "2025-12-13 vs 2003-03-21"), (KEY_B, "2025-12-13 vs 2011-02-20"),
                     ((E2024, NEW), "2025-12-13 vs 2024-03-21"), ((E2003, E2020), "2003-03-21 vs 2020-07-11")):
        v = boot.get(key)
        if v is None:
            continue
        v = v[np.isfinite(v)]
        if len(v):
            ax.hist(v, bins=40, alpha=0.5, label=lab)
    ax.set_xlabel("D at ILAR, element bootstrap (resampling the array's elements)")
    ax.set_ylabel("replicates")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=110)
    plt.close(fig)


# ----------------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", default="/data")
    ap.add_argument("--out-dir", default="/out")
    ap.add_argument("--nboot", type=int, default=500)
    ap.add_argument("--no-yka", action="store_true")
    args = ap.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)
    plots = os.path.join(args.out_dir, "plots")
    os.makedirs(plots, exist_ok=True)
    here = os.path.dirname(os.path.abspath(__file__))
    evs = read_events(os.path.join(here, "z004_events.csv"))
    edges = read_edges(os.path.join(here, "z004_edges.csv"))
    pub = read_published(os.path.join(here, "z004_published_D.csv"))
    s = dict(version=VERSION, started=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()), args=vars(args))

    # ---- ILAR
    proc = os.path.join(args.data_dir, PROC_TAG)
    data, missing = {}, []
    for m in evs:
        p = os.path.join(proc, m + ".npz")
        if os.path.exists(p):
            data[m] = load_ilar(p)
        else:
            missing.append(m)
    s["ilar_missing"] = missing
    log("ILAR beams loaded:", len(data), "missing:", missing)
    rows = []
    if len(data) == len(evs):
        try:
            frame = ilar_frame(list(evs), data, edges)
            s["ilar_frame"] = dict(root=frame["root"], onset=frame["onset"], dfw=frame["dfw"], refw=frame["refw"],
                                   pick_minus_pred_s=frame["pick_minus_pred_s"], n_members=len(frame["members"]))
            beams = {m: data[m]["beam"] for m in frame["members"]}
            D0 = all_D(beams, frame, evs)
            diffs = [abs(D0[p] - pub[p]["D"]) for p in D0 if p in pub and math.isfinite(D0[p])]
            s["R0_reproduction"] = dict(n=len(diffs), max_abs_diff=max(diffs) if diffs else None,
                                        reproduced=bool(diffs) and max(diffs) < 1e-3 and len(diffs) == len(pub))
            log("R0:", json.dumps(jsonable(s["R0_reproduction"])))
            common = common_elements(data, frame["members"])
            Dc = all_D(element_beams(data, frame["members"], common), frame, evs)
            alt = {k: all_D(beams, frame, evs, dfw=window(frame["onset"], w)) for k, w in ALT_WINDOWS.items()}
            log("bootstrap over", len(common), "common elements,", args.nboot, "replicates")
            common, boot = element_bootstrap(data, frame, evs, nboot=args.nboot)
            loo = leave_one_out(data, frame, evs, common)
            s["common_elements"] = common
            s["C1_element_bootstrap"] = grade_c1(boot, loo)
            s["C1_in_sample_contrast"] = dict(
                frac_2003_2011_gt_2003_2020=frac_greater(boot[(E2003, E2011)], boot[(E2003, E2020)]),
                frac_2003_2011_gt_2003_2024=frac_greater(boot[(E2003, E2011)], boot[(E2003, E2024)]),
                frac_2003_new_gt_2024_new=frac_greater(boot[KEY_A], boot[(E2024, NEW)]),
                frac_2003_new_gt_2008_new=frac_greater(boot[KEY_A], boot[(E2008, NEW)]))
            s["C2_windows"] = {k: dict(D_2003_new=v[KEY_A], D_2011_new=v[KEY_B],
                                       same_direction=bool(v[KEY_A] > v[KEY_B]))
                               for k, v in alt.items()}
            s["C2_windows_all_same_direction"] = all(v["same_direction"] for v in s["C2_windows"].values())
            for a, b in ordered_pairs(frame["members"], evs):
                p = (a, b)
                lv = [d[p] for d in loo.values() if math.isfinite(d[p])]
                bs = summarize(boot[p])
                rows.append(dict(ev1=a, ev2=b, date1=evs[a]["time"][:10], date2=evs[b]["time"][:10],
                                 D_published=pub.get(p, {}).get("D"), ok_published=pub.get(p, {}).get("ok"),
                                 D_reproduced=D0[p], D_common_elements=Dc[p], boot_median=bs["median"],
                                 boot_p05=bs["p05"], boot_p95=bs["p95"],
                                 loo_min=min(lv) if lv else float("nan"), loo_max=max(lv) if lv else float("nan"),
                                 D_W_early=alt["W_early"][p], D_W_late=alt["W_late"][p], D_W_long=alt["W_long"][p],
                                 dM=model.delta_M(evs[a]["t"], evs[b]["t"]), dB=model.delta_B(evs[a]["t"], evs[b]["t"])))
            write_csv(os.path.join(args.out_dir, "ilar_z004_matrix.csv"), rows)
            try:
                plot_matrix(os.path.join(plots, "ilar_z004_D_matrix.png"), frame["members"], evs, D0,
                            "M03-24 (Z004) at ILAR: PKIKP dissimilarity D, icr-m0324-v1 beams")
                plot_matrix(os.path.join(plots, "ilar_z004_D_matrix_W_long.png"), frame["members"], evs,
                            alt["W_long"], "M03-24 (Z004) at ILAR: D in the long DF+coda window (-0.3..6 s)")
                plot_boot(os.path.join(plots, "ilar_z004_bootstrap.png"), boot)
            except Exception:
                traceback.print_exc()
        except Exception as ex:
            traceback.print_exc()
            s["ilar_error"] = repr(ex)[:300]
    else:
        s["ilar_error"] = "missing cached ILAR beams; nothing measured at ILAR"

    # ---- YKA
    if not args.no_yka:
        try:
            run_yka(evs, list(evs), args.data_dir, args.out_dir, plots, s)
        except Exception as ex:
            traceback.print_exc()
            s["yka"] = dict(status="exception", error=repr(ex)[:300])
    s["finished"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    s["minutes"] = round((time.time() - T0) / 60.0, 1)
    with open(os.path.join(args.out_dir, "m0324_check_summary.json"), "w") as f:
        json.dump(jsonable(s), f, indent=1)
    log("done", json.dumps(jsonable({k: s.get(k) for k in ("R0_reproduction", "C1_element_bootstrap",
                                                           "C2_windows_all_same_direction")})))
    if "yka" in s:
        log("YKA:", json.dumps(jsonable(s["yka"]))[:2000])


if __name__ == "__main__":
    main()

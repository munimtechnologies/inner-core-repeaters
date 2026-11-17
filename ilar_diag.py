#!/usr/bin/env python3
"""Robustness diagnostics for the icr-ilar-v1 pair table (results/ilar-v1/pairs.csv).

Part of inner-core-repeaters by MV2, Munim, Inc. (made with AI agents). Licence: Apache-2.0.

These checks were written AFTER predictions.json of icr-ilar-v1 had been read, so they are
exploratory robustness checks, not pre-registered tests. They ask whether the in-sample
preference for the backtracking ordering (P1: rho_B = 0.19, p = 0.03; multiplet bootstrap
88%, below the pre-registered 90%) survives obvious confounds:

  E  ILAR instrument epoch: the digitiser/response changed on 2018-09-25. Same-epoch vs
     cross-epoch pairs, and P1 restricted to pairs with both events before the break.
  S  signal-to-noise: Spearman(D, min DF SNR), and partial rank correlations of D with
     delta_B and delta_M after removing the SNR rank.
  J  leave-one-multiplet-out: range of rho_B - rho_M; multiplets whose removal flips it.
  Q  without the 2021 South Sandwich sequence (events 2021.6-2022.0), which supplies many
     of the short-separation pairs.
  W  within-multiplet sign test: in each multiplet with >= 3 usable pairs, is Kendall's tau
     of D with delta_B larger than with delta_M?
  T  turn-date bootstrap: the P5 grid search repeated on 300 multiplet bootstrap samples,
     giving a spread for the best turn date T and the rho surface at T = 2008.5, ratio 2.5.

Usage: python ilar_diag.py --pairs <pairs.csv or URL> --out-dir /out
"""
import argparse
import csv
import json
import math
import os
import time
import urllib.request

import numpy as np

import model

PAIRS_URL = ("https://raw.githubusercontent.com/munimtechnologies/inner-core-repeaters/main/"
             "results/ilar-v1/pairs.csv")
EPOCH_2018 = 2018.73  # 2018-09-25
SEQ_2021 = (2021.6, 2022.0)


def load_pairs(src, out_dir, retries=12, wait_s=60):
    if os.path.exists(src):
        path = src
    else:
        path = os.path.join(out_dir, "pairs_input.csv")
        err = None
        for k in range(retries):
            try:
                with urllib.request.urlopen(src, timeout=60) as r, open(path, "wb") as f:
                    f.write(r.read())
                err = None
                break
            except Exception as ex:  # the publish commit may land after the job starts
                err = ex
                time.sleep(wait_s)
        if err is not None:
            raise SystemExit("could not fetch %s: %r" % (src, err))
    rows = []
    with open(path, newline="") as f:
        for r in csv.DictReader(f):
            try:
                d = dict(mult=r["mult"], ev1=r["ev1"], ev2=r["ev2"], t1=float(r["t1"]), t2=float(r["t2"]),
                         D=float(r["D"]), Dx=float(r["Dx"]), cc_ref=float(r["cc_ref"]),
                         snr=min(float(r["snr_df1"]), float(r["snr_df2"])),
                         ok=str(r["ok"]).strip() == "True")
            except (KeyError, ValueError):
                continue
            if all(math.isfinite(d[k]) for k in ("t1", "t2", "D", "cc_ref", "snr")):
                rows.append(d)
    return rows


def ranks(x):
    from scipy.stats import rankdata

    return rankdata(np.asarray(x, dtype=float))


def partial_spearman(x, y, z):
    """Spearman correlation of x and y after removing a linear dependence of both ranks on z's rank."""
    rx, ry, rz = ranks(x), ranks(y), ranks(z)
    A = np.column_stack([np.ones_like(rz), rz])
    ex = rx - A @ np.linalg.lstsq(A, rx, rcond=None)[0]
    ey = ry - A @ np.linalg.lstsq(A, ry, rcond=None)[0]
    if np.std(ex) == 0 or np.std(ey) == 0:
        return float("nan")
    return float(np.corrcoef(ex, ey)[0, 1])


def rho_pair(rows, key="D"):
    (rb, pb, n), (rm, pm, _) = model._rhos(rows, key=key)
    return dict(n=n, rho_B=rb, p_B=pb, rho_M=rm, p_M=pm,
                diff=(rb - rm) if math.isfinite(rb) and math.isfinite(rm) else float("nan"))


def epoch_class(r):
    a, b = r["t1"] >= EPOCH_2018, r["t2"] >= EPOCH_2018
    return "post" if (a and b) else ("pre" if not (a or b) else "cross")


def diag_E(rows):
    out = {}
    for c in ("pre", "cross", "post"):
        rr = [r for r in rows if epoch_class(r) == c]
        out[c] = dict(n=len(rr), median_D=float(np.median([r["D"] for r in rr])) if rr else None,
                      median_cc_ref=float(np.median([r["cc_ref"] for r in rr])) if rr else None,
                      median_dM=float(np.median([model.delta_M(r["t1"], r["t2"]) for r in rr])) if rr else None)
    pre = [r for r in rows if epoch_class(r) == "pre"]
    out["P1_pre_epoch_only"] = model.grade_set(pre, nboot=500) if len(pre) >= 8 else dict(n=len(pre))
    return out


def diag_S(rows):
    D = [r["D"] for r in rows]
    snr = [r["snr"] for r in rows]
    dB = [model.delta_B(r["t1"], r["t2"]) for r in rows]
    dM = [model.delta_M(r["t1"], r["t2"]) for r in rows]
    return dict(rho_D_snr=model.spearman(snr, D)[0], rho_dB_snr=model.spearman(snr, dB)[0],
                rho_dM_snr=model.spearman(snr, dM)[0],
                partial_rho_B=partial_spearman(dB, D, snr), partial_rho_M=partial_spearman(dM, D, snr))


def diag_J(rows):
    full = rho_pair(rows)
    mults = sorted({r["mult"] for r in rows})
    vals, flips = [], []
    for m in mults:
        rr = [r for r in rows if r["mult"] != m]
        d = rho_pair(rr)["diff"]
        if math.isfinite(d):
            vals.append((m, d))
            if math.isfinite(full["diff"]) and (d > 0) != (full["diff"] > 0):
                flips.append(m)
    v = [d for _, d in vals]
    return dict(full=full, min_diff=min(v) if v else None, max_diff=max(v) if v else None,
                most_influential=sorted(vals, key=lambda x: abs(x[1] - full["diff"]))[-5:], flips=flips)


def diag_Q(rows):
    def in_seq(t):
        return SEQ_2021[0] <= t < SEQ_2021[1]
    rr = [r for r in rows if not (in_seq(r["t1"]) or in_seq(r["t2"]))]
    res = model.grade_set(rr, nboot=500) if len(rr) >= 8 else dict(n=len(rr))
    res["P4_without_2021"] = model.p4(rr)
    return res


def diag_W(rows, min_pairs=3):
    from scipy.stats import kendalltau

    by = {}
    for r in rows:
        by.setdefault(r["mult"], []).append(r)
    nB = nM = ties = 0
    per = []
    for m, rr in sorted(by.items()):
        if len(rr) < min_pairs:
            continue
        D = [r["D"] for r in rr]
        tb = kendalltau([model.delta_B(r["t1"], r["t2"]) for r in rr], D)[0]
        tm = kendalltau([model.delta_M(r["t1"], r["t2"]) for r in rr], D)[0]
        if not (math.isfinite(tb) and math.isfinite(tm)) or abs(tb - tm) < 1e-9:
            ties += 1
            continue
        if tb > tm:
            nB += 1
        else:
            nM += 1
        per.append(dict(mult=m, n=len(rr), tau_B=float(tb), tau_M=float(tm)))
    from scipy.stats import binomtest

    p = binomtest(nB, nB + nM, 0.5, alternative="greater").pvalue if nB + nM > 0 else float("nan")
    return dict(n_B_better=nB, n_M_better=nM, n_tied=ties, p_one_sided=float(p), per_multiplet=per)


def diag_T(rows, nboot=300, seed=1):
    Ts = list(np.arange(1998.0, 2020.01, 0.5)) + [2030.0]
    ratios = (1.0, 1.5, 2.0, 2.5, 3.0, 4.0)
    base = model.fit_turn(rows, Ts=Ts, ratios=ratios)
    surface = []
    D = [r["D"] for r in rows]
    for T in Ts:
        vals = [model.spearman([model.delta_B(r["t1"], r["t2"], T, q) for r in rows], D)[0] for q in ratios]
        vals = [v for v in vals if math.isfinite(v)]
        surface.append(dict(T=float(T), best_rho_over_ratio=float(max(vals)) if vals else float("nan")))
    groups = {}
    for r in rows:
        groups.setdefault(r["mult"], []).append(r)
    keys = list(groups)
    rng = np.random.default_rng(seed)
    bestT = []
    for _ in range(nboot):
        rr = [x for k in rng.integers(0, len(keys), len(keys)) for x in groups[keys[k]]]
        f = model.fit_turn(rr, Ts=Ts, ratios=ratios)
        if f["best_T"] is not None:
            bestT.append(f["best_T"])
    bt = np.array(bestT)
    return dict(best=base, surface=surface, boot_n=len(bt),
                boot_T_p05=float(np.percentile(bt, 5)) if len(bt) else None,
                boot_T_median=float(np.median(bt)) if len(bt) else None,
                boot_T_p95=float(np.percentile(bt, 95)) if len(bt) else None,
                frac_T_2006_2011=float(np.mean((bt >= 2006.5) & (bt <= 2010.5))) if len(bt) else None,
                frac_no_turn=float(np.mean(bt >= 2029)) if len(bt) else None)


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


def run(rows, nboot_T=300):
    ok = [r for r in rows if r["ok"]]
    ins = [r for r in ok if r["t2"] < model.NEW_CUT]
    return dict(n_rows=len(rows), n_ok=len(ok), n_in=len(ins),
                E_epoch=diag_E(ins), S_snr=diag_S(ins), J_jackknife=diag_J(ins), Q_without_2021=diag_Q(ins),
                W_within_multiplet=diag_W(ins), T_turn_bootstrap=diag_T(ins, nboot=nboot_T))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pairs", default=PAIRS_URL)
    ap.add_argument("--out-dir", default="/out")
    args = ap.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)
    rows = load_pairs(args.pairs, args.out_dir)
    res = jsonable(run(rows))
    with open(os.path.join(args.out_dir, "ilar_diag.json"), "w") as f:
        json.dump(res, f, indent=1)
    print(json.dumps(res, indent=1)[:8000])


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
# Copyright 2026 MV2, Munim, Inc.
# Licensed under the Apache License, Version 2.0 (Apache-2.0). See LICENSE.
# Made with AI agents working in MV2's name.
"""mkar21.py - pre-registration 21: pre-registration 20's Chile -> Makanchi (MKAR)
backtracking-clock test, run unchanged on the M >= 4.5 event set.

Plan (published before this job ran): findings/2026-11-20-prereg21-mkar-m45-plan.md.

1. Download mkar20.py 0.1.0 at its pinned commit (it fetches its own pinned helpers).
2. Change only the magnitude floor (5.0 -> 4.5) and the plan mkar20 waits for; run mkar20.main().
3. Add the pre-registered extras: test N (pairs with at least one M < 5.0 event), test O (both
   events M >= 5.0), the long-span positive control L (verified before pairs with dM >= 0.5 yr),
   verified counts by kind and the forecasts -> mkar21_summary.json.

Standard library only at top level.

Usage: python mkar21.py --data-dir /data --out-dir /out --wait-min 120 --time-budget-min 420
"""

import argparse
import copy
import hashlib
import json
import math
import os
import sys
import time
import urllib.request

VERSION = "0.1.0"
RAW = "https://raw.githubusercontent.com/munimtechnologies/inner-core-repeaters/"
M20 = ("mkar20.py", "a7f1ceb06c422297020938e08805f2d5d81ce929", 'VERSION = "0.1.0"')
PLAN = "findings/2026-11-20-prereg21-mkar-m45-plan.md"
PLAN_TOKEN = b"pre-registration 21"
MIN_MAG = 4.5
OLD_MIN_MAG = 5.0
LONG_SPAN = 0.5
PC_MIN, PC_RHO = 15, 0.10
EPS = 1e-9
NAN = float("nan")


def isf(x):
    try:
        return math.isfinite(float(x))
    except (TypeError, ValueError):
        return False


def log(*a):
    print(time.strftime("%H:%M:%S"), "[mkar21]", *a, flush=True)


# ----------------------------------------------------------------------------- patching mkar20
class Capture:
    """Wraps mkar20.write_csv: writes as before and keeps a copy of the rows per file name."""

    def __init__(self, orig):
        self.orig = orig
        self.rows = {}

    def __call__(self, path, rows):
        self.orig(path, rows)
        self.rows[os.path.basename(path)] = copy.deepcopy(list(rows))


def patch(m20, min_mag=MIN_MAG):
    """Apply P21's only changes to an imported mkar20 module. Returns the Capture."""
    orig_cp = m20.candidate_pairs

    def candidate_pairs(ev, **kw):
        kw.setdefault("min_mag", min_mag)
        return orig_cp(ev, **kw)

    m20.candidate_pairs = candidate_pairs
    m20.PAIR["min_mag"] = min_mag
    m20.PLAN = PLAN
    m20.PLAN_TOKEN = PLAN_TOKEN
    cap = Capture(m20.write_csv)
    m20.write_csv = cap
    return cap


# ----------------------------------------------------------------------------- pure helpers
def scores_from_rows(rows):
    """mkar20_scores rows (ev1, ev2, S, n_st, ...) -> {(ev1, ev2): dict(S, n_st)}."""
    out = {}
    for r in rows:
        if not isf(r.get("S")):
            continue
        out[(str(r["ev1"]), str(r["ev2"]))] = dict(S=float(r["S"]), n_st=int(float(r["n_st"])),
                                                   D_late_median=r.get("D_late_median"))
    return out


def mags_of(ev_rows):
    return {str(r["id"]): float(r["mag"]) for r in ev_rows if isf(r.get("mag"))}


def is_new(row, mags, old=OLD_MIN_MAG):
    """True if at least one event of the pair is below P20's floor (so P20 could not contain it)."""
    return any(isf(mags.get(e)) and mags[e] < old - EPS for e in (row["ev1"], row["ev2"]))


def by_kind(rows):
    return {k: sum(1 for r in rows if r.get("kind") == k) for k in ("straddle", "before", "after")}


def long_span_pc(V, spearman, family_boot, pct, span=LONG_SPAN):
    """Positive control on verified before pairs with dM >= span."""
    B = [r for r in V if r.get("kind") == "before" and isf(r.get("dM")) and r["dM"] >= span - EPS]

    def stat(sub):
        return spearman([r["dM"] for r in sub], [r["D"] for r in sub])

    rho = stat(B)
    bo = [x for x in family_boot(B, stat) if isf(x)]
    p05, p95 = pct(bo, 0.05), pct(bo, 0.95)
    return dict(span_yr=span, n_pairs=len(B), n_families=len({r["mult"] for r in B}), rho=rho, p05=p05, p95=p95,
                passed=bool(len(B) >= PC_MIN and isf(rho) and rho >= PC_RHO and isf(p05) and p05 > 0))


def forecasts(primary_decision, n_decision, l_passed, n_verified, calib_passed):
    untest = ("untestable", "untestable-calibration")
    return dict(H1_primary_testable=primary_decision not in untest,
                H2_primary_decision=primary_decision,
                H3_new_pairs_testable=n_decision not in untest,
                H4_long_span_control_passes=bool(l_passed),
                H5_verified_ge_80=bool(n_verified >= 80),
                H6_calibration_passes=bool(calib_passed))


def replication(primary_decision, n_decision):
    if primary_decision in ("B", "M"):
        return "replicated" if n_decision == primary_decision else "not replicated in new pairs"
    return "n/a"


# ----------------------------------------------------------------------------- I/O
def fetch_bytes(url, tries=6, wait=20, timeout=90):
    last = None
    for k in range(tries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "inner-core-repeaters (MV2, Munim, Inc.)"})
            with urllib.request.urlopen(req, timeout=timeout) as r:
                return r.read()
        except Exception as ex:  # noqa: BLE001
            last = ex
            time.sleep(wait * (k + 1))
    raise RuntimeError("fetch failed %s: %r" % (url, last))


def ensure_m20(here):
    name, sha, must = M20
    p = os.path.join(here, name)
    if not os.path.exists(p):
        body = fetch_bytes(RAW + sha + "/" + name)
        if must.encode() not in body:
            raise RuntimeError("unexpected content in %s at %s" % (name, sha))
        with open(p, "wb") as f:
            f.write(body)
    with open(p, "rb") as f:
        return dict(file=name, commit=sha, sha256=hashlib.sha256(f.read()).hexdigest())


# ----------------------------------------------------------------------------- main
def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", default="/data")
    ap.add_argument("--out-dir", default="/out")
    ap.add_argument("--wait-min", type=float, default=0.0)
    ap.add_argument("--time-budget-min", type=float, default=420.0)
    a = ap.parse_args(argv)
    here = os.path.dirname(os.path.abspath(__file__))
    if here not in sys.path:
        sys.path.insert(0, here)
    info = ensure_m20(here)
    log("mkar20", info)
    import mkar20 as m20

    cap = patch(m20)
    rc = m20.main(["--data-dir", a.data_dir, "--out-dir", a.out_dir, "--wait-min", str(a.wait_min),
                   "--time-budget-min", str(a.time_budget_min)])
    import mantlep3 as m3

    with open(os.path.join(a.out_dir, "mkar20_summary.json"), encoding="utf-8") as f:
        base = json.load(f)
    frows = cap.rows.get("mkar20_pairs.csv", [])
    ev_rows = cap.rows.get("mkar20_events.csv", [])
    scores = scores_from_rows(cap.rows.get("mkar20_scores.csv", []))
    mags = mags_of(ev_rows)
    new = [r for r in frows if is_new(r, mags)]
    old = [r for r in frows if not is_new(r, mags)]
    g = dict(N_new_pairs=m3.clock_test(new, scores, "D"), O_both_m5=m3.clock_test(old, scores, "D"))
    V = m20.verified(frows, scores)
    g["L_long_span_control"] = long_span_pc(V, m20.spearman, m20.family_boot, m20.pct)
    prim = base.get("grade", {}).get("primary", {})
    pd = base.get("grade", {}).get("decision", "missing")
    g["replication"] = replication(pd, g["N_new_pairs"]["decision"])
    s = dict(version=VERSION, mkar20=info, min_mag=MIN_MAG, plan=PLAN, mkar20_main_rc=rc,
             plan_public_at_start=base.get("plan_public_at_start"), primary_decision=pd,
             n_rows=len(frows), n_rows_new=len(new), n_rows_both_m5=len(old),
             n_events_with_mag=len(mags), n_events_below_5=sum(1 for m in mags.values() if m < OLD_MIN_MAG - EPS),
             verified_by_kind=by_kind(V), n_verified=len(V), extras=g)
    s["forecasts"] = forecasts(pd, g["N_new_pairs"]["decision"], g["L_long_span_control"]["passed"], len(V),
                               prim.get("calibration", {}).get("passed", False))
    with open(os.path.join(a.out_dir, "mkar21_summary.json"), "w", encoding="utf-8") as f:
        json.dump(m20.jsonable(s), f, indent=1)
    log("primary", pd, "N", g["N_new_pairs"]["decision"], "O", g["O_both_m5"]["decision"],
        "L", g["L_long_span_control"]["passed"], json.dumps(m20.jsonable(s["forecasts"])))
    return rc


if __name__ == "__main__":
    sys.exit(main())

#!/usr/bin/env python3
"""ykab_grade.py: Rule A bookkeeping and hypothesis H4 for YKA stage B.

Part of inner-core-repeaters by MV2, Munim, Inc. (made with AI agents). Licence: Apache-2.0.

What it grades, all fixed before this file was written:
  * Rule A of findings/2026-11-18-yka-pick-amendment.md: a multiplet's stack pick is valid if
    |pick_minus_pred_s| <= 2.0 s; primary grading uses valid-pick multiplets only; if > 25% of
    the multiplets with >= 2 YKA members are mis-picked, stage B is re-run with the constrained
    picker (icr-ykab-v2), which becomes primary.
  * H4 of findings/2026-11-17-yka-hardware-prereg.md: the hardware floor (median D_hw of
    icr-ykagen-v1) is below the median YKA D of same-generation ok stage-B pairs less than
    2.0 years apart.

Reported, not graded: the v1 run regraded on valid-pick multiplets only ("v1, valid picks"),
which Rule A asks for alongside the v2 run.
"""
import csv
import json
import math
import os
import time

VALID_PICK_S = 2.0
RULE_A_MAX_FRAC = 0.25
# icr-ykagen-v1/out/ykagen_summary.json, "D_median" (16 usable events, H1-H3 hold)
D_HW_MEDIAN = 0.011288166121190768
NEAR_YEARS = 2.0
V1_RAW = "https://raw.githubusercontent.com/munimtechnologies/inner-core-repeaters/main/results/ykab-v1/"


def fnum(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return float("nan")


def median(xs):
    xs = sorted(x for x in xs if math.isfinite(x))
    n = len(xs)
    if n == 0:
        return float("nan")
    return xs[n // 2] if n % 2 else 0.5 * (xs[n // 2 - 1] + xs[n // 2])


def pick_valid(v, lim=VALID_PICK_S):
    v = fnum(v)
    return math.isfinite(v) and abs(v) <= lim


# ----------------------------------------------------------------------------- tables
def read_csv(path):
    if not path or not os.path.exists(path) or os.path.getsize(path) == 0:
        return []
    with open(path, newline="") as f:
        return list(csv.DictReader(f))


def read_pairs(path):
    rows = []
    for r in read_csv(path):
        rr = dict(r)
        for k in ("t1", "t2", "D", "Dx", "Ds", "D_c13", "D_c310", "D_full", "cc_ref", "snr_dir1", "snr_dir2",
                  "snr_coda1", "snr_coda2"):
            rr[k] = fnum(r.get(k))
        if not math.isfinite(rr["Ds"]):
            rr["Ds"] = None
        rr["ok"] = str(r.get("ok")) == "True"
        rows.append(rr)
    return rows


# ----------------------------------------------------------------------------- Rule A
def rule_a(mult_rows):
    """Pick validity over multiplets with >= 2 YKA members."""
    two = []
    for r in mult_rows:
        try:
            n_yka = int(float(r.get("n_yka") or 0))
        except ValueError:
            n_yka = 0
        if n_yka >= 2:
            two.append(r)
    picked = [r for r in two if math.isfinite(fnum(r.get("pick_minus_pred_s")))]
    invalid = [r["mult"] for r in picked if not pick_valid(r.get("pick_minus_pred_s"))]
    valid = sorted(r["mult"] for r in picked if pick_valid(r.get("pick_minus_pred_s")))
    frac = len(invalid) / len(two) if two else float("nan")
    picks = [fnum(r.get("pick_minus_pred_s")) for r in picked]
    return dict(n_multiplets_ge2_yka=len(two), n_with_pick=len(picked), n_invalid=len(invalid),
                n_valid=len(valid), frac_invalid_of_ge2=frac,
                frac_invalid_of_picked=(len(invalid) / len(picked)) if picked else float("nan"),
                triggers_v2=bool(two) and frac > RULE_A_MAX_FRAC,
                n_pick_before_minus4s=sum(1 for p in picks if p < -4.0),
                n_pick_minus4_to_minus2s=sum(1 for p in picks if -4.0 <= p < -2.0),
                n_pick_within_2s=sum(1 for p in picks if abs(p) <= 2.0),
                n_pick_after_2s=sum(1 for p in picks if p > 2.0),
                invalid=sorted(invalid), valid=valid)


# ----------------------------------------------------------------------------- H4
def h4(rows, d_hw=D_HW_MEDIAN, near=NEAR_YEARS, key="D"):
    sel = [r for r in rows if r["ok"] and r.get("gen1") == r.get("gen2") and r.get("gen1")
           and math.isfinite(r["t1"]) and math.isfinite(r["t2"]) and abs(r["t2"] - r["t1"]) < near
           and math.isfinite(r.get(key, float("nan")))]
    m = median([r[key] for r in sel])
    testable = len(sel) >= 3
    return dict(n_pairs=len(sel), n_multiplets=len({r["mult"] for r in sel}),
                n_legacy=sum(1 for r in sel if r["gen1"] == "legacy"),
                n_current=sum(1 for r in sel if r["gen1"] == "current"),
                median_D_same_gen_near=m, median_D_hw=d_hw, ratio=(m / d_hw) if d_hw > 0 else float("nan"),
                testable=testable, verdict=("holds" if m > d_hw else "fails") if testable else "untestable")


# ----------------------------------------------------------------------------- regrade
def regrade(rows, keep, model=None, yka_b=None, nboot=1000):
    """Grade a pairs table restricted to the multiplets in `keep` (model.evaluate + the F1-Y vote)."""
    sub = [r for r in rows if r["mult"] in keep]
    out = dict(n_rows=len(sub), n_ok=sum(r["ok"] for r in sub), n_multiplets=len({r["mult"] for r in sub}))
    if yka_b is not None:
        cases = yka_b.cases_from_pairs(sub, "D")
        out["F1_Y"] = yka_b.grade_vote(cases, nboot=nboot)
        out["F1x_Y"] = yka_b.grade_vote(yka_b.cases_from_pairs(sub, "Dx"), nboot=nboot)
        out["F1_Y_D_c310"] = yka_b.grade_vote(yka_b.cases_from_pairs(sub, "D_c310"), nboot=300)
        out["cases"] = cases
    if model is not None:
        try:
            out["predictions"] = model.evaluate(sub, nboot=nboot)
        except Exception as ex:  # report, never stop the job
            out["predictions"] = dict(error=repr(ex)[:300])
    return out


def fetch(url, path, must_contain, tries=8, wait_s=60):
    import requests

    for attempt in range(tries):
        try:
            r = requests.get(url, timeout=90)
            if r.status_code == 200 and must_contain in r.text:
                with open(path, "w", encoding="utf-8", newline="") as f:
                    f.write(r.text)
                return True
        except Exception:
            pass
        time.sleep(wait_s)
    return False


def clean(o):
    if isinstance(o, dict):
        return {str(k): clean(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [clean(v) for v in o]
    if isinstance(o, float) and not math.isfinite(o):
        return None
    if hasattr(o, "item") and not isinstance(o, (str, bytes)):
        try:
            return clean(o.item())
        except Exception:
            return str(o)
    return o


def run(out_dir, model=None, yka_b=None, nboot=1000, fetch_v1=True):
    """Grade this run (v2) and, if it can be fetched, the published v1 run."""
    res = dict(rule=dict(valid_pick_s=VALID_PICK_S, rule_a_max_frac=RULE_A_MAX_FRAC, d_hw_median=D_HW_MEDIAN,
                         near_years=NEAR_YEARS))
    m2 = read_csv(os.path.join(out_dir, "yka_multiplets.csv"))
    p2 = read_pairs(os.path.join(out_dir, "yka_pairs.csv"))
    ra2 = rule_a(m2)
    res["v2_rule_a"] = ra2
    res["H4_v2_primary"] = h4(p2)
    res["H4_v2_D_c310"] = h4(p2, key="D_c310")
    v2v = regrade(p2, set(ra2["valid"]), model=model, yka_b=yka_b, nboot=nboot)
    res["v2_valid_only"] = {k: v for k, v in v2v.items() if k != "cases"}
    if fetch_v1:
        d = os.path.join(out_dir, "v1")
        os.makedirs(d, exist_ok=True)
        ok1 = fetch(V1_RAW + "yka_multiplets.csv", os.path.join(d, "yka_multiplets.csv"), "pick_minus_pred_s")
        ok2 = ok1 and fetch(V1_RAW + "yka_pairs.csv", os.path.join(d, "yka_pairs.csv"), "array,mult")
        if ok1 and ok2:
            m1 = read_csv(os.path.join(d, "yka_multiplets.csv"))
            p1 = read_pairs(os.path.join(d, "yka_pairs.csv"))
            ra1 = rule_a(m1)
            res["v1_rule_a"] = ra1
            res["H4_v1_as_written"] = h4(p1)
            res["H4_v1_valid_picks"] = h4([r for r in p1 if r["mult"] in set(ra1["valid"])])
            v1v = regrade(p1, set(ra1["valid"]), model=model, yka_b=yka_b, nboot=nboot)
            res["v1_valid_only"] = {k: v for k, v in v1v.items() if k != "cases"}
        else:
            res["v1"] = "could not fetch results/ykab-v1 from the public repo"
    with open(os.path.join(out_dir, "ykab_grade.json"), "w") as f:
        json.dump(clean(res), f, indent=1)
    return res

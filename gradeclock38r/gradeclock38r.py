#!/usr/bin/env python3
"""gradeclock38r.py: robustness of pre-registration 38 arm W (does it all hang on one multiplet or on YKA's bad-clock years?).

Part of inner-core-repeaters by MV2, Munim, Inc. (mv2@munimtech.com). Made with AI agents working
in MV2's name; no human expert has checked it. Licence: Apache-2.0.

Written AFTER the arm W results were known (findings/2026-11-27-prereg38-W-results.md), so these are
labelled robustness checks, not a new confirmatory test. Their decision rules were fixed and committed
(findings/2026-11-27-prereg38r-robustness-plan.md) before this code ran.

It re-uses gradeclock38.py v0.1.0 unchanged (the committed copy is imported from this folder or, if
absent, downloaded from the repository) and changes only which pairs enter:

  R1  drop multiplet A (the largest family, which carries most of the arm-W gain); W1 and W2 per array.
  R2  YKA only: pairs with any event inside the YKA timing-error eras (Ackerley 2026, SRL:
      2014-12-01..2020-04-08 and 2021-12-01..2022-04-01) are treated as ungraded; W1, W2, W3.
  R3  transfer: the best fine-grid clocks fitted on multiplet A alone, scored on all other multiplets.
      Graded percentile = share of fine-grid clocks scoring strictly better on the others than the
      (upper) median score, on the others, of A's tied best clocks.
      transfers if <= 0.10; does-not-transfer if > 0.25; else mixed.

The Wang et al. (2024) table is read at run time and never written to the output folder.
"""
import argparse
import importlib
import json
import os
import sys
import urllib.request

VERSION = "0.1.0"
BASE_URL = ("https://raw.githubusercontent.com/munimtechnologies/inner-core-repeaters/main/"
            "gradeclock38/gradeclock38.py")
BASE_VERSION = "0.1.0"
ERAS = (("2014-12-01", "2020-04-08"), ("2021-12-01", "2022-04-01"))   # YKA timing-error eras
UA = "inner-core-repeaters/gradeclock38r (MV2, Munim, Inc.; mv2@munimtech.com)"


def load_base(cache_dir="/tmp/gc38"):
    """Import gradeclock38 from this folder, else download the committed copy."""
    here = os.path.dirname(os.path.abspath(__file__))
    for d in (here, cache_dir):
        if os.path.exists(os.path.join(d, "gradeclock38.py")) and d not in sys.path:
            sys.path.insert(0, d)
    try:
        g = importlib.import_module("gradeclock38")
    except ImportError:
        os.makedirs(cache_dir, exist_ok=True)
        req = urllib.request.Request(BASE_URL, headers={"User-Agent": UA})
        with urllib.request.urlopen(req, timeout=120) as r:
            body = r.read()
        with open(os.path.join(cache_dir, "gradeclock38.py"), "wb") as f:
            f.write(body)
        if cache_dir not in sys.path:
            sys.path.insert(0, cache_dir)
        importlib.invalidate_caches()
        g = importlib.import_module("gradeclock38")
    if g.VERSION != BASE_VERSION:
        raise RuntimeError("gradeclock38 version %s, expected %s" % (g.VERSION, BASE_VERSION))
    return g


def in_eras(t, g):
    return any(g.decimal_year(a) <= t <= g.decimal_year(b) for a, b in ERAS)


def drop_multiplet(rows, m):
    return [r for r in rows if r["mult"] != m]


def mask_yka_eras(rows, g):
    """Copy of rows with the YKA grade set to -2 (no data) when either event lies in an era."""
    out = []
    for r in rows:
        r2 = dict(r)
        if in_eras(r["t1"], g) or in_eras(r["t2"], g):
            r2["YKA"] = -2
        out.append(r2)
    return out


def grade_r3(pct):
    if pct <= 0.10:
        return "transfers"
    return "does-not-transfer" if pct > 0.25 else "mixed"


def transfer(rows, g, array, source="A", turns=None, ratios=None):
    """Best clocks on multiplet `source` alone, scored on every other multiplet."""
    times = g.event_times(rows)
    clocks = g.grid_clocks(turns or g.FINE_TURNS, ratios or g.WIDE_RATIOS)
    spec_src = g.build([r for r in rows if r["mult"] == source], array)
    spec_rest = g.build([r for r in rows if r["mult"] != source], array)
    if not spec_src["comps"] or not spec_rest["comps"]:
        return dict(grade="untestable", n_src=len(spec_src["comps"]), n_rest=len(spec_rest["comps"]))
    s_src = g.scan(spec_src, times, clocks)
    s_rest = g.scan(spec_rest, times, clocks)
    best = g.best_set(clocks, s_src)
    idx = {c: i for i, c in enumerate(clocks)}
    vals = sorted(s_rest[idx[c]] for c in best)
    med = vals[len(vals) // 2]                       # upper median over the source's tied best clocks
    pct_med = g.frac_below(med, s_rest)              # graded quantity
    pct_best = g.frac_below(vals[0], s_rest)         # report only (most favourable member of the set)
    wang = g.score(spec_rest, times, g.NAMED["wang_2024_T2008.5_r2.5"])
    return dict(n_src=len(spec_src["comps"]), n_rest=len(spec_rest["comps"]),
                src_best=g.describe_set(best), src_s_min=round(min(s_src), 6),
                rest_s_median_at_src_best=round(med, 6), rest_s_min_at_src_best=round(vals[0], 6),
                rest_s_noturn=round(s_rest[0], 6), rest_s_min=round(min(s_rest), 6),
                rest_pct_median=round(pct_med, 6), rest_pct_best=round(pct_best, 6),
                rest_s_wang=round(wang, 6), rest_wang_pct=round(g.frac_below(wang, s_rest), 6),
                grade=grade_r3(pct_med))


def strip(res):
    """Keep the graded numbers; drop the long lists."""
    out = dict(checks=res["checks"], W4=res.get("W4"), arrays={})
    for a, d in res["arrays"].items():
        out["arrays"][a] = {k: v for k, v in d.items() if k not in ("best_test",)}
        if "best_fine" in d:
            bf = dict(d["best_fine"])
            bf.pop("first", None)
            out["arrays"][a]["best_fine"] = bf
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--table", default="")
    ap.add_argument("--data-dir", default="/data")
    ap.add_argument("--out-dir", default="/out")
    ap.add_argument("--shuffles", type=int, default=500)
    args = ap.parse_args()
    g = load_base()
    os.makedirs(args.out_dir, exist_ok=True)
    if args.table:
        rows = g.parse_table2_csv(args.table)
    else:
        os.makedirs(args.data_dir, exist_ok=True)
        rows = g.fetch_table2(args.data_dir)

    full, _, _ = g.analyse(rows, n_shuffles=0, fine=False)          # sanity: reproduce arm W counts
    checks = dict(n_rows=len(rows), n_comp_YKA=full["arrays"]["YKA"]["n_comp"],
                  n_comp_ILAR=full["arrays"]["ILAR"]["n_comp"],
                  s_star_test_YKA=full["arrays"]["YKA"]["s_star_test"],
                  s_star_test_ILAR=full["arrays"]["ILAR"]["s_star_test"])
    expected = dict(n_rows=143, n_comp_YKA=193, n_comp_ILAR=141,
                    s_star_test_YKA=0.124352, s_star_test_ILAR=0.120567)
    checks_ok = all(abs(checks[k] - v) < 1e-6 for k, v in expected.items())

    r1, _, _ = g.analyse(drop_multiplet(rows, "A"), n_shuffles=args.shuffles)
    masked = mask_yka_eras(rows, g)
    r2, _, _ = g.analyse(masked, n_shuffles=args.shuffles)
    n_masked = sum(1 for a, b in zip(rows, masked) if a["YKA"] in g.GRADES and b["YKA"] not in g.GRADES)
    r3 = {a: transfer(rows, g, a) for a in g.ARRAYS}

    res = dict(version=VERSION, base_version=g.VERSION, n_shuffles=args.shuffles, seed=g.SEED,
               checks=checks, expected=expected, checks_ok=checks_ok,
               R1_drop_A=dict(strip(r1), grades={a: dict(W1=r1["arrays"][a].get("W1", "untestable"),
                                                         W2=r1["arrays"][a].get("W2", "untestable"))
                                                  for a in g.ARRAYS}),
               R2_yka_eras=dict(n_yka_pairs_masked=n_masked, YKA=strip(r2)["arrays"]["YKA"],
                                grades=dict(W1=r2["arrays"]["YKA"].get("W1", "untestable"),
                                            W2=r2["arrays"]["YKA"].get("W2", "untestable"),
                                            W3=r2["arrays"]["YKA"].get("W3", "untestable"))),
               R3_transfer_from_A=r3)
    with open(os.path.join(args.out_dir, "gradeclock38r_summary.json"), "w") as f:
        json.dump(res, f, indent=1)
    print(json.dumps(dict(checks_ok=checks_ok, R1=res["R1_drop_A"]["grades"], R2=res["R2_yka_eras"]["grades"],
                          R2_masked=n_masked, R3={a: r3[a]["grade"] for a in r3})))
    return 0


if __name__ == "__main__":
    sys.exit(main())

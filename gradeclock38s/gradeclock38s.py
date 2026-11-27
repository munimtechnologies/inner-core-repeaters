#!/usr/bin/env python3
"""gradeclock38s.py: does arm W's turn signal survive the mantle-P source gate? (pre-registration 38, arm WS)

Part of inner-core-repeaters by MV2, Munim, Inc. (mv2@munimtech.com). Made with AI agents working
in MV2's name; no human expert has checked it. Licence: Apache-2.0.

Rules: findings/2026-11-27-prereg38s-source-gated-plan.md, written and committed AFTER arm W was graded
but BEFORE any result of the source census wsc39 (pre-registration 39) existed. It is a labelled
follow-up of a graded result, not a new confirmatory test.

It re-uses gradeclock38.py v0.1.0 unchanged (same grids, 500 shuffles, seed 38, W1-W4 rules) and
changes only which of Wang et al. (2024)'s Table 2 pairs enter:

  flag of a pair (from wsc39 v0.2.0's wsc39_summary.json):
    odd-member   either event is an odd member of any multiplet
    S            wsc39 decided the pair itself S (not one source)
    N / inconclusive / untestable   wsc39's decision for the pair
    yardstick    both events are M03-24 reference events (not graded by wsc39; relocated repeats)
    untested     no wsc39 decision (e.g. Table 2 members that Table 3 omits)

  WS1 (primary): drop every pair flagged odd-member or S; grade W1-W4 exactly as arm W.
  WS2 (strict):  keep only pairs flagged N or yardstick; grade W1-W4 (under 30 comparisons at an
                 array: untestable there, as in arm W).
  Per array, for WS1 and WS2: survives if W1 structure and W2 turn-needed; lost if W1 no-structure
  or W2 not-needed; otherwise weakened. Headline from WS1: survives at both arrays -> "survives";
  lost at both -> "lost"; else "mixed". If WS1 drops no pair, WS1 is arm W itself: "unchanged".

Rows are dropped, not masked, so a dropped event's date also leaves the shuffle pool.
The Wang et al. table is read at run time and never written to the output folder.
"""
import argparse
import importlib
import json
import os
import sys
import urllib.request

VERSION = "0.1.0"
RAW = "https://raw.githubusercontent.com/munimtechnologies/inner-core-repeaters/main/"
BASE_URL = RAW + "gradeclock38/gradeclock38.py"
BASE_VERSION = "0.1.0"
WSC_URL = RAW + "results/wsc39-v2/wsc39_summary.json"
WSC_VERSION = "0.2.0"
UA = "inner-core-repeaters/gradeclock38s (MV2, Munim, Inc.; mv2@munimtech.com)"

# Wang et al. Table 1 numbers of the M03-24 reference events (usp000bt5x, usp000fz4d, usp000hute,
# us7000am0e, us6000h2wz); the other three M03-24 events are not in Wang et al.'s table.
REF_WANG = frozenset(("718", "734", "742", "848", "854"))
DROP_WS1 = ("odd-member", "S")
KEEP_WS2 = ("N", "yardstick")
ARM_W_EXPECTED = dict(n_rows=143, n_comp_YKA=193, n_comp_ILAR=141,
                      s_star_test_YKA=0.124352, s_star_test_ILAR=0.120567)
ARM_W_GRADES = {"YKA": ("structure", "turn-needed"), "ILAR": ("structure", "turn-needed")}


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


def load_wsc(src):
    if src.startswith("http"):
        req = urllib.request.Request(src, headers={"User-Agent": UA})
        with urllib.request.urlopen(req, timeout=120) as r:
            return json.loads(r.read().decode("utf-8"))
    with open(src, encoding="utf-8") as f:
        return json.load(f)


def wsc_gate(wsc):
    """True if the census is the amended, finished version 0.2.0 with 734 in multiplet A."""
    if wsc.get("version") != WSC_VERSION or not wsc.get("finished"):
        return False
    pairs_a = (wsc.get("multiplets", {}).get("A", {}) or {}).get("pairs", {}) or {}
    return any("734" in k.split("-") for k in pairs_a)


def wsc_flags(wsc):
    """-> (odd: set of Wang numbers as str, decisions: {frozenset({str, str}): decision})"""
    odd, dec = set(), {}
    for m in (wsc.get("multiplets") or {}).values():
        for w in m.get("odd_members") or []:
            odd.add(str(w))
        for pk, d in (m.get("pairs") or {}).items():
            a, b = pk.split("-")
            dec[frozenset((str(a), str(b)))] = d.get("decision")
    return odd, dec


def flag(row, odd, dec):
    e1, e2 = str(row["e1"]), str(row["e2"])
    if e1 in odd or e2 in odd:
        return "odd-member"
    d = dec.get(frozenset((e1, e2)))
    if d is None:
        return "yardstick" if (e1 in REF_WANG and e2 in REF_WANG) else "untested"
    return d


def select(rows, flags, drop=None, keep=None):
    out = []
    for r, f in zip(rows, flags):
        if drop is not None and f in drop:
            continue
        if keep is not None and f not in keep:
            continue
        out.append(r)
    return out


def read_array(grades):
    """grades: (W1, W2) -> survives / lost / weakened / untestable"""
    w1, w2 = grades
    if w1 == "untestable" or w2 == "untestable":
        return "untestable"
    if w1 == "structure" and w2 == "turn-needed":
        return "survives"
    if w1 == "no-structure" or w2 == "not-needed":
        return "lost"
    return "weakened"


def headline(per_array):
    v = list(per_array.values())
    if all(x == "survives" for x in v):
        return "survives"
    if all(x == "lost" for x in v):
        return "lost"
    return "mixed"


def grades_of(res, g):
    out = {}
    for a in g.ARRAYS:
        d = res["arrays"].get(a, {})
        out[a] = dict(W1=d.get("W1", "untestable"), W2=d.get("W2", "untestable"),
                      W3=d.get("W3", "untestable"), n_pairs=d.get("n_pairs"), n_comp=d.get("n_comp"),
                      s_noturn=d.get("s_noturn"), s_star_test=d.get("s_star_test"),
                      s_star_fine=d.get("s_star_fine"), p_structure=d.get("p_structure"),
                      p_turn=d.get("p_turn"), improvement=d.get("improvement"),
                      best_fine=({k: v for k, v in d["best_fine"].items() if k != "first"}
                                 if "best_fine" in d else None),
                      wang=(d.get("named") or {}).get("wang_2024_T2008.5_r2.5"))
    out["W4"] = res.get("W4", {}).get("grade", "untestable")
    out["W4_detail"] = res.get("W4")
    return out


def flag_table(rows, flags, g):
    """Counts of flags among pairs that are graded 0-2 at each array (report only)."""
    t = {}
    for r, f in zip(rows, flags):
        for a in g.ARRAYS:
            if r[a] in g.GRADES and r["mult"] not in g.EXCLUDED_MULTIPLETS and r["mult"] != "-":
                t.setdefault(a, {}).setdefault(f, {}).setdefault(str(r[a]), 0)
                t[a][f][str(r[a])] += 1
    return t


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--table", default="")
    ap.add_argument("--wsc", default=WSC_URL, help="wsc39_summary.json path or URL")
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
    wsc = load_wsc(args.wsc)

    full, _, _ = g.analyse(rows, n_shuffles=0, fine=False)
    checks = dict(n_rows=len(rows), n_comp_YKA=full["arrays"]["YKA"]["n_comp"],
                  n_comp_ILAR=full["arrays"]["ILAR"]["n_comp"],
                  s_star_test_YKA=full["arrays"]["YKA"]["s_star_test"],
                  s_star_test_ILAR=full["arrays"]["ILAR"]["s_star_test"])
    checks_ok = all(abs(checks[k] - v) < 1e-6 for k, v in ARM_W_EXPECTED.items())
    gate_ok = wsc_gate(wsc)
    odd, dec = wsc_flags(wsc)
    flags = [flag(r, odd, dec) for r in rows]
    res = dict(version=VERSION, base_version=g.VERSION, wsc_version=wsc.get("version"),
               wsc_headline=(wsc.get("census") or {}).get("headline"), checks=checks,
               expected=ARM_W_EXPECTED, checks_ok=checks_ok, wsc_gate_ok=gate_ok,
               odd_members=sorted(odd), flags_report_only=flag_table(rows, flags, g),
               pair_flags={r["pair"]: f for r, f in zip(rows, flags)})
    if checks_ok and gate_ok:
        ws1_rows = select(rows, flags, drop=DROP_WS1)
        if len(ws1_rows) == len(rows):
            res["WS1"] = dict(status="unchanged", n_rows=len(rows))
            res["headline"] = "unchanged"
        else:
            r1, _, _ = g.analyse(ws1_rows, n_shuffles=args.shuffles)
            gr = grades_of(r1, g)
            read = {a: read_array((gr[a]["W1"], gr[a]["W2"])) for a in g.ARRAYS}
            res["WS1"] = dict(n_rows=len(ws1_rows), grades=gr, reading=read)
            res["headline"] = headline(read)
        ws2_rows = select(rows, flags, keep=KEEP_WS2)
        r2, _, _ = g.analyse(ws2_rows, n_shuffles=args.shuffles)
        gr2 = grades_of(r2, g)
        res["WS2"] = dict(n_rows=len(ws2_rows), grades=gr2,
                          reading={a: read_array((gr2[a]["W1"], gr2[a]["W2"])) for a in g.ARRAYS})
    else:
        res["headline"] = "not graded (gate failed)"
    with open(os.path.join(args.out_dir, "gradeclock38s_summary.json"), "w") as f:
        json.dump(res, f, indent=1)
    print(json.dumps(dict(checks_ok=checks_ok, wsc_gate_ok=gate_ok, headline=res["headline"],
                          WS1=(res.get("WS1") or {}).get("reading"),
                          WS2=(res.get("WS2") or {}).get("reading"))))
    return 0


if __name__ == "__main__":
    sys.exit(main())

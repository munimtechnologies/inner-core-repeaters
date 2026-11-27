#!/usr/bin/env python3
"""gradez.py: do Zhang & Wen (2026)'s own change / no-change calls order on one inner-core rotation clock?

Part of inner-core-repeaters by MV2, Munim, Inc. (mv2@munimtech.com). Made with AI agents working
in MV2's name; no human expert has checked it. Licence: Apache-2.0.
Pre-registration 38, arm Z (findings/2026-11-27-prereg38-armZ-plan.md), written and committed
before the grades were transcribed or read.

Arm W (gradeclock38.py) found that Wang, Vidale et al. (2024)'s by-eye grades order on one
one-turn rotation clock far better than date-shuffled grades. Zhang & Wen (2026, JGR Solid Earth,
https://doi.org/10.1029/2026JB034595, "ZW26") read the same kind of data as episodic change,
quiescence and recovery of the inner-core surface and say no backtracking model fits. Arm Z runs
the identical frozen statistic on ZW26's own pairwise calls, and compares the best clocks and the
shared pairs with arm W and with Wang's table.

Statistic (as arm W): within each multiplet, every two pairs with different grades form a
comparison; it is discordant if the less similar pair is closer on the clock (ties 1/2).
    s(clock) = discordant comparisons / all comparisons      (0 = perfect, ~0.5 = no relation)

Input: a transcription CSV with columns multiplet,array,date1,date2,region,label,source
(label S = similar / no change, C = changed; anything else is ignored). Regions of one pair at
one array merge to 0 (all S), 2 (all C) or 1 (mixed). Standard library only, except that Wang's
workbook (for the report on shared pairs, Z5) is read with openpyxl if it is available.
"""
import argparse
import calendar
import csv
import datetime
import json
import os
import random
import sys

VERSION = "0.1.0"

# ---- pre-registered constants (prereg 38, arm Z); grid and null identical to arm W except the seed
ARRAYS = ("YKA", "ILAR")
ARRAY_ALIASES = {"YKA": "YKA", "YKAR": "YKA", "ILAR": "ILAR"}
WIDE_RATIOS = [round(0.1 * 1.05 ** k, 4) for k in range(109)]          # 0.1 .. ~19.4
FINE_TURNS = [round(2000.0 + 0.05 * k, 2) for k in range(401)]         # 2000.00 .. 2020.00
TEST_TURNS = [round(2000.0 + 0.25 * k, 2) for k in range(81)]          # 2000.00 .. 2020.00
TEST_RATIOS = WIDE_RATIOS[::4]                                          # 28 values
N_SHUFFLES = 500
SEED = 40
TIE = 1e-9
MIN_COMPARISONS = 30
MIN_SHARED = 5
WANG_CLOCK = (2008.5, 2.5)
ODD_EVENT_F = ("F", "2012-12-03")      # prereg 35 G2: not a mantle-P repeat (report-only variant)
ARMW_URL = ("https://raw.githubusercontent.com/munimtechnologies/inner-core-repeaters/main/"
            "results/gradeclock38-v1/gradeclock38_summary.json")
XLSX_URLS = [
    "https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-024-07536-4/MediaObjects/"
    "41586_2024_7536_MOESM1_ESM.xlsx",
    "https://pmc.ncbi.nlm.nih.gov/articles/instance/11236701/bin/41586_2024_7536_MOESM1_ESM.xlsx",
]
UA = "inner-core-repeaters/gradez (MV2, Munim, Inc.; mv2@munimtech.com)"


# ----------------------------------------------------------------------------- input
def iso_date(s):
    s = str(s).strip().replace("/", "-")
    parts = s[:10].split("-")
    return datetime.date(int(parts[0]), int(parts[1]), int(parts[2])).isoformat()


def decimal_year(date_str):
    d = datetime.date.fromisoformat(iso_date(date_str))
    n = 366 if calendar.isleap(d.year) else 365
    return d.year + (d.timetuple().tm_yday - 1) / float(n)


def merge_grade(labels):
    """Region labels of one pair at one array -> 0 (all S), 2 (all C), 1 (mixed), None (no usable label)."""
    labs = [x for x in labels if x in ("S", "C")]
    if not labs:
        return None
    if all(x == "S" for x in labs):
        return 0
    if all(x == "C" for x in labs):
        return 2
    return 1


def read_transcription(records):
    """records: list of dicts with multiplet,array,date1,date2,region,label -> merged pair grades.
    Returns {array: [(mult, d1, d2, grade)]} with d1 < d2 (ISO dates)."""
    acc = {}
    for r in records:
        arr = ARRAY_ALIASES.get(str(r.get("array", "")).strip().upper())
        if arr is None:
            continue
        lab = str(r.get("label", "")).strip().upper()
        m = str(r.get("multiplet", "")).strip()
        d1, d2 = sorted((iso_date(r["date1"]), iso_date(r["date2"])))
        if d1 == d2:
            continue
        acc.setdefault((arr, m, d1, d2), []).append(lab)
    out = {a: [] for a in ARRAYS}
    for (arr, m, d1, d2) in sorted(acc):
        g = merge_grade(acc[(arr, m, d1, d2)])
        if g is not None:
            out[arr].append((m, d1, d2, g))
    return out


def load_csv(path):
    with open(path, newline="", encoding="utf-8-sig") as f:
        return list(csv.DictReader(f))


def build(pairs, drop=None):
    """pairs: [(mult, d1, d2, grade)] -> spec with event keys 'mult|date'; comps within multiplet."""
    spec_pairs, groups = [], {}
    for m, d1, d2, g in pairs:
        if drop and ((m, d1) == drop or (m, d2) == drop):
            continue
        groups.setdefault(m, []).append(len(spec_pairs))
        spec_pairs.append((m + "|" + d1, m + "|" + d2, g, m))
    comps = []
    for m in sorted(groups):
        idx = groups[m]
        for a in idx:
            for b in idx:
                if spec_pairs[a][2] < spec_pairs[b][2]:
                    comps.append((a, b))
    times, mults = {}, {}
    for e1, e2, _, m in spec_pairs:
        for e in (e1, e2):
            times[e] = decimal_year(e.split("|", 1)[1])
            mults.setdefault(m, set()).add(e)
    return dict(pairs=spec_pairs, comps=comps, times=times, mults=mults)


# ----------------------------------------------------------------------------- clocks (as gradeclock38)
def equiv(t, turn, ratio):
    return t if t <= turn else turn - (t - turn) / ratio


def position(t, clock):
    return t if clock is None else equiv(t, clock[0], clock[1])


def grid_clocks(turns, ratios):
    return [None] + [(T, r) for T in turns for r in ratios]


def discordance(spec, times, clock):
    d = [abs(position(times[a], clock) - position(times[b], clock)) for a, b, _, _ in spec["pairs"]]
    tot = 0.0
    for i, j in spec["comps"]:
        if d[i] > d[j] + TIE:
            tot += 1.0
        elif abs(d[i] - d[j]) <= TIE:
            tot += 0.5
    return tot


def score(spec, times, clock):
    n = len(spec["comps"])
    return discordance(spec, times, clock) / n if n else None


def scan(spec, times, clocks):
    n = len(spec["comps"])
    return [discordance(spec, times, c) / n for c in clocks]


def shuffle_times(times, mults, rng):
    new = dict(times)
    for m in sorted(mults):
        evs = sorted(mults[m])
        ts = [times[e] for e in evs]
        rng.shuffle(ts)
        for e, t in zip(evs, ts):
            new[e] = t
    return new


# ----------------------------------------------------------------------------- grading rules
def p_le(real, nulls):
    return (1 + sum(1 for x in nulls if x <= real + TIE)) / (1.0 + len(nulls))


def p_ge(real, nulls):
    return (1 + sum(1 for x in nulls if x >= real - TIE)) / (1.0 + len(nulls))


def frac_below(value, scores):
    return sum(1 for x in scores if x < value - TIE) / float(len(scores))


def grade_z1(p):
    if p <= 0.05:
        return "structure"
    return "no-structure" if p > 0.20 else "weak"


def grade_z2(p, improvement):
    if improvement < 0.02 or p > 0.20:
        return "not-needed"
    return "turn-needed" if p <= 0.05 else "weak"


def grade_z3(frac):
    return "top5" if frac <= 0.05 else "not-top5"


def grade_z4(pct):
    if pct <= 0.10:
        return "agrees"
    return "disagrees" if pct > 0.25 else "mixed"


def grade_z5(n, n_agree):
    if n < MIN_SHARED:
        return "untestable"
    share = n_agree / float(n)
    if share >= 0.8:
        return "agree"
    return "disagree" if share <= 0.5 else "mixed"


def headline(z1, z2, z4):
    """Pre-registered reading per array."""
    if z1 == "untestable":
        return "untestable"
    if z1 == "no-structure":
        return "no-clock"
    if z1 == "structure" and z2 == "turn-needed":
        return {"agrees": "same-clock", "disagrees": "other-clock"}.get(z4, "turn-unclear")
    if z1 == "structure" and z2 == "not-needed":
        return "calendar-only"
    return "weak"


def best_set(clocks, scores):
    m = min(scores)
    return [c for c, s in zip(clocks, scores) if s <= m + TIE]


def describe_set(clocks):
    turns = [c[0] for c in clocks if c is not None]
    ratios = [c[1] for c in clocks if c is not None]
    return dict(n=len(clocks), includes_noturn=None in clocks,
                turn_range=[min(turns), max(turns)] if turns else None,
                ratio_range=[min(ratios), max(ratios)] if ratios else None,
                first=[list(c) if c else "noturn" for c in clocks[:40]])


def turn_profile(clocks, scores):
    prof = {}
    for c, s in zip(clocks, scores):
        if c is None:
            continue
        rec = prof.setdefault(c[0], [9.0, []])
        if s < rec[0] - TIE:
            prof[c[0]] = [s, [c[1]]]
        elif abs(s - rec[0]) <= TIE:
            rec[1].append(c[1])
    return [dict(turn=T, s_min=round(v[0], 6), r_lo=min(v[1]), r_hi=max(v[1])) for T, v in sorted(prof.items())]


def best_percentile(clocks, spec, times, fine_scores):
    """Lowest percentile, among fine-grid scores, reached by any clock in `clocks`."""
    vals = [frac_below(score(spec, times, None if c in (None, "noturn") else tuple(c)), fine_scores)
            for c in clocks]
    return min(vals) if vals else None


# ----------------------------------------------------------------------------- Wang table (Z5)
def wang_pairs(records):
    """Supplementary Data Table 2 rows -> {(d1, d2): {"YKA": g, "ILAR": g}} with ISO dates, d1 < d2."""
    out, started = {}, False
    for rec in records:
        if not rec:
            continue
        head = str(rec[0]).strip()
        if head.startswith("#"):
            started = True
            continue
        if not started or len(rec) < 11 or not (head.startswith("P") and head[1:].isdigit()):
            continue
        try:
            d1, d2 = sorted((iso_date(rec[6]), iso_date(rec[7])))
            out[(d1, d2)] = {"YKA": int(float(str(rec[2]).strip())), "ILAR": int(float(str(rec[3]).strip()))}
        except (ValueError, IndexError):
            continue
    return out


def fetch_wang(data_dir):
    import urllib.request

    path = os.path.join(data_dir, "wang2024_supp.xlsx")
    if not os.path.exists(path):
        for u in XLSX_URLS:
            try:
                req = urllib.request.Request(u, headers={"User-Agent": UA})
                with urllib.request.urlopen(req, timeout=180) as r:
                    body = r.read()
                if body[:2] == b"PK":
                    with open(path, "wb") as f:
                        f.write(body)
                    break
            except Exception as ex:
                print("download failed", u, repr(ex)[:200], flush=True)
    if not os.path.exists(path):
        return None
    import openpyxl

    wb = openpyxl.load_workbook(path, data_only=True, read_only=True)
    ws = wb.worksheets[1]
    rows = []
    for r in ws.iter_rows(values_only=True):
        rows.append(["" if v is None else (v.isoformat() if hasattr(v, "isoformat") else str(v)) for v in r])
    return wang_pairs(rows)


def shared_agreement(zpairs, wang, array):
    """Pairs graded 0 or 2 by both ZW26 and Wang at `array` (same two dates): count and agreements."""
    n = agree = 0
    rows = []
    for m, d1, d2, g in zpairs:
        w = wang.get((d1, d2), {}).get(array)
        if w is None:
            continue
        rows.append(dict(array=array, multiplet=m, date1=d1, date2=d2, zw26=g, wang=w))
        if g in (0, 2) and w in (0, 2):
            n += 1
            agree += 1 if g == w else 0
    return n, agree, rows


def fetch_armw():
    import urllib.request

    req = urllib.request.Request(ARMW_URL, headers={"User-Agent": UA})
    with urllib.request.urlopen(req, timeout=120) as r:
        return json.loads(r.read().decode("utf-8"))


# ----------------------------------------------------------------------------- analysis
def analyse_array(spec, n_shuffles, seed, fine=True):
    if len(spec["comps"]) < MIN_COMPARISONS:
        return None
    times, mults = spec["times"], spec["mults"]
    test_clocks = grid_clocks(TEST_TURNS, TEST_RATIOS)
    fine_clocks = grid_clocks(FINE_TURNS, WIDE_RATIOS) if fine else test_clocks
    ts = scan(spec, times, test_clocks)
    fs = scan(spec, times, fine_clocks)
    s_star, s_none = min(ts), ts[0]
    rng = random.Random(seed)
    stars, imps = [], []
    for _ in range(n_shuffles):
        sh = shuffle_times(times, mults, rng)
        x = scan(spec, sh, test_clocks)
        stars.append(min(x))
        imps.append(x[0] - min(x))
    p1, p2 = p_le(s_star, stars), p_ge(s_none - s_star, imps)
    s_wang = score(spec, times, WANG_CLOCK)
    w3 = frac_below(s_wang, fs)
    return dict(n_pairs=len(spec["pairs"]), n_comp=len(spec["comps"]),
                s_noturn=round(s_none, 6), s_star_test=round(s_star, 6), s_star_fine=round(min(fs), 6),
                improvement=round(s_none - s_star, 6), p_structure=round(p1, 6), p_turn=round(p2, 6),
                null_s_star_pcts={q: round(sorted(stars)[int(q / 100.0 * (len(stars) - 1))], 6)
                                  for q in (5, 50, 95)} if stars else None,
                s_wang=round(s_wang, 6), wang_frac_fine_below=round(w3, 6),
                Z1=grade_z1(p1), Z2=grade_z2(p2, s_none - s_star), Z3=grade_z3(w3),
                best_fine=describe_set(best_set(fine_clocks, fs)),
                _fine=fs, _fine_clocks=fine_clocks)


def analyse(zpairs, armw_best, wang=None, expected=None, n_shuffles=N_SHUFFLES, seed=SEED, fine=True):
    res = dict(version=VERSION, seed=seed, n_shuffles=n_shuffles, arrays={}, report_only={})
    checks = {}
    specs = {a: build(zpairs[a]) for a in ARRAYS}
    for a in ARRAYS:
        checks["n_pairs_" + a] = len(specs[a]["pairs"])
        checks["n_comp_" + a] = len(specs[a]["comps"])
    res["checks"] = checks
    res["expected"] = expected or {}
    res["checks_ok"] = bool(expected) and all(checks.get(k) == v for k, v in expected.items())
    profiles, full = {}, {}
    for a in ARRAYS:
        R = analyse_array(specs[a], n_shuffles, seed, fine)
        if R is None:
            res["arrays"][a] = dict(n_pairs=len(specs[a]["pairs"]), n_comp=len(specs[a]["comps"]),
                                    Z1="untestable", headline="untestable")
            continue
        fs, fc = R.pop("_fine"), R.pop("_fine_clocks")
        full[a] = (fs, fc)
        pct = best_percentile(armw_best.get(a, []), specs[a], specs[a]["times"], fs) if armw_best.get(a) else None
        R["armW_best_pct"] = None if pct is None else round(pct, 6)
        R["Z4"] = "untestable" if pct is None else grade_z4(pct)
        R["headline"] = headline(R["Z1"], R["Z2"], R["Z4"])
        res["arrays"][a] = R
        profiles[a] = turn_profile(fc, fs)
        # report-only: without the odd 2012-12-03 member of multiplet F
        sd = build(zpairs[a], drop=ODD_EVENT_F)
        res["report_only"]["no_F2012_" + a] = dict(
            n_comp=len(sd["comps"]),
            s_noturn=round(score(sd, sd["times"], None), 6) if sd["comps"] else None,
            s_wang=round(score(sd, sd["times"], WANG_CLOCK), 6) if sd["comps"] else None)
    if len(full) == 2:
        bests = {a: best_set(full[a][1], full[a][0]) for a in ARRAYS}
        py = best_percentile(bests["YKA"], specs["ILAR"], specs["ILAR"]["times"], full["ILAR"][0])
        pi = best_percentile(bests["ILAR"], specs["YKA"], specs["YKA"]["times"], full["YKA"][0])
        res["report_only"]["cross_array"] = dict(YKA_best_in_ILAR_pct=round(py, 6), ILAR_best_in_YKA_pct=round(pi, 6))
    z5 = {}
    rows = []
    for a in ARRAYS:
        if wang is None:
            z5[a] = dict(grade="untestable", reason="Wang table unavailable")
            continue
        n, ag, rr = shared_agreement(zpairs[a], wang, a)
        rows += rr
        z5[a] = dict(n_shared_binary=n, n_agree=ag, n_shared_any=len(rr), grade=grade_z5(n, ag))
    res["Z5"] = z5
    return res, profiles, rows


def write_csv(path, rows):
    keys = []
    for r in rows:
        for k in r:
            if k not in keys:
                keys.append(k)
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys or ["empty"])
        w.writeheader()
        w.writerows(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--grades", required=True, help="ZW26 transcription CSV")
    ap.add_argument("--expect", default="", help="JSON (or a JSON file) of hand-tallied counts, e.g. {\"n_comp_YKA\": 40}")
    ap.add_argument("--data-dir", default="/data")
    ap.add_argument("--out-dir", default="/out")
    ap.add_argument("--shuffles", type=int, default=N_SHUFFLES)
    args = ap.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)
    os.makedirs(args.data_dir, exist_ok=True)
    zpairs = read_transcription(load_csv(args.grades))
    if args.expect and os.path.exists(args.expect):
        with open(args.expect) as f:
            expected = json.load(f)
    else:
        expected = json.loads(args.expect) if args.expect else {}
    try:
        aw = fetch_armw()
        armw_best = {a: aw["arrays"][a]["best_fine"]["first"] for a in ARRAYS
                     if a in aw.get("arrays", {}) and "best_fine" in aw["arrays"][a]}
    except Exception as ex:
        print("arm W summary unavailable", repr(ex)[:200], flush=True)
        armw_best = {}
    try:
        wang = fetch_wang(args.data_dir)
    except Exception as ex:
        print("Wang table unavailable", repr(ex)[:200], flush=True)
        wang = None
    res, profiles, rows = analyse(zpairs, armw_best, wang=wang, expected=expected, n_shuffles=args.shuffles)
    with open(os.path.join(args.out_dir, "gradez_summary.json"), "w") as f:
        json.dump(res, f, indent=1)
    write_csv(os.path.join(args.out_dir, "gradez_turn_profile.csv"),
              [dict(r, array=a) for a in profiles for r in profiles[a]])
    write_csv(os.path.join(args.out_dir, "gradez_shared_pairs.csv"), rows)
    print(json.dumps(dict(checks_ok=res["checks_ok"], checks=res["checks"], Z5=res["Z5"],
                          arrays={a: {k: v for k, v in d.items() if k != "best_fine"}
                                  for a, d in res["arrays"].items()}))[:6000])
    return 0


if __name__ == "__main__":
    sys.exit(main())

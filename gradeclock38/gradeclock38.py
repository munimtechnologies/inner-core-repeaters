#!/usr/bin/env python3
"""gradeclock38.py: do Wang et al. (2024)'s own hand grades single out one inner-core rotation clock?

Part of inner-core-repeaters by MV2, Munim, Inc. (mv2@munimtech.com). Made with AI agents working
in MV2's name; no human expert has checked it. Licence: Apache-2.0.
Pre-registration 38, arm W (findings/2026-11-27-prereg38-gradeclock-plan.md), written and committed
before this code ran on the real table.

Wang, Vidale et al. 2024 (Nature 631:340, https://doi.org/10.1038/s41586-024-07536-4) graded 143
South Sandwich Islands repeating-earthquake pairs by eye at the YKA and ILAR arrays
(Supplementary Data Table 2: 0 similar, 1 somewhat similar, 2 different, -1 noisy, -2 no data) and
read the pattern as the inner core moving forward until about 2008.5 and then back about 2.5 times
more slowly ("backtracking").

If the PKIKP coda is set only by the inner core's angle, then on the right clock the grade of a pair
can only rise with the distance between the two events' positions. For each clock we count, within
each multiplet, the pairs of pairs that break this (a less similar pair closer on the clock than a
more similar one; ties count one half) and divide by the number of comparisons:

    s(clock) = discordant comparisons / all comparisons       (0 = perfect, ~0.5 = no relation)

Clocks: one-turn clocks equiv(t, T, r) (forward until turn T, then back r times slower; as in
clocks34.py / m34.py) plus the no-turn clock (position = date). The null shuffles event dates within
each multiplet (fixed seed) and repeats the whole search, so the flexibility of the clock family is
paid for.

Standard library only, except that the publisher's workbook is read with openpyxl when no CSV is
given. The table itself is never written to the output folder.
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

XLSX_URLS = [
    "https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-024-07536-4/MediaObjects/"
    "41586_2024_7536_MOESM1_ESM.xlsx",
    "https://pmc.ncbi.nlm.nih.gov/articles/instance/11236701/bin/41586_2024_7536_MOESM1_ESM.xlsx",
]
UA = "inner-core-repeaters/gradeclock38 (MV2, Munim, Inc.; mv2@munimtech.com)"

# ---- pre-registered constants (prereg 38, arm W)
ARRAYS = ("YKA", "ILAR")
GRADES = (0, 1, 2)
EXCLUDED_MULTIPLETS = ("I", "M")          # Wang Table 3: "not good repeaters", "too noisy"
WIDE_RATIOS = [round(0.1 * 1.05 ** k, 4) for k in range(109)]          # 0.1 .. ~19.4 (clocks34)
FINE_TURNS = [round(2000.0 + 0.05 * k, 2) for k in range(401)]         # 2000.00 .. 2020.00
TEST_TURNS = [round(2000.0 + 0.25 * k, 2) for k in range(81)]          # 2000.00 .. 2020.00
TEST_RATIOS = WIDE_RATIOS[::4]                                          # 28 values
N_SHUFFLES = 500
SEED = 38
TIE = 1e-9
MIN_COMPARISONS = 30
NAMED = {
    "wang_2024_T2008.5_r2.5": (2008.5, 2.5),
    "yang_song_2023_T2009.0_r1": (2009.0, 1.0),
    "yang_song_2023_T2009.5_r1": (2009.5, 1.0),
}
# hand tallies made before the run (findings/2026-11-27-prereg38-gradeclock-plan.md)
EXPECTED = {"n_rows": 143, "n_comp_YKA": 193, "n_comp_ILAR": 141,
            "A_YKA_wang_discordant": 4.0, "A_YKA_comparisons": 72}


# ----------------------------------------------------------------------------- table
def decimal_year(date_str):
    d = datetime.date.fromisoformat(str(date_str).strip()[:10])
    n = 366 if calendar.isleap(d.year) else 365
    return d.year + (d.timetuple().tm_yday - 1) / float(n)


def _event_id(v):
    s = str(v).strip()
    try:
        return str(int(float(s)))
    except ValueError:
        return s


def parse_rows(records):
    """records: iterable of lists of cell strings (one per spreadsheet row)."""
    rows, started = [], False
    for rec in records:
        if not rec:
            continue
        head = str(rec[0]).strip()
        if head.startswith("#"):
            started = True
            continue
        if not started or len(rec) < 11:
            continue
        if not (head.startswith("P") and head[1:].isdigit()):
            continue
        rows.append(dict(pair=head, YKA=int(float(str(rec[2]).strip())), ILAR=int(float(str(rec[3]).strip())),
                         e1=_event_id(rec[4]), e2=_event_id(rec[5]),
                         t1=decimal_year(rec[6]), t2=decimal_year(rec[7]),
                         mult=str(rec[10]).strip()))
    return rows


def parse_table2_csv(path):
    with open(path, newline="", encoding="utf-8-sig", errors="replace") as f:
        return parse_rows(list(csv.reader(f)))


def _cell(v):
    if v is None:
        return ""
    if isinstance(v, (datetime.datetime, datetime.date)):
        return v.isoformat()
    return str(v)


def fetch_table2(data_dir):
    """Download the publisher's workbook to data_dir (never to /out) and return sheet 2's rows."""
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
                print("not an xlsx from", u, flush=True)
            except Exception as ex:  # network hiccup: try the next mirror
                print("download failed", u, repr(ex)[:200], flush=True)
    if not os.path.exists(path):
        raise RuntimeError("workbook unavailable")
    import openpyxl

    wb = openpyxl.load_workbook(path, data_only=True, read_only=True)
    ws = wb.worksheets[1]
    return parse_rows([[_cell(v) for v in r] for r in ws.iter_rows(values_only=True)])


def event_times(rows):
    times = {}
    for r in rows:
        for e, t in ((r["e1"], r["t1"]), (r["e2"], r["t2"])):
            if e in times and abs(times[e] - t) > 1e-6:
                raise ValueError("event %s has two dates" % e)
            times[e] = t
    return times


def multiplet_events(rows):
    out = {}
    for r in rows:
        m = r["mult"]
        if m == "-" or m in EXCLUDED_MULTIPLETS:
            continue
        out.setdefault(m, set()).update((r["e1"], r["e2"]))
    return out


def build(rows, array, pooled=False):
    """Pairs graded 0-2 at `array` and the comparisons (i, j) with grade_i < grade_j.
    Within-multiplet by default; pooled=True compares all graded pairs (singletons included)."""
    pairs, groups = [], {}
    for r in rows:
        g = r[array]
        if g not in GRADES:
            continue
        m = r["mult"]
        if m in EXCLUDED_MULTIPLETS:
            continue
        if not pooled and m == "-":
            continue
        key = "all" if pooled else m
        groups.setdefault(key, []).append(len(pairs))
        pairs.append((r["e1"], r["e2"], g, key))
    comps = []
    for key in sorted(groups):
        idx = groups[key]
        for a in idx:
            for b in idx:
                if pairs[a][2] < pairs[b][2]:
                    comps.append((a, b))
    return dict(pairs=pairs, comps=comps)


# ----------------------------------------------------------------------------- clocks
def equiv(t, turn, ratio):
    """Effective position on a one-turn clock (clocks34.py / m34.py)."""
    return t if t <= turn else turn - (t - turn) / ratio


def position(t, clock):
    return t if clock is None else equiv(t, clock[0], clock[1])


def grid_clocks(turns, ratios):
    return [None] + [(T, r) for T in turns for r in ratios]


def distances(spec, times, clock):
    return [abs(position(times[a], clock) - position(times[b], clock)) for a, b, _, _ in spec["pairs"]]


def discordance(spec, times, clock):
    d = distances(spec, times, clock)
    tot = 0.0
    for i, j in spec["comps"]:
        di, dj = d[i], d[j]
        if di > dj + TIE:
            tot += 1.0
        elif abs(di - dj) <= TIE:
            tot += 0.5
    return tot


def discordance_groups(spec, times, clock):
    """{group: [discordant, comparisons]}"""
    d = distances(spec, times, clock)
    out = {}
    for i, j in spec["comps"]:
        key = spec["pairs"][i][3]
        rec = out.setdefault(key, [0.0, 0])
        rec[1] += 1
        if d[i] > d[j] + TIE:
            rec[0] += 1.0
        elif abs(d[i] - d[j]) <= TIE:
            rec[0] += 0.5
    return out


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


# ----------------------------------------------------------------------------- statistics
def p_le(real, nulls):
    return (1 + sum(1 for x in nulls if x <= real + TIE)) / (1.0 + len(nulls))


def p_ge(real, nulls):
    return (1 + sum(1 for x in nulls if x >= real - TIE)) / (1.0 + len(nulls))


def frac_below(value, scores):
    return sum(1 for x in scores if x < value - TIE) / float(len(scores))


def grade_w1(p):
    if p <= 0.05:
        return "structure"
    return "no-structure" if p > 0.20 else "weak"


def grade_w2(p, improvement):
    if improvement < 0.02 or p > 0.20:
        return "not-needed"
    return "turn-needed" if p <= 0.05 else "weak"


def grade_w3(frac):
    return "top5" if frac <= 0.05 else "not-top5"


def grade_w4(pct_a_in_b, pct_b_in_a):
    if pct_a_in_b <= 0.10 and pct_b_in_a <= 0.10:
        return "consistent"
    if pct_a_in_b > 0.25 and pct_b_in_a > 0.25:
        return "inconsistent"
    return "mixed"


def best_set(clocks, scores):
    m = min(scores)
    return [c for c, s in zip(clocks, scores) if s <= m + TIE]


def cross_percentile(best_clocks_a, clocks, scores_b):
    """Best (lowest) percentile, in array b's grid scores, of any clock in array a's best set."""
    by = dict(zip([c if c is not None else "noturn" for c in clocks], scores_b))
    return min(frac_below(by[c if c is not None else "noturn"], scores_b) for c in best_clocks_a)


def describe_set(clocks):
    turns = [c[0] for c in clocks if c is not None]
    ratios = [c[1] for c in clocks if c is not None]
    return dict(n=len(clocks), includes_noturn=None in clocks,
                turn_range=[min(turns), max(turns)] if turns else None,
                ratio_range=[min(ratios), max(ratios)] if ratios else None,
                first=[list(c) if c else "noturn" for c in clocks[:40]])


def turn_profile(clocks, scores):
    """Per turn date: best score over ratios and the ratios that reach it."""
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


# ----------------------------------------------------------------------------- whole analysis
def analyse(rows, n_shuffles=N_SHUFFLES, seed=SEED, fine=True):
    times = event_times(rows)
    mults = multiplet_events(rows)
    specs = {a: build(rows, a) for a in ARRAYS}
    pooled = {a: build(rows, a, pooled=True) for a in ARRAYS}
    test_clocks = grid_clocks(TEST_TURNS, TEST_RATIOS)
    fine_clocks = grid_clocks(FINE_TURNS, WIDE_RATIOS) if fine else test_clocks

    checks = dict(n_rows=len(rows), n_comp_YKA=len(specs["YKA"]["comps"]),
                  n_comp_ILAR=len(specs["ILAR"]["comps"]))
    ga = discordance_groups(specs["YKA"], times, NAMED["wang_2024_T2008.5_r2.5"]).get("A", [None, 0])
    checks["A_YKA_wang_discordant"], checks["A_YKA_comparisons"] = ga[0], ga[1]
    checks_ok = all(checks[k] == v for k, v in EXPECTED.items())

    real = {}
    for a in ARRAYS:
        if len(specs[a]["comps"]) < MIN_COMPARISONS:
            real[a] = None
            continue
        ts = scan(specs[a], times, test_clocks)
        fs = scan(specs[a], times, fine_clocks)
        real[a] = dict(test=ts, fine=fs, s_noturn=ts[0], s_star=min(ts), improvement=ts[0] - min(ts))

    rng = random.Random(seed)
    null = {a: [] for a in ARRAYS}
    for k in range(n_shuffles):
        sh = shuffle_times(times, mults, rng)
        for a in ARRAYS:
            if real[a] is None:
                continue
            ts = scan(specs[a], sh, test_clocks)
            null[a].append((ts[0], min(ts)))

    res = dict(version=VERSION, n_shuffles=n_shuffles, seed=seed, checks=checks, expected=EXPECTED,
               checks_ok=checks_ok, n_test_clocks=len(test_clocks), n_fine_clocks=len(fine_clocks),
               arrays={})
    best = {}
    for a in ARRAYS:
        if real[a] is None:
            res["arrays"][a] = dict(grade="untestable", n_comp=len(specs[a]["comps"]))
            continue
        R = real[a]
        stars = [x[1] for x in null[a]]
        imps = [x[0] - x[1] for x in null[a]]
        noturns = [x[0] for x in null[a]]
        p1 = p_le(R["s_star"], stars)
        p2 = p_ge(R["improvement"], imps)
        p0 = p_le(R["s_noturn"], noturns)
        named = {}
        for name, c in NAMED.items():
            s = score(specs[a], times, c)
            named[name] = dict(s=round(s, 6), frac_fine_below=round(frac_below(s, R["fine"]), 6))
        best[a] = best_set(fine_clocks, R["fine"])
        groups_best = discordance_groups(specs[a], times, best[a][0])
        groups_wang = discordance_groups(specs[a], times, NAMED["wang_2024_T2008.5_r2.5"])
        groups_none = discordance_groups(specs[a], times, None)
        res["arrays"][a] = dict(
            n_pairs=len(specs[a]["pairs"]), n_comp=len(specs[a]["comps"]),
            s_noturn=round(R["s_noturn"], 6), s_star_test=round(R["s_star"], 6),
            improvement=round(R["improvement"], 6), s_star_fine=round(min(R["fine"]), 6),
            p_structure=round(p1, 6), p_turn=round(p2, 6), p_noturn_vs_null=round(p0, 6),
            null_s_star_pcts={q: round(sorted(stars)[int(q / 100.0 * (len(stars) - 1))], 6)
                              for q in (5, 50, 95)} if stars else None,
            best_test=describe_set(best_set(test_clocks, R["test"])),
            best_fine=describe_set(best[a]), named=named,
            W1=grade_w1(p1), W2=grade_w2(p2, R["improvement"]),
            W3=grade_w3(named["wang_2024_T2008.5_r2.5"]["frac_fine_below"]),
            per_multiplet={m: dict(n=groups_none[m][1], disc_noturn=groups_none[m][0],
                                   disc_wang=groups_wang.get(m, [None])[0],
                                   disc_best=groups_best.get(m, [None])[0]) for m in sorted(groups_none)},
            pooled_report_only=dict(n_comp=len(pooled[a]["comps"]),
                                    s_noturn=round(score(pooled[a], times, None), 6),
                                    s_star_test=round(min(scan(pooled[a], times, test_clocks)), 6),
                                    s_wang=round(score(pooled[a], times, NAMED["wang_2024_T2008.5_r2.5"]), 6)))
    if real["YKA"] is not None and real["ILAR"] is not None:
        py = cross_percentile(best["YKA"], fine_clocks, real["ILAR"]["fine"])
        pi = cross_percentile(best["ILAR"], fine_clocks, real["YKA"]["fine"])
        res["W4"] = dict(YKA_best_in_ILAR_pct=round(py, 6), ILAR_best_in_YKA_pct=round(pi, 6),
                         grade=grade_w4(py, pi))
    else:
        res["W4"] = dict(grade="untestable")
    profiles = {a: turn_profile(fine_clocks, real[a]["fine"]) for a in ARRAYS if real[a] is not None}
    return res, null, profiles


def write_csv(path, rows):
    keys = []
    for r in rows:
        for k in r:
            if k not in keys:
                keys.append(k)
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        w.writerows(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--table", default="", help="CSV dump of Supplementary Data Table 2 (else downloaded)")
    ap.add_argument("--data-dir", default="/data")
    ap.add_argument("--out-dir", default="/out")
    ap.add_argument("--shuffles", type=int, default=N_SHUFFLES)
    args = ap.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)
    if args.table:
        rows = parse_table2_csv(args.table)
    else:
        os.makedirs(args.data_dir, exist_ok=True)
        rows = fetch_table2(args.data_dir)
    res, null, profiles = analyse(rows, n_shuffles=args.shuffles)
    with open(os.path.join(args.out_dir, "gradeclock38_summary.json"), "w") as f:
        json.dump(res, f, indent=1)
    write_csv(os.path.join(args.out_dir, "gradeclock38_null.csv"),
              [dict(k=k, array=a, s_noturn=round(x[0], 6), s_star=round(x[1], 6))
               for a in ARRAYS for k, x in enumerate(null[a])])
    write_csv(os.path.join(args.out_dir, "gradeclock38_turn_profile.csv"),
              [dict(r, array=a) for a in profiles for r in profiles[a]])
    short = dict(checks_ok=res["checks_ok"], checks=res["checks"], W4=res["W4"],
                 arrays={a: {k: v for k, v in d.items() if k not in ("per_multiplet", "best_test", "best_fine")}
                         for a, d in res["arrays"].items()})
    print(json.dumps(short)[:6000])
    return 0


if __name__ == "__main__":
    sys.exit(main())

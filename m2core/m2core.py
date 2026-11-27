#!/usr/bin/env python3
"""m2core.py: does one inner-core rotation clock order every South Sandwich multiplet, at both arrays?

Part of inner-core-repeaters by MV2, Munim, Inc. (mv2@munimtech.com). Made with AI agents working
in MV2's name; no human expert has checked it. Licence: Apache-2.0.
Pre-registration 38, arm M (findings/2026-11-27-prereg38-armM-plan.md), frozen and committed before
any M2 waveform was fetched. This file is the decision core: it takes measured pair dissimilarities
D = 1 - cc (one row per event pair, multiplet and array) and applies the frozen rules. The
measurement front end (m2measure.py, pinned to m0324_check v0.8.0 + yka_repick 0.9.2) writes the
input CSV; it is written and frozen separately, before it runs.

Rigid rotation (Wang, Vidale et al. 2024, Nature 631:340, https://doi.org/10.1038/s41586-024-07536-4)
makes the inner core's angle one function of time. Then:
  M2a  one clock (turn date T, back-speed ratio r) must order the coda states of every multiplet at
       both arrays at once: sorted by clock position, each state occupies one contiguous run.
  M2b  a pair that re-matches at ILAR must also re-match at YKA (ILAR decorrelates more than YKA:
       Wang, Pang & Vidale 2026, EPSL 692:120274), outside the YKA timing-error eras.

Input CSV columns: array (YKA|ILAR), mult (Wang Table 3 letter), ev1, ev2, t1, t2 (decimal years),
D (1 - cc, late-coda window), ok (true/false: passed the pinned pipeline's quality gates).
Standard library only.
"""
import argparse
import calendar
import csv
import datetime
import json
import math
import os
import random
import sys

VERSION = "0.1.0"

# ---- pre-registered constants (prereg 38, arm M)
ARRAYS = ("YKA", "ILAR")
EXCLUDED_MULTIPLETS = ("I", "M", "-", "")     # Wang Table 3: I "not good repeaters", M "too noisy"
CUT_YKA = 0.098                               # P34 BETWEEN_MIN: same if D <= cut
CUT_ILAR_DEFAULT = 0.098
ILAR_CUT_MAX_DT = 2.0                         # ILAR cut = 75th pct (nearest rank) of usable
ILAR_CUT_Q = 0.75                             #   within-multiplet ILAR pairs < 2 yr apart,
ILAR_CUT_MIN_N = 8                            #   if at least 8 such pairs; else the default
HW_ERAS_ISO = (("2014-12-01", "2020-04-08"),  # YKA quantised element timing errors (Ackerley 2026)
               ("2021-12-01", "2022-04-01"))  # YKA B4/R4 non-quantised errors (Ackerley 2026)
SPLIT = 2009.0                                # a testable unit has members on both sides
MIN_MEMBERS = 3
MIN_STATES = 2
MIN_MULTIPLETS = 3
WIDE_RATIOS = [round(0.1 * 1.05 ** k, 4) for k in range(109)]          # 0.1 .. ~19.4 (clocks34)
FINE_TURNS = [round(2000.0 + 0.05 * k, 2) for k in range(401)]         # 2000.00 .. 2020.00
TEST_TURNS = [round(2000.0 + 0.25 * k, 2) for k in range(81)]          # 2000.00 .. 2020.00
TEST_RATIOS = WIDE_RATIOS[::4]                                          # 28 values
N_SHUFFLES = 1000
SEED = 3813
P_HOLDS = 0.05
CLOSURE_BROKEN_MIN = 2
CLOSURE_BROKEN_FRAC = 0.25
CLOSURE_HOLDS_MIN = 5
TIE = 1e-9
NAMED = {
    "wang_2024_T2008.5_r2.5": (2008.5, 2.5),
    "yang_song_2023_T2009.0_r1": (2009.0, 1.0),
    "yang_song_2023_T2009.5_r1": (2009.5, 1.0),
}


# ----------------------------------------------------------------------------- input
def decimal_year(date_str):
    d = datetime.date.fromisoformat(str(date_str).strip()[:10])
    n = 366 if calendar.isleap(d.year) else 365
    return d.year + (d.timetuple().tm_yday - 1) / float(n)


HW_ERAS = tuple((decimal_year(a), decimal_year(b)) for a, b in HW_ERAS_ISO)


def in_hw(t):
    return any(a <= t < b for a, b in HW_ERAS)


def truthy(x):
    return str(x).strip().lower() in ("1", "true", "yes", "y")


def fnum(x):
    try:
        v = float(x)
    except (TypeError, ValueError):
        return float("nan")
    return v


def read_pairs(path):
    with open(path, newline="", encoding="utf-8-sig") as f:
        return clean_pairs(list(csv.DictReader(f)))


def clean_pairs(records):
    out = []
    for r in records:
        a = str(r.get("array", "")).strip().upper()
        m = str(r.get("mult", "")).strip()
        e1, e2 = str(r.get("ev1", "")).strip(), str(r.get("ev2", "")).strip()
        if a not in ARRAYS or m in EXCLUDED_MULTIPLETS or not e1 or not e2 or e1 == e2:
            continue
        out.append(dict(array=a, mult=m, ev1=e1, ev2=e2, t1=fnum(r["t1"]), t2=fnum(r["t2"]),
                        D=fnum(r.get("D")), ok=truthy(r.get("ok", ""))))
    return out


def usable(p):
    return p["ok"] and math.isfinite(p["D"])


def event_times(pairs):
    times = {}
    for p in pairs:
        for e, t in ((p["ev1"], p["t1"]), (p["ev2"], p["t2"])):
            if e in times and abs(times[e] - t) > 1e-6:
                raise ValueError("event %s has two dates" % e)
            times[e] = t
    return times


def multiplet_events(pairs):
    out = {}
    for p in pairs:
        out.setdefault(p["mult"], set()).update((p["ev1"], p["ev2"]))
    return out


def ilar_cut(pairs):
    vals = sorted(p["D"] for p in pairs
                  if p["array"] == "ILAR" and usable(p) and abs(p["t1"] - p["t2"]) < ILAR_CUT_MAX_DT)
    if len(vals) < ILAR_CUT_MIN_N:
        return CUT_ILAR_DEFAULT, len(vals)
    k = int(math.ceil(ILAR_CUT_Q * len(vals))) - 1
    return vals[k], len(vals)


# ----------------------------------------------------------------------------- states
def unit_data(pairs, array, mult, exclude_hw=True, drop=None):
    """Usable pairs of one multiplet at one array -> ({event: t}, {(a, b): D})."""
    times, D = {}, {}
    for p in pairs:
        if p["array"] != array or p["mult"] != mult or not usable(p):
            continue
        if drop is not None and drop in (p["ev1"], p["ev2"]):
            continue
        if exclude_hw and array == "YKA" and (in_hw(p["t1"]) or in_hw(p["t2"])):
            continue
        key = tuple(sorted((p["ev1"], p["ev2"])))
        D.setdefault(key, p["D"])
        times[p["ev1"]], times[p["ev2"]] = p["t1"], p["t2"]
    return times, D


def complete_linkage(members, D, cut):
    """Agglomerative clustering; linkage = max D over the MEASURED cross pairs (none measured:
    infinite). Merge while the smallest linkage is <= cut; ties go to the first pair in sorted
    order. Returns {event: state}, states numbered by their first member in sorted order."""
    clusters = [[m] for m in sorted(members)]
    while len(clusters) > 1:
        best = None
        for i in range(len(clusters)):
            for j in range(i + 1, len(clusters)):
                link, seen = 0.0, False
                for a in clusters[i]:
                    for b in clusters[j]:
                        d = D.get((a, b) if a < b else (b, a))
                        if d is not None:
                            seen = True
                            link = max(link, d)
                if not seen:
                    link = float("inf")
                if best is None or link < best[0] - TIE:
                    best = (link, i, j)
        if best[0] > cut + TIE:
            break
        _, i, j = best
        clusters[i] = sorted(clusters[i] + clusters[j])
        del clusters[j]
    clusters.sort(key=lambda c: c[0])
    return {m: k for k, c in enumerate(clusters) for m in c}


def build_units(pairs, cuts, exclude_hw=True, drop=None, arrays=ARRAYS, skip_mults=()):
    mults = sorted(set(p["mult"] for p in pairs))
    units = []
    for a in arrays:
        for m in mults:
            if m in skip_mults:
                continue
            times, D = unit_data(pairs, a, m, exclude_hw, drop)
            if not times:
                continue
            states = complete_linkage(list(times), D, cuts[a])
            ts = list(times.values())
            n_states = len(set(states.values()))
            testable = (len(times) >= MIN_MEMBERS and n_states >= MIN_STATES
                        and min(ts) < SPLIT <= max(ts))
            units.append(dict(array=a, mult=m, items=sorted((e, states[e]) for e in times),
                              n_members=len(times), n_states=n_states, testable=testable))
    return units


# ----------------------------------------------------------------------------- clocks
def equiv(t, turn, ratio):
    """Effective position on a one-turn clock (clocks34.py / m34.py)."""
    return t if t <= turn else turn - (t - turn) / ratio


def position(t, clock):
    return t if clock is None else equiv(t, clock[0], clock[1])


def grid_clocks(turns, ratios):
    return [None] + [(T, r) for T in turns for r in ratios]


def ordered(items, times, clock):
    """True if, sorted by clock position (ties by date), each state is one contiguous run."""
    seq = sorted((position(times[e], clock), times[e], s) for e, s in items)
    runs, prev = 1, seq[0][2]
    for x in seq[1:]:
        if x[2] != prev:
            runs += 1
            prev = x[2]
    return runs == len(set(s for _, s in items))


def intersection(units, times, clocks):
    return [c for c in clocks if all(ordered(u["items"], times, c) for u in units)]


def any_clock(units, times, clocks):
    for c in clocks:
        if all(ordered(u["items"], times, c) for u in units):
            return True
    return False


def shuffle_times(times, mult_events, rng):
    new = dict(times)
    for m in sorted(mult_events):
        evs = sorted(mult_events[m])
        ts = [times[e] for e in evs]
        rng.shuffle(ts)
        for e, t in zip(evs, ts):
            new[e] = t
    return new


# ----------------------------------------------------------------------------- grades
def grade_m2a(n_mult, i_test_nonempty, noturn_in_i, p_consist, i_fine_empty, drop_robust):
    if n_mult < MIN_MULTIPLETS:
        return "untestable"
    if i_fine_empty and drop_robust:
        return "broken"
    if i_test_nonempty and noturn_in_i:
        return "monotone"
    if i_test_nonempty and p_consist <= P_HOLDS:
        return "holds"
    return "mixed"


def grade_m2b(n, v):
    if v >= CLOSURE_BROKEN_MIN and v >= CLOSURE_BROKEN_FRAC * n:
        return "broken"
    if n >= CLOSURE_HOLDS_MIN and v == 0:
        return "holds"
    if n < CLOSURE_HOLDS_MIN:
        return "untestable"
    return "mixed"


def closure_test(pairs, cuts):
    idx = {}
    for p in pairs:
        if usable(p):
            key = (p["mult"],) + tuple(sorted((p["ev1"], p["ev2"])))
            idx.setdefault(key, {}).setdefault(p["array"], p)
    rows = []
    n = v = hw_n = hw_v = rev = 0
    for key in sorted(idx):
        rec = idx[key]
        if "YKA" not in rec or "ILAR" not in rec:
            continue
        py, pi = rec["YKA"], rec["ILAR"]
        hw = in_hw(py["t1"]) or in_hw(py["t2"])
        ilar_same = pi["D"] <= cuts["ILAR"] + TIE
        yka_same = py["D"] <= CUT_YKA + TIE
        if ilar_same:
            if hw:
                hw_n += 1
                hw_v += 0 if yka_same else 1
            else:
                n += 1
                v += 0 if yka_same else 1
        elif yka_same and not hw:
            rev += 1
        rows.append(dict(mult=key[0], ev1=key[1], ev2=key[2], t1=round(py["t1"], 4), t2=round(py["t2"], 4),
                         D_YKA=round(py["D"], 6), D_ILAR=round(pi["D"], 6), ilar_same=ilar_same,
                         yka_same=yka_same, hw_era=hw,
                         violation=ilar_same and not yka_same and not hw))
    return dict(n_ilar_same=n, violations=v, hw_ilar_same=hw_n, hw_violations=hw_v,
                yka_same_ilar_changed=rev, grade=grade_m2b(n, v)), rows


# ----------------------------------------------------------------------------- whole analysis
def evaluate(pairs, cuts, times, test_clocks, fine_clocks, exclude_hw=True, arrays=ARRAYS, skip_mults=()):
    units = build_units(pairs, cuts, exclude_hw, arrays=arrays, skip_mults=skip_mults)
    tu = [u for u in units if u["testable"]]
    n_mult = len(set(u["mult"] for u in tu))
    I_test = intersection(tu, times, test_clocks) if tu else []
    fine_nonempty = any_clock(tu, times, fine_clocks) if tu else False
    named = {k: (all(ordered(u["items"], times, c) for u in tu) if tu else None) for k, c in NAMED.items()}
    return dict(units=units, testable=tu, n_mult=n_mult, I_test=I_test, fine_nonempty=fine_nonempty,
                noturn_in=(None in I_test), named=named)


def describe_clocks(clocks):
    turns = [c[0] for c in clocks if c is not None]
    ratios = [c[1] for c in clocks if c is not None]
    return dict(n=len(clocks), includes_noturn=None in clocks,
                turn_range=[min(turns), max(turns)] if turns else None,
                ratio_range=[min(ratios), max(ratios)] if ratios else None)


def drop_one(pairs, cuts, times, fine_clocks, events, exclude_hw=True):
    out = []
    for e in events:
        tu = [u for u in build_units(pairs, cuts, exclude_hw, drop=e) if u["testable"]]
        if len(set(u["mult"] for u in tu)) < MIN_MULTIPLETS:
            out.append(dict(dropped=e, result="untestable"))
        else:
            out.append(dict(dropped=e, result="nonempty" if any_clock(tu, times, fine_clocks) else "empty"))
    return out


def analyse(pairs, n_shuffles=N_SHUFFLES, seed=SEED, test_clocks=None, fine_clocks=None, variants=True):
    test_clocks = test_clocks or grid_clocks(TEST_TURNS, TEST_RATIOS)
    fine_clocks = fine_clocks or grid_clocks(FINE_TURNS, WIDE_RATIOS)
    times = event_times(pairs)
    mev = multiplet_events(pairs)
    seen, overlap = set(), []
    for m in sorted(mev):
        overlap += sorted(mev[m] & seen)
        seen |= mev[m]
    cut_i, n_short = ilar_cut(pairs)
    cuts = {"YKA": CUT_YKA, "ILAR": cut_i}
    checks = dict(n_rows=len(pairs), n_usable=sum(1 for p in pairs if usable(p)),
                  n_usable_YKA=sum(1 for p in pairs if usable(p) and p["array"] == "YKA"),
                  n_usable_ILAR=sum(1 for p in pairs if usable(p) and p["array"] == "ILAR"),
                  multiplet_overlap=overlap)
    checks_ok = checks["n_usable"] > 0 and not overlap

    ev = evaluate(pairs, cuts, times, test_clocks, fine_clocks)
    tu = ev["testable"]
    unit_rows = []
    for u in ev["units"]:
        row = dict(array=u["array"], mult=u["mult"], n_members=u["n_members"], n_states=u["n_states"],
                   testable=u["testable"],
                   states=" ".join("%s:%d@%.3f" % (e, s, times[e]) for e, s in u["items"]))
        if u["testable"]:
            row["share_test"] = round(len(intersection([u], times, test_clocks)) / float(len(test_clocks)), 6)
            row["share_fine"] = round(len(intersection([u], times, fine_clocks)) / float(len(fine_clocks)), 6)
            row["noturn_ordered"] = ordered(u["items"], times, None)
            for k, c in NAMED.items():
                row[k] = ordered(u["items"], times, c)
        unit_rows.append(row)

    drops, robust = [], None
    if ev["n_mult"] >= MIN_MULTIPLETS and not ev["fine_nonempty"]:
        evs = sorted(set(e for u in tu for e, _ in u["items"]))
        drops = drop_one(pairs, cuts, times, fine_clocks, evs)
        robust = all(d["result"] == "empty" for d in drops)

    rng = random.Random(seed)
    null = []
    if tu:
        for k in range(n_shuffles):
            sh = shuffle_times(times, mev, rng)
            I = intersection(tu, sh, test_clocks)
            null.append(dict(k=k, n_clocks=len(I), noturn_in=(None in I)))
    n_nonempty = sum(1 for x in null if x["n_clocks"] > 0)
    p_consist = (1 + n_nonempty) / (1.0 + len(null)) if null else None

    grade = grade_m2a(ev["n_mult"], len(ev["I_test"]) > 0, ev["noturn_in"],
                      p_consist if p_consist is not None else 1.0, not ev["fine_nonempty"], bool(robust))
    closure, closure_rows = closure_test(pairs, cuts)

    res = dict(version=VERSION, checks=checks, checks_ok=checks_ok, cuts=cuts, n_short_ilar_pairs=n_short,
               n_test_clocks=len(test_clocks), n_fine_clocks=len(fine_clocks), n_shuffles=len(null), seed=seed,
               M2a=dict(grade=grade, n_testable_units=len(tu), n_testable_multiplets=ev["n_mult"],
                        units_empty_fine=[(u["array"], u["mult"]) for u, r in zip(ev["units"], unit_rows)
                                          if u["testable"] and r["share_fine"] == 0],
                        I_test=describe_clocks(ev["I_test"]), I_test_share=round(len(ev["I_test"]) / float(len(test_clocks)), 6),
                        I_fine_nonempty=ev["fine_nonempty"], noturn_in_I=ev["noturn_in"],
                        p_consist=round(p_consist, 6) if p_consist is not None else None,
                        null_nonempty=n_nonempty, drop_one_robust=robust,
                        drop_one_breakers=[d["dropped"] for d in drops if d["result"] != "empty"],
                        named_in_I=ev["named"]),
               M2b=closure)
    if variants:
        rep = {}
        specs = {"yka_only": dict(arrays=("YKA",)), "ilar_only": dict(arrays=("ILAR",)),
                 "with_hw_era": dict(exclude_hw=False), "without_A": dict(skip_mults=("A",))}
        for name, kw in specs.items():
            v = evaluate(pairs, cuts, times, test_clocks, fine_clocks, **kw)
            rep[name] = dict(n_testable_multiplets=v["n_mult"], I_test_share=round(len(v["I_test"]) / float(len(test_clocks)), 6),
                             I_fine_nonempty=v["fine_nonempty"], noturn_in_I=v["noturn_in"], named_in_I=v["named"])
        for f in (0.8, 1.2):
            c2 = {a: cuts[a] * f for a in cuts}
            v = evaluate(pairs, c2, times, test_clocks, fine_clocks)
            rep["cut_x%.1f" % f] = dict(n_testable_multiplets=v["n_mult"], I_test_share=round(len(v["I_test"]) / float(len(test_clocks)), 6),
                                        I_fine_nonempty=v["fine_nonempty"], noturn_in_I=v["noturn_in"], named_in_I=v["named"])
        res["report_only_variants"] = rep
    return res, unit_rows, null, drops, closure_rows


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
    ap.add_argument("--pairs", required=True, help="m2measure output: array,mult,ev1,ev2,t1,t2,D,ok")
    ap.add_argument("--out-dir", default="/out")
    ap.add_argument("--shuffles", type=int, default=N_SHUFFLES)
    args = ap.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)
    pairs = read_pairs(args.pairs)
    res, units, null, drops, closure_rows = analyse(pairs, n_shuffles=args.shuffles)
    with open(os.path.join(args.out_dir, "m2_summary.json"), "w") as f:
        json.dump(res, f, indent=1)
    write_csv(os.path.join(args.out_dir, "m2_units.csv"), units)
    write_csv(os.path.join(args.out_dir, "m2_null.csv"), null)
    write_csv(os.path.join(args.out_dir, "m2_drop_one.csv"), drops)
    write_csv(os.path.join(args.out_dir, "m2_closure_pairs.csv"), closure_rows)
    print(json.dumps(dict(checks_ok=res["checks_ok"], M2a=res["M2a"], M2b=res["M2b"]))[:6000])
    return 0


if __name__ == "__main__":
    sys.exit(main())

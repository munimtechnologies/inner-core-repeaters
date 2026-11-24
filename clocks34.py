#!/usr/bin/env python3
"""clocks34.py: which inner-core rotation clocks can order the YKA coda states of multiplet M03-24?

Part of inner-core-repeaters by MV2, Munim, Inc. (made with AI agents). Licence: Apache-2.0.
Amendment 1 to pre-registration 34 (findings/2026-11-24-prereg34-amendment1-clock-families.md),
written and published before job m34-mirror-v1 produced any result.

If the late PKIKP coda state at YKA (X or Y, optionally intermediate I) is set only by the inner
core's angular position, and the structure has one boundary between the X and Y regions, then the
members' positions must sort as X < I < Y (or the mirror image). This module asks which clock
families allow that, for the old eight members (exploratory: their labels are already known) and,
after m34-mirror-v1, for old plus new members (pre-registered):

  one-turn clocks   forward until `turn`, then back `ratio` times slower (m34's `equiv`)
  symmetric clocks  one turn, equal speed both ways (ratio 1; also any long-period sinusoid)
  periodic clocks   position = A sin(2 pi t / P + c) with free phase c and no drift

Standard library only.
"""
import argparse
import calendar
import csv
import json
import math
import os
import sys

VERSION = "0.1.0"

X_GROUP = ("usp000bt5x", "us7000am0e", "us6000h2wz", "us6000mk65", "us6000runq")
Y_GROUP = ("usp000fz4d", "usp000hute", "usc000t4wj")
OLD_TIMES = {
    "usp000bt5x": "2003-03-21T10:47:47.620Z", "usp000fz4d": "2008-02-17T19:32:09.100Z",
    "usp000hute": "2011-02-20T17:37:37.320Z", "usc000t4wj": "2014-12-09T19:23:36.400Z",
    "us7000am0e": "2020-07-11T18:58:06.501Z", "us6000h2wz": "2022-03-07T23:25:40.208Z",
    "us6000mk65": "2024-03-21T02:47:21.923Z", "us6000runq": "2025-12-13T14:04:33.132Z"}
ORDER = {"X": 0, "I": 1, "Y": 2}
STATE_CODE = {"X": "X", "Y": "Y", "intermediate": "I"}

# ---- pre-registered grids (amendment 1)
WANG_TURN, WANG_RATIO = 2008.5, 2.5
P_MIN, P_MAX, P_STEP = 2.0, 100.0, 0.01
PERIODS = [round(P_MIN + P_STEP * k, 2) for k in range(int(round((P_MAX - P_MIN) / P_STEP)) + 1)]
WIDE_TURNS = [round(2003.0 + 0.05 * k, 2) for k in range(361)]          # 2003.00 .. 2021.00
WIDE_RATIOS = [round(0.1 * 1.05 ** k, 4) for k in range(109)]            # 0.1 .. ~19.4
SYM_TURNS = [round(2003.0 + 0.001 * k, 3) for k in range(18001)]         # 2003.000 .. 2021.000
BANDS = [(2.0, 5.0), (5.0, 7.0), (7.0, 10.0), (10.0, 15.0), (15.0, 18.0), (18.0, 30.0), (30.0, 100.01)]
NAMED_PERIOD_BANDS = {"wang_vidale_2022_6yr": (5.5, 6.5), "tkalcic_20_30yr": (20.0, 30.0)}


def decimal_year(iso):
    """Same convention as z004_events.csv and m34.py."""
    s = iso.strip().rstrip("Z")
    date, hms = s.split("T")
    y, mo, d = (int(v) for v in date.split("-"))
    h, mi, sec = hms.split(":")
    t0 = calendar.timegm((y, 1, 1, 0, 0, 0, 0, 0, 0))
    t = calendar.timegm((y, mo, d, int(h), int(mi), 0, 0, 0, 0)) + float(sec)
    days = 366 if calendar.isleap(y) else 365
    return y + (t - t0) / 86400.0 / days


def old_members():
    m = {e: (decimal_year(OLD_TIMES[e]), "X") for e in X_GROUP}
    m.update({e: (decimal_year(OLD_TIMES[e]), "Y") for e in Y_GROUP})
    return m


# ----------------------------------------------------------------------------- one-turn clocks
def equiv(t, turn, ratio):
    """Effective date (position) on a one-turn clock, as in m34.py."""
    return t if t <= turn else turn - (t - turn) / ratio


def ordered_ok(values):
    """values: list of (position, rank). True if positions sort strictly by rank (X < I < Y)."""
    for k in (1, 2):
        lo = [v for v, r in values if r < k]
        hi = [v for v, r in values if r >= k]
        if lo and hi and not max(lo) < min(hi):
            return False
    return True


def one_turn_ok(members, turn, ratio):
    return ordered_ok([(equiv(t, turn, ratio), ORDER[c]) for t, c in members.values()])


def intervals(xs, flags):
    """Merge consecutive grid points with flag True into [first, last] intervals."""
    out, start, prev = [], None, None
    for x, f in zip(xs, flags):
        if f and start is None:
            start = x
        if not f and start is not None:
            out.append([start, prev])
            start = None
        prev = x
    if start is not None:
        out.append([start, prev])
    return out


def symmetric_window(members, turns=SYM_TURNS):
    return intervals(turns, [one_turn_ok(members, T, 1.0) for T in turns])


def wide_grid(members, turns=WIDE_TURNS, ratios=WIDE_RATIOS):
    rows, n = [], 0
    for T in turns:
        ok = [r for r in ratios if one_turn_ok(members, T, r)]
        n += len(ok)
        rows.append(dict(turn=T, n_ok=len(ok), ratio_min=min(ok) if ok else None,
                         ratio_max=max(ok) if ok else None))
    return n, rows


# ----------------------------------------------------------------------------- periodic clocks
def phases(members, P):
    return [((t - 2000.0) % P, ORDER[c]) for t, c in members.values()]


def circ(a, b, P):
    d = abs(a - b) % P
    return min(d, P - d)


def contiguous_y(pts):
    """X/Y only: True if the Y phases are consecutive around the circle (no X between them)."""
    s = sorted(pts)
    codes = [r for _, r in s]
    if all(c == 2 for c in codes) or all(c != 2 for c in codes):
        return True
    starts = sum(1 for i in range(len(codes)) if codes[i] == 2 and codes[i - 1] != 2)
    # a tie in phase between an X and a Y cannot be separated
    for i in range(len(s)):
        if s[i][0] == s[i - 1][0] and s[i][1] != s[i - 1][1]:
            return False
    return starts == 1


def nested_ok(pts, P):
    """Exact test for nested arcs with a common centre (Y inside, I in the ring, X outside).
    The truth of every inequality d(c,a) < d(c,b) can only change where c is the midpoint of a and b
    or its antipode, so testing one centre inside each cell between those points is exact."""
    crit = set()
    for i in range(len(pts)):
        for j in range(i + 1, len(pts)):
            m = ((pts[i][0] + pts[j][0]) / 2.0) % P
            crit.add(m)
            crit.add((m + P / 2.0) % P)
    crit = sorted(crit)
    if not crit:
        return True
    tests = [(crit[i] + crit[i + 1]) / 2.0 for i in range(len(crit) - 1)]
    tests.append(((crit[-1] + crit[0] + P) / 2.0) % P)
    for c in tests:
        if ordered_ok([(-circ(p, c, P), r) for p, r in pts]):   # closer to c = higher rank
            return True
    return False


def periodic_ok(members, P):
    pts = phases(members, P)
    if all(r != 1 for _, r in pts):
        return contiguous_y(pts)
    return nested_ok(pts, P)


def periodic_scan(members, periods=PERIODS):
    return [periodic_ok(members, P) for P in periods]


def continuum_from(periods, flags):
    """Smallest grid period from which every larger grid period is compatible (None if the last isn't)."""
    if not flags or not flags[-1]:
        return None
    k = len(flags) - 1
    while k > 0 and flags[k - 1]:
        k -= 1
    return periods[k]


def band_fractions(periods, flags, bands=BANDS):
    out = {}
    for lo, hi in bands:
        sel = [f for P, f in zip(periods, flags) if lo <= P < hi]
        out["%g-%g" % (lo, min(hi, 100.0))] = round(sum(sel) / float(len(sel)), 4) if sel else None
    return out


def chance_fraction(n_x, n_y):
    """Share of circular arrangements of n_x X and n_y Y in which the Y are consecutive."""
    n = n_x + n_y
    if n_y == 0 or n_x == 0:
        return 1.0
    return n / float(math.comb(n, n_y))


# ----------------------------------------------------------------------------- whole analysis
def analyse(members):
    flags = periodic_scan(members)
    n_wide, wide_rows = wide_grid(members)
    named = {}
    for name, (lo, hi) in NAMED_PERIOD_BANDS.items():
        sel = [f for P, f in zip(PERIODS, flags) if lo <= P <= hi]
        named[name] = dict(band=[lo, hi], fraction_compatible=round(sum(sel) / float(len(sel)), 4))
    named["wang_vidale_2022_P6.00"] = periodic_ok(members, 6.0)
    named["wang_2024_turn2008.5_ratio2.5"] = one_turn_ok(members, WANG_TURN, WANG_RATIO)
    named["yang_song_2023_symmetric_turn2009.0"] = one_turn_ok(members, 2009.0, 1.0)
    named["yang_song_2023_symmetric_turn2009.5"] = one_turn_ok(members, 2009.5, 1.0)
    nx = sum(1 for _, c in members.values() if c == "X")
    ny = sum(1 for _, c in members.values() if c == "Y")
    return dict(n_members=len(members), n_periods=len(PERIODS), n_periodic_ok=sum(flags),
                periodic_intervals=intervals(PERIODS, flags),
                continuum_from=continuum_from(PERIODS, flags),
                band_fractions=band_fractions(PERIODS, flags),
                chance_fraction_xy=round(chance_fraction(nx, ny), 6),
                symmetric_window=symmetric_window(members),
                wide_n_clocks=len(WIDE_TURNS) * len(WIDE_RATIOS), wide_n_ok=n_wide,
                named=named), flags, wide_rows


def drop_one(members, keys):
    """Per dropped key: (periodic n_ok, wide n_ok, symmetric window)."""
    out = {}
    for k in keys:
        sub = {e: v for e, v in members.items() if e != k}
        out[k] = dict(periodic_n_ok=sum(periodic_scan(sub)), wide_n_ok=wide_grid(sub)[0],
                      symmetric_window=symmetric_window(sub), continuum_from=continuum_from(
                          PERIODS, periodic_scan(sub)))
    return out


def load_new(path):
    """Confirmed new members with a usable state from m34_candidates.csv (pre-registration 34)."""
    new = {}
    with open(path, newline="") as f:
        for r in csv.DictReader(f):
            if str(r.get("used_in_mirror", "")).strip() not in ("True", "true", "1"):
                continue
            st = (r.get("state") or "").strip()
            if st in STATE_CODE:
                new[r["id"]] = (float(r["t"]), STATE_CODE[st])
    return new


def grade(old_res, all_res, n_new, drop):
    """Pre-registered amendment-1 labels (only when n_new >= 1)."""
    if n_new == 0:
        return dict(R3="untestable", R4="untestable", R5="untestable")
    r3 = "compatible" if all_res["n_periodic_ok"] else (
        "excluded-fragile" if any(v["periodic_n_ok"] for v in drop.values()) else "excluded")
    r4 = "compatible" if all_res["wide_n_ok"] else (
        "contradicted-fragile" if any(v["wide_n_ok"] for v in drop.values()) else "contradicted")
    ow, aw = old_res["symmetric_window"], all_res["symmetric_window"]
    if not aw:
        r5 = "closed"
    elif aw == ow:
        r5 = "unchanged"
    else:
        r5 = "narrowed"
    return dict(R3=r3, R4=r4, R5=r5,
                continuum_from_old=old_res["continuum_from"], continuum_from_all=all_res["continuum_from"],
                continuum_rose=(all_res["continuum_from"] or 1e9) > (old_res["continuum_from"] or 1e9))


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
    ap.add_argument("--candidates", default="", help="m34_candidates.csv from m34-mirror-v1 (optional)")
    ap.add_argument("--out-dir", default="/out")
    args = ap.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)
    old = old_members()
    old_res, old_flags, old_wide = analyse(old)
    s = dict(version=VERSION, old8=old_res, old8_drop_one=drop_one(old, list(old)))
    rows = [dict(P=P, ok_old8=f) for P, f in zip(PERIODS, old_flags)]
    wide = [dict(r, set="old8") for r in old_wide]
    if args.candidates and os.path.exists(args.candidates):
        new = load_new(args.candidates)
        allm = dict(old)
        allm.update(new)
        all_res, all_flags, all_wide = analyse(allm)
        drop = drop_one(allm, list(new))
        s.update(new_members={k: list(v) for k, v in new.items()}, all=all_res, all_drop_one_new=drop,
                 grade=grade(old_res, all_res, len(new), drop))
        for r, f in zip(rows, all_flags):
            r["ok_all"] = f
        wide += [dict(r, set="all") for r in all_wide]
    write_csv(os.path.join(args.out_dir, "clocks34_periods.csv"), rows)
    write_csv(os.path.join(args.out_dir, "clocks34_oneturn.csv"), wide)
    with open(os.path.join(args.out_dir, "clocks34_summary.json"), "w") as f:
        json.dump(s, f, indent=1)
    print(json.dumps({k: v for k, v in s.items() if k != "old8_drop_one"})[:6000])
    return 0


if __name__ == "__main__":
    sys.exit(main())

#!/usr/bin/env python3
"""wsy1b.py: can the M03-24 YKA coda states be ordered by Wu, Song & Yang's (2026) rotation rates?

Part of inner-core-repeaters by MV2, Munim, Inc. (made with AI agents). Licence: Apache-2.0.
Amendment 1b to pre-registration 34 (findings/2026-11-24-prereg34-amendment1b-wsy-rates.md),
written before job m34-mirror-v1 produced any result. Exploratory: the eight labels are known.

Clock: the inner core moves forward (toward the Y region) at a constant rate f (deg/yr) from
2003-03-21 to the turn T, then back at a constant rate b. With one X/Y boundary and X before Y,
the eight states fit iff
  C1  f (T - t_2003) > b (t_2014 - T)      (Dec 2014 still past Mar 2003)
  C2  b (t_2020 - T) > f (T - t_2008)      (Jul 2020 back past Feb 2008)
(C1 is automatic for T >= t_2014, C2 for T <= t_2008.) Only the ratio f/b matters, which is
clocks34's return `ratio`. Standard library only.
"""
import calendar
import csv
import sys

VERSION = "0.1.0"

X_TIMES = ("2003-03-21T10:47:47.620Z", "2020-07-11T18:58:06.501Z", "2022-03-07T23:25:40.208Z",
           "2024-03-21T02:47:21.923Z", "2025-12-13T14:04:33.132Z")
Y_TIMES = ("2008-02-17T19:32:09.100Z", "2011-02-20T17:37:37.320Z", "2014-12-09T19:23:36.400Z")

# rates read from Wu, Song & Yang (2026, GRL 53:e2025GL119295), Section 4 and Discussion
WSY_B_MEAN, WSY_B_MAX = 0.01, 0.03     # westward (backward) rate, "average" and "maximum", deg/yr
WSY_PEAK = 0.10                        # forward rate from the early 1980s, deg/yr
TURNS = [2004.0, 2005.0, 2006.0, 2007.0, 2008.0, 2009.0, 2010.0, 2011.0, 2012.0]
NEAR_ZERO = (0.01, 0.02, 0.03)         # candidate readings of "near-zero" forward rate, deg/yr


def decimal_year(iso):
    """Same convention as clocks34.py / m34.py."""
    s = iso.strip().rstrip("Z")
    date, hms = s.split("T")
    y, mo, d = (int(v) for v in date.split("-"))
    h, mi, sec = hms.split(":")
    t0 = calendar.timegm((y, 1, 1, 0, 0, 0, 0, 0, 0))
    t = calendar.timegm((y, mo, d, int(h), int(mi), 0, 0, 0, 0)) + float(sec)
    days = 366 if calendar.isleap(y) else 365
    return y + (t - t0) / 86400.0 / days


T03 = decimal_year(X_TIMES[0])
T08 = decimal_year(Y_TIMES[0])
T14 = decimal_year(Y_TIMES[2])
T20 = decimal_year(X_TIMES[1])


def position(t, turn, f, b):
    """Angle relative to the 2003 member (deg) on the forward-then-back clock."""
    if t <= turn:
        return f * (t - T03)
    return f * (turn - T03) - b * (t - turn)


def all_eight_ok(turn, f, b):
    """Direct check on all eight members: every X strictly below every Y."""
    xs = [position(decimal_year(s), turn, f, b) for s in X_TIMES]
    ys = [position(decimal_year(s), turn, f, b) for s in Y_TIMES]
    return max(xs) < min(ys)


def f_min_c1(turn, b):
    """Smallest forward rate that satisfies C1 (0 if C1 holds for any f > 0)."""
    if turn >= T14:
        return 0.0
    if turn <= T03:
        return float("inf")
    return b * (T14 - turn) / (turn - T03)


def f_max_c2(turn, b):
    """Largest forward rate that satisfies C2 (inf if C2 holds for any f)."""
    if turn <= T08:
        return float("inf")
    return b * (T20 - turn) / (turn - T08)


def earliest_turn(f, b):
    """Earliest turn at which forward rate f and backward rate b satisfy C1 (closed form)."""
    # f (T - T03) = b (T14 - T)  ->  T = (f T03 + b T14) / (f + b)
    return (f * T03 + b * T14) / (f + b)


def table():
    rows = []
    for T in TURNS:
        for b in (WSY_B_MEAN, WSY_B_MAX):
            lo, hi = f_min_c1(T, b), f_max_c2(T, b)
            rows.append(dict(turn=T, b_deg_per_yr=b, f_min_c1=round(lo, 5),
                             f_max_c2="" if hi == float("inf") else round(hi, 5),
                             f_min_over_wsy_peak=round(lo / WSY_PEAK, 4)))
    return rows


def earliest_rows():
    rows = []
    for b in (WSY_B_MEAN, WSY_B_MAX):
        for f in NEAR_ZERO:
            rows.append(dict(b_deg_per_yr=b, f_deg_per_yr=f, earliest_turn=round(earliest_turn(f, b), 4)))
    return rows


def write_csv(path, rows):
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)


def main(argv):
    out = argv[1] if len(argv) > 1 else "."
    write_csv(out + "/wsy1b_turn_table.csv", table())
    write_csv(out + "/wsy1b_earliest_turn.csv", earliest_rows())
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))

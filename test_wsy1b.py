"""Tests for wsy1b.py (amendment 1b of pre-registration 34). Expected values were worked by hand.
Part of inner-core-repeaters by MV2, Munim, Inc. (made with AI agents). Licence: Apache-2.0."""
import csv
import os
import unittest

import wsy1b as w

HERE = os.path.dirname(os.path.abspath(__file__))
CSV_DIRS = [os.path.join(HERE, "results", "wsy1b"), os.path.join(HERE, "results")]

# hand-computed table: (turn, b, f_min_c1, f_max_c2 or None)
EXPECTED = [
    (2004.0, 0.01, 0.13983, None), (2004.0, 0.03, 0.41949, None),
    (2005.0, 0.01, 0.05577, None), (2005.0, 0.03, 0.16730, None),
    (2006.0, 0.01, 0.03213, None), (2006.0, 0.03, 0.09639, None),
    (2007.0, 0.01, 0.02099, None), (2007.0, 0.03, 0.06297, None),
    (2008.0, 0.01, 0.01451, None), (2008.0, 0.03, 0.04353, None),
    (2009.0, 0.01, 0.01027, 0.13259), (2009.0, 0.03, 0.03081, 0.39777),
    (2010.0, 0.01, 0.00728, 0.05631), (2010.0, 0.03, 0.02185, 0.16894),
    (2011.0, 0.01, 0.00506, 0.03320), (2011.0, 0.03, 0.01519, 0.09960),
    (2012.0, 0.01, 0.00335, 0.02204), (2012.0, 0.03, 0.01004, 0.06611),
]
EXPECTED_EARLIEST = [(0.01, 0.01, 2009.0784), (0.01, 0.02, 2007.1248), (0.01, 0.03, 2006.1481),
                     (0.03, 0.01, 2012.0088), (0.03, 0.02, 2010.2506), (0.03, 0.03, 2009.0784)]


def close(a, b):
    return abs(a - b) <= max(2e-5, 1e-3 * abs(b))


class TestDates(unittest.TestCase):
    def test_decimal_years(self):
        self.assertAlmostEqual(w.T03, 2003.217671, places=5)
        self.assertAlmostEqual(w.T08, 2008.130639, places=5)
        self.assertAlmostEqual(w.T14, 2014.939200, places=5)
        self.assertAlmostEqual(w.T20, 2020.526750, places=5)


class TestConditions(unittest.TestCase):
    def test_turn_2010(self):
        self.assertAlmostEqual(w.f_min_c1(2010.0, 0.01), 0.0072825, places=6)
        self.assertAlmostEqual(w.f_max_c2(2010.0, 0.01), 0.056312, places=5)

    def test_matches_amendment1(self):
        # amendment 1: turn 2005 needs return > 5.58 times slower; turn 2012 ratio 0.335-2.204
        self.assertAlmostEqual(w.f_min_c1(2005.0, 1.0), 5.5765, places=3)
        self.assertAlmostEqual(w.f_min_c1(2012.0, 1.0), 0.33467, places=4)
        self.assertAlmostEqual(w.f_max_c2(2012.0, 1.0), 2.20366, places=4)

    def test_edges(self):
        self.assertEqual(w.f_min_c1(2015.0, 0.01), 0.0)
        self.assertEqual(w.f_min_c1(2003.0, 0.01), float("inf"))
        self.assertEqual(w.f_max_c2(2008.0, 0.01), float("inf"))

    def test_symmetric_earliest_turn(self):
        self.assertAlmostEqual(w.earliest_turn(0.01, 0.01), (w.T03 + w.T14) / 2.0, places=9)
        self.assertAlmostEqual(w.earliest_turn(0.02, 0.01), 2007.12485, places=4)

    def test_wsy_readings(self):
        self.assertTrue(w.all_eight_ok(2010.0, 0.02, 0.01))     # spline turn, modest forward rate
        self.assertFalse(w.all_eight_ok(2005.0, 0.02, 0.01))    # "recent 20 years" turn
        self.assertTrue(w.all_eight_ok(2005.0, 0.06, 0.01))     # needs > 0.0558 deg/yr in 2003-05
        self.assertFalse(w.all_eight_ok(2010.0, 0.01, 0.03))    # needs >= 0.0218 at the 0.03 max
        self.assertFalse(w.all_eight_ok(2010.0, 0.06, 0.01))    # C2 fails above 0.0563

    def test_direct_check_equals_c1_and_c2(self):
        n = 0
        for k in range(36):
            T = 2003.6 + 0.37 * k
            for f in (0.003, 0.0071, 0.013, 0.027, 0.051, 0.093, 0.17):
                for b in (0.01, 0.03):
                    c1 = f > w.f_min_c1(T, b)
                    c2 = f < w.f_max_c2(T, b)
                    self.assertEqual(w.all_eight_ok(T, f, b), c1 and c2, (T, f, b))
                    n += 1
        self.assertEqual(n, 36 * 7 * 2)


class TestTables(unittest.TestCase):
    def test_table_matches_hand_values(self):
        rows = w.table()
        self.assertEqual(len(rows), len(EXPECTED))
        for r, (T, b, lo, hi) in zip(rows, EXPECTED):
            self.assertEqual((r["turn"], r["b_deg_per_yr"]), (T, b))
            self.assertTrue(close(r["f_min_c1"], lo), (T, b, r["f_min_c1"], lo))
            if hi is None:
                self.assertEqual(r["f_max_c2"], "")
            else:
                self.assertTrue(close(r["f_max_c2"], hi), (T, b, r["f_max_c2"], hi))

    def test_earliest_matches_hand_values(self):
        for r, (b, f, T) in zip(w.earliest_rows(), EXPECTED_EARLIEST):
            self.assertEqual((r["b_deg_per_yr"], r["f_deg_per_yr"]), (b, f))
            self.assertTrue(abs(r["earliest_turn"] - T) < 2e-4, (b, f, r["earliest_turn"], T))

    def test_committed_csv_if_present(self):
        for d in CSV_DIRS:
            p = os.path.join(d, "wsy1b_turn_table.csv")
            if not os.path.exists(p):
                continue
            with open(p, newline="") as fh:
                got = list(csv.DictReader(fh))
            want = w.table()
            self.assertEqual(len(got), len(want))
            for g, r in zip(got, want):
                self.assertTrue(close(float(g["f_min_c1"]), r["f_min_c1"]))
                if r["f_max_c2"] == "":
                    self.assertEqual(g["f_max_c2"], "")
                else:
                    self.assertTrue(close(float(g["f_max_c2"]), r["f_max_c2"]))


if __name__ == "__main__":
    unittest.main()

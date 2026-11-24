"""Tests for clocks34.py (amendment 1 to pre-registration 34). Standard library only.
Every expected value below was worked out by hand before the code was run (see the amendment).
Part of inner-core-repeaters by MV2, Munim, Inc. (made with AI agents). Licence: Apache-2.0."""
import csv
import os
import tempfile
import unittest

import clocks34 as c


class TestDates(unittest.TestCase):
    def test_decimal_years(self):
        m = c.old_members()
        self.assertAlmostEqual(m["usp000bt5x"][0], 2003.21767, places=4)
        self.assertAlmostEqual(m["usp000fz4d"][0], 2008.13064, places=4)
        self.assertAlmostEqual(m["usp000hute"][0], 2011.13900, places=4)
        self.assertAlmostEqual(m["usc000t4wj"][0], 2014.93920, places=4)
        self.assertAlmostEqual(m["us7000am0e"][0], 2020.52675, places=4)
        self.assertAlmostEqual(m["us6000runq"][0], 2025.94955, places=4)
        self.assertEqual(sorted(v for _, v in m.values()).count("X"), 5)

    def test_grids(self):
        self.assertEqual(len(c.PERIODS), 9801)
        self.assertEqual((c.PERIODS[0], c.PERIODS[-1]), (2.0, 100.0))
        self.assertEqual(len(c.WIDE_TURNS), 361)
        self.assertEqual((c.WIDE_TURNS[0], c.WIDE_TURNS[-1]), (2003.0, 2021.0))
        self.assertEqual(len(c.WIDE_RATIOS), 109)
        self.assertEqual(c.WIDE_RATIOS[0], 0.1)
        self.assertTrue(19.0 < c.WIDE_RATIOS[-1] < 20.0)
        self.assertEqual((c.SYM_TURNS[0], c.SYM_TURNS[-1], len(c.SYM_TURNS)), (2003.0, 2021.0, 18001))


class TestOneTurn(unittest.TestCase):
    def setUp(self):
        self.m = c.old_members()

    def test_wang(self):
        self.assertTrue(c.one_turn_ok(self.m, 2008.5, 2.5))

    def test_turn_2005_needs_ratio_above_5_577(self):
        self.assertFalse(c.one_turn_ok(self.m, 2005.0, 5.5))
        self.assertTrue(c.one_turn_ok(self.m, 2005.0, 5.7))

    def test_turn_2012_ratio_window(self):
        self.assertFalse(c.one_turn_ok(self.m, 2012.0, 0.3))
        self.assertTrue(c.one_turn_ok(self.m, 2012.0, 0.4))
        self.assertTrue(c.one_turn_ok(self.m, 2012.0, 2.2))
        self.assertFalse(c.one_turn_ok(self.m, 2012.0, 2.25))

    def test_turn_before_2003_impossible(self):
        for r in (0.2, 1.0, 3.0, 10.0):
            self.assertFalse(c.one_turn_ok(self.m, 2002.0, r))

    def test_symmetric_edges(self):
        self.assertFalse(c.one_turn_ok(self.m, 2009.0, 1.0))
        self.assertFalse(c.one_turn_ok(self.m, 2009.05, 1.0))
        self.assertTrue(c.one_turn_ok(self.m, 2009.10, 1.0))
        self.assertTrue(c.one_turn_ok(self.m, 2014.30, 1.0))
        self.assertFalse(c.one_turn_ok(self.m, 2014.36, 1.0))

    def test_symmetric_window_coarse(self):
        turns = [round(2009.0 + 0.01 * k, 2) for k in range(551)]
        self.assertEqual(c.symmetric_window(self.m, turns), [[2009.08, 2014.32]])


class TestPeriodic(unittest.TestCase):
    def setUp(self):
        self.m = c.old_members()

    def test_excluded_periods(self):
        for P in (6.0, 10.0, 12.0, 17.0, 17.5, 17.8):
            self.assertFalse(c.periodic_ok(self.m, P), P)

    def test_allowed_periods(self):
        for P in (17.85, 18.0, 20.0, 40.0, 70.0):
            self.assertTrue(c.periodic_ok(self.m, P), P)

    def test_nested_agrees_with_contiguity(self):
        for P in (6.0, 10.0, 12.0, 17.8, 17.85, 18.0, 20.0, 40.0, 70.0, 3.3, 7.7, 25.0):
            pts = c.phases(self.m, P)
            self.assertEqual(c.nested_ok(pts, P), c.contiguous_y(pts), P)

    def test_intermediate_cases(self):
        good = {"a": (2010.0, "Y"), "b": (2009.0, "I"), "c": (2011.0, "I"), "d": (2005.0, "X"),
                "e": (2015.0, "X")}
        self.assertTrue(c.periodic_ok(good, 100.0))
        bad = {"a": (2010.0, "Y"), "b": (2011.0, "X"), "c": (2012.0, "Y"), "d": (2030.0, "I"),
               "e": (2060.0, "X")}
        self.assertFalse(c.periodic_ok(bad, 100.0))

    def test_chance_fraction(self):
        self.assertAlmostEqual(c.chance_fraction(5, 3), 8.0 / 56.0, places=9)


class TestHelpers(unittest.TestCase):
    def test_intervals(self):
        self.assertEqual(c.intervals([1, 2, 3, 4, 5], [False, True, True, False, True]), [[2, 3], [5, 5]])
        self.assertEqual(c.intervals([1, 2], [False, False]), [])

    def test_continuum(self):
        self.assertEqual(c.continuum_from([1, 2, 3, 4], [True, False, True, True]), 3)
        self.assertIsNone(c.continuum_from([1, 2], [True, False]))
        self.assertEqual(c.continuum_from([1, 2], [True, True]), 1)

    def test_load_new(self):
        d = tempfile.mkdtemp()
        p = os.path.join(d, "m34_candidates.csv")
        with open(p, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=["id", "t", "state", "used_in_mirror"])
            w.writeheader()
            w.writerow(dict(id="a", t="2005.5", state="Y", used_in_mirror="True"))
            w.writerow(dict(id="b", t="2016.0", state="intermediate", used_in_mirror="True"))
            w.writerow(dict(id="c", t="2017.0", state="X", used_in_mirror="False"))
            w.writerow(dict(id="d", t="2018.0", state="new-state", used_in_mirror="True"))
        self.assertEqual(c.load_new(p), {"a": (2005.5, "Y"), "b": (2016.0, "I")})

    def test_grade_untestable(self):
        self.assertEqual(c.grade({}, {}, 0, {})["R4"], "untestable")


if __name__ == "__main__":
    unittest.main()

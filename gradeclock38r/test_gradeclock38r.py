"""Tests for gradeclock38r.py (robustness of pre-registration 38 arm W). Expected values derived by hand.
Part of inner-core-repeaters by MV2, Munim, Inc. (made with AI agents). Licence: Apache-2.0."""
import unittest

import gradeclock38r as r

g = r.load_base()


def row(pair, e1, t1, e2, t2, yka, ilar=-2, mult="X"):
    return dict(pair=pair, YKA=yka, ILAR=ilar, e1=e1, e2=e2, t1=t1, t2=t2, mult=mult)


# multiplet A: backtracking-like (a3 in 2012 matches a1 in 2000, not a2 in 2004)
A = [row("Q1", "a1", 2000.0, "a2", 2004.0, 2, mult="A"),
     row("Q2", "a1", 2000.0, "a3", 2012.0, 0, mult="A"),
     row("Q3", "a2", 2004.0, "a3", 2012.0, 2, mult="A")]
# multiplet X: calendar-like (only the 4-year pair is similar)
X = [row("P1", "e1", 2000.0, "e2", 2004.0, 0),
     row("P2", "e1", 2000.0, "e3", 2012.0, 2),
     row("P3", "e2", 2004.0, "e3", 2012.0, 2)]
TURNS, RATIOS = [2006.0, 2008.0], [1.0, 2.0]


class TestRobust(unittest.TestCase):
    def test_base_version(self):
        self.assertEqual(g.VERSION, "0.1.0")

    def test_eras(self):
        self.assertTrue(r.in_eras(2016.0, g))
        self.assertTrue(r.in_eras(2022.1, g))
        self.assertFalse(r.in_eras(2021.0, g))
        self.assertFalse(r.in_eras(2014.9, g))      # 2014-12-01 = 2014 + 334/365 = 2014.915
        self.assertTrue(r.in_eras(2014.92, g))

    def test_mask(self):
        rows = [row("P1", "e1", 2010.0, "e2", 2016.0, 1, ilar=0), row("P2", "e1", 2010.0, "e3", 2012.0, 0)]
        m = r.mask_yka_eras(rows, g)
        self.assertEqual([x["YKA"] for x in m], [-2, 0])
        self.assertEqual(m[0]["ILAR"], 0)
        self.assertEqual(rows[0]["YKA"], 1)          # original untouched

    def test_drop(self):
        self.assertEqual([x["pair"] for x in r.drop_multiplet(A + X, "A")], ["P1", "P2", "P3"])

    def test_grade_r3(self):
        self.assertEqual(r.grade_r3(0.10), "transfers")
        self.assertEqual(r.grade_r3(0.25), "mixed")
        self.assertEqual(r.grade_r3(0.2501), "does-not-transfer")

    def test_transfer(self):
        # A on clocks [none, (2006,1), (2006,2), (2008,1), (2008,2)] scores [1, 0, .5, .75, 1]: best (2006,1)
        # X scores [0, .75, 1, .75, .5]; at (2006,1) X = 0.75; 2 of 5 clocks are lower -> 0.4
        # Wang's clock puts e3 at 2007.1: d = 4, 7.1, 3.1 -> X = 0.5; 1 of 5 lower -> 0.2
        t = r.transfer(A + X, g, "YKA", turns=TURNS, ratios=RATIOS)
        self.assertEqual((t["n_src"], t["n_rest"]), (2, 2))
        self.assertEqual(t["src_s_min"], 0.0)
        self.assertEqual(t["src_best"]["n"], 1)
        self.assertEqual(t["rest_s_median_at_src_best"], 0.75)
        self.assertEqual(t["rest_s_noturn"], 0.0)
        self.assertEqual(t["rest_pct_median"], 0.4)
        self.assertEqual(t["rest_pct_best"], 0.4)
        self.assertEqual(t["rest_s_wang"], 0.5)
        self.assertEqual(t["rest_wang_pct"], 0.2)
        self.assertEqual(t["grade"], "does-not-transfer")

    def test_transfer_same_pattern(self):
        a_as_x = [dict(x, mult="X", e1=x["e1"].replace("a", "b"), e2=x["e2"].replace("a", "b")) for x in A]
        t = r.transfer(A + a_as_x, g, "YKA", turns=TURNS, ratios=RATIOS)
        self.assertEqual(t["rest_pct_median"], 0.0)
        self.assertEqual(t["grade"], "transfers")

    def test_untestable(self):
        self.assertEqual(r.transfer(X, g, "YKA", turns=TURNS, ratios=RATIOS)["grade"], "untestable")


if __name__ == "__main__":
    unittest.main()

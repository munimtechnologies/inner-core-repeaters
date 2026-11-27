"""Tests for gradez.py (prereg 38 arm Z). Expected values derived by hand before any run.
MV2, Munim, Inc.; made with AI agents. Apache-2.0."""
import os
import sys
import unittest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gradez as z  # noqa: E402


def rec(m, arr, d1, d2, lab, region="all"):
    return dict(multiplet=m, array=arr, date1=d1, date2=d2, region=region, label=lab, source="test")


TRIPLE = [rec("X", "YKAR", "2000-01-01", "2005-01-01", "C"),
          rec("X", "YKA", "2010-01-01", "2000-01-01", "S"),      # reversed dates on purpose
          rec("X", "YKA", "2005/01/01", "2010/01/01", "C")]


class TestInput(unittest.TestCase):
    def test_merge(self):
        self.assertEqual(z.merge_grade(["S", "S"]), 0)
        self.assertEqual(z.merge_grade(["C"]), 2)
        self.assertEqual(z.merge_grade(["S", "C"]), 1)
        self.assertEqual(z.merge_grade(["?", "S"]), 0)
        self.assertIsNone(z.merge_grade(["?"]))
        self.assertIsNone(z.merge_grade([]))

    def test_transcription(self):
        recs = TRIPLE + [rec("X", "ILAR", "2000-01-01", "2005-01-01", "S", "east"),
                         rec("X", "ILAR", "2000-01-01", "2005-01-01", "C", "west"),
                         rec("X", "ILAR", "2000-01-01", "2000-01-01", "C"),
                         rec("X", "PFO", "2000-01-01", "2005-01-01", "C")]
        out = z.read_transcription(recs)
        self.assertEqual(out["YKA"], [("X", "2000-01-01", "2005-01-01", 2),
                                      ("X", "2000-01-01", "2010-01-01", 0),
                                      ("X", "2005-01-01", "2010-01-01", 2)])
        self.assertEqual(out["ILAR"], [("X", "2000-01-01", "2005-01-01", 1)])

    def test_decimal_year(self):
        self.assertEqual(z.decimal_year("2000-01-01"), 2000.0)
        self.assertAlmostEqual(z.decimal_year("2001-07-02"), 2001 + 182 / 365.0)


class TestStatistic(unittest.TestCase):
    def setUp(self):
        self.spec = z.build(z.read_transcription(TRIPLE)["YKA"])

    def test_comparisons(self):
        # pair 1 (2000-2010, grade 0) against pairs 0 and 2 (grade 2)
        self.assertEqual(self.spec["comps"], [(1, 0), (1, 2)])

    def test_scores(self):
        t = self.spec["times"]
        self.assertEqual(z.score(self.spec, t, None), 1.0)          # d = 5, 10, 5
        self.assertEqual(z.score(self.spec, t, (2005.0, 1.0)), 0.0)  # d = 5, 0, 5
        self.assertEqual(z.score(self.spec, t, (2005.0, 2.0)), 0.25)  # d = 5, 2.5, 2.5 (one tie)

    def test_drop(self):
        sd = z.build(z.read_transcription(TRIPLE)["YKA"], drop=("X", "2005-01-01"))
        self.assertEqual(len(sd["pairs"]), 1)
        self.assertEqual(sd["comps"], [])

    def test_grades(self):
        self.assertEqual(z.grade_z1(0.05), "structure")
        self.assertEqual(z.grade_z1(0.2), "weak")
        self.assertEqual(z.grade_z1(0.21), "no-structure")
        self.assertEqual(z.grade_z2(0.01, 0.019), "not-needed")
        self.assertEqual(z.grade_z2(0.05, 0.02), "turn-needed")
        self.assertEqual(z.grade_z2(0.1, 0.1), "weak")
        self.assertEqual(z.grade_z2(0.3, 0.1), "not-needed")
        self.assertEqual(z.grade_z3(0.05), "top5")
        self.assertEqual(z.grade_z3(0.051), "not-top5")
        self.assertEqual(z.grade_z4(0.10), "agrees")
        self.assertEqual(z.grade_z4(0.25), "mixed")
        self.assertEqual(z.grade_z4(0.26), "disagrees")
        self.assertEqual(z.grade_z5(4, 4), "untestable")
        self.assertEqual(z.grade_z5(5, 4), "agree")
        self.assertEqual(z.grade_z5(5, 3), "mixed")
        self.assertEqual(z.grade_z5(6, 3), "disagree")

    def test_headline(self):
        self.assertEqual(z.headline("structure", "turn-needed", "agrees"), "same-clock")
        self.assertEqual(z.headline("structure", "turn-needed", "disagrees"), "other-clock")
        self.assertEqual(z.headline("structure", "turn-needed", "mixed"), "turn-unclear")
        self.assertEqual(z.headline("structure", "not-needed", "agrees"), "calendar-only")
        self.assertEqual(z.headline("no-structure", "not-needed", "agrees"), "no-clock")
        self.assertEqual(z.headline("weak", "weak", "mixed"), "weak")
        self.assertEqual(z.headline("untestable", "", ""), "untestable")


class TestWang(unittest.TestCase):
    def test_parse(self):
        rows = [["Table 2"], ["# pair"],
                ["P1", "x", "0", "2", "101", "102", "2003-03-21T00:00:00", "2008-02-17", "", "", "A"],
                ["note"]]
        self.assertEqual(z.wang_pairs(rows), {("2003-03-21", "2008-02-17"): {"YKA": 0, "ILAR": 2}})

    def test_shared(self):
        zp = [("F", "2005-06-10", "2012-12-03", 2), ("F", "2005-06-10", "2022-03-07", 0),
              ("A", "2001-01-01", "2002-01-01", 1), ("A", "2001-01-01", "2003-01-01", 0)]
        w = {("2005-06-10", "2012-12-03"): {"YKA": 0, "ILAR": 2},
             ("2005-06-10", "2022-03-07"): {"YKA": 0, "ILAR": 0},
             ("2001-01-01", "2002-01-01"): {"YKA": 0, "ILAR": 0}}
        n, ag, rows = z.shared_agreement(zp, w, "YKA")
        self.assertEqual((n, ag, len(rows)), (2, 1, 3))
        n, ag, rows = z.shared_agreement(zp, w, "ILAR")
        self.assertEqual((n, ag, len(rows)), (2, 2, 3))


class TestAnalysis(unittest.TestCase):
    def test_untestable_and_checks(self):
        zp = z.read_transcription(TRIPLE)
        exp = {"n_pairs_YKA": 3, "n_comp_YKA": 2, "n_pairs_ILAR": 0, "n_comp_ILAR": 0}
        res, prof, rows = z.analyse(zp, {}, wang=None, expected=exp, n_shuffles=0)
        self.assertTrue(res["checks_ok"])
        self.assertEqual(res["arrays"]["YKA"]["Z1"], "untestable")
        self.assertEqual(res["Z5"]["YKA"]["grade"], "untestable")
        res2, _, _ = z.analyse(zp, {}, wang=None, expected={"n_comp_YKA": 3}, n_shuffles=0)
        self.assertFalse(res2["checks_ok"])
        res3, _, _ = z.analyse(zp, {}, wang=None, expected={}, n_shuffles=0)
        self.assertFalse(res3["checks_ok"])

    def test_calendar_ordered(self):
        # 9 events 2000..2008 (1 Jan): consecutive pairs similar (8), pairs >= 4 yr apart changed (15)
        recs = []
        for i in range(9):
            for j in range(i + 1, 9):
                if j - i == 1:
                    recs.append(rec("Q", "YKA", "%d-01-01" % (2000 + i), "%d-01-01" % (2000 + j), "S"))
                elif j - i >= 4:
                    recs.append(rec("Q", "YKA", "%d-01-01" % (2000 + i), "%d-01-01" % (2000 + j), "C"))
        spec = z.build(z.read_transcription(recs)["YKA"])
        self.assertEqual(len(spec["pairs"]), 23)
        R = z.analyse_array(spec, n_shuffles=3, seed=1, fine=False)
        self.assertEqual(R["n_comp"], 120)
        self.assertEqual(R["s_noturn"], 0.0)
        self.assertEqual(R["s_star_test"], 0.0)
        self.assertEqual(R["improvement"], 0.0)
        self.assertEqual(R["Z2"], "not-needed")
        self.assertIn(R["Z1"], ("structure", "weak", "no-structure"))
        armw = [[2009.2, 2.2705], "noturn"]
        pct = z.best_percentile(armw, spec, spec["times"], R["_fine"])
        self.assertEqual(pct, 0.0)   # the no-turn clock scores 0, nothing is strictly below it


if __name__ == "__main__":
    unittest.main()

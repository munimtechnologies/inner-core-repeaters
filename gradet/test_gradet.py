# Copyright 2026 MV2, Munim, Inc.
# Licensed under the Apache License, Version 2.0 (Apache-2.0). See LICENSE.
# Made with AI agents working in MV2's name.
"""Hand-derived tests for gradet.py 0.1.0 (standard library only)."""

import math
import os
import sys
import unittest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gradet as g  # noqa: E402


def row(mult="B", array="YKA", ev1="x1", ev2="x2", t1=2005.0, t2=2006.0, D=0.01, ok="True"):
    return dict(array=array, mult=mult, ev1=ev1, ev2=ev2, t1=str(t1), t2=str(t2), D=str(D), ok=ok)


class TestGradeT(unittest.TestCase):
    def test_version(self):
        self.assertEqual(g.VERSION, "0.1.0")

    def test_decimal_year(self):
        # m2measure 0.1.0 convention: 79 days + 10:47:47.62 over 365 days
        self.assertAlmostEqual(g.decimal_year("2003-03-21T10:47:47.620Z"), 2003.2176708, places=6)
        # 2022-03-07: day 66 -> 65/365
        self.assertAlmostEqual(g.decimal_year("2022-03-07T00:00:00Z"), 2022 + 65 / 365.0, places=9)

    def test_eras(self):
        self.assertTrue(g.in_era(g.decimal_year("2022-03-07T00:00:00Z")))
        self.assertTrue(g.in_era(2017.0))
        self.assertFalse(g.in_era(g.decimal_year("2020-07-11T00:00:00Z")))   # 2020.5246 > 2020.268
        self.assertFalse(g.in_era(g.decimal_year("2014-11-30T00:00:00Z")))
        self.assertFalse(g.in_era(2021.5))

    def test_nearest_rank(self):
        self.assertEqual(g.nearest_rank(list(range(1, 21)), 0.95), 19)
        self.assertEqual(g.nearest_rank(list(range(1, 11)), 0.50), 5)
        self.assertEqual(g.nearest_rank(list(range(10, 0, -1)), 0.95), 10)
        self.assertTrue(math.isnan(g.nearest_rank([], 0.5)))

    def test_reference_filter(self):
        rows = [row(),                                   # kept
                row(mult="A"),                           # multiplet A excluded
                row(array="ILAR"),                       # ILAR excluded
                row(ok="False"),                         # not ok
                row(t1=2005.0, t2=2007.0),               # dt = 2.0 -> excluded (< 2 required)
                row(t1=2016.0, t2=2017.0),               # era
                row(t1=2021.5, t2=2022.0),               # second era (2022.0 inside)
                row(D="nan"),                            # not finite
                row(ev1="y1", ev2="y2", t1=2010.0, t2=2011.9, D=0.03)]   # kept
        ref = g.reference(rows)
        self.assertEqual([r["D"] for r in ref], [0.01, 0.03])

    def test_grade_beyond(self):
        ref = [0.001 * i for i in range(1, 21)]          # q50 = 0.010, q95 = 0.019
        res = g.grade(ref)
        self.assertAlmostEqual(res["q50"], 0.010, places=12)
        self.assertAlmostEqual(res["q95"], 0.019, places=12)
        v = {p["pair"]: p["verdict"] for p in res["pairs"]}
        self.assertEqual(v["2008-02-17/2011-02-20"], "within")    # 0.017939 <= 0.019
        self.assertEqual(v["2003-03-21/2020-07-11"], "above")     # 0.022813
        self.assertEqual(v["2020-07-11/2022-03-07"], "above")     # 0.020867
        self.assertEqual(res["n_above"], 2)
        self.assertEqual(res["headline"], "beyond-scatter")
        p = [p for p in res["pairs"] if p["pair"] == "2008-02-17/2011-02-20"][0]
        self.assertAlmostEqual(p["frac_ref_ge"], 3 / 20.0)        # 0.018, 0.019, 0.020
        era = {p["pair"]: p["era_pair"] for p in res["pairs"]}
        self.assertEqual(era, {"2008-02-17/2011-02-20": False, "2003-03-21/2020-07-11": False,
                               "2020-07-11/2022-03-07": True})
        self.assertEqual(res["report_only"]["frac_ref_above_cut"]["armM_cut"], 0.0)
        self.assertEqual(res["report_only"]["key_pair_2020_2022"], ["above"])

    def test_grade_within(self):
        res = g.grade([0.01 * i for i in range(1, 11)])  # q95 = 0.10
        self.assertEqual(res["n_above"], 0)
        self.assertEqual(res["headline"], "within-scatter")

    def test_grade_untestable(self):
        res = g.grade([0.01] * 7)
        self.assertEqual(res["headline"], "untestable")

    def test_check_published(self):
        rows = [dict(ev1="usp000hute", ev2="usp000fz4d", D_c310="0.017939"),
                dict(ev1="usp000bt5x", ev2="us7000am0e", D_c310="0.02281300"),
                dict(ev1="us7000am0e", ev2="us6000h2wz", D_c310="0.020867"),
                dict(ev1="usp000bt5x", ev2="us6000h2wz", D_c310="0.014711"),
                dict(ev1="usp000bt5x", ev2="us6000mk65", D_c310="0.013100")]   # off by 1e-4
        out = g.check_published(rows)
        self.assertTrue(out["2008-02-17/2011-02-20"]["ok"])
        self.assertTrue(out["2020-07-11/2022-03-07"]["ok"])
        self.assertFalse(out["2003-03-21/2024-03-21"]["ok"])


if __name__ == "__main__":
    unittest.main()

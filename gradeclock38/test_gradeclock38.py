"""Tests for gradeclock38.py (pre-registration 38, arm W). Expected values derived by hand.
Part of inner-core-repeaters by MV2, Munim, Inc. (made with AI agents). Licence: Apache-2.0."""
import os
import random
import tempfile
import unittest

import gradeclock38 as g


def row(pair, e1, t1, e2, t2, yka, ilar=-2, mult="X"):
    return dict(pair=pair, YKA=yka, ILAR=ilar, e1=e1, e2=e2, t1=t1, t2=t2, mult=mult)


# toy backtracking multiplet X: e3 (2012) matches e1 (2000), not e2 (2004)
TOY = [row("P1", "e1", 2000.0, "e2", 2004.0, 2),
       row("P2", "e1", 2000.0, "e3", 2012.0, 0),
       row("P3", "e2", 2004.0, "e3", 2012.0, 2)]
EXTRA = [row("P4", "f1", 2001.0, "f2", 2002.0, 1, mult="Y"),
         row("P5", "h1", 2000.0, "h2", 2020.0, 0, mult="-"),
         row("P6", "i1", 2000.0, "i2", 2001.0, 0, mult="I")]

CSV_TEXT = (
    "Supplementary Data Table 2  List of 142 event pairs,,,,,,,,,,,\n"
    ",,,,,,,,,,,\n"
    "#   ,Q ,YM,IM ,Events #1 ,Events #2 ,  date 1 ,        date 2 ,              lat ,          lon,    Multiplet,reference\n"
    "P21,1,  1,-2,712,718,1999-07-01T00:00:00,2003-03-21T00:00:00,-59.05,-25.718,A,PK\n"
    "P84,1,  2,-2,801,727,1995-08-08T00:00:00,2005-03-24T00:00:00,-55.985,-27.078,P,S2 \n"
    ",,,,,,,,,,,\n"
    "Key to columns,,,,,,,,,,,\n"
    ",Not all numbers are present,,,,,,,,,,\n")


class TestTable(unittest.TestCase):
    def test_decimal_year(self):
        self.assertAlmostEqual(g.decimal_year("1999-07-01T00:00:00"), 1999 + 181 / 365.0)
        self.assertAlmostEqual(g.decimal_year("2008-02-17T00:00:00"), 2008 + 47 / 366.0)

    def test_parse_csv(self):
        fd, path = tempfile.mkstemp(suffix=".csv")
        with os.fdopen(fd, "w", newline="") as f:
            f.write(CSV_TEXT)
        try:
            rows = g.parse_table2_csv(path)
        finally:
            os.remove(path)
        self.assertEqual([r["pair"] for r in rows], ["P21", "P84"])
        self.assertEqual((rows[0]["YKA"], rows[0]["ILAR"], rows[0]["e1"], rows[0]["e2"]), (1, -2, "712", "718"))
        self.assertAlmostEqual(rows[0]["t2"], 2003 + 79 / 365.0)
        self.assertEqual(rows[1]["mult"], "P")
        self.assertEqual(rows[1]["YKA"], 2)

    def test_event_id_normalised(self):
        self.assertEqual(g._event_id(" 701.0"), "701")
        self.assertEqual(g._event_id(701), "701")


class TestScore(unittest.TestCase):
    def setUp(self):
        self.times = g.event_times(TOY + EXTRA)
        self.spec = g.build(TOY, "YKA")

    def test_equiv(self):
        self.assertEqual(g.equiv(2004.0, 2006.0, 2.0), 2004.0)
        self.assertEqual(g.equiv(2012.0, 2006.0, 2.0), 2003.0)
        self.assertEqual(g.position(2012.0, None), 2012.0)

    def test_comparisons(self):
        self.assertEqual(sorted(self.spec["comps"]), [(1, 0), (1, 2)])
        self.assertEqual(len(g.build(TOY, "ILAR")["comps"]), 0)

    def test_noturn(self):          # d = 4, 12, 8: the similar pair is the farthest
        self.assertEqual(g.discordance(self.spec, self.times, None), 2.0)
        self.assertEqual(g.score(self.spec, self.times, None), 1.0)

    def test_symmetric_turn(self):  # e3 -> 2000: d = 4, 0, 4
        self.assertEqual(g.score(self.spec, self.times, (2006.0, 1.0)), 0.0)

    def test_slow_return(self):     # e3 -> 2003: d = 4, 3, 1
        self.assertEqual(g.score(self.spec, self.times, (2006.0, 2.0)), 0.5)

    def test_tie(self):             # e3 -> 2002: d = 4, 2, 2 (tie counts one half)
        self.assertEqual(g.score(self.spec, self.times, (2006.0, 1.5)), 0.25)

    def test_late_turn(self):       # e3 -> 2006: d = 4, 6, 2
        self.assertEqual(g.score(self.spec, self.times, (2008.0, 2.0)), 1.0)

    def test_groups_and_exclusions(self):
        rows = TOY + EXTRA
        within = g.build(rows, "YKA")
        self.assertEqual(len(within["comps"]), 2)           # Y has one pair, '-' and I excluded
        pooled = g.build(rows, "YKA", pooled=True)
        self.assertEqual(len(pooled["comps"]), 8)           # P2,P5 vs P1,P3,P4; P4 vs P1,P3
        self.assertEqual(g.discordance(pooled, self.times, None), 6.0)
        grp = g.discordance_groups(within, self.times, None)
        self.assertEqual(grp, {"X": [2.0, 2]})

    def test_scan_order(self):
        clocks = g.grid_clocks([2006.0, 2008.0], [1.0, 2.0])
        self.assertEqual(clocks, [None, (2006.0, 1.0), (2006.0, 2.0), (2008.0, 1.0), (2008.0, 2.0)])
        # (2008, 1): e3 -> 2004: d = 4, 4, 0 -> P2 vs P1 tie 0.5, P2 vs P3 1 -> 0.75
        self.assertEqual(g.scan(self.spec, self.times, clocks), [1.0, 0.0, 0.5, 0.75, 1.0])
        self.assertEqual(g.best_set(clocks, g.scan(self.spec, self.times, clocks)), [(2006.0, 1.0)])

    def test_turn_profile(self):
        clocks = g.grid_clocks([2006.0, 2008.0], [1.0, 2.0])
        prof = g.turn_profile(clocks, g.scan(self.spec, self.times, clocks))
        self.assertEqual(prof, [dict(turn=2006.0, s_min=0.0, r_lo=1.0, r_hi=1.0),
                                dict(turn=2008.0, s_min=0.75, r_lo=1.0, r_hi=1.0)])

    def test_shuffle_keeps_dates(self):
        rows = TOY + EXTRA
        times = g.event_times(rows)
        mults = g.multiplet_events(rows)
        self.assertEqual(sorted(mults), ["X", "Y"])
        sh = g.shuffle_times(times, mults, random.Random(1))
        self.assertEqual(sorted(sh[e] for e in ("e1", "e2", "e3")), [2000.0, 2004.0, 2012.0])
        self.assertEqual(sorted(sh[e] for e in ("f1", "f2")), [2001.0, 2002.0])
        self.assertEqual((sh["h1"], sh["h2"], sh["i1"]), (2000.0, 2020.0, 2000.0))


class TestStats(unittest.TestCase):
    def test_p(self):
        nulls = [0.05, 0.1, 0.2, 0.3]
        self.assertAlmostEqual(g.p_le(0.1, nulls), 0.6)
        self.assertAlmostEqual(g.p_ge(0.2, nulls), 0.6)
        self.assertAlmostEqual(g.p_le(0.0, nulls), 0.2)

    def test_frac_below(self):
        self.assertEqual(g.frac_below(0.3, [0.1, 0.3, 0.5, 0.2]), 0.5)

    def test_cross_percentile(self):
        clocks = [None, (2006.0, 1.0), (2008.0, 2.0)]
        self.assertAlmostEqual(g.cross_percentile([(2008.0, 2.0)], clocks, [0.2, 0.5, 0.1]), 0.0)
        self.assertAlmostEqual(g.cross_percentile([None, (2006.0, 1.0)], clocks, [0.2, 0.5, 0.1]), 1 / 3.0)

    def test_grades(self):
        self.assertEqual(g.grade_w1(0.05), "structure")
        self.assertEqual(g.grade_w1(0.2), "weak")
        self.assertEqual(g.grade_w1(0.2001), "no-structure")
        self.assertEqual(g.grade_w2(0.01, 0.019), "not-needed")
        self.assertEqual(g.grade_w2(0.05, 0.02), "turn-needed")
        self.assertEqual(g.grade_w2(0.1, 0.05), "weak")
        self.assertEqual(g.grade_w2(0.25, 0.3), "not-needed")
        self.assertEqual(g.grade_w3(0.05), "top5")
        self.assertEqual(g.grade_w3(0.0501), "not-top5")
        self.assertEqual(g.grade_w4(0.1, 0.05), "consistent")
        self.assertEqual(g.grade_w4(0.3, 0.26), "inconsistent")
        self.assertEqual(g.grade_w4(0.05, 0.3), "mixed")

    def test_grids(self):
        self.assertEqual(len(g.TEST_RATIOS), 28)
        self.assertEqual((g.TEST_TURNS[0], g.TEST_TURNS[-1], len(g.TEST_TURNS)), (2000.0, 2020.0, 81))
        self.assertEqual((g.FINE_TURNS[0], g.FINE_TURNS[-1], len(g.FINE_TURNS)), (2000.0, 2020.0, 401))
        self.assertIn(2008.5, g.FINE_TURNS)


if __name__ == "__main__":
    unittest.main()

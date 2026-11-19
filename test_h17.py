#!/usr/bin/env python3
"""Offline tests for h17.py (inner-core-repeaters, MV2, Munim, Inc.; made with AI agents). Apache-2.0."""
import math
import random
import unittest

import h17


def _row(mult, a, b, t1, t2, D, g1="new", g2="new", D310=None):
    return dict(mult=mult, ev1=a, ev2=b, key=tuple(sorted((a, b))), t1=t1, t2=t2, D=D,
                D_c13=D, D_c310=D if D310 is None else D310, dB=h17.d_B(t1, t2), dM=h17.d_M(t1, t2),
                cross=h17.is_cross(t1, t2), m0324=False, gen1=g1, gen2=g2, gclass=h17.gen_class(g1, g2))


class TestCopied(unittest.TestCase):
    def test_rank_ties(self):
        r = list(h17.rankdata([3.0, 1.0, 3.0, 2.0]))
        self.assertEqual(r[1], 0.0)
        self.assertEqual(r[3], 1.0)
        self.assertAlmostEqual(r[0], 2.5)
        self.assertAlmostEqual(r[2], 2.5)

    def test_spearman(self):
        self.assertAlmostEqual(h17.spearman([1, 2, 3, 4, 5], [2, 4, 6, 8, 10]), 1.0, places=9)
        self.assertAlmostEqual(h17.spearman([1, 2, 3, 4, 5], [5, 4, 3, 2, 1]), -1.0, places=9)
        self.assertTrue(math.isnan(h17.spearman([1, 2, 3], [1, 2, 3])))
        self.assertTrue(math.isnan(h17.spearman([1, 2, float("nan"), 4, 5], [1, 2, 3, 4, 5])))

    def test_clock_times(self):
        self.assertAlmostEqual(h17.equiv(2018.5), 2004.5, places=9)
        self.assertAlmostEqual(h17.d_B(2000.0, 2018.5), 4.5, places=9)
        self.assertTrue(h17.is_cross(2005.0, 2010.0))
        self.assertFalse(h17.is_cross(2009.0, 2010.0))

    def test_decide13(self):
        self.assertEqual(h17.decide13(0.3, 0.1, 0.95, 0.05, 0.4, 0.2), "B")
        self.assertEqual(h17.decide13(0.1, 0.3, 0.05, 0.95, 0.2, 0.4), "M")
        self.assertEqual(h17.decide13(0.02, 0.03, 0.5, 0.5, 0.1, 0.15), "none")
        self.assertEqual(h17.decide13(0.12, 0.15, 0.3, 0.7, 0.2, 0.25), "inconclusive")
        self.assertEqual(h17.decide13(float("nan"), 0.1, 0.5, 0.5, 0.1, 0.1), "untestable")


class TestNew(unittest.TestCase):
    def test_gen_class(self):
        self.assertEqual(h17.gen_class("legacy", "legacy"), "same-legacy")
        self.assertEqual(h17.gen_class("new", "new"), "same-new")
        self.assertEqual(h17.gen_class("legacy", "new"), "cross")
        self.assertEqual(h17.gen_class("", "new"), "unknown")

    def test_read_yka(self):
        text = ("mult,ev1,ev2,t1,t2,gen1,gen2,ok,D,D_c13,D_c310\n"
                "Y1,a,b,2001.0,2002.0,legacy,legacy,True,0.1,0.05,0.2\n"
                "Y1,a,c,2001.0,2016.0,legacy,new,True,0.3,0.1,0.4\n"
                "Y1,b,c,2002.0,2016.0,legacy,new,False,0.3,0.1,0.4\n")
        rows = h17.read_yka(text)
        self.assertEqual(len(rows), 2)
        self.assertEqual(rows[0]["gclass"], "same-legacy")
        self.assertEqual(rows[1]["gclass"], "cross")
        self.assertTrue(rows[1]["cross"])
        self.assertAlmostEqual(rows[1]["dM"], 15.0, places=9)

    def test_read_scores_and_verified(self):
        text = ("array,mult,ev1,ev2,S,n_st\n"
                "yka,Y1,b,a,0.2,3\n"
                "yka,Y1,a,c,0.5,2\n"
                "yka,Y1,c,d,0.1,1\n"
                "yka,Y1,d,e,,0\n")
        sc = h17.read_scores(text)
        self.assertEqual(set(sc), {("a", "b"), ("a", "c"), ("c", "d")})
        rows = [_row("Y1", "a", "b", 2001, 2002, 0.1), _row("Y1", "a", "c", 2001, 2003, 0.2),
                _row("Y1", "c", "d", 2003, 2004, 0.3), _row("Y1", "d", "e", 2004, 2005, 0.3)]
        scored, V = h17.verified(rows, sc)
        self.assertEqual(len(scored), 2)
        self.assertEqual([r["key"] for r in V], [("a", "b")])

    def test_decide17(self):
        self.assertEqual(h17.decide17(0.2, 0.05, 0.3, 100, 10), "time-signal")
        self.assertEqual(h17.decide17(0.05, -0.05, 0.15, 100, 10), "no-time-signal")
        self.assertEqual(h17.decide17(0.05, -0.05, 0.25, 100, 10), "inconclusive")
        self.assertEqual(h17.decide17(0.2, -0.01, 0.3, 100, 10), "inconclusive")
        self.assertEqual(h17.decide17(0.2, 0.05, 0.3, 39, 10), "untestable")
        self.assertEqual(h17.decide17(0.2, 0.05, 0.3, 100, 3), "untestable")

    def _synthetic(self, hardware_only):
        rng = random.Random(5)
        rows = []
        for f in range(8):
            times = [2001 + rng.random() * 11 for _ in range(4)] + [2015.5 + rng.random() * 10 for _ in range(4)]
            gens = ["legacy"] * 4 + ["new"] * 4
            for i in range(8):
                for j in range(i + 1, 8):
                    g = h17.gen_class(gens[i], gens[j])
                    dm = abs(times[j] - times[i])
                    if hardware_only:   # no time trend inside a generation (slightly falling, by design)
                        D = 0.25 + 0.01 * rng.random() if g == "cross" else 0.05 - 0.001 * dm + 1e-5 * rng.random()
                    else:
                        D = 0.02 + 0.02 * dm + 0.01 * rng.random()
                    rows.append(_row("F%d" % f, "e%d_%d" % (f, i), "e%d_%d" % (f, j), times[i], times[j], D,
                                     gens[i], gens[j]))
        return rows

    def test_time_signal_detects_trend(self):
        V = [r for r in self._synthetic(False) if r["gclass"].startswith("same-")]
        res = h17.time_signal(V, "D", nboot=200)
        self.assertEqual(res["decision"], "time-signal")

    def test_hardware_only_has_no_same_gen_signal(self):
        rows = self._synthetic(True)
        V = [r for r in rows if r["gclass"].startswith("same-")]
        res = h17.time_signal(V, "D", nboot=200)
        self.assertEqual(res["decision"], "no-time-signal")
        allres = h17.time_signal(rows, "D", nboot=200)
        self.assertEqual(allres["decision"], "time-signal")

    def test_bins_table(self):
        rows = [_row("F", "a", "b", 2001, 2002, 0.1, "legacy", "legacy"),
                _row("F", "a", "c", 2001, 2016, 0.5, "legacy", "new"),
                _row("F", "b", "c", 2002, 2016, 0.3, "legacy", "new")]
        tb = h17.bins_table(rows)
        self.assertEqual(tb[0]["n_same"], 1)
        self.assertEqual(tb[3]["n_cross_gen"], 1)
        self.assertEqual(tb[4]["n_cross_gen"], 1)
        self.assertAlmostEqual(tb[0]["median_same"], 0.1)
        self.assertTrue(math.isnan(tb[1]["median_same"]))

    def test_locks(self):
        fams = [dict(station="IU.QSPA", mult="I:M000", status="ok", onset_minus_pred_s="4.25",
                     ref_onset_minus_pred_s="3.25"),
                dict(station="IU.QSPA", mult="I:M001", status="ok", onset_minus_pred_s="2.55",
                     ref_onset_minus_pred_s="2.35"),
                dict(station="IU.QSPA", mult="I:M013", status="too_few", onset_minus_pred_s="",
                     ref_onset_minus_pred_s="")]
        self.assertEqual(h17.locked_family_stations(fams), {("IU.QSPA", "I:M000")})

    def test_source_scores(self):
        pr = [dict(station="X", mult="F", ev1="a%d" % i, ev2="b%d" % i, D_late=str(i / 10.0)) for i in range(20)]
        pr += [dict(station="Y", mult="F", ev1="b%d" % i, ev2="a%d" % i, D_late=str(i / 10.0)) for i in range(20)]
        pr += [dict(station="Z", mult="F", ev1="a0", ev2="b0", D_late="0.1")]
        sc, used = h17.source_scores(pr)
        self.assertEqual(used, ["X", "Y"])
        self.assertEqual(sc[("a0", "b0")]["n_st"], 2)
        self.assertAlmostEqual(sc[("a0", "b0")]["S"], 1.0 / 20.0)
        self.assertAlmostEqual(sc[("a19", "b19")]["S"], 1.0)

    def test_analyse_runs(self):
        rows = self._synthetic(False)
        scores = {r["key"]: dict(S=0.1, n_st=3) for r in rows}
        res, V = h17.analyse(rows, scores, nboot=50)
        self.assertEqual(len(V), len(rows))
        self.assertIn(res["headline"], ("time-signal", "no-time-signal", "inconclusive", "untestable"))
        self.assertFalse(res["reproduction"]["reproduced"])
        self.assertIn("same-legacy", res["verified_by_generation"])


if __name__ == "__main__":
    unittest.main()

#!/usr/bin/env python3
"""Offline tests for hw18.py (inner-core-repeaters, MV2, Munim, Inc.; made with AI agents). Apache-2.0.

No network, no obspy: only the pure functions are tested here.
"""
import math
import random
import unittest

import hw18


def _row(mult, a, b, t1, t2, D, g1="new", g2="new"):
    return dict(mult=mult, ev1=a, ev2=b, key=tuple(sorted((a, b))), t1=t1, t2=t2, D=D, D_c13=D, D_c310=D,
                D_full=D, dB=hw18.d_B(t1, t2), dM=hw18.d_M(t1, t2), cross=hw18.is_cross(t1, t2), m0324=False,
                gen1=g1, gen2=g2, gclass=hw18.gen_class(g1, g2), pol_flip=False)


def _synthetic(kind, h=0.2):
    """8 families, 4 legacy (2001-2012) + 4 new (2015.5-2025.5) events each.

    kind 'hardware': no time trend inside a generation (slightly falling, by design, so the
    sign cannot flip); cross-generation pairs carry an extra hardware decorrelation h.
    kind 'trend': D grows with time apart for every pair; no hardware term.
    """
    rng = random.Random(5)
    rows = []
    for f in range(8):
        times = [2001 + rng.random() * 11 for _ in range(4)] + [2015.5 + rng.random() * 10 for _ in range(4)]
        gens = ["legacy"] * 4 + ["new"] * 4
        for i in range(8):
            for j in range(i + 1, 8):
                g = hw18.gen_class(gens[i], gens[j])
                dm = abs(times[j] - times[i])
                if kind == "hardware":
                    base = 0.05 - 0.001 * dm + 1e-5 * rng.random()
                    D = 1.0 - (1.0 - base) * (1.0 - h) if g == "cross" else base
                else:
                    D = 0.02 + 0.02 * dm + 0.01 * rng.random()
                rows.append(_row("F%d" % f, "e%d_%d" % (f, i), "e%d_%d" % (f, j), times[i], times[j], D,
                                 gens[i], gens[j]))
    return rows


class TestPure(unittest.TestCase):
    def test_choose_providers(self):
        sets = {"nrcan": {"legacy": {"B1", "B2", "R1"}, "new": set()},
                "earthscope": {"legacy": set(), "new": {"B1", "B2", "R1", "R2"}}}
        self.assertEqual(hw18.choose_providers(sets), {"legacy": "nrcan", "new": "earthscope"})
        tie = {"nrcan": {"legacy": {"B1"}, "new": {"B1"}}, "earthscope": {"legacy": {"B2"}, "new": {"B2"}}}
        self.assertEqual(hw18.choose_providers(tie), {"legacy": "earthscope", "new": "earthscope"})
        none = {"nrcan": {"legacy": {"B1"}, "new": set()}}
        self.assertEqual(hw18.choose_providers(none), {"legacy": "nrcan", "new": None})

    def test_correct_D(self):
        self.assertAlmostEqual(hw18.correct_D(0.3, 0.1), 1.0 - 0.7 / 0.9, places=12)
        self.assertAlmostEqual(hw18.correct_D(0.3, 0.1, "add"), 0.2, places=12)
        self.assertAlmostEqual(hw18.correct_D(0.3, 0.0), 0.3, places=12)
        self.assertTrue(math.isnan(hw18.correct_D(0.3, 1.0)))
        self.assertTrue(math.isnan(hw18.correct_D(float("nan"), 0.1)))
        # exact inverse of the product rule
        base, h = 0.07, 0.15
        self.assertAlmostEqual(hw18.correct_D(1.0 - (1.0 - base) * (1.0 - h), h), base, places=12)

    def test_corrected_rows_only_cross(self):
        rows = [_row("F", "a", "b", 2001, 2002, 0.1, "legacy", "legacy"),
                _row("F", "a", "c", 2001, 2016, 0.5, "legacy", "new")]
        out = hw18.corrected_rows(rows, 0.2)
        self.assertAlmostEqual(out[0]["D"], 0.1, places=12)
        self.assertAlmostEqual(out[1]["D"], 1.0 - 0.5 / 0.8, places=12)
        self.assertAlmostEqual(rows[1]["D"], 0.5, places=12)  # input untouched

    def test_hw_stats(self):
        a = hw18.hw_stats([0.01, 0.02, 0.03, 0.04, 0.05], nboot=300)
        b = hw18.hw_stats([0.01, 0.02, 0.03, 0.04, 0.05], nboot=300)
        self.assertEqual(a, b)
        self.assertEqual(a["n"], 5)
        self.assertAlmostEqual(a["median"], 0.03, places=12)
        self.assertLessEqual(a["p05"], a["median"])
        self.assertGreaterEqual(a["p95"], a["median"])
        self.assertEqual(hw18.hw_stats([float("nan"), 0.2])["n"], 1)
        self.assertEqual(hw18.hw_stats([])["n"], 0)

    def test_decide18(self):
        self.assertEqual(hw18.decide18(False, 10, 0.2, 0.05, 0.2, 0.3), "not-reproduced")
        self.assertEqual(hw18.decide18(True, 3, 0.2, 0.05, 0.2, 0.3), "untestable")
        self.assertEqual(hw18.decide18(True, 10, 0.2, 0.05, 0.2, 0.3), "trend-survives-hardware")
        self.assertEqual(hw18.decide18(True, 10, 0.02, -0.05, 0.05, 0.15), "hardware-explains")
        self.assertEqual(hw18.decide18(True, 10, 0.12, -0.01, 0.13, 0.25), "inconclusive")
        self.assertEqual(hw18.decide18(True, 10, 0.05, -0.05, 0.08, 0.25), "inconclusive")
        self.assertEqual(hw18.decide18(True, 10, float("nan"), 0.0, 0.1, 0.1), "untestable")

    def test_pair_ok(self):
        self.assertTrue(hw18.pair_ok(0.95, 0.05, (6.0, 7.0), (3.0, 2.5)))
        self.assertFalse(hw18.pair_ok(0.85, 0.05, (6.0, 7.0), (3.0, 2.5)))
        self.assertFalse(hw18.pair_ok(0.95, 0.05, (4.0, 7.0), (3.0, 2.5)))
        self.assertFalse(hw18.pair_ok(0.95, 0.05, (6.0, 7.0), (1.5, 2.5)))
        self.assertFalse(hw18.pair_ok(0.95, float("nan"), (6.0, 7.0), (3.0, 2.5)))

    def test_parse_catalog(self):
        text = ("time,latitude,longitude,depth,mag,magType,id,type\n"
                "2013-07-15T14:03:39.880Z,-60.857,-25.07,11,7.3,mww,usb000ief9,earthquake\n"
                "2013-07-03T20:16:39.880Z,-58.376,-26.33,146.8,5.3,mb,usb000i6aa,earthquake\n"
                "2013-07-04T00:00:00.000Z,-58.0,-26.0,,4.6,mb,xx1,earthquake\n"
                "2013-07-05T00:00:00.000Z,-58.0,-26.0,10,4.6,mb,xx2,explosion\n")
        evs = hw18.parse_catalog(text, lambda s: 2013.0 + int(s[8:10]) / 365.0)
        self.assertEqual([e["id"] for e in evs], ["usb000i6aa", "xx1", "usb000ief9"])
        self.assertAlmostEqual(evs[1]["depth"], 10.0)
        self.assertEqual(hw18.parse_catalog("", lambda s: 0.0), [])

    def test_read_yka_polarity(self):
        text = ("mult,ev1,ev2,t1,t2,gen1,gen2,polarity_flipped,ok,D,D_c13,D_c310,D_full\n"
                "Y1,a,b,2001.0,2002.0,legacy,legacy,False,True,0.1,0.05,0.2,0.1\n"
                "Y1,a,c,2001.0,2016.0,legacy,new,True,True,0.3,0.1,0.4,0.3\n")
        rows = hw18.read_yka(text)
        self.assertEqual([r["pol_flip"] for r in rows], [False, True])
        self.assertEqual(rows[1]["gclass"], "cross")
        self.assertAlmostEqual(rows[1]["D_full"], 0.3)

    def test_bins_excess(self):
        rows = [_row("F", "a", "b", 2001, 2002, 0.1, "legacy", "legacy"),
                _row("F", "c", "d", 2001, 2002.5, 0.3, "legacy", "new")]
        tb = hw18.bins_table(rows)
        self.assertAlmostEqual(tb[0]["excess"], 0.2, places=12)
        self.assertTrue(math.isnan(tb[1]["excess"]))


class TestGrade(unittest.TestCase):
    def test_hardware_only_is_explained(self):
        V = _synthetic("hardware", h=0.2)
        raw = hw18.grade(V, hw18.hw_stats([0.0] * 6, nboot=50), True, nboot=200)
        self.assertEqual(raw["decision"], "trend-survives-hardware")  # uncorrected: looks like a trend
        res = hw18.grade(V, hw18.hw_stats([0.2] * 6, nboot=50), True, nboot=200)
        self.assertEqual(res["decision"], "hardware-explains")
        self.assertLess(res["central"]["rho_M"], 0.0)

    def test_real_trend_survives(self):
        V = _synthetic("trend")
        res = hw18.grade(V, hw18.hw_stats([0.005, 0.01, 0.015, 0.01, 0.012], nboot=100), True, nboot=200)
        self.assertEqual(res["decision"], "trend-survives-hardware")
        self.assertGreater(res["conservative"]["rho_M"], 0.5)

    def test_too_few_events(self):
        V = _synthetic("trend")
        res = hw18.grade(V, hw18.hw_stats([0.01, 0.02, 0.03], nboot=50), True, nboot=50)
        self.assertEqual(res["decision"], "untestable")
        res = hw18.grade(V, hw18.hw_stats([0.01] * 5, nboot=50), False, nboot=50)
        self.assertEqual(res["decision"], "not-reproduced")


if __name__ == "__main__":
    unittest.main()

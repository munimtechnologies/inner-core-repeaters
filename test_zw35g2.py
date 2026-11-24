# Copyright 2026 MV2, Munim, Inc.
# Licensed under the Apache License, Version 2.0 (Apache-2.0). See LICENSE.
# Made with AI agents working in MV2's name.
"""Tests for zw35g2.py (standard library only; expected values worked out by hand)."""
import unittest

import zw35g2 as z

E = list(z.OLD8)
A, B, C = z.F2005, z.F2012, z.F2022


def ref_D(n_gated=5):
    """Pairs among the first n_gated old events get D = 0.01, 0.02, ... in pair order."""
    D, k = {}, 1
    for i in range(n_gated):
        for j in range(i + 1, n_gated):
            D[z.key(E[i], E[j])] = 0.01 * k
            k += 1
    gated = {e: (i < n_gated) for i, e in enumerate(E)}
    return D, gated


class TestScore(unittest.TestCase):
    def test_q_below_hit_level(self):
        D, g = ref_D(5)                       # 10 reference pairs: 0.01 .. 0.10
        g.update({A: True, C: True})
        D[z.key(A, C)] = 0.095                # 9 of 10 below -> q 0.9 < 0.93
        r = z.station_pair_score(D, g, E, A, C)
        self.assertTrue(r["counted"])
        self.assertEqual(r["n_ref"], 10)
        self.assertAlmostEqual(r["q"], 0.9)
        self.assertFalse(r["hit"])
        self.assertAlmostEqual(r["ref_median"], 0.055)
        self.assertAlmostEqual(r["ref_max"], 0.10)

    def test_q_above_all_is_hit(self):
        D, g = ref_D(5)
        g.update({A: True, C: True})
        D[z.key(A, C)] = 0.2
        r = z.station_pair_score(D, g, E, A, C)
        self.assertAlmostEqual(r["q"], 1.0)
        self.assertTrue(r["hit"])

    def test_small_pair_D(self):
        D, g = ref_D(5)
        g.update({A: True, C: True})
        D[z.key(A, C)] = 0.005
        r = z.station_pair_score(D, g, E, A, C)
        self.assertAlmostEqual(r["q"], 0.0)
        self.assertFalse(r["hit"])

    def test_too_few_reference_pairs(self):
        D, g = ref_D(4)                       # 6 pairs < 10
        g.update({A: True, C: True})
        D[z.key(A, C)] = 0.5
        r = z.station_pair_score(D, g, E, A, C)
        self.assertFalse(r["counted"])
        self.assertEqual(r["reason"], "too_few_reference_pairs")
        self.assertAlmostEqual(r["d_pair"], 0.5)

    def test_f_event_not_gated(self):
        D, g = ref_D(6)
        g.update({A: True, C: False})
        D[z.key(A, C)] = 0.5
        r = z.station_pair_score(D, g, E, A, C)
        self.assertFalse(r["counted"])
        self.assertEqual(r["reason"], "f_event_not_gated")
        self.assertEqual(r["n_ref"], 15)

    def test_missing_pair(self):
        D, g = ref_D(5)
        g.update({A: True, B: True})
        r = z.station_pair_score(D, g, E, A, B)
        self.assertFalse(r["counted"])
        self.assertEqual(r["reason"], "no_pair_D")


class TestDecide(unittest.TestCase):
    def test_rule(self):
        self.assertEqual(z.decide([]), ("untestable", 0))
        self.assertEqual(z.decide([True, False]), ("untestable", 1))
        self.assertEqual(z.decide([True, True, False]), ("S", 2))
        self.assertEqual(z.decide([False, False, False]), ("N", 0))
        self.assertEqual(z.decide([True, False, False]), ("N", 1))
        self.assertEqual(z.decide([True, True, False, False, False]), ("inconclusive", 2))
        self.assertEqual(z.decide([True] * 2 + [False] * 8), ("N", 2))
        self.assertEqual(z.decide([True] * 5 + [False] * 5), ("S", 5))
        self.assertEqual(z.decide([True] * 3 + [False] * 7), ("inconclusive", 3))


class TestApply(unittest.TestCase):
    def test_mapping(self):
        self.assertEqual(z.apply_g2("violated", {"2005-2012": "N", "2005-2022": "N"}),
                         ("violated", "source-checked"))
        self.assertEqual(z.apply_g2("holds", {"2005-2012": "S", "2005-2022": "N"}),
                         ("untestable (source)", "different source"))
        self.assertEqual(z.apply_g2("violated", {"2005-2012": "N", "2005-2022": "S"}),
                         ("untestable (source)", "different source"))
        self.assertEqual(z.apply_g2("violated", {"2005-2012": "N", "2005-2022": "inconclusive"}),
                         ("violated", "source check incomplete"))
        self.assertEqual(z.apply_g2(None, {}), ("pending zw33-v1", "source check incomplete"))
        # the report-only pair never changes the grade
        self.assertEqual(z.apply_g2("holds", {"2005-2012": "N", "2005-2022": "N", "2012-2022": "S"}),
                         ("holds", "source-checked"))


class TestGradePairs(unittest.TestCase):
    def test_grade(self):
        sc = {}
        for st in ("IU.AAA", "IU.BBB", "GT.CCC"):
            sc[st] = {"2005-2012": dict(counted=True, q=1.0, hit=True),
                      "2005-2022": dict(counted=True, q=0.5, hit=False),
                      "2012-2022": dict(counted=False, q=None, hit=None)}
        g = z.grade_pairs(sc)
        self.assertEqual(g["2005-2012"]["decision"], "S")
        self.assertEqual(g["2005-2012"]["hits"], 3)
        self.assertEqual(g["2005-2022"]["decision"], "N")
        self.assertAlmostEqual(g["2005-2022"]["median_q"], 0.5)
        self.assertEqual(g["2012-2022"]["decision"], "untestable")
        self.assertIsNone(g["2012-2022"]["median_q"])
        self.assertTrue(g["2005-2012"]["graded"])
        self.assertFalse(g["2012-2022"]["graded"])
        self.assertEqual(g["2005-2012"]["stations"], ["GT.CCC", "IU.AAA", "IU.BBB"])


class TestConstants(unittest.TestCase):
    def test_constants(self):
        self.assertEqual(len(z.OLD8), 8)
        self.assertEqual(z.label((A, C)), "2005-2022")
        self.assertEqual(z.label((A, B)), "2005-2012")
        self.assertEqual([e[1][:10] for e in z.F_EVENTS], ["2005-06-10", "2012-12-03", "2022-05-19"])
        self.assertTrue(z.MANTLEP2_COMMIT.startswith("b4712f0e"))
        self.assertEqual(z.HIT_Q, 0.93)
        self.assertEqual(z.MIN_REF_PAIRS, 10)
        self.assertIsNone(z.read_p35("/nonexistent/summary.json"))


if __name__ == "__main__":
    unittest.main()

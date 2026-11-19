"""Tests for m27.py (stdlib only; no network). Part of inner-core-repeaters by MV2, Munim, Inc.
(made with AI agents). Licence: Apache-2.0."""
import math
import unittest

import m27 as M

NEW = "us7000ty3q"
B = "usp000a1fh"   # 2000-10-08
MP = "us6000lfpq"  # 2023-10-15


def prow(a, b, d, ok=True):
    # values written as yka.analyze_multiplet writes them via csv ("True"/"False", numbers as text)
    return dict(ev1=a, ev2=b, D=str(d), D_c310=str(d), D_c13=str(d), ok="True" if ok else "False")


def times():
    return {e[0]: M.decimal_year(e[1]) for e in M.Z006}


class TestTime(unittest.TestCase):
    def test_decimal_year_matches_catalog(self):
        # values from results/m0324-v1/m0324_catalog.csv (column t)
        self.assertAlmostEqual(M.decimal_year("2000-10-08T14:14:45.800Z"), 2000.7693813815524, places=7)
        self.assertAlmostEqual(M.decimal_year("2026-11-13T19:04:07.703Z"), 2026.8679302290398, places=7)
        self.assertAlmostEqual(M.decimal_year("2023-10-15T08:58:34.423Z"), 2023.7873260534943, places=7)

    def test_clock_picks(self):
        t = times()
        tn = t[NEW]
        self.assertTrue(2001.1 < M.equiv(tn) < 2001.2)
        # B partner is the 2000 event, M partner the 2023 event, as pre-registered
        b = min((e for e in t if e != NEW), key=lambda e: abs(M.equiv(tn) - M.equiv(t[e])))
        m = min((e for e in t if e != NEW), key=lambda e: abs(tn - t[e]))
        self.assertEqual(b, M.FAMILIES["Z006"]["b_partner"])
        self.assertEqual(m, M.FAMILIES["Z006"]["m_partner"])

    def test_family_lists(self):
        for fam, spec in M.FAMILIES.items():
            ids = [e[0] for e in spec["events"]]
            self.assertIn(spec["new"], ids)
            self.assertEqual(len(ids), len(set(ids)))
            self.assertEqual(ids[-1], spec["new"])          # newest last
            ts = [M.decimal_year(e[1]) for e in spec["events"]]
            self.assertEqual(ts, sorted(ts))


class TestContrast(unittest.TestCase):
    def test_labels(self):
        self.assertEqual(M.key_contrast({B: 0.05, MP: 0.10}, B, MP)["label"], "B")
        self.assertEqual(M.key_contrast({B: 0.10, MP: 0.05}, B, MP)["label"], "M")
        self.assertEqual(M.key_contrast({B: 0.10, MP: 0.12}, B, MP)["label"], "inconclusive")
        self.assertEqual(M.key_contrast({B: 0.5, MP: 0.75}, B, MP)["label"], "B")        # exactly 1.5
        self.assertEqual(M.key_contrast({MP: 0.15}, B, MP)["label"], "untestable")
        self.assertEqual(M.key_contrast({B: float("nan"), MP: 0.15}, B, MP)["label"], "untestable")
        r = M.key_contrast({B: 0.0, MP: 0.1}, B, MP)
        self.assertEqual(r["label"], "B")
        self.assertIsNone(r["K"])

    def test_from_rows_ok_filter(self):
        rows = [prow(B, NEW, 0.04), prow(MP, NEW, 0.12, ok=False), prow(B, MP, 0.5)]
        pv = M.partner_values(rows, NEW, "D", True)
        self.assertEqual(pv, {B: 0.04})
        self.assertEqual(M.key_contrast(pv, B, MP)["label"], "untestable")
        pv2 = M.partner_values(rows, NEW, "D", False)
        self.assertEqual(M.key_contrast(pv2, B, MP)["label"], "B")
        self.assertAlmostEqual(M.key_contrast(pv2, B, MP)["K"], 3.0)

    def test_empty_values_skipped(self):
        rows = [dict(ev1=B, ev2=NEW, D="", ok="True"), prow(MP, NEW, 0.1)]
        self.assertEqual(M.partner_values(rows, NEW, "D"), {MP: 0.1})


class TestClock(unittest.TestCase):
    def test_backtracking_truth(self):
        t = times()
        tn = t[NEW]
        pv = {e: 0.01 + abs(M.equiv(tn) - M.equiv(t[e])) for e in t if e != NEW}
        r = M.clock_test(pv, tn, t)
        self.assertEqual(r["n"], 6)
        self.assertAlmostEqual(r["rho_B"], 1.0)
        self.assertEqual(r["best_match"], B)
        self.assertEqual(r["b_pick"], B)
        self.assertEqual(r["m_pick"], MP)

    def test_monotone_truth(self):
        t = times()
        tn = t[NEW]
        pv = {e: 0.01 + abs(tn - t[e]) for e in t if e != NEW}
        r = M.clock_test(pv, tn, t)
        self.assertAlmostEqual(r["rho_M"], 1.0)
        self.assertEqual(r["best_match"], MP)
        self.assertLess(r["rho_B"], 0.5)

    def test_too_few(self):
        t = times()
        r = M.clock_test({B: 0.1, MP: 0.2}, t[NEW], t)
        self.assertTrue(math.isnan(r["rho_B"]) and math.isnan(r["rho_M"]))
        self.assertEqual(r["best_match"], B)


class TestSource(unittest.TestCase):
    def test_ratios(self):
        g = {NEW: True, B: True, MP: True}
        D = {tuple(sorted((NEW, B))): 0.2, tuple(sorted((NEW, MP))): 0.3}
        self.assertAlmostEqual(M.source_ratios(D, g, NEW, B, MP), 1.5)
        g[B] = False
        self.assertTrue(math.isnan(M.source_ratios(D, g, NEW, B, MP)))

    def test_headline(self):
        self.assertEqual(M.headline("B", [1.0, 1.1]), "B-supported")
        self.assertEqual(M.headline("B", [1.6, 2.0]), "source-confounded:B")
        self.assertEqual(M.headline("M", [1.0, 0.9]), "M-supported")
        self.assertEqual(M.headline("M", [0.5, 0.6]), "source-confounded:M")
        self.assertEqual(M.headline("B", [1.0]), "source-unverified:B")
        self.assertEqual(M.headline("B", [1.0, float("nan")]), "source-unverified:B")
        self.assertEqual(M.headline("inconclusive", [1.0, 1.0]), "inconclusive")
        self.assertEqual(M.headline("untestable", [1.0, 1.0]), "untestable")

    def test_membership_score(self):
        olds = [e[0] for e in M.Z006 if e[0] != NEW]
        D = {}
        for i, a in enumerate(olds):
            for b in olds[i + 1:]:
                D[tuple(sorted((a, b)))] = 0.10
            D[tuple(sorted((a, NEW)))] = 0.5
        g = {e: True for e in olds + [NEW]}
        r = M.source_score(D, g, olds, NEW)
        self.assertTrue(r["counted"])
        self.assertEqual(r["q"], 1.0)
        g[NEW] = False
        self.assertFalse(M.source_score(D, g, olds, NEW)["counted"])


class TestHelpers(unittest.TestCase):
    def test_median_spearman_json(self):
        self.assertEqual(M.median([3, "", 1, float("nan"), 2]), 2.0)
        self.assertAlmostEqual(M.spearman([1, 2, 3, 4], [4, 3, 2, 1]), -1.0)
        self.assertEqual(M.jsonable({"a": float("inf"), "b": (1, 2)}), {"a": None, "b": [1, 2]})


if __name__ == "__main__":
    unittest.main()

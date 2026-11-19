"""Tests for m26.py (stdlib only; no network). Part of inner-core-repeaters by MV2, Munim, Inc.
(made with AI agents). Licence: Apache-2.0."""
import math
import unittest

import m26 as M


def prow(a, b, d310, ok=True, D=0.005, c13=0.004):
    return dict(ev1=a, ev2=b, D_c310=d310, D=D, D_c13=c13, ok="True" if ok else "False")


class TestTime(unittest.TestCase):
    def test_decimal_year_matches_z004_convention(self):
        # values from z004_events.csv / yka_z004_pairs.csv
        self.assertAlmostEqual(M.decimal_year("2003-03-21T10:47:47.620Z"), 2003.2176708403094, places=7)
        self.assertAlmostEqual(M.decimal_year("2025-12-13T14:04:33.132Z"), 2025.9495520399544, places=7)

    def test_new_event_time_and_clock(self):
        t = M.event_row(M.NEW_EVENT)["t"]
        self.assertTrue(2026.86 < t < 2026.88)
        self.assertAlmostEqual(M.equiv(2025.9495520399544), 2001.52, places=2)
        self.assertTrue(2001.1 < M.equiv(t) < 2001.2)
        self.assertEqual(M.equiv(2003.0), 2003.0)


class TestState(unittest.TestCase):
    def test_classify(self):
        self.assertEqual(M.classify_state(0.01, 0.12, 5, 3), "X")
        self.assertEqual(M.classify_state(0.12, 0.02, 5, 3), "Y")
        self.assertEqual(M.classify_state(0.09, 0.08, 5, 3), "new-state")
        self.assertEqual(M.classify_state(0.04, 0.06, 5, 3), "ambiguous")
        self.assertEqual(M.classify_state(0.03, 0.08, 5, 3), "ambiguous")   # ratio < 3
        self.assertEqual(M.classify_state(0.01, 0.12, 2, 3), "untestable")
        self.assertEqual(M.classify_state(0.01, 0.12, 5, 1), "untestable")
        self.assertEqual(M.classify_state(float("nan"), 0.12, 5, 3), "untestable")

    def test_state_from_rows(self):
        n = M.NEW_ID
        rows = [prow(e, n, 0.01 + 0.001 * k) for k, e in enumerate(M.X_GROUP)]
        rows += [prow(e, n, 0.12) for e in M.Y_GROUP]
        rows.append(prow("us6000ittv", n, 0.5))           # excluded event never counts
        rows.append(prow(M.X_GROUP[0], M.X_GROUP[1], 0.9))  # pair without the new event
        r = M.state_test(rows, n)
        self.assertEqual((r["nX"], r["nY"]), (5, 3))
        self.assertAlmostEqual(r["dX"], 0.012)
        self.assertEqual(r["label"], "X")

    def test_ok_filter(self):
        n = M.NEW_ID
        rows = [prow(e, n, 0.01, ok=(k < 2)) for k, e in enumerate(M.X_GROUP)]
        rows += [prow(e, n, 0.12) for e in M.Y_GROUP]
        self.assertEqual(M.state_test(rows, n)["label"], "untestable")
        self.assertEqual(M.state_test(rows, n, ok_only=False)["label"], "X")

    def test_backtracking_order(self):
        n = M.NEW_ID
        times = {"usp000bt5x": 2003.2177, "us7000am0e": 2020.5267, "us6000h2wz": 2022.1808,
                 "us6000mk65": 2024.2189, "us6000runq": 2025.9496}
        t_new = M.event_row(M.NEW_EVENT)["t"]
        # D grows with effective-time distance -> rho = +1
        rows = [prow(e, n, 0.001 + abs(M.equiv(t_new) - M.equiv(t))) for e, t in times.items()]
        r = M.backtracking_order(rows, n, t_new, times)
        self.assertEqual(r["n"], 5)
        self.assertAlmostEqual(r["rho"], 1.0)


class TestSource(unittest.TestCase):
    def make(self, new_level):
        old = list(M.OLD8)
        D = {}
        for i, a in enumerate(old):
            for b in old[i + 1:]:
                D[tuple(sorted((a, b)))] = 0.10 + 0.001 * (old.index(a) + old.index(b))
            D[tuple(sorted((a, M.NEW_ID)))] = new_level
        gated = {e: True for e in old + [M.NEW_ID]}
        return D, gated

    def test_same_source(self):
        D, g = self.make(0.09)
        r = M.source_score(D, g, list(M.OLD8), M.NEW_ID)
        self.assertTrue(r["counted"])
        self.assertEqual(r["q"], 0.0)

    def test_different_source(self):
        D, g = self.make(0.5)
        r = M.source_score(D, g, list(M.OLD8), M.NEW_ID)
        self.assertEqual(r["q"], 1.0)

    def test_not_counted(self):
        D, g = self.make(0.09)
        g[M.NEW_ID] = False
        self.assertFalse(M.source_score(D, g, list(M.OLD8), M.NEW_ID)["counted"])
        D, g = self.make(0.09)
        for e in M.OLD8[:4]:
            g[e] = False
        self.assertFalse(M.source_score(D, g, list(M.OLD8), M.NEW_ID)["counted"])

    def test_decide(self):
        self.assertEqual(M.decide_source([0.1, 0.5]), "untestable")
        self.assertEqual(M.decide_source([0.1, 0.5, 0.75]), "same")
        self.assertEqual(M.decide_source([1.0, 1.0, 0.2]), "different")
        self.assertEqual(M.decide_source([0.8, 0.9, 1.0, 1.0]), "inconclusive")
        self.assertEqual(M.decide_source([0.1, float("nan"), 0.2]), "untestable")

    def test_headline(self):
        self.assertEqual(M.headline("same", "X"), "X-confirmed")
        self.assertEqual(M.headline("same", "Y"), "change-detected")
        self.assertEqual(M.headline("same", "new-state"), "change-detected")
        self.assertEqual(M.headline("same", "ambiguous"), "ambiguous")
        self.assertEqual(M.headline("different", "X"), "not-a-repeat")
        self.assertEqual(M.headline("untestable", "X"), "source-unverified:X")
        self.assertEqual(M.headline("same", "untestable"), "untestable")


class TestHelpers(unittest.TestCase):
    def test_median_and_truthy(self):
        self.assertEqual(M.median([3, "", 1, float("nan"), 2]), 2.0)
        self.assertTrue(math.isnan(M.median([])))
        self.assertTrue(M.truthy("True") and M.truthy(True))
        self.assertFalse(M.truthy("False") or M.truthy(""))

    def test_max_abs_diff(self):
        a = [dict(ev1="x", ev2="y", D="0.10"), dict(ev1="y", ev2="z", D="0.20")]
        b = [dict(ev1="y", ev2="x", D="0.13"), dict(ev1="q", ev2="z", D="0.0")]
        r = M.max_abs_diff(a, b, "D")
        self.assertEqual(r["n"], 1)
        self.assertAlmostEqual(r["max_abs_diff"], 0.03)
        sa = [dict(station="S", ev1="x", ev2="y", D_late="0.1")]
        sb = [dict(station="T", ev1="x", ev2="y", D_late="0.2")]
        self.assertEqual(M.max_abs_diff(sa, sb, "D_late", stations=True)["n"], 0)

    def test_spearman(self):
        self.assertAlmostEqual(M.spearman([1, 2, 3, 4], [10, 20, 30, 40]), 1.0)
        self.assertAlmostEqual(M.spearman([1, 2, 3, 4], [4, 3, 2, 1]), -1.0)
        self.assertTrue(math.isnan(M.spearman([1, 2], [1, 2])))

    def test_jsonable(self):
        self.assertEqual(M.jsonable({"a": float("nan"), "b": [1.5, (2, 3)]}), {"a": None, "b": [1.5, [2, 3]]})


if __name__ == "__main__":
    unittest.main()

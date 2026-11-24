"""Tests for m34.py (stdlib only; no network). Part of inner-core-repeaters by MV2, Munim, Inc.
(made with AI agents). Licence: Apache-2.0."""
import math
import unittest

import m34 as M


def prow(a, b, d310, ok=True):
    return dict(ev1=a, ev2=b, D_c310=d310, D=0.005, D_c13=0.004, ok="True" if ok else "False")


class TestTime(unittest.TestCase):
    def test_decimal_years(self):
        self.assertAlmostEqual(M.decimal_year("2003-03-21T10:47:47.620Z"), 2003.2176708403094, places=7)
        self.assertAlmostEqual(M.decimal_year("2008-02-17T19:32:09.100Z"), 2008.130639, places=5)
        self.assertAlmostEqual(M.decimal_year("2014-12-09T19:23:36.400Z"), 2014.939200, places=5)
        self.assertAlmostEqual(M.decimal_year("2020-07-11T18:58:06.501Z"), 2020.526750, places=5)

    def test_equiv_and_boundary(self):
        self.assertEqual(M.equiv(2003.0), 2003.0)
        self.assertAlmostEqual(M.equiv(2014.939200), 2005.92432, places=4)
        self.assertAlmostEqual(M.equiv(2018.0, 2009.0, 3.0), 2006.0, places=9)
        lo, hi = M.wang_boundary()
        self.assertAlmostEqual(lo, 2003.6893, places=3)
        self.assertAlmostEqual(hi, 2005.9243, places=3)

    def test_wang_prediction(self):
        self.assertEqual(M.wang_prediction(2003.37), "X")
        self.assertEqual(M.wang_prediction(2006.74), "Y")
        self.assertEqual(M.wang_prediction(2016.0), "free")     # equiv 2005.5
        self.assertEqual(M.wang_prediction(2021.0), "X")        # equiv 2003.5
        self.assertEqual(M.wang_prediction(2013.0), "Y")        # equiv 2006.7

    def test_candidates(self):
        self.assertEqual(len(M.CANDIDATES), 31)
        ids = [c[0] for c in M.CANDIDATES]
        self.assertEqual(len(set(ids)), 31)
        self.assertFalse(set(ids) & set(M.OLD8))
        ts = [M.cand_dict(c)["t"] for c in M.CANDIDATES]
        self.assertEqual(sum(1 for t in ts if t < 2011.0), 11)
        self.assertTrue(all(2003.22 < t < 2008.13 or 2014.94 < t < 2020.52 for t in ts))


class TestState(unittest.TestCase):
    def test_classify(self):
        self.assertEqual(M.classify34(0.01, 0.12, 5, 3), "X")
        self.assertEqual(M.classify34(0.12, 0.02, 5, 3), "Y")
        self.assertEqual(M.classify34(0.05, 0.16, 5, 3), "X")
        self.assertEqual(M.classify34(0.06, 0.15, 5, 3), "ambiguous")      # 0.15 < 3 * 0.06
        self.assertEqual(M.classify34(0.03, 0.06, 5, 3), "intermediate")
        self.assertEqual(M.classify34(0.12, 0.11, 5, 3), "new-state")
        self.assertEqual(M.classify34(0.01, 0.12, 2, 3), "untestable")
        self.assertEqual(M.classify34(0.01, 0.12, 5, 1), "untestable")
        self.assertEqual(M.classify34(float("nan"), 0.12, 5, 3), "untestable")

    def test_state_from_rows(self):
        n = "cand"
        rows = [prow(e, n, 0.01 + 0.001 * k) for k, e in enumerate(M.X_GROUP)]
        rows += [prow(e, n, 0.12) for e in M.Y_GROUP]
        rows.append(prow("us6000ittv", n, 0.5))
        rows.append(prow(M.X_GROUP[0], M.X_GROUP[1], 0.9))
        r = M.state_test(rows, n)
        self.assertEqual((r["nX"], r["nY"]), (5, 3))
        self.assertAlmostEqual(r["dX"], 0.012)
        self.assertAlmostEqual(r["s"], 0.012 / 0.132)
        self.assertEqual(r["label"], "X")

    def test_ok_filter(self):
        n = "cand"
        rows = [prow(e, n, 0.01, ok=(k < 2)) for k, e in enumerate(M.X_GROUP)]
        rows += [prow(e, n, 0.12) for e in M.Y_GROUP]
        self.assertEqual(M.state_test(rows, n)["label"], "untestable")
        self.assertEqual(M.state_test(rows, n, ok_only=False)["label"], "X")


class TestSource(unittest.TestCase):
    def make(self, new_level, new="cand"):
        old = list(M.OLD8)
        D = {}
        for i, a in enumerate(old):
            for b in old[i + 1:]:
                D[tuple(sorted((a, b)))] = 0.10 + 0.001 * (old.index(a) + old.index(b))
            D[tuple(sorted((a, new)))] = new_level
        return D, {e: True for e in old + [new]}

    def test_same_and_different(self):
        D, g = self.make(0.09)
        self.assertEqual(M.source_score(D, g, list(M.OLD8), "cand")["q"], 0.0)
        D, g = self.make(0.5)
        self.assertEqual(M.source_score(D, g, list(M.OLD8), "cand")["q"], 1.0)

    def test_not_counted(self):
        D, g = self.make(0.09)
        g["cand"] = False
        self.assertFalse(M.source_score(D, g, list(M.OLD8), "cand")["counted"])

    def test_decide(self):
        self.assertEqual(M.decide_source([0.1, 0.5]), "untestable")
        self.assertEqual(M.decide_source([0.1, 0.5, 0.75]), "same")
        self.assertEqual(M.decide_source([1.0, 1.0, 0.2]), "different")
        self.assertEqual(M.decide_source([0.8, 0.9, 1.0, 1.0]), "inconclusive")


class TestMirror(unittest.TestCase):
    def test_old_members_consistent(self):
        old = M.old_members()
        self.assertEqual(M.violations(old, M.WANG_TURN, M.WANG_RATIO), [])
        comp = M.compatible_clocks(old)
        self.assertIn((2008.5, 2.5), comp)
        self.assertNotIn((2007.0, 2.0), comp)      # 2014-12 Y maps to 2003.03, before the 2003.22 X

    def test_untestable(self):
        g = M.mirror_tests({})
        self.assertEqual(g["R1"]["label"], "untestable")
        self.assertEqual(g["R2"]["label"], "untestable")
        self.assertGreater(g["baseline_compatible"], 0)

    def test_agreeing_member(self):
        g = M.mirror_tests({"a": (2006.74, "Y")})
        self.assertEqual(g["R1"]["label"], "consistent")
        self.assertEqual(g["R2"]["label"], "compatible")
        self.assertTrue(g["R2"]["wang_in_set"])

    def test_single_wang_violation(self):
        new = {"a": (2006.74, "X")}
        g = M.mirror_tests(new)
        self.assertEqual(g["R1"]["label"], "violated-fragile")
        self.assertEqual(g["R1"]["violations"], [("a", "usc000t4wj")])
        self.assertEqual(g["R2"]["label"], "compatible")
        self.assertFalse(g["R2"]["wang_in_set"])
        allm = M.old_members()
        allm.update(new)
        self.assertIn((2009.0, 3.0), M.compatible_clocks(allm))

    def test_return_flip_flop(self):
        g = M.mirror_tests({"a": (2016.0, "X"), "b": (2018.0, "Y")})
        self.assertEqual(g["R1"]["label"], "violated-fragile")
        self.assertEqual(g["R2"]["label"], "contradicted-fragile")
        self.assertEqual(g["R2"]["n_compatible"], 0)

    def test_robust_flip_flops(self):
        g = M.mirror_tests({"a": (2016.0, "X"), "b": (2017.0, "Y"), "c": (2018.0, "X"), "d": (2019.0, "Y")})
        self.assertEqual(g["R1"]["label"], "violated")
        self.assertEqual(g["R2"]["label"], "contradicted")

    def test_intermediate_order(self):
        # intermediate between the blocks on the Wang clock: fine; after the Y block: violation
        self.assertEqual(M.mirror_tests({"i": (2016.0, "I")})["R1"]["label"], "consistent")
        self.assertEqual(M.mirror_tests({"i": (2003.37, "I")})["R1"]["label"], "violated-fragile")

    def test_brackets(self):
        b = M.brackets({})
        self.assertEqual(b["return"]["label"], "clean")
        self.assertAlmostEqual(b["return"]["last_before"], 2014.9392, places=3)
        self.assertAlmostEqual(b["return"]["first_after"], 2020.5268, places=3)
        self.assertAlmostEqual(b["outbound"]["width_yr"], 4.913, places=2)
        b = M.brackets({"a": (2016.0, "X"), "b": (2018.0, "Y")})
        self.assertEqual(b["return"]["label"], "overlap")
        self.assertEqual(b["return"]["n_out_of_order"], 2)


class TestGrade(unittest.TestCase):
    def test_grade_empty(self):
        table, g = M.grade({}, {})
        self.assertEqual(len(table), 31)
        self.assertEqual(g["n_confirmed"], 0)
        self.assertEqual(g["R1"]["label"], "untestable")
        self.assertTrue(all(r["state"] == "untestable" for r in table))

    def test_grade_one_member(self):
        cid = "usp000eu1v"   # 2006-09-29: the Wang clock requires Y
        rows = [prow(e, cid, 0.012) for e in M.X_GROUP] + [prow(e, cid, 0.13) for e in M.Y_GROUP]
        y = {cid: dict(status="ran", state=M.state_test(rows, cid))}
        m = {cid: dict(decision="same", n_stations=5, median_q=0.4)}
        table, g = M.grade(y, m)
        row = [r for r in table if r["id"] == cid][0]
        self.assertEqual(row["state"], "X")
        self.assertEqual(row["wang_prediction"], "Y")
        self.assertEqual(row["wang_check"], "disagrees")
        self.assertEqual(g["n_used"], 1)
        self.assertEqual(g["R1"]["label"], "violated-fragile")


class TestHelpers(unittest.TestCase):
    def test_median_truthy_json(self):
        self.assertEqual(M.median([3, "", 1, float("nan"), 2]), 2.0)
        self.assertTrue(math.isnan(M.median([])))
        self.assertTrue(M.truthy("True") and not M.truthy("False"))
        self.assertEqual(M.jsonable({"a": float("nan"), "b": [1.5, (2, 3)]}), {"a": None, "b": [1.5, [2, 3]]})

    def test_max_abs_diff(self):
        a = [dict(ev1="x", ev2="y", D="0.10")]
        b = [dict(ev1="y", ev2="x", D="0.13")]
        self.assertAlmostEqual(M.max_abs_diff(a, b, "D")["max_abs_diff"], 0.03)


if __name__ == "__main__":
    unittest.main()

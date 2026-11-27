"""Tests for m2core.py (pre-registration 38, arm M). Expected values derived by hand
(findings/2026-11-27-prereg38-armM-plan.md, "Hand checks").
Part of inner-core-repeaters by MV2, Munim, Inc. (made with AI agents). Licence: Apache-2.0."""
import os
import tempfile
import unittest

import m2core as m

T = {"a": 2003.0, "b": 2006.0, "c": 2011.0, "d": 2013.0}
# pattern "mirror": {a, d} vs {b, c} -> needs a turn (ordered at T 2009 r 1, not by date order)
# pattern "date":   {a, b} vs {c, d} -> ordered by date order (no turn)
GROUPS = {"mirror": ({"a", "d"}, {"b", "c"}), "date": ({"a", "b"}, {"c", "d"})}


def pr(array, mult, e1, t1, e2, t2, D, ok=True):
    return dict(array=array, mult=mult, ev1=e1, ev2=e2, t1=t1, t2=t2, D=D, ok=ok)


def multiplet(mult, pattern, array="ILAR"):
    g1, _ = GROUPS[pattern]
    out = []
    keys = sorted(T)
    for i, x in enumerate(keys):
        for y in keys[i + 1:]:
            same = (x in g1) == (y in g1)
            out.append(pr(array, mult, mult + x, T[x], mult + y, T[y], 0.05 if same else 0.2))
    return out


class TestBasics(unittest.TestCase):
    def test_hw(self):
        self.assertAlmostEqual(m.HW_ERAS[0][0], 2014 + 334 / 365.0)
        self.assertAlmostEqual(m.HW_ERAS[0][1], 2020 + 98 / 366.0)
        self.assertFalse(m.in_hw(2014.9))
        self.assertTrue(m.in_hw(2016.0))
        self.assertFalse(m.in_hw(2020.3))
        self.assertTrue(m.in_hw(2022.0))
        self.assertFalse(m.in_hw(2022.3))

    def test_clean_and_read(self):
        fd, path = tempfile.mkstemp(suffix=".csv")
        with os.fdopen(fd, "w") as f:
            f.write("array,mult,ev1,ev2,t1,t2,D,ok\n"
                    "YKA,A,e1,e2,2003.0,2010.0,0.05,True\n"
                    "ilar,A,e1,e2,2003.0,2010.0,0.30,false\n"
                    "YKA,I,f1,f2,2003.0,2010.0,0.05,True\n"
                    "YKA,-,g1,g2,2003.0,2010.0,0.05,True\n"
                    "YKA,B,h1,h2,2003.0,2010.0,,True\n")
        try:
            rows = m.read_pairs(path)
        finally:
            os.remove(path)
        self.assertEqual(len(rows), 3)
        self.assertEqual(rows[1]["array"], "ILAR")
        self.assertFalse(rows[1]["ok"])
        self.assertTrue(m.usable(rows[0]))
        self.assertFalse(m.usable(rows[2]))      # D blank -> nan

    def test_ilar_cut(self):
        vals = [0.12, 0.05, 0.11, 0.06, 0.10, 0.07, 0.09, 0.08]
        rows = [pr("ILAR", "A", "x%d" % i, 2003.0, "y%d" % i, 2004.0, v) for i, v in enumerate(vals)]
        rows.append(pr("ILAR", "A", "z1", 2003.0, "z2", 2008.0, 0.5))      # 5 yr apart: not used
        rows.append(pr("YKA", "A", "z3", 2003.0, "z4", 2003.5, 0.5))       # YKA: not used
        self.assertEqual(m.ilar_cut(rows), (0.10, 8))                      # nearest rank ceil(6)-1 = 5
        self.assertEqual(m.ilar_cut(rows[1:]), (m.CUT_ILAR_DEFAULT, 7))


class TestStates(unittest.TestCase):
    def test_mirror_states(self):
        times, D = m.unit_data(multiplet("A", "mirror"), "ILAR", "A")
        st = m.complete_linkage(list(times), D, 0.098)
        self.assertEqual(st, {"Aa": 0, "Ad": 0, "Ab": 1, "Ac": 1})

    def test_chain_complete_linkage(self):
        D = {("e1", "e2"): 0.05, ("e2", "e3"): 0.05, ("e1", "e3"): 0.15}
        self.assertEqual(m.complete_linkage(["e3", "e1", "e2"], D, 0.098), {"e1": 0, "e2": 0, "e3": 1})

    def test_missing_pair_blocks_only_unmeasured(self):
        rows = [r for r in multiplet("A", "mirror") if {r["ev1"], r["ev2"]} != {"Aa", "Ad"}]
        times, D = m.unit_data(rows, "ILAR", "A")
        st = m.complete_linkage(list(times), D, 0.098)
        self.assertEqual(len(set(st.values())), 3)
        self.assertEqual(st["Ab"], st["Ac"])

    def test_hw_exclusion_at_yka_only(self):
        rows = [pr("YKA", "A", "e1", 2003.0, "e2", 2016.0, 0.05), pr("YKA", "A", "e1", 2003.0, "e3", 2011.0, 0.05)]
        times, _ = m.unit_data(rows, "YKA", "A")
        self.assertEqual(sorted(times), ["e1", "e3"])
        times, _ = m.unit_data(rows, "YKA", "A", exclude_hw=False)
        self.assertEqual(sorted(times), ["e1", "e2", "e3"])
        rows_i = [dict(r, array="ILAR") for r in rows]
        times, _ = m.unit_data(rows_i, "ILAR", "A")
        self.assertEqual(sorted(times), ["e1", "e2", "e3"])

    def test_testable(self):
        u = m.build_units(multiplet("A", "mirror"), {"YKA": 0.098, "ILAR": 0.098})
        self.assertEqual(len(u), 1)
        self.assertTrue(u[0]["testable"])
        self.assertEqual((u[0]["n_members"], u[0]["n_states"]), (4, 2))
        rows = [r for r in multiplet("A", "mirror") if r["t1"] < 2009 and r["t2"] < 2009]
        self.assertEqual(m.build_units(rows, {"YKA": 0.098, "ILAR": 0.098}), [
            dict(array="ILAR", mult="A", items=[("Aa", 0), ("Ab", 1)], n_members=2, n_states=2, testable=False)])


class TestClocks(unittest.TestCase):
    def test_ordered_mirror(self):
        items = [("a", 0), ("b", 1), ("c", 1), ("d", 0)]
        self.assertFalse(m.ordered(items, T, None))            # a b c d -> X Y Y X
        self.assertTrue(m.ordered(items, T, (2009.0, 1.0)))    # a 2003, d 2005, b 2006, c 2007
        self.assertFalse(m.ordered(items, T, (2008.5, 2.5)))   # a 2003, b 2006, d 2006.7, c 2007.5

    def test_ordered_date(self):
        items = [("a", 0), ("b", 0), ("c", 1), ("d", 1)]
        self.assertTrue(m.ordered(items, T, None))
        self.assertFalse(m.ordered(items, T, (2009.0, 1.0)))   # a X, d Y, b X, c Y

    def test_grids(self):
        self.assertEqual(len(m.grid_clocks(m.TEST_TURNS, m.TEST_RATIOS)), 2269)
        self.assertEqual(len(m.grid_clocks(m.FINE_TURNS, m.WIDE_RATIOS)), 43710)
        self.assertTrue(set(m.TEST_TURNS) <= set(m.FINE_TURNS))
        self.assertIn(0.5003, m.WIDE_RATIOS)

    def test_mirror_and_date_never_together(self):
        units = [dict(items=[("a", 0), ("b", 1), ("c", 1), ("d", 0)]),
                 dict(items=[("a", 0), ("b", 0), ("c", 1), ("d", 1)])]
        self.assertEqual(m.intersection(units, T, m.grid_clocks(m.FINE_TURNS, m.WIDE_RATIOS)), [])


class TestGrades(unittest.TestCase):
    def test_m2a(self):
        self.assertEqual(m.grade_m2a(2, True, False, 0.01, False, None), "untestable")
        self.assertEqual(m.grade_m2a(3, False, False, 0.5, True, True), "broken")
        self.assertEqual(m.grade_m2a(3, False, False, 0.5, True, False), "mixed")
        self.assertEqual(m.grade_m2a(3, True, True, 0.01, False, None), "monotone")
        self.assertEqual(m.grade_m2a(3, True, False, 0.05, False, None), "holds")
        self.assertEqual(m.grade_m2a(3, True, False, 0.06, False, None), "mixed")
        self.assertEqual(m.grade_m2a(3, False, False, 0.01, False, None), "mixed")

    def test_m2b(self):
        self.assertEqual(m.grade_m2b(3, 2), "broken")
        self.assertEqual(m.grade_m2b(8, 2), "broken")         # 2 >= 0.25 * 8
        self.assertEqual(m.grade_m2b(9, 2), "mixed")          # 2 < 0.25 * 9
        self.assertEqual(m.grade_m2b(10, 2), "mixed")
        self.assertEqual(m.grade_m2b(5, 0), "holds")
        self.assertEqual(m.grade_m2b(4, 0), "untestable")
        self.assertEqual(m.grade_m2b(4, 1), "untestable")

    def test_closure(self):
        t = {"e1": 2003.0, "e2": 2004.0, "e3": 2010.0, "e4": 2016.0, "e5": 2012.0}
        spec = [("e1", "e2", 0.05, 0.05), ("e1", "e3", 0.2, 0.05), ("e2", "e3", 0.2, 0.05),
                ("e1", "e4", 0.2, 0.05), ("e1", "e5", 0.05, 0.2)]
        rows = []
        for a, b, dy, di in spec:
            rows.append(pr("YKA", "Z", a, t[a], b, t[b], dy))
            rows.append(pr("ILAR", "Z", b, t[b], a, t[a], di))          # reversed order on purpose
        res, out = m.closure_test(rows, {"YKA": 0.098, "ILAR": 0.098})
        self.assertEqual((res["n_ilar_same"], res["violations"]), (3, 2))
        self.assertEqual((res["hw_ilar_same"], res["hw_violations"]), (1, 1))
        self.assertEqual(res["yka_same_ilar_changed"], 1)
        self.assertEqual(res["grade"], "broken")
        self.assertEqual(sum(1 for r in out if r["violation"]), 2)


class TestAnalyse(unittest.TestCase):
    def test_mixed_not_robust(self):
        rows = multiplet("A", "mirror") + multiplet("B", "date") + multiplet("C", "mirror")
        res, units, null, drops, _ = m.analyse(rows, n_shuffles=5, variants=False)
        a = res["M2a"]
        self.assertTrue(res["checks_ok"])
        self.assertEqual(res["cuts"]["ILAR"], m.CUT_ILAR_DEFAULT)
        self.assertEqual(a["n_testable_multiplets"], 3)
        self.assertFalse(a["I_fine_nonempty"])
        self.assertEqual(a["I_test_share"], 0.0)
        self.assertFalse(a["drop_one_robust"])
        self.assertIn("Bc", a["drop_one_breakers"])        # B without c: T 2009 r 0.5003 orders all
        self.assertEqual(a["grade"], "mixed")
        self.assertEqual(len(null), 5)
        self.assertEqual(len(drops), 12)
        self.assertEqual(res["M2b"]["grade"], "untestable")

    def test_monotone(self):
        rows = multiplet("A", "date") + multiplet("B", "date") + multiplet("C", "date")
        res, units, null, drops, _ = m.analyse(rows, n_shuffles=3)
        a = res["M2a"]
        self.assertTrue(a["noturn_in_I"])
        self.assertEqual(a["grade"], "monotone")
        self.assertEqual(drops, [])
        self.assertIn("without_A", res["report_only_variants"])
        self.assertEqual(res["report_only_variants"]["without_A"]["n_testable_multiplets"], 2)
        self.assertEqual(res["report_only_variants"]["yka_only"]["n_testable_multiplets"], 0)

    def test_untestable(self):
        rows = multiplet("A", "mirror") + multiplet("B", "date")
        res, _, _, _, _ = m.analyse(rows, n_shuffles=2, variants=False)
        self.assertEqual(res["M2a"]["grade"], "untestable")


if __name__ == "__main__":
    unittest.main()

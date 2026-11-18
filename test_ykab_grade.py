"""Tests for ykab_grade.py and the v0.10.1 additions to yka_b.py (stdlib only).

MV2, Munim, Inc. (made with AI agents). Apache-2.0.
"""
import csv
import json
import os
import tempfile
import unittest

import yka_b
import ykab_grade as g


def prow(mult, ev1, t1, ev2, t2, D, gen1="legacy", gen2="legacy", ok=True):
    nan = float("nan")
    return dict(mult=mult, ev1=ev1, ev2=ev2, t1=t1, t2=t2, D=D, Dx=D, Ds=None, D_c13=D, D_c310=D, D_full=D,
                cc_ref=0.95, snr_dir1=10.0, snr_dir2=10.0, snr_coda1=3.0, snr_coda2=3.0, ok=ok, gen1=gen1,
                gen2=gen2, ignored=nan)


def mrow(mult, n_yka, pick):
    return dict(mult=mult, n_yka=str(n_yka), pick_minus_pred_s="" if pick is None else str(pick))


class TestRuleA(unittest.TestCase):
    def test_counts_and_trigger(self):
        rows = [mrow("Y0", 5, -5.9), mrow("Y1", 3, -1.5), mrow("Y2", 2, 2.0), mrow("Y3", 2, -2.05),
                mrow("Y4", 1, None), mrow("Y5", 2, None)]
        r = g.rule_a(rows)
        self.assertEqual(r["n_multiplets_ge2_yka"], 5)
        self.assertEqual(r["n_with_pick"], 4)
        self.assertEqual(r["invalid"], ["Y0", "Y3"])
        self.assertEqual(r["valid"], ["Y1", "Y2"])
        self.assertAlmostEqual(r["frac_invalid_of_ge2"], 2 / 5)
        self.assertTrue(r["triggers_v2"])
        self.assertEqual(r["n_pick_before_minus4s"], 1)

    def test_no_trigger_at_quarter(self):
        rows = [mrow("Y%d" % k, 2, -1.0) for k in range(3)] + [mrow("Y9", 2, -6.0)]
        self.assertFalse(g.rule_a(rows)["triggers_v2"])  # exactly 25% is not "more than 25%"

    def test_pick_valid(self):
        self.assertTrue(g.pick_valid("-2.0"))
        self.assertFalse(g.pick_valid("-2.001"))
        self.assertFalse(g.pick_valid(""))


class TestH4(unittest.TestCase):
    def test_selection(self):
        rows = [prow("Y0", "a", 2001.0, "b", 2002.5, 0.10),            # near, same gen: in
                prow("Y0", "a", 2001.0, "c", 2004.0, 0.90),            # 3 yr apart: out
                prow("Y1", "d", 2021.0, "e", 2021.5, 0.20, "current", "current"),  # in
                prow("Y1", "d", 2013.9, "f", 2014.2, 0.70, "legacy", "current"),   # cross-gen: out
                prow("Y2", "g", 2021.0, "h", 2021.2, 0.30, "current", "current"),  # in
                prow("Y2", "g", 2021.0, "i", 2021.1, 0.05, "current", "current", ok=False)]  # not ok
        r = g.h4(rows, d_hw=0.011)
        self.assertEqual(r["n_pairs"], 3)
        self.assertEqual(r["n_legacy"], 1)
        self.assertEqual(r["n_current"], 2)
        self.assertAlmostEqual(r["median_D_same_gen_near"], 0.20)
        self.assertEqual(r["verdict"], "holds")

    def test_fails_and_untestable(self):
        rows = [prow("Y0", "a%d" % k, 2001.0, "b%d" % k, 2001.5, 0.005) for k in range(3)]
        self.assertEqual(g.h4(rows, d_hw=0.011)["verdict"], "fails")
        self.assertEqual(g.h4(rows[:2], d_hw=0.011)["verdict"], "untestable")

    def test_median(self):
        self.assertAlmostEqual(g.median([3.0, 1.0, 2.0, float("nan")]), 2.0)
        self.assertAlmostEqual(g.median([1.0, 2.0, 3.0, 4.0]), 2.5)


class TestRegrade(unittest.TestCase):
    def test_vote_restricted_to_valid(self):
        rows = []
        for k in range(7):
            m = "Y%d" % k
            rows.append(prow(m, "a%d" % k, 2001.5, "n%d" % k, 2026.0, 0.05))
            rows.append(prow(m, "b%d" % k, 2007.5, "n%d" % k, 2026.0, 0.30))
        out = g.regrade(rows, {"Y0", "Y1", "Y2", "Y3", "Y4", "Y5", "Y6"}, yka_b=yka_b, nboot=100)
        self.assertEqual(out["F1_Y"]["verdict"], "supported")
        out = g.regrade(rows, {"Y0", "Y1"}, yka_b=yka_b, nboot=100)
        self.assertEqual(out["F1_Y"]["verdict"], "untestable")
        self.assertEqual(out["n_multiplets"], 2)


class TestRunOffline(unittest.TestCase):
    def test_run_writes_json(self):
        with tempfile.TemporaryDirectory() as d:
            with open(os.path.join(d, "yka_multiplets.csv"), "w", newline="") as f:
                w = csv.DictWriter(f, fieldnames=["mult", "n_yka", "status", "pick_minus_pred_s"])
                w.writeheader()
                w.writerow(dict(mult="Y0", n_yka=3, status="ok", pick_minus_pred_s=-1.2))
                w.writerow(dict(mult="Y1", n_yka=2, status="ok", pick_minus_pred_s=-5.9))
            rows = [prow("Y0", "a", 2001.0, "b", 2001.5, 0.1), prow("Y0", "a", 2001.0, "c", 2001.6, 0.2),
                    prow("Y0", "b", 2001.5, "c", 2001.6, 0.3), prow("Y1", "x", 2001.0, "y", 2001.2, 0.9)]
            with open(os.path.join(d, "yka_pairs.csv"), "w", newline="") as f:
                w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
                w.writeheader()
                for r in rows:
                    w.writerow(r)
            res = g.run(d, yka_b=yka_b, nboot=50, fetch_v1=False)
            self.assertTrue(os.path.exists(os.path.join(d, "ykab_grade.json")))
            with open(os.path.join(d, "ykab_grade.json")) as f:
                back = json.load(f)
        self.assertEqual(res["v2_rule_a"]["invalid"], ["Y1"])
        self.assertEqual(back["H4_v2_primary"]["n_pairs"], 4)
        self.assertEqual(back["v2_valid_only"]["n_multiplets"], 1)


class TestPickRecords(unittest.TestCase):
    def _mults(self, d, rows):
        p = os.path.join(d, "m.csv")
        with open(p, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=["mult", "status", "pick_minus_pred_s"])
            w.writeheader()
            for r in rows:
                w.writerow(r)
        return p

    def test_consistent_match(self):
        with tempfile.TemporaryDirectory() as d:
            p = self._mults(d, [dict(mult="Y0", status="ok", pick_minus_pred_s="-1.1"),
                                dict(mult="Y1", status="fewer-than-2-at-YKA", pick_minus_pred_s=""),
                                dict(mult="Y2", status="no-DF-pick", pick_minus_pred_s=""),
                                dict(mult="Y3", status="ok", pick_minus_pred_s="0.4")])
            recs = [dict(pick_constrained_s=-1.1, pick_original_s=-5.9), dict(pick_constrained_s=None),
                    dict(pick_constrained_s=0.4, pick_original_s=0.4)]
            diag = yka_b.match_pick_records(p, recs)
            self.assertTrue(diag["consistent"])
            self.assertEqual([r["mult"] for r in diag["rows"]], ["Y0", "Y2", "Y3"])
            self.assertEqual(diag["rows"][0]["pick_original_s"], -5.9)
            out = os.path.join(d, "diag.csv")
            yka_b.write_pick_diag(out, diag)
            with open(out, newline="") as f:
                self.assertEqual(len(list(csv.DictReader(f))), 3)

    def test_inconsistent_keeps_raw(self):
        with tempfile.TemporaryDirectory() as d:
            p = self._mults(d, [dict(mult="Y0", status="ok", pick_minus_pred_s="-1.1")])
            diag = yka_b.match_pick_records(p, [dict(pick_constrained_s=0.9)])
            self.assertFalse(diag["consistent"])
            self.assertEqual(diag["rows"], [])
            self.assertEqual(len(diag["raw"]), 1)


if __name__ == "__main__":
    unittest.main()

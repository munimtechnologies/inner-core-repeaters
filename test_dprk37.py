# Copyright 2026 MV2, Munim, Inc.
# Licensed under the Apache License, Version 2.0 (Apache-2.0). See LICENSE.
# Made with AI agents working in MV2's name.
"""Stdlib-only tests for dprk37.py (pre-registration 37)."""

import math
import unittest

import dprk37 as d


def row(**kw):
    r = {"status": "ok", "covers": True, "channel": "00.BHZ", "epoch_start": "2010-01-01T00:00:00",
         "snr_df_0.7-2.0": 6.0, "snr_ab_0.7-2.0": 8.0}
    r.update(kw)
    return r


class Geometry(unittest.TestCase):
    def test_gc_point_east(self):
        la, lo = d.gc_point(0.0, 0.0, 90.0, 90.0)
        self.assertAlmostEqual(la, 0.0, places=6)
        self.assertAlmostEqual(lo, 90.0, places=6)

    def test_round_trip(self):
        la, lo = d.gc_point(41.3, 129.05, 37.0, 160.0)
        self.assertAlmostEqual(d.gc_deg(41.3, 129.05, la, lo), 160.0, places=6)
        self.assertAlmostEqual(d.azimuth(41.3, 129.05, la, lo), 37.0, places=6)

    def test_azimuth_cardinal(self):
        self.assertAlmostEqual(d.azimuth(0.0, 0.0, 0.0, 10.0), 90.0, places=9)
        self.assertAlmostEqual(d.azimuth(0.0, 0.0, 10.0, 0.0), 0.0, places=9)

    def test_gc_zero(self):
        self.assertEqual(d.gc_deg(41.3, 129.05, 41.3, 129.05), 0.0)

    def test_icb_crossings(self):
        pts = [(0.0, 0.0), (30.0, 2889.0), (60.0, 5153.5), (90.0, 6000.0), (120.0, 5153.5),
               (150.0, 0.0)]
        self.assertEqual(d.icb_crossings(pts), (60.0, 120.0))
        self.assertIsNone(d.icb_crossings([(0.0, 0.0), (10.0, 2889.0)]))

    def test_is_patch(self):
        self.assertTrue(d.is_patch([(0.0, 0.0)], [(0.0, 9.9)]))
        self.assertFalse(d.is_patch([(0.0, 0.0)], [(0.0, 10.1)]))


class Numbers(unittest.TestCase):
    def test_parabolic(self):
        f = [-(x - 0.3) ** 2 for x in (-1.0, 0.0, 1.0)]
        delta, peak = d.parabolic_peak(*f)
        self.assertAlmostEqual(delta, 0.3, places=9)
        self.assertAlmostEqual(peak, 0.0, places=9)
        self.assertEqual(d.parabolic_peak(1.0, 1.0, 1.0), (0.0, 1.0))

    def test_expected_cc(self):
        self.assertAlmostEqual(d.expected_cc(2.0, 2.0), 0.75, places=12)
        self.assertEqual(d.expected_cc(0.0, 5.0), 0.0)
        self.assertEqual(d.signal_fraction(1.0), 0.0)
        self.assertEqual(d.signal_fraction(""), 0.0)

    def test_excess_x(self):
        e_df = d.expected_cc(10.0, 10.0)
        e_ab = d.expected_cc(20.0, 20.0)
        self.assertAlmostEqual(d.excess_x(e_df, e_ab, 10.0, 10.0, 20.0, 20.0), 0.0, places=12)
        self.assertAlmostEqual(d.excess_x(e_df - 0.09, e_ab, 10.0, 10.0, 20.0, 20.0), 0.09,
                               places=9)

    def test_quantile(self):
        self.assertEqual(d.quantile([1, 2, 3, 4, 5], 0.5), 3)
        self.assertAlmostEqual(d.quantile([1, 2, 3, 4], 0.25), 1.75, places=12)
        self.assertAlmostEqual(d.quantile([5, 1], 0.975), 4.9, places=12)
        self.assertTrue(math.isnan(d.quantile([], 0.5)))

    def test_fit_cos(self):
        rows = []
        for k in range(8):
            baz = 45.0 * k
            rows.append((baz, 0.03 + 0.01 * math.cos(math.radians(baz))
                         - 0.02 * math.sin(math.radians(baz))))
        c, a, b = d.fit_cos(rows)
        self.assertAlmostEqual(c, 0.03, places=9)
        self.assertAlmostEqual(a, 0.01, places=9)
        self.assertAlmostEqual(b, -0.02, places=9)

    def test_fit_degenerate(self):
        self.assertIsNone(d.fit_cos([(10.0, 0.1)] * 6))
        self.assertIsNone(d.fit_cos([(0.0, 1.0), (90.0, 2.0)]))

    def test_bootstrap_deterministic(self):
        items = [{"X": v} for v in (0.0, 0.1, 0.2, 0.3)]
        b1 = d.bootstrap(items, d.stat_median_x, n=50)
        b2 = d.bootstrap(items, d.stat_median_x, n=50)
        self.assertEqual(b1, b2)
        self.assertEqual(len(b1), 50)
        self.assertTrue(all(0.0 <= v <= 0.3 for v in b1))

    def test_perm_separated(self):
        # 10 vs 10 well-separated values: only splits that keep both middle pairs reach d
        # (70 of C(20,10) = 184,756), so p is far below 0.05.
        hi = [1.0 + 0.01 * k for k in range(10)]
        lo = [0.01 * k for k in range(10)]
        dd, p = d.perm_p(hi, lo, n=2000)
        self.assertAlmostEqual(dd, 1.0, places=9)
        self.assertLess(p, 0.01)

    def test_perm_four_tied(self):
        # Four tied values per group (P37b's minimum): median(A) - median(B) >= 1 whenever A
        # holds at least three of the four 1s, so the exact p is (4*4 + 1)/C(8,4) = 17/70.
        # This documents that P37b cannot reach p <= 0.05 with such groups.
        dd, p = d.perm_p([1.0] * 4, [0.0] * 4, n=2000)
        self.assertEqual(dd, 1.0)
        self.assertAlmostEqual(p, 17.0 / 70.0, delta=0.04)

    def test_perm_null(self):
        dd, p = d.perm_p([0.5] * 4, [0.5] * 4, n=100)
        self.assertEqual(dd, 0.0)
        self.assertEqual(p, 1.0)


class Rules(unittest.TestCase):
    def test_version(self):
        self.assertEqual(d.VERSION, "0.1.1")

    def test_labels(self):
        self.assertEqual(d.label_w(0.05, 0.02, 0.08), "changed")
        self.assertEqual(d.label_w(0.01, -0.01, 0.025), "same")
        self.assertEqual(d.label_w(0.03, 0.005, 0.06), "inconclusive")
        self.assertEqual(d.label_t(0.03, 0.01, 0.05), "shift")
        self.assertEqual(d.label_t(-0.025, -0.04, -0.01), "shift")
        self.assertEqual(d.label_t(0.005, -0.01, 0.015), "no shift")
        self.assertEqual(d.label_t(0.01, -0.03, 0.04), "inconclusive")
        self.assertEqual(d.headline("changed", "no shift"), "DF changed")
        self.assertEqual(d.headline("same", "shift"), "DF changed")
        self.assertEqual(d.headline("same", "no shift"), "no change at this level")
        self.assertEqual(d.headline("same", "inconclusive"), "inconclusive")
        self.assertEqual(d.label_patch(0.04, 0.01), "patch excess")
        self.assertEqual(d.label_patch(0.005, 0.9), "no patch excess")
        self.assertEqual(d.label_patch(0.02, 0.01), "inconclusive")
        self.assertEqual(d.label_patch(0.04, 0.2), "inconclusive")

    def test_gate(self):
        bl = "0.7-2.0"
        self.assertEqual(d.gate(row(), row(), bl, True), "")
        self.assertEqual(d.gate(row(), row(covers="True"), bl, True), "")
        self.assertEqual(d.gate(row(), row(), bl, False), "windows not clear")
        self.assertEqual(d.gate(row(), row(**{"snr_df_0.7-2.0": 4.9}), bl, True), "snr")
        self.assertEqual(d.gate(row(), row(**{"snr_ab_0.7-2.0": ""}), bl, True), "snr")
        self.assertEqual(d.gate(row(), row(channel="10.BHZ"), bl, True),
                         "channel or epoch differs")
        self.assertEqual(d.gate(row(), row(epoch_start="2016-03-01T00:00:00"), bl, True),
                         "channel or epoch differs")
        self.assertEqual(d.gate(row(), row(status="gap"), bl, True), "record not ok")
        self.assertEqual(d.gate(row(covers=False), row(), bl, True),
                         "record does not cover window")
        self.assertEqual(d.gate(None, row(), bl, True), "record not ok")

    def test_all_pairs(self):
        p = d.all_pairs()
        self.assertEqual(len(p), 15)
        self.assertEqual(p[0], ("G1", "us10004bnm", "us10006n8a"))
        self.assertEqual(p[1], ("G2", "usc000f5t0", "us10004bnm"))
        self.assertEqual(len({(a, b) for _, a, b in p}), 15)
        self.assertEqual(sum(1 for lab, _, _ in p if lab.startswith("R:")), 13)

    def test_pick_branch_p(self):
        arr = [("PKIKP", 1190.0, 1.9), ("PKP", 1200.5, 4.1)]
        self.assertEqual(d.pick_branch_p(arr, 1200.5), 4.1)
        self.assertEqual(d.pick_branch_p(arr, 1190.0), 1.9)
        self.assertIsNone(d.pick_branch_p(arr, None))

    def test_summarise_untestable(self):
        meas = [{"station": "N.S%d" % k, "X": 0.0, "dt": 0.0, "baz": 40.0 * k} for k in range(7)]
        self.assertEqual(d.summarise_pair(meas, {})["headline"], "untestable")

    def test_summarise_same(self):
        xs = (-0.004, -0.002, 0.0, 0.002, 0.004, 0.001) * 2
        meas, flags = [], {}
        for k, x in enumerate(xs):
            baz = 30.0 * k
            st = "N.S%d" % k
            meas.append({"station": st, "X": x, "baz": baz,
                         "dt": 0.001 * math.cos(math.radians(baz))})
            flags[st] = k < 6
        out = d.summarise_pair(meas, flags)
        self.assertEqual(out["W"], "same")
        self.assertEqual(out["T"], "no shift")
        self.assertEqual(out["headline"], "no change at this level")
        self.assertEqual(out["P37b"], "no patch excess")

    def test_summarise_changed(self):
        meas = [{"station": "N.S%d" % k, "X": 0.1, "dt": 0.05, "baz": 36.0 * k}
                for k in range(10)]
        out = d.summarise_pair(meas, {})
        self.assertEqual(out["W"], "changed")
        self.assertEqual(out["T"], "shift")
        self.assertEqual(out["headline"], "DF changed")
        self.assertEqual(out["P37b"], "untestable")


if __name__ == "__main__":
    unittest.main()

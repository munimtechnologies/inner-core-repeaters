# Copyright 2026 MV2, Munim, Inc.
# Licensed under the Apache License, Version 2.0 (Apache-2.0). See LICENSE.
# Made with AI agents working in MV2's name.
"""Synthetic tests for ks32.py (no network, numpy only)."""

import math
import unittest

import numpy as np

import ks32

SR = 20.0
TB = 0.0
N = 6000  # 300 s
NAN = float("nan")


def pulse(t_center, sigma=3.0, amp=1.0, tb=TB):
    t = tb + np.arange(N) / SR
    return amp * np.exp(-0.5 * ((t - t_center) / sigma) ** 2)


def m(d, cc=0.95, edge=False, snr_a=10.0, snr_b=10.0):
    return dict(d=d, cc=cc, edge=edge, snr_a=snr_a, snr_b=snr_b)


class TestConstants(unittest.TestCase):
    def test_plan_constants(self):
        P = ks32.PARAMS
        self.assertEqual(P["band_long"], (0.01, 0.05))
        self.assertEqual(P["band_short"], (0.1, 0.5))
        self.assertEqual(P["advance_s"], -0.5)
        self.assertEqual(P["min_stations"], 3)
        self.assertEqual(P["cc_long"], 0.90)
        self.assertEqual(P["cc_short"], 0.80)
        self.assertEqual(P["snr_min"], 3.0)
        self.assertEqual(P["maxlag_short"], 1.5)
        self.assertEqual(ks32.ZHOU_WIN["COLA"], (1411.0, 1463.0))
        self.assertEqual(ks32.ZHOU_WIN["YAK"], (1429.0, 1471.0))
        self.assertEqual(ks32.PRIMARY_PAIR, ("E1997", "E2018"))
        self.assertEqual(ks32.REFS, ["P1", "PP", "S1", "SS"])


class TestPhaseChange(unittest.TestCase):
    def test_recovers_subsample_delay(self):
        # A: arrival at its prediction (100 s). B: predicted 103.3 s, observed 0.47 s late.
        a = pulse(100.0)
        b = pulse(103.3 + 0.47)
        d, cc, edge = ks32.phase_change(a, b, SR, TB, 100.0, 103.3, 25.0, 25.0, 8.0)
        self.assertFalse(edge)
        self.assertGreater(cc, 0.999)
        self.assertAlmostEqual(d, 0.47, delta=0.01)

    def test_negative_delay(self):
        a = pulse(120.0)
        b = pulse(118.0 - 0.9)   # predicted 118.0, observed 0.9 s early
        d, cc, edge = ks32.phase_change(a, b, SR, TB, 120.0, 118.0, 25.0, 25.0, 8.0)
        self.assertFalse(edge)
        self.assertAlmostEqual(d, -0.9, delta=0.01)

    def test_origin_error_cancels_in_differential(self):
        # Two phases 60 s apart; B has a +1.2 s origin error on both, SKS a further -0.8 s.
        a = pulse(100.0) + pulse(160.0)
        b = pulse(100.0 + 1.2 - 0.8) + pulse(160.0 + 1.2)
        d_sks, _, e1 = ks32.phase_change(a, b, SR, TB, 100.0, 100.0, 20.0, 20.0, 8.0)
        d_ref, _, e2 = ks32.phase_change(a, b, SR, TB, 160.0, 160.0, 20.0, 20.0, 8.0)
        self.assertFalse(e1)
        self.assertFalse(e2)
        self.assertAlmostEqual(d_sks, 0.4, delta=0.01)
        self.assertAlmostEqual(d_ref, 1.2, delta=0.01)
        diff = ks32.differentials({"SKS": d_sks, "S1": d_ref}, refs=["S1"])
        self.assertAlmostEqual(diff["S1"], -0.8, delta=0.02)

    def test_window_off_the_end_is_edge(self):
        a = pulse(5.0)
        b = pulse(5.0)
        d, cc, edge = ks32.phase_change(a, b, SR, TB, 5.0, 5.0, 25.0, 25.0, 8.0)
        self.assertTrue(edge)
        self.assertTrue(math.isnan(d))


class TestLagInWindow(unittest.TestCase):
    def test_absolute_window_and_prediction_shift(self):
        # Records start 1300 s after origin. A's SKS at 1440 s; B predicted 1442 s,
        # observed 1441 s (1.0 s early). Zhou-style window 1411-1463 s on A.
        tb = 1300.0
        a = pulse(1440.0, sigma=5.0, tb=tb)
        b = pulse(1441.0, sigma=5.0, tb=tb)
        tot, cc, edge = ks32.lag_in_window(a, b, SR, tb, 1411.0, 1463.0, 2.0, 8.0)
        self.assertFalse(edge)
        self.assertGreater(cc, 0.999)
        self.assertAlmostEqual(tot, 1.0, delta=0.01)
        self.assertAlmostEqual(tot - 2.0, -1.0, delta=0.01)

    def test_extra_centre_shift(self):
        # Narrower pulses, short 7-s window, search +/-1.5 s around a given centre.
        tb = 1300.0
        a = pulse(1440.0, sigma=2.0, tb=tb)
        b = pulse(1441.0, sigma=2.0, tb=tb)
        tot, _cc, edge = ks32.lag_in_window(a, b, SR, tb, 1437.0, 1444.0, 2.0, 1.5)
        self.assertFalse(edge)
        self.assertAlmostEqual(tot, 1.0, delta=0.01)
        # centre 0.0 s, true 1.0 s, maxlag 1.5 s -> found
        tot2, _cc2, edge2 = ks32.lag_in_window(a, b, SR, tb, 1437.0, 1444.0, 0.0, 1.5)
        self.assertFalse(edge2)
        self.assertAlmostEqual(tot2, 1.0, delta=0.01)
        # centre 3.0 s, maxlag 1.5 s: the true 1.0 s is 2.0 s away -> edge
        _t3, _c3, edge3 = ks32.lag_in_window(a, b, SR, tb, 1437.0, 1444.0, 3.0, 1.5)
        self.assertTrue(edge3)


class TestPeakAndWindows(unittest.TestCase):
    def test_peak_time(self):
        tb = 1300.0
        x = pulse(1452.3, sigma=1.0, amp=-2.0, tb=tb) + pulse(1470.0, sigma=1.0, amp=5.0, tb=tb)
        t = ks32.peak_time(x, SR, tb, 1450.0, 5.0)
        self.assertAlmostEqual(t, 1452.3, places=6)

    def test_zhou_window(self):
        self.assertEqual(ks32.zhou_window("COLA", 1440.0, 1442.0), (1413.0, 1465.0, "table-s2"))
        self.assertEqual(ks32.zhou_window("COLA", 1411.0, 1411.0), (1411.0, 1463.0, "table-s2"))
        self.assertEqual(ks32.zhou_window("COLA", 1470.0, 1471.0), (1461.0, 1511.0, "fallback"))
        self.assertEqual(ks32.zhou_window("XYZ", 1400.0, 1400.0), (1390.0, 1440.0, "fallback"))


class TestStationValue(unittest.TestCase):
    def test_median_over_passing_refs(self):
        refs = {"P1": m(0.1), "PP": m(0.2, cc=0.92), "S1": m(0.5, cc=0.85), "SS": m(NAN)}
        v = ks32.station_value(m(-0.9, cc=0.97), refs, 0.90, 3.0)
        self.assertEqual(v["status"], "ok")
        self.assertEqual(v["refs"], ["P1", "PP"])
        self.assertEqual(v["n_ref"], 2)
        self.assertAlmostEqual(v["x"], -1.05, places=9)
        self.assertAlmostEqual(v["c"], 0.15, places=9)

    def test_gates(self):
        refs = {"P1": m(0.1)}
        self.assertEqual(ks32.station_value(m(-0.9, cc=0.85), refs, 0.90, 3.0)["status"],
                         "sks-gate")
        self.assertEqual(ks32.station_value(m(-0.9, snr_b=2.0), refs, 0.90, 3.0)["status"],
                         "sks-gate")
        self.assertEqual(ks32.station_value(m(-0.9, snr_b=NAN), refs, 0.90, 3.0)["status"],
                         "sks-gate")
        self.assertEqual(ks32.station_value(m(-0.9, edge=True), refs, 0.90, 3.0)["status"],
                         "sks-gate")
        bad = {"P1": m(0.1, cc=0.5), "PP": m(0.1, snr_a=1.0)}
        v = ks32.station_value(m(-0.9), bad, 0.90, 3.0)
        self.assertEqual(v["status"], "no-ref")
        self.assertTrue(math.isnan(v["x"]))


class TestGrades(unittest.TestCase):
    L = ks32.T1_LABELS

    def test_grade_count(self):
        g = ks32.grade_count
        self.assertEqual(g({"A": -1.0, "B": -0.8, "C": -0.6, "D": NAN}, -0.5, 3, self.L),
                         "reproduced")
        self.assertEqual(g({"A": -1.0, "B": -0.8, "C": -0.1, "D": 0.2}, -0.5, 3, self.L),
                         "partial")
        self.assertEqual(g({"A": -1.0, "B": -0.1, "C": 0.0, "D": NAN}, -0.5, 3, self.L),
                         "not-reproduced")
        self.assertEqual(g({"A": -1.0, "B": NAN, "C": NAN, "D": -0.9}, -0.5, 3, self.L),
                         "untestable")
        self.assertEqual(g({"A": -0.5, "B": -0.5, "C": -0.5}, -0.5, 3, self.L), "reproduced")
        self.assertEqual(g({"A": -0.1, "B": 0.1, "C": 0.0, "D": -0.2}, -0.5, 3,
                           ks32.T2_LABELS), "absent")

    def test_joint(self):
        j = ks32.grade_joint
        self.assertEqual(j("reproduced", "present"), "direct-SKS-advance")
        self.assertEqual(j("reproduced", "absent"), "window-content")
        self.assertEqual(j("reproduced", "mixed"), "mixed")
        self.assertEqual(j("reproduced", "untestable"), "mixed")
        self.assertEqual(j("partial", "present"), "mixed")
        self.assertEqual(j("not-reproduced", "present"), "not-reproduced")
        self.assertEqual(j("untestable", "absent"), "untestable")


class TestDifferentials(unittest.TestCase):
    def test_skips_missing_and_nan(self):
        out = ks32.differentials({"SKS": 0.5, "P1": 0.2, "PP": NAN}, refs=["P1", "PP", "SS"])
        self.assertEqual(set(out), {"P1"})
        self.assertAlmostEqual(out["P1"], 0.3)

    def test_no_target(self):
        self.assertEqual(ks32.differentials({"P1": 0.2}), {})


class TestSnr(unittest.TestCase):
    def test_snr_of_clean_pulse_with_noise(self):
        rng = np.random.default_rng(32)
        x = pulse(150.0, amp=10.0) + rng.normal(0.0, 1.0, N)
        s = ks32.snr_at(x, SR, TB, 150.0, (-60.0, -10.0), 10.0, 40.0)
        self.assertGreater(s, 5.0)
        self.assertLess(s, 20.0)

    def test_snr_win_exact(self):
        x = np.zeros(N)
        x[1000:1200] = 1.0     # 50-60 s: noise rms 1
        x[2000] = -4.0         # 100 s: peak 4
        self.assertAlmostEqual(ks32.snr_win(x, SR, TB, 50.0, 60.0, 95.0, 105.0), 4.0)


if __name__ == "__main__":
    unittest.main()

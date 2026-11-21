# Copyright 2026 MV2, Munim, Inc.
# Licensed under the Apache License, Version 2.0 (Apache-2.0). See LICENSE.
# Made with AI agents working in MV2's name.
"""Unit tests for sdvwidth25.py (stdlib + numpy only; no network)."""

import math
import unittest

import numpy as np

import sdvwidth25 as W

SR = 200.0


def gabor(t, tc, sig, f0):
    u = np.asarray(t, float) - tc
    return np.exp(-0.5 * (u / sig) ** 2) * np.cos(2 * np.pi * f0 * u)


class TestOnset(unittest.TestCase):
    def setUp(self):
        self.t = np.arange(0.0, 20.0, 1.0 / SR)
        self.x = gabor(self.t, 10.0, 0.5, 3.0)
        self.env = W.envelope(self.x)

    def test_envelope_peak_near_one(self):
        self.assertAlmostEqual(float(np.max(self.env)), 1.0, delta=0.02)

    def test_onset_gaussian(self):
        # exp(-u^2 / (2 sig^2)) = 0.2 at u = sig * sqrt(2 ln 5) = 0.897 s before the peak
        on, ip = W.onset_index(self.env, int(8 * SR), int(12 * SR), 0.2, int(3 * SR))
        self.assertTrue(abs(ip - 2000) <= 2)
        self.assertAlmostEqual(on / SR, 10.0 - 0.5 * math.sqrt(2 * math.log(5)), delta=0.02)

    def test_lower_fraction_is_earlier(self):
        o1, _ = W.onset_index(self.env, int(8 * SR), int(12 * SR), 0.1, int(3 * SR))
        o3, _ = W.onset_index(self.env, int(8 * SR), int(12 * SR), 0.3, int(3 * SR))
        self.assertLess(o1, o3)

    def test_onset_not_found_when_scan_too_short(self):
        on, _ = W.onset_index(self.env, int(8 * SR), int(12 * SR), 0.2, int(0.5 * SR))
        self.assertTrue(math.isnan(on))


class TestStretch(unittest.TestCase):
    def setUp(self):
        self.t = np.arange(0.0, 20.0, 1.0 / SR)
        self.a = gabor(self.t, 10.0, 0.3, 3.0)
        self.i1, self.i2, self.ic = int(8.5 * SR), int(11.5 * SR), int(10.0 * SR)

    def test_identical_pulses(self):
        f = W.stretch_fit(self.a, self.a.copy(), self.i1, self.i2, self.ic, 0.7, 1.4, 0.005, 60)
        self.assertLessEqual(abs(f["s"] - 1.0), 0.006)
        self.assertGreater(f["cc"], 0.999)
        self.assertGreater(f["cc_s1"], 0.999)

    def test_broader_and_shifted_pulse(self):
        # b(t) = a((t - 10 - 0.06) / 1.2 + 10): 20% broader, 60 ms later
        b = gabor((self.t - 10.0 - 0.06) / 1.2 + 10.0, 10.0, 0.3, 3.0)
        f = W.stretch_fit(self.a, b, self.i1, self.i2, self.ic, 0.7, 1.4, 0.005, 60)
        self.assertLessEqual(abs(f["s"] - 1.2), 0.011)
        self.assertGreater(f["cc"], 0.995)
        self.assertLess(f["cc_s1"], f["cc"])

    def test_narrower_pulse(self):
        b = gabor((self.t - 10.0) / 0.9 + 10.0, 10.0, 0.3, 3.0)
        f = W.stretch_fit(self.a, b, self.i1, self.i2, self.ic, 0.7, 1.4, 0.005, 60)
        self.assertLessEqual(abs(f["s"] - 0.9), 0.011)

    def test_window_too_close_to_edge(self):
        f = W.stretch_fit(self.a, self.a, 10, 400, 200, 0.7, 1.4, 0.005, 60)
        self.assertTrue(math.isnan(f["s"]))


class TestLagAndBoot(unittest.TestCase):
    def test_window_lag_sign(self):
        t = np.arange(0.0, 20.0, 1.0 / SR)
        a = gabor(t, 10.0, 0.3, 3.0)
        b = gabor(t, 10.05, 0.3, 3.0)   # 10 samples later
        k, cc, edge = W.window_lag(a, b, int(9 * SR), int(11 * SR), 60)
        self.assertAlmostEqual(k, 10.0, delta=0.2)
        self.assertGreater(cc, 0.999)
        self.assertFalse(edge)

    def test_boot_median_brackets(self):
        med, lo, hi = W.boot_median([1.0, 1.02, 0.98, 1.05, 1.01, float("nan")], 500, 1)
        self.assertAlmostEqual(med, 1.01, places=9)
        self.assertLessEqual(lo, med)
        self.assertGreaterEqual(hi, med)

    def test_boot_median_empty(self):
        self.assertTrue(math.isnan(W.boot_median([], 10, 1)[0]))


class TestGrades(unittest.TestCase):
    def test_onset(self):
        self.assertEqual(W.grade_onset(-0.16, True), "onset-advance")
        self.assertEqual(W.grade_onset(0.03, True), "no-onset-advance")
        self.assertEqual(W.grade_onset(-0.07, True), "partial")
        self.assertEqual(W.grade_onset(-0.2, False), "untestable")
        self.assertEqual(W.grade_onset(float("nan"), True), "untestable")

    def test_width(self):
        self.assertEqual(W.grade_width(1.10, True), "broader")
        self.assertEqual(W.grade_width(0.90, True), "narrower")
        self.assertEqual(W.grade_width(1.01, True), "same")
        self.assertEqual(W.grade_width(1.04, True), "slight")
        self.assertEqual(W.grade_width(1.10, False), "untestable")

    def test_source(self):
        self.assertEqual(W.grade_source("same", 1.0, 1.0, 1.0), "not-needed")
        self.assertEqual(W.grade_source("broader", 1.10, 1.09, 1.08), "source")
        self.assertEqual(W.grade_source("broader", 1.10, 1.00, 1.01), "path-specific")
        self.assertEqual(W.grade_source("broader", 1.10, 1.06, 1.02), "unclear")
        self.assertEqual(W.grade_source("broader", 1.10, float("nan"), 1.0), "untestable")
        self.assertEqual(W.grade_source("untestable", 1.10, 1.0, 1.0), "untestable")

    def test_headline(self):
        self.assertEqual(W.headline("onset-advance", "source"), "onset-advance-from-source-width")
        self.assertEqual(W.headline("onset-advance", "path-specific"),
                         "onset-advance-path-specific")
        self.assertEqual(W.headline("onset-advance", "not-needed"), "onset-advance-unexplained")
        self.assertEqual(W.headline("no-onset-advance", "source"), "no-advance-in-either-reading")
        self.assertEqual(W.headline("partial", "unclear"), "partial-onset-advance")
        self.assertEqual(W.headline("untestable", "source"), "untestable")

    def test_pinned_and_paths(self):
        self.assertEqual(len(W.SDV23_COMMIT), 40)
        self.assertTrue(W.PLAN_PATH.startswith("findings/") and "prereg25" in W.PLAN_PATH)
        self.assertEqual(W.PULSES[0], "SKS_R")


if __name__ == "__main__":
    unittest.main()

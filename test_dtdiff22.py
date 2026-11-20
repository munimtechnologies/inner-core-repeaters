# Copyright 2026 MV2, Munim, Inc.
# Licensed under the Apache License, Version 2.0 (Apache-2.0). See LICENSE.
# Made with AI agents working in MV2's name.
"""Unit tests for dtdiff22.py (standard library + numpy only)."""

import unittest

import numpy as np

import dtdiff22 as d

FS = 20.0
N = 3000
C0 = 1200  # ak135 PKIKP at 60 s


def gabor(t0, amp, f=1.5, w=0.3):
    t = np.arange(N) / FS
    return amp * np.exp(-((t - t0) / w) ** 2) * np.cos(2 * np.pi * f * (t - t0))


def noise(seed):
    return np.random.default_rng(seed).normal(0.0, 1.0, N)


class TestClocks(unittest.TestCase):
    def test_backtracking_clock(self):
        self.assertAlmostEqual(d.g_B(2008.5), 0.0)
        self.assertAlmostEqual(d.g_B(2003.5), -5.0)
        self.assertAlmostEqual(d.g_B(2018.5), -4.0)

    def test_pause_and_linear(self):
        self.assertAlmostEqual(d.g_P(2015.0), 0.0)
        self.assertAlmostEqual(d.g_P(2004.0), -5.0)
        self.assertAlmostEqual(d.g_L(2010.5), 2.0)

    def test_fit_origin(self):
        a, r2 = d.fit_origin([1.0, 2.0, 3.0], [2.0, 4.0, 6.0])
        self.assertAlmostEqual(a, 2.0)
        self.assertAlmostEqual(r2, 1.0)

    def test_pct(self):
        self.assertAlmostEqual(d.pct([4, 1, 3, 2], 0.5), 2.5)
        self.assertAlmostEqual(d.pct([1, 2, 3, 4, 5], 0.0), 1.0)


def synthetic_rows(fn):
    rows, i = [], 0
    for f in range(12):
        t1 = 2001.0 + 0.5 * f
        for t2, kind in ((2012.0 + 0.7 * f, "straddle"), (2016.0 + 0.5 * f, "straddle"), (t1 + 1.5, "before")):
            r = dict(t1=t1, t2=t2, kind=kind, mult="F%02d" % f)
            r["ddt_s"] = fn(r, i)
            rows.append(r)
            i += 1
    return rows


class TestGrade(unittest.TestCase):
    def test_backtracking_signal_is_B(self):
        rows = synthetic_rows(lambda r, i: 0.01 * d.dg(d.g_B, r) + 0.002 * (-1) ** i)
        g = d.grade(rows, nboot=200)
        self.assertEqual(g["decision"], "B")
        self.assertGreater(g["models"]["B"]["a"], 0.009)
        self.assertLess(g["models"]["B"]["a"], 0.011)

    def test_zero_signal_is_none(self):
        g = d.grade(synthetic_rows(lambda r, i: 0.0), nboot=200)
        self.assertEqual(g["decision"], "none")

    def test_tiny_noise_slopes_small(self):
        g = d.grade(synthetic_rows(lambda r, i: 0.0005 * (-1) ** i), nboot=200)
        self.assertNotEqual(g["decision"], "untestable")
        for m in g["models"].values():
            self.assertLess(abs(m["a"]), d.NONE_RATE)

    def test_few_pairs_untestable(self):
        rows = synthetic_rows(lambda r, i: 0.01 * d.dg(d.g_B, r))[:5]
        self.assertEqual(d.grade(rows, nboot=50)["decision"], "untestable")


class TestSelect(unittest.TestCase):
    def test_select_pairs(self):
        rows = [
            dict(ev1="a", ev2="b", t1="2001.0", t2="2012.0", kind="straddle", mult="F1", S="0.3", n_st="2",
                 mkar_ok="True", dB="1", dM="11"),
            dict(ev1="c", ev2="d", t1="2015.0", t2="2004.0", kind="straddle", mult="F2", S="0.4", n_st="5",
                 mkar_ok="True", dB="1", dM="11"),
            dict(ev1="e", ev2="f", t1="2001.0", t2="2002.0", kind="before", mult="F3", S="", n_st="0",
                 mkar_ok="True", dB="1", dM="1"),
            dict(ev1="g", ev2="h", t1="2001.0", t2="2012.0", kind="straddle", mult="F4", S="0.1", n_st="4",
                 mkar_ok="False", dB="1", dM="11"),
        ]
        out = d.select_pairs(rows)
        self.assertEqual(len(out), 3)
        self.assertTrue(out[0]["verified"])
        self.assertFalse(out[1]["verified"])
        self.assertFalse(out[2]["verified"])
        self.assertEqual((out[1]["ev1"], out[1]["ev2"]), ("d", "c"))
        self.assertLess(out[1]["t1"], out[1]["t2"])


class TestMeasure(unittest.TestCase):
    def test_recovers_df_shift(self):
        dl = 2.0
        xa = noise(1) + gabor(60.0, 20.0) + gabor(60.0 + dl, 40.0)
        xb = noise(2) + gabor(60.8, 20.0) + gabor(60.7 + dl, 40.0)
        m = d.measure_pair(xa, xb, FS, C0, dl)
        self.assertEqual(m["status"], "ok")
        self.assertEqual(m["onset_class"], "df-onset")
        self.assertLess(abs(m["lag_late_s"] - 0.7), 0.05)
        self.assertLess(abs(m["ddt_s"] - 0.1), 0.03)

    def test_missing_df(self):
        dl = 2.0
        xa = noise(3) + gabor(60.0 + dl, 40.0)
        xb = noise(4) + gabor(60.5 + dl, 40.0)
        self.assertEqual(d.measure_pair(xa, xb, FS, C0, dl)["status"], "df-not-onset")
        self.assertNotEqual(d.measure_pair(xa, xb, FS, C0, dl, late_onset_ok=True)["status"], "ok")

    def test_guards(self):
        x = noise(5)
        self.assertEqual(d.measure_pair(None, x, FS, C0, 2.0)["status"], "no-beam")
        self.assertEqual(d.measure_pair(x, x, FS, C0, 0.5)["status"], "dl-out-of-range")
        self.assertEqual(d.measure_pair(x, x, FS, C0, float("nan"))["status"], "dl-out-of-range")


if __name__ == "__main__":
    unittest.main()

# Copyright 2026 MV2, Munim, Inc.
# Licensed under the Apache License, Version 2.0 (Apache-2.0). See LICENSE.
# Made with AI agents working in MV2's name.
"""Unit tests for sdv23.py (stdlib + numpy only)."""

import math
import unittest

import numpy as np

import sdv23 as m


def gauss(n, c, w):
    i = np.arange(n, dtype=float)
    return np.exp(-((i - c) / w) ** 2)


def wavelet(n, c, w, f):
    i = np.arange(n, dtype=float)
    return np.cos(2 * math.pi * f * (i - c)) * np.exp(-((i - c) / w) ** 2)


def rt_to_ne(r, t, baz):
    b = math.radians(baz)
    n = -r * math.cos(b) + t * math.sin(b)
    e = -r * math.sin(b) - t * math.cos(b)
    return n, e


def ne_to_channels(n, e, az1, az2):
    a1, a2 = math.radians(az1), math.radians(az2)
    return (n * math.cos(a1) + e * math.sin(a1), n * math.cos(a2) + e * math.sin(a2))


class TestAngles(unittest.TestCase):
    def test_wrap180(self):
        self.assertAlmostEqual(float(m.wrap180(357.0)), -3.0, places=9)
        self.assertAlmostEqual(float(m.wrap180(-190.0)), 170.0, places=9)
        self.assertAlmostEqual(float(m.wrap180(343.0)), -17.0, places=9)

    def test_sdv_location(self):
        self.assertEqual(m.sdv_location("1996-01-20T07:16:37.84"), "--")
        self.assertEqual(m.sdv_location("1998-11-08T23:48:47.86"), "--")
        self.assertEqual(m.sdv_location("2000-07-07T15:46:44.56"), "00")
        self.assertEqual(m.sdv_location("2009-12-17T20:01:21.61"), "00")


class TestSignal(unittest.TestCase):
    def test_upsample_keeps_samples(self):
        n = 64
        x = np.sin(2 * math.pi * 3 * np.arange(n) / n)
        y = m.upsample(x, 5)
        self.assertEqual(len(y), 5 * n)
        self.assertTrue(np.allclose(y[::5], x, atol=1e-9))
        self.assertAlmostEqual(y[2], math.sin(2 * math.pi * 3 * (2 / 5) / n), places=6)

    def test_hilbert_cos_to_sin(self):
        n = 128
        i = np.arange(n)
        x = np.cos(2 * math.pi * 5 * i / n)
        self.assertTrue(np.allclose(m.hilbert(x), np.sin(2 * math.pi * 5 * i / n), atol=1e-9))

    def test_window_lag_positive_and_negative(self):
        a = gauss(400, 100.0, 8.0)
        for s in (7.3, -4.6):
            b = gauss(400, 100.0 + s, 8.0)
            lag, cc, edge = m.window_lag(a, b, 60, 140, 20)
            self.assertFalse(edge)
            self.assertAlmostEqual(lag, s, delta=0.15)
            self.assertGreater(cc, 0.99)

    def test_window_lag_edge(self):
        a = gauss(400, 100.0, 8.0)
        b = gauss(400, 125.0, 8.0)
        lag, cc, edge = m.window_lag(a, b, 60, 140, 20)
        self.assertTrue(edge)

    def test_window_lag_out_of_range(self):
        a = gauss(100, 50.0, 5.0)
        lag, cc, edge = m.window_lag(a, a, 5, 60, 20)
        self.assertTrue(edge)
        self.assertTrue(math.isnan(lag))

    def test_snr(self):
        x = np.zeros(200)
        x[:100] = np.tile([1.0, -1.0], 50)
        x[150] = 10.0
        self.assertAlmostEqual(m.snr(x, 0, 100, 140, 160), 10.0, places=9)


class TestRotation(unittest.TestCase):
    def test_rotate_to_ne_roundtrip(self):
        rng = np.random.default_rng(3)
        n, e = rng.normal(size=50), rng.normal(size=50)
        for az1, az2 in ((0.0, 90.0), (357.0, 87.0), (343.0, 73.0), (357.0, -273.0)):
            c1, c2 = ne_to_channels(n, e, az1, az2)
            n2, e2 = m.rotate_to_ne(c1, c2, az1, az2)
            self.assertTrue(np.allclose(n2, n, atol=1e-9))
            self.assertTrue(np.allclose(e2, e, atol=1e-9))

    def test_parallel_channels_refused(self):
        with self.assertRaises(ValueError):
            m.rotate_to_ne([1.0], [1.0], 10.0, 190.0)

    def test_ne_to_rt_convention(self):
        n, e = np.array([1.0]), np.array([2.0])
        r, t = m.ne_to_rt(n, e, 0.0)        # source due north: radial points south
        self.assertAlmostEqual(float(r[0]), -1.0, places=9)
        self.assertAlmostEqual(float(t[0]), -2.0, places=9)
        r, t = m.ne_to_rt(n, e, 90.0)       # source due east
        self.assertAlmostEqual(float(r[0]), -2.0, places=9)
        self.assertAlmostEqual(float(t[0]), 1.0, places=9)

    def test_rt_inverse_helper(self):
        rng = np.random.default_rng(5)
        r, t = rng.normal(size=20), rng.normal(size=20)
        n, e = rt_to_ne(r, t, 47.0)
        r2, t2 = m.ne_to_rt(n, e, 47.0)
        self.assertTrue(np.allclose(r2, r) and np.allclose(t2, t))

    def test_orient_rz_recovers_misorientation(self):
        nn = 2000
        r = wavelet(nn, 1000.0, 200.0, 0.02)
        t = 0.8 * wavelet(nn, 1300.0, 150.0, 0.035)
        z = m.hilbert(r)
        baz = 40.0
        n, e = rt_to_ne(r, t, baz)
        c1, c2 = ne_to_channels(n, e, 355.0, 85.0)
        a_hat, c = m.orient_rz(c1, c2, z, 90.0, baz, 0.0, 45.0, 0.1)
        self.assertLess(abs(float(m.wrap180(a_hat - 355.0))), 1.0)
        self.assertGreater(c, 0.9)

    def test_orient_pair_recovers_rotation(self):
        nn = 1500
        n = wavelet(nn, 600.0, 120.0, 0.03)
        e = 0.6 * wavelet(nn, 800.0, 150.0, 0.02)
        a1, a2 = ne_to_channels(n, e, 0.0, 90.0)
        b1, b2 = ne_to_channels(0.7 * n, 0.7 * e, 357.0, 87.0)
        for daz_b in (90.0, -270.0):
            d, cc = m.orient_pair(a1, a2, b1, b2, 90.0, daz_b, 30.0, 0.1)
            self.assertAlmostEqual(d, -3.0, delta=0.15)
            self.assertGreater(cc, 0.99)


class TestLocation(unittest.TestCase):
    def test_relocate_exact(self):
        az = np.arange(0.0, 360.0, 60.0)
        p = np.full(len(az), 0.06)
        dt0, dxn, dxe = 0.02, 0.5, -0.3
        lags = dt0 - p * (dxn * np.cos(np.radians(az)) + dxe * np.sin(np.radians(az)))
        r = m.relocate(lags, p, az)
        self.assertAlmostEqual(r["dt0"], dt0, places=6)
        self.assertAlmostEqual(r["dxN"], dxn, places=6)
        self.assertAlmostEqual(r["dxE"], dxe, places=6)
        self.assertLess(r["rms"], 1e-9)
        self.assertEqual(r["n"], 6)

    def test_relocate_needs_four(self):
        self.assertIsNone(m.relocate([0.1, 0.2, 0.3], [0.06] * 3, [0.0, 90.0, 180.0]))

    def test_loc_pred_sign(self):
        a = 30.0
        dxn, dxe = math.cos(math.radians(a)), math.sin(math.radians(a))
        self.assertAlmostEqual(m.loc_pred(dxn, dxe, a, 0.02), -0.02, places=9)
        self.assertAlmostEqual(m.loc_pred(dxn, dxe, a + 180.0, 0.02), 0.02, places=9)


class TestGrades(unittest.TestCase):
    G = dict(a=True, b=True)

    def test_primary(self):
        g = m.grade_primary
        self.assertEqual(g(self.G, -0.16, -0.16, -0.15), "change-survives")
        self.assertEqual(g(self.G, -0.16, -0.01, 0.0), "orientation-explains")
        self.assertEqual(g(self.G, -0.05, -0.05, -0.05), "not-reproduced")
        self.assertEqual(g(self.G, 0.16, 0.16, 0.16), "not-reproduced")
        self.assertEqual(g(self.G, -0.16, -0.16, -0.07), "partly-orientation")
        self.assertEqual(g(dict(a=True, b=False), -0.16, -0.16, -0.15), "untestable")
        self.assertEqual(g({}, -0.16, -0.16, -0.15), "untestable")
        self.assertEqual(g(self.G, -0.16, float("nan"), -0.15), "untestable")

    def test_history(self):
        self.assertEqual(m.grade_history(0.01, True), "stable")
        self.assertEqual(m.grade_history(-0.08, True), "changed")
        self.assertEqual(m.grade_history(0.045, True), "unclear")
        self.assertEqual(m.grade_history(0.01, False), "untestable")
        self.assertEqual(m.grade_history(float("nan"), True), "untestable")

    def test_rotation_physical(self):
        before = np.array([-0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3]) % 360.0
        after = (-3.0 + np.array([-0.2, -0.1, 0.0, 0.1, 0.2])) % 360.0
        r = m.grade_rotation(before, after, -3.0)
        self.assertEqual(r["grade"], "physical-rotation")
        self.assertAlmostEqual(r["diff"], -3.0, places=6)

    def test_rotation_none(self):
        before = np.array([-0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3]) % 360.0
        after = np.array([-0.2, -0.1, 0.0, 0.1, 0.2]) % 360.0
        self.assertEqual(m.grade_rotation(before, after, -3.0)["grade"], "no-physical-rotation")

    def test_rotation_too_few(self):
        self.assertEqual(m.grade_rotation([0.0] * 4, [0.0] * 5, -3.0)["grade"], "untestable")


class TestConstants(unittest.TestCase):
    def test_events_and_params(self):
        ids = [e[0] for e in m.EVENTS]
        self.assertIn(m.MASTER, ids)
        self.assertIn(m.PARTNER, ids)
        self.assertEqual(len(ids), len(set(ids)))
        self.assertGreater(m.PARAMS["reproduce_s"], m.PARAMS["explain_s"])
        self.assertEqual(m.META_STEP["00"], -3.0)

    def test_jsonable(self):
        self.assertIs(m.jsonable(np.bool_(True)), True)
        self.assertEqual(m.jsonable(np.int64(3)), 3)
        self.assertEqual(m.jsonable(np.arange(2)), [0, 1])


if __name__ == "__main__":
    unittest.main()

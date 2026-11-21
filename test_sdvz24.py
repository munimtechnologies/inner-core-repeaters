# Copyright 2026 MV2, Munim, Inc.
# Licensed under the Apache License, Version 2.0 (Apache-2.0). See LICENSE.
# Made with AI agents working in MV2's name.
"""Unit tests for sdvz24.py (stdlib + numpy only)."""

import math
import os
import tempfile
import unittest
import zipfile

import numpy as np

import sdv23 as base
import sdvz24 as m


def wavelet(n, c, w, f):
    i = np.arange(n, dtype=float)
    return np.cos(2 * math.pi * f * (i - c)) * np.exp(-((i - c) / w) ** 2)


def ne_to_channels(n, e, az1, az2):
    a1, a2 = math.radians(az1), math.radians(az2)
    return (n * math.cos(a1) + e * math.sin(a1), n * math.cos(a2) + e * math.sin(a2))


def field(nn, seed):
    rng = np.random.default_rng(seed)
    n = np.zeros(nn)
    e = np.zeros(nn)
    for _ in range(12):
        c = rng.uniform(50, nn - 50)
        n += rng.normal() * wavelet(nn, c, 8.0, rng.uniform(0.02, 0.06))
        e += rng.normal() * wavelet(nn, c + rng.uniform(-20, 20), 8.0, rng.uniform(0.02, 0.06))
    return n, e


class TestReframe(unittest.TestCase):
    def test_reframe_matches_to_rt(self):
        n, e = field(400, 1)
        for a, da, baz in ((0.0, -3.0, 312.4), (357.0, 3.0, 40.0), (10.0, 17.0, 180.0)):
            c1, c2 = ne_to_channels(n, e, 0.0, 90.0)
            r, t = base.to_rt(c1, c2, a, a + 90.0, baz)
            r2, t2 = m.reframe(r, t, da)
            r3, t3 = base.to_rt(c1, c2, a + da, a + da + 90.0, baz)
            self.assertTrue(np.allclose(r2, r3, atol=1e-9))
            self.assertTrue(np.allclose(t2, t3, atol=1e-9))

    def test_reframe_zero_and_inverse(self):
        rng = np.random.default_rng(2)
        r, t = rng.normal(size=30), rng.normal(size=30)
        r0, t0 = m.reframe(r, t, 0.0)
        self.assertTrue(np.allclose(r0, r) and np.allclose(t0, t))
        r1, t1 = m.reframe(*m.reframe(r, t, 7.5), -7.5)
        self.assertTrue(np.allclose(r1, r) and np.allclose(t1, t))


class TestFitFrame(unittest.TestCase):
    def setUp(self):
        self.nn, self.mlag = 900, 20
        self.N, self.E = field(self.nn, 7)
        self.baz = 312.4

    def theirs(self, assumed_vs_true, k0, reflect=False):
        """Their R/T made with a frame rotated by assumed_vs_true from the geographic
        N/E, cut so that ours[m + k0 + i] matches theirs[i]."""
        n = self.nn - 2 * self.mlag
        r, t = base.to_rt(self.N, self.E, assumed_vs_true, assumed_vs_true + 90.0, self.baz)
        sl = slice(self.mlag + k0, self.mlag + k0 + n)
        return r[sl], (-t[sl] if reflect else t[sl])

    def test_recovers_angle_and_lag(self):
        for a_true, k0 in ((0.0, 0), (3.0, 4), (-17.0, -6), (357.0, 11)):
            rp, tp = self.theirs(a_true, k0)
            fr = m.fit_frame(self.N, self.E, rp, tp, self.baz, self.mlag)
            self.assertLess(abs(float(base.wrap180(fr["a"] - a_true))), 0.01)
            self.assertEqual(fr["lag"], k0)
            self.assertGreater(fr["c"], 0.9999)
            self.assertEqual(fr["conv"], "proper")

    def test_reflected_convention(self):
        rp, tp = self.theirs(3.0, 2, reflect=True)
        fr = m.fit_frame(self.N, self.E, rp, tp, self.baz, self.mlag)
        self.assertEqual(fr["conv"], "reflected")
        self.assertGreater(fr["c"], 0.9999)
        self.assertEqual(fr["lag"], 2)

    def test_noise_keeps_angle_close(self):
        rng = np.random.default_rng(11)
        rp, tp = self.theirs(-3.0, 0)
        s = 0.05 * np.std(rp)
        fr = m.fit_frame(self.N, self.E, rp + s * rng.normal(size=len(rp)),
                         tp + s * rng.normal(size=len(tp)), self.baz, self.mlag)
        self.assertLess(abs(float(base.wrap180(fr["a"] + 3.0))), 0.5)
        self.assertGreater(fr["c"], 0.95)

    def test_bad_lengths(self):
        with self.assertRaises(ValueError):
            m.fit_frame(self.N[:-1], self.E, self.N[:100], self.E[:100], 0.0, self.mlag)


class TestGrid(unittest.TestCase):
    def test_on_grid_identity_and_nan(self):
        x = np.arange(100, dtype=float)
        y = m.on_grid(x, 10.0, 5.0, 6.0, 8.0, 10.0)      # samples 10..29
        self.assertEqual(len(y), 20)
        self.assertTrue(np.allclose(y, np.arange(10, 30)))
        y = m.on_grid(x, 10.0, 5.0, 4.0, 5.5, 10.0)
        self.assertTrue(np.isnan(y[0]) and np.isfinite(y[-1]))

    def test_on_grid_upsamples(self):
        t = np.arange(200) / 20.0
        x = np.sin(2 * math.pi * 0.5 * t)
        y = m.on_grid(x, 20.0, 0.0, 2.0, 6.0, 200.0)
        tg = 2.0 + np.arange(len(y)) / 200.0
        self.assertEqual(len(y), 800)
        self.assertLess(float(np.max(np.abs(y - np.sin(2 * math.pi * 0.5 * tg)))), 0.02)

    def test_covered(self):
        x = np.ones(100)
        x[90:] = np.nan
        self.assertTrue(m.covered(x, 0.0, 10.0, 1.0, 8.0))
        self.assertFalse(m.covered(x, 0.0, 10.0, 1.0, 9.5))
        self.assertFalse(m.covered(x, 0.0, 10.0, -1.0, 5.0))


class TestGrades(unittest.TestCase):
    def test_classify_frame(self):
        self.assertEqual(m.classify_frame(0.4), "same-frame")
        self.assertEqual(m.classify_frame(-2.6), "metadata-frame")
        self.assertEqual(m.classify_frame(357.5), "metadata-frame")
        self.assertEqual(m.classify_frame(-17.0), "other")
        self.assertEqual(m.classify_frame(float("nan")), "untestable")

    def test_combined(self):
        self.assertEqual(m.grade_combined(-3.0, -2.6)[0], "consistent")
        g, eps = m.grade_combined(0.0, -3.0)
        self.assertEqual(g, "misrotated")
        self.assertAlmostEqual(eps, 3.0, places=9)
        self.assertEqual(m.grade_combined(0.0, -1.5)[0], "unclear")
        self.assertEqual(m.grade_combined(0.0, float("nan"))[0], "unresolved")

    def test_reproduce(self):
        G = dict(a=True)
        self.assertEqual(m.grade_reproduce(G, -0.16), "reproduced")
        self.assertEqual(m.grade_reproduce(G, -0.05), "not-reproduced")
        self.assertEqual(m.grade_reproduce(G, 0.2), "not-reproduced")
        self.assertEqual(m.grade_reproduce(dict(a=False), -0.16), "untestable")
        self.assertEqual(m.grade_reproduce(G, float("nan")), "untestable")


class TestMeasure(unittest.TestCase):
    def make(self, sr, shift_all, sks_extra, n_s=100.0):
        nn = int(n_s * sr)
        w = 0.25 * sr
        t_sks, t_s = 40.0, 80.0
        r = (wavelet(nn, (t_sks + shift_all + sks_extra) * sr, w, 0.8 / sr)
             + 0.6 * wavelet(nn, (t_s + shift_all) * sr, w, 0.7 / sr))
        t = wavelet(nn, (t_s + 1.0 + shift_all) * sr, w, 0.6 / sr)
        return r, t, t_sks, t_s

    def test_windows_and_x(self):
        sr, g0 = 50.0, 0.0
        ra, ta, t_sks, t_s = self.make(sr, 0.0, 0.0)
        rb, tb, _, _ = self.make(sr, 1.3, -0.16)
        wins = m.pick_windows(ra, ta, sr, g0, t_sks, t_s)
        self.assertAlmostEqual(wins["SKS_R"][1], t_sks - 2.0, delta=0.05)
        self.assertAlmostEqual(wins["S_T"][1], t_s + 1.0 - 3.0, delta=0.05)
        res, shift = m.measure(ra, ta, rb, tb, sr, g0, wins)
        self.assertEqual(shift, 65)
        self.assertAlmostEqual(res["S_T"]["lag"], 1.3, delta=0.02)
        self.assertAlmostEqual(res["SKS_R"]["lag"] - res["S_T"]["lag"], -0.16, delta=0.02)
        self.assertFalse(res["SKS_R"]["edge"])

    def test_rotation_leaks_transverse(self):
        """A pure rotation of b changes the radial S window more than SKS when T is only
        at S time (sanity of the sensitivity loop)."""
        sr, g0 = 50.0, 0.0
        ra, ta, t_sks, t_s = self.make(sr, 0.0, 0.0)
        wins = m.pick_windows(ra, ta, sr, g0, t_sks, t_s)
        r0, _ = m.measure(ra, ta, ra, ta, sr, g0, wins, shift=0)
        rb, tb = m.reframe(ra, ta, 10.0)
        r1, _ = m.measure(ra, ta, rb, tb, sr, g0, wins, shift=0)
        self.assertAlmostEqual(r0["SKS_R"]["lag"], 0.0, delta=1e-3)
        self.assertAlmostEqual(r1["SKS_R"]["lag"], 0.0, delta=0.01)
        self.assertLess(r1["S_R_early"]["cc"], r0["S_R_early"]["cc"])

    def test_needed_range(self):
        wins = {"a": ("R", 10.0, 15.0), "b": ("T", 20.0, 33.0)}
        lo, hi = m.needed_range(wins)
        pad = m.PARAMS["coarse_maxlag_s"] + m.PARAMS["maxlag_s"] + 0.5
        self.assertAlmostEqual(lo, 10.0 - pad, places=9)
        self.assertAlmostEqual(hi, 33.0 + pad, places=9)


class TestFiles(unittest.TestCase):
    def test_set_key_and_comp(self):
        k1 = m.set_key("d/2009.351.IU.SDV.00.BHR.SAC", "BHR", 20.0, "2009-12-17T20:00:00")
        k2 = m.set_key("d/2009.351.IU.SDV.00.BHT.SAC", "BHT", 20.0, "2009-12-17T20:00:00")
        self.assertEqual(k1, k2)
        self.assertEqual(m.comp_of("BHT", "x"), "T")
        self.assertEqual(m.comp_of("", "a/AI-0009.2000.SDV.r"), "R")
        self.assertEqual(m.comp_of("", "a/README"), "")

    def test_safe_member(self):
        self.assertEqual(m.safe_member("a/b.sac"), os.path.normpath("a/b.sac"))
        self.assertIsNone(m.safe_member("../x"))
        self.assertIsNone(m.safe_member("a/../../x"))
        self.assertIsNone(m.safe_member("/etc/passwd"))

    def test_extract_zip_nested(self):
        with tempfile.TemporaryDirectory() as d:
            inner = os.path.join(d, "inner.zip")
            with zipfile.ZipFile(inner, "w") as z:
                z.writestr("x/a.sac", b"abc")
            outer = os.path.join(d, "outer.zip")
            with zipfile.ZipFile(outer, "w") as z:
                z.write(inner, "data/inner.zip")
                z.writestr("b.txt", b"hello")
                z.writestr("../evil.txt", b"no")
            files = m.extract_zip(outer, os.path.join(d, "out"), 10)
            names = sorted(os.path.basename(f) for f in files)
            self.assertEqual(names, ["a.sac", "b.txt"])
            self.assertFalse(os.path.exists(os.path.join(d, "evil.txt")))

    def test_md5(self):
        self.assertEqual(m.md5(b""), "d41d8cd98f00b204e9800998ecf8427e")


class TestConstants(unittest.TestCase):
    def test_constants(self):
        self.assertEqual(set(m.EVENTS), {m.MASTER, m.PARTNER})
        ev = {r[0]: r for r in base.EVENTS}
        for k, v in m.EVENTS.items():
            self.assertEqual(ev[k][1], v[0])
        self.assertEqual(len(m.ZENODO_MD5), 32)
        self.assertGreater(m.PARAMS["misrot_min"], m.PARAMS["consistent_tol"])
        self.assertEqual(m.PARAMS["meta_step"], base.META_STEP["00"])


if __name__ == "__main__":
    unittest.main()

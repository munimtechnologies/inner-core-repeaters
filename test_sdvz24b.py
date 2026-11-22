# Copyright 2026 MV2, Munim, Inc.
# Licensed under the Apache License, Version 2.0 (Apache-2.0). See LICENSE.
# Made with AI agents working in MV2's name.
"""Unit tests for sdvz24b.py (stdlib + numpy only; no network)."""

import math
import unittest

import numpy as np

import sdvz24b as m


def wavelet(n, c, w, f):
    i = np.arange(n, dtype=float)
    return np.cos(2 * math.pi * f * (i - c)) * np.exp(-((i - c) / w) ** 2)


class TestKeys(unittest.TestCase):
    def test_authors_file_names_split(self):
        k = lambda p, sr: m.set_key_v2(p, "", sr, "2009-12-17T20:11:43.948013Z")  # noqa: E731
        self.assertEqual(k("data_IU.SDV/2009.R.sac", 100.0), "data_IU.SDV/2009.?.sac")
        self.assertEqual(k("data_IU.SDV/2009.T.sac", 100.0), "data_IU.SDV/2009.?.sac")
        self.assertEqual(k("data_IU.SDV/2009.R.sac.bp", 100.0), "data_IU.SDV/2009.?.sac.bp")
        self.assertNotEqual(k("data_IU.SDV/2009.R.sac", 100.0),
                            k("data_IU.SDV/2009.R.sac.bp", 100.0))
        self.assertEqual(k("data_IU.SDV/2000.Z.sac.bp", 1000.0), "data_IU.SDV/2000.?.sac.bp")

    def test_channel_code_and_fallback(self):
        k1 = m.set_key_v2("d/2009.351.IU.SDV.00.BHR.SAC", "BHR", 20.0, "2009-12-17T20:00:00")
        k2 = m.set_key_v2("d/2009.351.IU.SDV.00.BHT.SAC", "BHT", 20.0, "2009-12-17T20:00:00")
        self.assertEqual(k1, k2)
        self.assertEqual(m.set_key_v2("d/README", "", 1.0, "2009-12-17T20:00:00.000"),
                         "d|1.000|2009-12-17T20:00:00.0")

    def test_suffix_and_bp_first(self):
        self.assertEqual(m.set_suffix("data_IU.SDV/2000.?.sac.bp"), ".sac.bp")
        self.assertEqual(m.set_suffix("data_IU.SDV/2000.?.sac"), ".sac")
        self.assertEqual(m.set_suffix("data_IU.SDV|100.000|x"), "")
        a, b = "data_IU.SDV/2000.?.sac", "data_IU.SDV/2000.?.sac.bp"
        self.assertEqual(sorted([a, b])[0], a)
        self.assertEqual(sorted([m.bp_first(a), m.bp_first(b)])[0], m.bp_first(b))
        self.assertEqual(m.set_suffix(m.bp_first(b)), ".sac.bp")
        self.assertEqual(m.bp_first(a), a)


class TestGrades(unittest.TestCase):
    def test_combined(self):
        self.assertEqual(m.grade_combined(0.4, 0.55)[0], "consistent")
        g, eps = m.grade_combined(-2.23, 0.55)
        self.assertEqual(g, "misrotated")
        self.assertAlmostEqual(eps, -2.78, places=6)
        self.assertEqual(m.grade_combined(-1.0, 0.55)[0], "unclear")
        self.assertEqual(m.grade_combined(177.85, 0.55)[0], "misrotated")
        self.assertEqual(m.grade_combined(float("nan"), 0.55)[0], "unresolved")

    def test_wrap(self):
        self.assertAlmostEqual(m.wrap180(359.0), -1.0)
        self.assertAlmostEqual(m.wrap180(-181.0), 179.0)

    def test_grades_of_like(self):
        S = dict(primary_set={"E2000a": "d/2000.?.sac", "E2009a": "d/2009.?.sac"},
                 frame=dict(grade="metadata-frame", rho=-2.5, gates={}),
                 measurement=dict(grade="not-reproduced", x0=0.03, gates={}))
        g = m.grades_of(S, ".sac")
        self.assertTrue(g["like_with_like"])
        self.assertEqual(g["combined_p23"], "misrotated")
        S["primary_set"]["E2009a"] = "d/2009.?.sac.bp"
        self.assertFalse(m.grades_of(S, ".sac")["like_with_like"])


class TestXcorr(unittest.TestCase):
    def test_shift_and_negation(self):
        n = 800
        a = wavelet(n, 400.0, 25.0, 0.02)
        b = -wavelet(n, 413.3, 25.0, 0.02)
        r = m.xcorr_signed(a, b, 60)
        self.assertLess(r["cc_abs"], -0.99)
        self.assertAlmostEqual(r["lag_abs"], 13.3, delta=0.2)
        self.assertGreater(r["cc_pos"], 0.0)
        self.assertGreater(abs(r["lag_pos"] - 13.3), 5.0)

    def test_identity(self):
        a = wavelet(500, 250.0, 20.0, 0.03)
        r = m.xcorr_signed(a, a.copy(), 30)
        self.assertAlmostEqual(r["lag_abs"], 0.0, delta=1e-6)
        self.assertAlmostEqual(r["cc_abs"], 1.0, places=9)

    def test_bad_lengths(self):
        with self.assertRaises(ValueError):
            m.xcorr_signed(np.ones(10), np.ones(11), 2)


class TestConstants(unittest.TestCase):
    def test_constants(self):
        self.assertEqual(len(m.CODE_COMMIT), 40)
        self.assertTrue(m.P23_G5)
        self.assertAlmostEqual(m.P23_DELTA, 0.55, places=6)
        self.assertEqual(set(m.MUST), {"sdv23.py", "sdvz24.py"})


if __name__ == "__main__":
    unittest.main()

# Copyright 2026 MV2, Munim, Inc.
# Licensed under the Apache License, Version 2.0 (Apache-2.0). See LICENSE.
# Made with AI agents working in MV2's name.
"""Synthetic tests for sksnbr26.py (no network, no ObsPy)."""

import math
import unittest

import numpy as np

import sksnbr26 as K


def gabor(n, c, sr, f=1.0, w=0.25):
    i = np.arange(n, dtype=float)
    t = (i - c) / sr
    return np.exp(-(t / w) ** 2) * np.cos(2 * math.pi * f * t)


class TestGeometry(unittest.TestCase):
    def test_dest_point_axes(self):
        lat, lon = K.dest_point(0.0, 0.0, 0.0, 10.0)
        self.assertAlmostEqual(lat, 10.0, places=6)
        self.assertAlmostEqual(lon, 0.0, places=6)
        lat, lon = K.dest_point(0.0, 0.0, 90.0, 10.0)
        self.assertAlmostEqual(lat, 0.0, places=6)
        self.assertAlmostEqual(lon, 10.0, places=6)

    def test_dest_point_roundtrip(self):
        for az, d in ((37.0, 25.0), (250.0, 80.0), (69.3, 9.0)):
            lat, lon = K.dest_point(51.4, 179.98, az, d)
            self.assertAlmostEqual(K.ang_dist(51.4, 179.98, lat, lon), d, places=5)

    def test_dateline_wrap(self):
        lat, lon = K.dest_point(51.4, 179.98, 90.0, 1.0)
        self.assertTrue(-180.0 <= lon < 180.0)
        self.assertLess(lon, -178.0)

    def test_arc_km(self):
        self.assertAlmostEqual(K.arc_km(0, 0, 0, 1), math.radians(1) * 3480.0, places=6)
        self.assertAlmostEqual(K.arc_km(10, 20, 10, 20), 0.0, places=9)

    def test_small_azimuth_offset(self):
        # entry points 9 deg from source, 4 deg apart in azimuth -> roughly 9*sin(4) deg apart
        e1 = K.dest_point(51.4, 179.98, 69.3, 9.0)
        e2 = K.dest_point(51.4, 179.98, 65.3, 9.0)
        km = K.arc_km(e1[0], e1[1], e2[0], e2[1])
        self.assertGreater(km, 30.0)
        self.assertLess(km, 45.0)

    def test_cmb_indices(self):
        d = [31.0, 35.0, 410.0, 2889.0, 4000.0, 5153.5, 4000.0, 2889.0, 660.0, 0.0]
        self.assertEqual(K.cmb_indices(d, 2889.0), (3, 7))
        self.assertEqual(K.cmb_indices([0.0, 100.0], 2889.0), (None, None))


class TestSelection(unittest.TestCase):
    def test_select_loc(self):
        self.assertEqual(K.select_loc({"10": {"BHZ", "BH1", "BH2"}, "00": {"BHZ", "BHN"}}), "10")
        self.assertEqual(K.select_loc({"": {"BHZ", "BHN", "BHE"},
                                       "00": {"BHZ", "BHN", "BHE"}}), "00")
        self.assertEqual(K.select_loc({"": {"BHZ", "BHN", "BHE"}}), "")
        self.assertIsNone(K.select_loc({"00": {"BHN", "BHE"}}))

    def test_classify(self):
        self.assertEqual(K.classify(40.0, 250.0, "VE"), "entry-neighbour")
        self.assertEqual(K.classify(150.0, 80.0, "IU"), "exit-neighbour")
        self.assertEqual(K.classify(400.0, 900.0, "IU"), "control")
        self.assertEqual(K.classify(400.0, 900.0, "XY"), "other")
        self.assertEqual(K.classify(200.0, 900.0, "IU"), "other")
        self.assertEqual(K.classify(float("nan"), float("nan"), "IU"), "other")


class TestLags(unittest.TestCase):
    def test_pair_lags_recovers_shifts(self):
        sr, n, c = 200.0, 6000, 2000
        A = {k: gabor(n, c, sr) for k in ("R", "T", "Z")}
        sh = 400
        delays = {"R": -24, "T": 6, "Z": 6}
        B = {k: gabor(n, c + sh + d, sr) for k, d in delays.items()}
        wi = {"SKS_R": ("R", c - 400, c + 600), "S_T": ("T", c - 600, c + 2000),
              "SKS_Z": ("Z", c - 400, c + 600)}
        res = K.pair_lags(A, B, wi, sh, 120)
        self.assertAlmostEqual(res["SKS_R"][0], sh - 24, delta=0.3)
        self.assertAlmostEqual(res["S_T"][0], sh + 6, delta=0.3)
        self.assertGreater(res["SKS_R"][1], 0.99)
        self.assertFalse(res["SKS_R"][2])

    def test_window_lag_edge(self):
        a = np.zeros(100)
        _, _, edge = K.window_lag(a, a, 0, 50, 10)
        self.assertTrue(edge)

    def test_snr(self):
        x = np.r_[0.1 * np.ones(100), np.zeros(50), 5.0 * np.ones(10)]
        self.assertAlmostEqual(K.snr(x, 0, 100, 150, 160), 50.0, places=6)


class TestStationX(unittest.TestCase):
    def setUp(self):
        self.lags = {"SKS_R": 1.00, "SKS_Z": 1.01, "S_T": 1.12, "P_Z": 1.05}
        self.locp = {"SKS_R": 0.002, "SKS_Z": 0.002, "S_T": 0.004, "P_Z": 0.001}

    def test_reference_s(self):
        cc = {"S_T": 0.95, "P_Z": 0.95}
        ed = {"S_T": False, "P_Z": False}
        o = K.station_x(self.lags, cc, ed, 1.0, self.locp)
        self.assertEqual(o["ref"], "S_T")
        self.assertAlmostEqual(o["x"], 1.00 - 1.12 - (0.002 - 0.004), places=9)
        self.assertAlmostEqual(o["xz"], 1.01 - 1.12 - (0.002 - 0.004), places=9)

    def test_reference_p_then_abs(self):
        o = K.station_x(self.lags, {"S_T": 0.5, "P_Z": 0.93}, {"S_T": False, "P_Z": False},
                        1.0, self.locp)
        self.assertEqual(o["ref"], "P_Z")
        self.assertAlmostEqual(o["x"], 1.00 - 1.05 - (0.002 - 0.001), places=9)
        o = K.station_x(self.lags, {"S_T": 0.95, "P_Z": 0.95}, {"S_T": True, "P_Z": True},
                        0.95, self.locp)
        self.assertEqual(o["ref"], "abs")
        self.assertAlmostEqual(o["x"], 1.00 - 0.95 - 0.002, places=9)

    def test_loc_pred(self):
        self.assertAlmostEqual(K.loc_pred(1.0, 0.0, 0.0, 0.1), -0.1, places=12)
        self.assertAlmostEqual(K.loc_pred(0.0, 2.0, 90.0, 0.05), -0.1, places=12)


class TestGrades(unittest.TestCase):
    def test_grade_station(self):
        self.assertEqual(K.grade_station(0.03, 0.0, True), "stable")
        self.assertEqual(K.grade_station(-0.12, -0.06, True), "changed")
        self.assertEqual(K.grade_station(-0.12, 0.06, True), "unclear")
        self.assertEqual(K.grade_station(-0.12, -0.02, True), "unclear")
        self.assertEqual(K.grade_station(0.06, 0.06, True), "unclear")
        self.assertEqual(K.grade_station(-0.2, -0.2, False), "unusable")
        self.assertEqual(K.grade_station(float("nan"), 0.0, True), "unusable")

    def test_grade_group(self):
        gates = {"G0": True, "G1": True}
        quiet = [("stable", 0.01)] * 5 + [("unclear", 0.09)]
        noisy = [("stable", 0.01)] * 3 + [("changed", 0.2)] * 2
        g = K.grade_group(gates, [("changed", -0.15)], quiet)
        self.assertEqual(g["grade"], "shared-change")
        self.assertAlmostEqual(g["ctrl_false_rate"], 1 / 6, places=9)
        self.assertEqual(K.grade_group(gates, [("changed", -0.15)], noisy)["grade"], "unclear")
        self.assertEqual(K.grade_group(gates, [("changed", -0.15)], quiet[:3])["grade"],
                         "unclear")
        self.assertEqual(K.grade_group(gates, [("stable", 0.01), ("stable", -0.03)],
                                       [])["grade"], "no-shared-change")
        self.assertEqual(K.grade_group(gates, [("stable", 0.01), ("unclear", 0.06)],
                                       quiet)["grade"], "unclear")
        self.assertEqual(K.grade_group(gates, [("unusable", None)], quiet)["grade"],
                         "untestable")
        self.assertEqual(K.grade_group({"G0": False}, [("stable", 0.0)], quiet)["grade"],
                         "untestable")


if __name__ == "__main__":
    unittest.main()

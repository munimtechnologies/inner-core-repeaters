# Copyright 2026 MV2, Munim, Inc.
# Licensed under the Apache License, Version 2.0 (Apache-2.0). See LICENSE.
# Made with AI agents working in MV2's name.
"""Synthetic tests for sdvexit27.py (no network, no ObsPy)."""

import math
import unittest

import numpy as np

import sdvexit27 as X


def gabor(n, c, sr, f=1.0, w=0.25):
    i = np.arange(n, dtype=float)
    t = (i - c) / sr
    return np.exp(-(t / w) ** 2) * np.cos(2 * math.pi * f * t)


def ev(eid, iso, lat, lon, depth=30.0, mag=5.8):
    return dict(id=eid, iso=iso, lat=lat, lon=lon, depth=depth, mag=mag, magtype="mww",
                ts=X.iso_epoch(iso))


CSV = (
    "time,latitude,longitude,depth,mag,magType,nst,gap,dmin,rms,net,id,updated,place,type\n"
    "2000-07-07T15:46:44.560Z,51.411,179.978,31,5.9,mwc,,,,1.1,us,usp0009vx1,"
    "2016-11-09T21:33:00.000Z,\"Rat Islands, Aleutian Islands, Alaska\",earthquake\n"
    "2001-01-01T00:00:00.000Z,10.0,20.0,10,5.5,mb,,,,1.0,us,usx1,"
    "2016-11-09T21:33:00.000Z,\"somewhere, far\",explosion\n"
    "2009-12-17T20:01:21.610Z,51.398,179.963,35,,mwc,,,,1.0,us,usnomag,"
    "2016-11-09T21:33:00.000Z,\"no magnitude\",earthquake\n"
    "2012-03-04T05:06:07Z,52.0,-178.0,40.5,5.4,mww,,,,1.0,us,usx2,"
    "2016-11-09T21:33:00.000Z,place,earthquake\n"
)


class TestSphere(unittest.TestCase):
    def test_dest_point_roundtrip(self):
        for az, d in ((37.0, 25.0), (250.0, 80.0), (329.3, 9.0)):
            lat, lon = X.dest_point(8.88, -70.63, az, d)
            self.assertAlmostEqual(X.ang_dist(8.88, -70.63, lat, lon), d, places=5)
            self.assertAlmostEqual(X.azimuth(8.88, -70.63, lat, lon), az, places=4)

    def test_azimuth_axes(self):
        self.assertAlmostEqual(X.azimuth(0, 0, 0, 10), 90.0, places=6)
        self.assertAlmostEqual(X.azimuth(0, 0, 10, 0), 0.0, places=6)
        self.assertAlmostEqual(X.azimuth(0, 0, 0, -10), 270.0, places=6)

    def test_wrap180(self):
        self.assertAlmostEqual(X.wrap180(350.0), -10.0)
        self.assertAlmostEqual(X.wrap180(-190.0), 170.0)

    def test_hav_vec_matches_scalar_and_dateline(self):
        d = X.hav_km_vec(51.4, 179.98, [51.4, 51.4], [-179.98, 179.98])
        self.assertAlmostEqual(d[0], math.radians(X.ang_dist(51.4, 179.98, 51.4, -179.98))
                               * X.R_EARTH, places=6)
        self.assertLess(d[0], 3.0)
        self.assertAlmostEqual(d[1], 0.0, places=9)

    def test_arc_km_and_cmb(self):
        self.assertAlmostEqual(X.arc_km(0, 0, 0, 1), math.radians(1) * 3480.0, places=6)
        d = [31.0, 410.0, 2889.0, 5153.5, 2889.0, 0.0]
        self.assertEqual(X.cmb_indices(d, 2889.0), (2, 4))
        self.assertEqual(X.cmb_indices([0.0], 2889.0), (None, None))


class TestCatalogue(unittest.TestCase):
    def test_iso_epoch(self):
        self.assertAlmostEqual(X.iso_epoch("1970-01-01T00:00:10.5Z"), 10.5)
        self.assertEqual(X.iso_epoch("2000-01-01"), 946684800)
        self.assertAlmostEqual(X.iso_epoch("2000-01-01T01:02:03"), 946684800 + 3723)

    def test_parse_usgs_csv(self):
        evs = X.parse_usgs_csv(CSV)
        self.assertEqual([e["id"] for e in evs], ["usp0009vx1", "usx2"])
        e = evs[0]
        self.assertEqual(e["iso"], "2000-07-07T15:46:44.560")
        self.assertAlmostEqual(e["lat"], 51.411)
        self.assertAlmostEqual(e["mag"], 5.9)
        self.assertAlmostEqual(e["ts"], X.iso_epoch("2000-07-07T15:46:44.56"))
        self.assertAlmostEqual(evs[1]["depth"], 40.5)

    def test_find_pairs(self):
        A = ev("a", "2000-07-07T15:46:44", 51.40, 179.98)
        post = [
            ev("ok", "2010-01-01T00:00:00", 51.43, 179.98),          # 3.3 km, 9.5 yr
            ev("short", "2004-01-01T00:00:00", 51.41, 179.98),       # span < 5 yr
            ev("far", "2012-01-01T00:00:00", 51.70, 179.98),         # 33 km
            ev("dl", "2015-01-01T00:00:00", 51.40, -179.99),         # across the dateline
            ev("bigm", "2015-01-01T00:00:00", 51.40, 179.99, mag=6.7),  # dmag 0.9
            ev("deep", "2015-01-01T00:00:00", 51.40, 179.99, depth=70.0),
        ]
        ps = X.find_pairs([A], post)
        self.assertEqual(sorted(p["B"]["id"] for p in ps), ["dl", "ok"])
        p = [q for q in ps if q["B"]["id"] == "ok"][0]
        self.assertAlmostEqual(p["sep_km"], 3.336, delta=0.05)
        self.assertGreater(p["span_yr"], 9.4)
        self.assertEqual(X.find_pairs([], post), [])

    def test_greedy_unique_and_doublet(self):
        a1 = ev("a1", "2000-01-01", 0, 0)
        a2 = ev("a2", "2001-01-01", 0, 0)
        b1 = ev("b1", "2010-01-01", 0, 0)
        b2 = ev("b2", "2011-01-01", 0, 0)
        ps = [dict(A=a1, B=b1, sep_km=5.0, mag_min=5.5),
              dict(A=a1, B=b2, sep_km=2.0, mag_min=5.5),
              dict(A=a2, B=b2, sep_km=1.0, mag_min=5.5),
              dict(A=a2, B=b1, sep_km=3.0, mag_min=5.5)]
        g = X.greedy_unique(ps)
        self.assertEqual([(p["A"]["id"], p["B"]["id"]) for p in g], [("a2", "b2"), ("a1", "b1")])
        self.assertTrue(X.is_doublet(ps[0], a1["ts"] + 10, b1["ts"] - 10))
        self.assertFalse(X.is_doublet(ps[0], a1["ts"] + 100, b1["ts"]))


class TestClasses(unittest.TestCase):
    def test_classify_pair(self):
        self.assertEqual(X.classify_pair(50.0, 60.0), "near")
        self.assertEqual(X.classify_pair(800.0, 60.0), "exit")
        self.assertEqual(X.classify_pair(50.0, 400.0), "other")
        self.assertEqual(X.classify_pair(400.0, 900.0), "control")
        self.assertEqual(X.classify_pair(200.0, 900.0), "other")
        self.assertEqual(X.classify_pair(float("nan"), 50.0), "other")


class TestLags(unittest.TestCase):
    def test_pair_lags_recovers_shifts(self):
        sr, n, c = 200.0, 6000, 2000
        A = {k: gabor(n, c, sr) for k in ("R", "T", "Z")}
        sh = 400
        delays = {"R": -24, "T": 6, "Z": 6}
        B = {k: gabor(n, c + sh + d, sr) for k, d in delays.items()}
        wi = {"SKS_R": ("R", c - 400, c + 600), "S_T": ("T", c - 600, c + 2000)}
        res = X.pair_lags(A, B, wi, sh, 120)
        self.assertAlmostEqual(res["SKS_R"][0], sh - 24, delta=0.3)
        self.assertAlmostEqual(res["S_T"][0], sh + 6, delta=0.3)
        self.assertGreater(res["SKS_R"][1], 0.99)
        self.assertFalse(res["SKS_R"][2])

    def test_window_lag_edge_and_snr(self):
        _, _, edge = X.window_lag(np.zeros(100), np.zeros(100), 0, 50, 10)
        self.assertTrue(edge)
        x = np.r_[0.1 * np.ones(100), np.zeros(50), 5.0 * np.ones(10)]
        self.assertAlmostEqual(X.snr(x, 0, 100, 150, 160), 50.0, places=6)

    def test_station_x(self):
        lags = {"SKS_R": 1.00, "SKS_Z": 1.01, "S_T": 1.12, "P_Z": 1.05}
        locp = {"SKS_R": 0.002, "SKS_Z": 0.002, "S_T": 0.004, "P_Z": 0.001}
        o = X.station_x(lags, {"S_T": 0.95, "P_Z": 0.95}, {"S_T": False, "P_Z": False}, 1.0,
                        locp)
        self.assertEqual(o["ref"], "S_T")
        self.assertAlmostEqual(o["x"], 1.00 - 1.12 + 0.002, places=9)
        o = X.station_x(lags, {"S_T": 0.5, "P_Z": 0.5}, {"S_T": False, "P_Z": False}, 0.95,
                        locp)
        self.assertEqual(o["ref"], "abs")
        self.assertAlmostEqual(o["x"], 1.00 - 0.95 - 0.002, places=9)
        self.assertAlmostEqual(X.loc_pred(0.0, 2.0, 90.0, 0.05), -0.1, places=12)


class TestGrades(unittest.TestCase):
    def test_grade_pair(self):
        self.assertEqual(X.grade_pair(0.03, 0.0, True), "stable")
        self.assertEqual(X.grade_pair(-0.12, -0.06, True), "changed")
        self.assertEqual(X.grade_pair(-0.12, 0.06, True), "unclear")
        self.assertEqual(X.grade_pair(0.06, 0.06, True), "unclear")
        self.assertEqual(X.grade_pair(-0.2, -0.2, False), "unusable")

    def test_grade_class(self):
        G = {"G0": True, "G1": True}
        quiet = [("stable", 0.01), ("stable", -0.02), ("unclear", 0.05)]
        noisy = [("stable", 0.01), ("stable", 0.0), ("changed", 0.2)]
        two = [("changed", -0.12), ("changed", -0.10), ("stable", 0.01)]
        self.assertEqual(X.grade_class(G, two, quiet)["grade"], "change")
        self.assertEqual(X.grade_class(G, two, quiet[:2])["grade"], "unclear")
        self.assertEqual(X.grade_class(G, two, noisy)["grade"], "unclear")
        mixed = [("changed", -0.12), ("changed", 0.10), ("stable", 0.01)]
        self.assertEqual(X.grade_class(G, mixed, quiet)["grade"], "unclear")
        self.assertEqual(X.grade_class(G, [("changed", -0.1)], quiet)["grade"], "change")
        one_of_three = [("changed", -0.12), ("stable", 0.0), ("stable", 0.01)]
        self.assertEqual(X.grade_class(G, one_of_three, quiet)["grade"], "unclear")
        st = [("stable", 0.0), ("stable", 0.02), ("unclear", 0.06)]
        r = X.grade_class(G, st, [])
        self.assertEqual(r["grade"], "stable")
        self.assertAlmostEqual(r["median_x"], 0.02)
        self.assertEqual(X.grade_class(G, [("stable", 0.0), ("unclear", 0.06),
                                           ("unclear", -0.05)], quiet)["grade"], "unclear")
        self.assertEqual(X.grade_class(G, [("unusable", None)], quiet)["grade"], "untestable")
        self.assertEqual(X.grade_class({"G0": False}, st, quiet)["grade"], "untestable")
        self.assertAlmostEqual(X.grade_class(G, two, noisy)["ctrl_false_rate"], 1 / 3)


if __name__ == "__main__":
    unittest.main()

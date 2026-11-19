"""Tests for geom.py (pure helpers only; no network, no ObsPy). Apache-2.0.
Part of inner-core-repeaters by MV2, Munim, Inc. (made with AI agents)."""
import math
import unittest

import geom as g


class TestSphere(unittest.TestCase):
    def test_gc(self):
        self.assertAlmostEqual(g.gc_deg(0, 0, 0, 90), 90.0, places=9)
        # asin near 1 loses precision: allow 1e-4 deg at the antipode
        self.assertAlmostEqual(g.gc_deg(10, 20, -10, -160), 180.0, delta=1e-4)

    def test_azimuth(self):
        self.assertAlmostEqual(g.azimuth_deg(0, 0, 0, 10), 90.0, places=9)
        self.assertAlmostEqual(g.azimuth_deg(0, 0, 10, 0), 0.0, places=9)
        self.assertAlmostEqual(g.azimuth_deg(0, 0, -10, 0), 180.0, places=9)

    def test_destination_round_trip(self):
        for az, ang in ((30.0, 40.0), (200.0, 150.0), (355.0, 10.0)):
            lat, lon = g.destination(g.CENTROID[0], g.CENTROID[1], az, ang)
            self.assertAlmostEqual(g.gc_deg(g.CENTROID[0], g.CENTROID[1], lat, lon), ang, places=6)
            self.assertAlmostEqual(g.azimuth_deg(g.CENTROID[0], g.CENTROID[1], lat, lon), az, places=5)

    def test_destination_antipode(self):
        lat, lon = g.destination(g.CENTROID[0], g.CENTROID[1], 77.0, 180.0)
        self.assertAlmostEqual(lat, g.ANTIPODE[0], places=6)
        self.assertLess(g.gc_deg(lat, lon, g.ANTIPODE[0], g.ANTIPODE[1]), 1e-6)

    def test_icb_km(self):
        r = 1217.5
        self.assertAlmostEqual(g.icb_km((0, 0), (0, 1), r), r * math.pi / 180.0, places=6)
        d, lab = g.nearest_ref((0, 0), [("a", (0, 10)), ("b", (0, 2))], r)
        self.assertEqual(lab, "b")
        self.assertAlmostEqual(d, 2 * r * math.pi / 180.0, places=6)


class TestRules(unittest.TestCase):
    def test_classify(self):
        self.assertEqual(g.classify(400.0, 600.0), "separate")
        self.assertEqual(g.classify(100.0, 600.0), "shared-entry")
        self.assertEqual(g.classify(400.0, 100.0), "shared-exit")
        self.assertEqual(g.classify(100.0, 100.0), "shared-both")
        self.assertEqual(g.classify(float("nan"), 100.0), "unknown")
        self.assertEqual(g.classify(200.0, 600.0, thr=150.0), "separate")

    def test_f_separate(self):
        cl = {"A": "separate", "B": "shared-entry", "C": "separate"}
        self.assertAlmostEqual(g.f_separate(cl, {"A": 30, "B": 10, "C": 0}), 0.75)
        self.assertAlmostEqual(g.f_separate(cl, {"A": 1, "B": 1, "C": 1}), 2.0 / 3.0)
        self.assertTrue(math.isnan(g.f_separate(cl, {"Z": 5})))
        self.assertTrue(math.isnan(g.f_separate(cl, {})))

    def test_reading(self):
        self.assertEqual(g.reading("B-supported", 0.8), "R1")
        self.assertEqual(g.reading("B-supported", 0.75), "R1")
        self.assertEqual(g.reading("B-supported", 0.5), "R2")
        self.assertEqual(g.reading("B-supported", float("nan")), "R2")
        self.assertEqual(g.reading("M-supported", 0.9), "R3-separate")
        self.assertEqual(g.reading("M-supported", 0.1), "R3-shared")
        self.assertEqual(g.reading("no-time-signal", 0.1), "R4")
        for h in ("untestable", "inconclusive", "source-confounded:B"):
            self.assertEqual(g.reading(h, 0.9), "none")
        for k in ("R1", "R2", "R3-separate", "R3-shared", "R4", "none"):
            self.assertIn(k, g.READINGS)


class TestPath(unittest.TestCase):
    def test_ic_metrics_synthetic(self):
        icb = 5153.5
        # source at the surface, down to the ICB at 10 deg, bottom 200 km below it at 20 deg,
        # back up to the ICB at 30 deg, then to the surface at 40 deg
        d = [math.radians(x) for x in (0.0, 10.0, 20.0, 30.0, 40.0)]
        z = [0.0, icb, icb + 200.0, icb, 0.0]
        m = g.ic_metrics(d, z, icb, (0.0, 0.0), 90.0)
        self.assertAlmostEqual(m["entry_ang"], 10.0, places=9)
        self.assertAlmostEqual(m["exit_ang"], 30.0, places=9)
        self.assertAlmostEqual(m["bottom_ang"], 20.0, places=9)
        self.assertAlmostEqual(m["entry_lat"], 0.0, places=6)
        self.assertAlmostEqual(m["entry_lon"], 10.0, places=6)
        self.assertAlmostEqual(m["exit_lon"], 30.0, places=6)
        self.assertAlmostEqual(m["bottom_below_icb_km"], 200.0, places=9)
        self.assertAlmostEqual(m["upper_ic_frac"], 1.0, places=9)
        r_icb, r_b = g.R_EARTH - icb, g.R_EARTH - icb - 200.0
        a, b = math.radians(10.0), math.radians(20.0)
        seg = math.hypot(r_b * math.cos(b) - r_icb * math.cos(a), r_b * math.sin(b) - r_icb * math.sin(a))
        self.assertAlmostEqual(m["ic_path_km"], 2 * seg, places=6)
        # on the equator every point is (r) from the spin axis
        self.assertGreater(m["mean_axis_dist_km"], r_b - 1e-6)
        self.assertLess(m["mean_axis_dist_km"], r_icb + 1e-6)
        self.assertAlmostEqual(m["entry_km_per_deg_rot"], r_icb * math.pi / 180.0, places=6)

    def test_ic_metrics_deep_share(self):
        icb = 5153.5
        d = [math.radians(x) for x in (0.0, 10.0, 20.0, 30.0, 40.0)]
        z = [0.0, icb, icb + 700.0, icb, 0.0]
        m = g.ic_metrics(d, z, icb, (0.0, 0.0), 90.0)
        self.assertEqual(m["upper_ic_frac"], 0.0)   # both segment midpoints are 350 km below the ICB

    def test_ic_metrics_no_core(self):
        self.assertIsNone(g.ic_metrics([0.0, 0.1], [0.0, 100.0], 5153.5, (0.0, 0.0), 0.0))


class TestStations(unittest.TestCase):
    TEXT = (
        "#Network | Station | Location | Channel | Latitude | Longitude | Elevation | Depth | Azimuth | Dip | "
        "SensorDescription | Scale | ScaleFreq | ScaleUnits | SampleRate | StartTime | EndTime\n"
        "IU|PET|00|BHZ|53.02|158.65|110|0|0|-90|STS-1|1|1|M/S|20|1993-01-01T00:00:00.0000|2010-01-01T00:00:00.0000\n"
        "IU|PET|10|BHZ|53.02|158.65|110|0|0|-90|STS-2|1|1|M/S|40|2005-01-01T00:00:00.0000|\n"
        "II|FAR|00|BHZ|0.0|0.0|0|0|0|-90|X|1|1|M/S|20|2000-01-01T00:00:00|\n"
    )

    def test_parse(self):
        st = g.parse_station_text(self.TEXT)
        self.assertEqual(set(st), {("IU", "PET"), ("II", "FAR")})
        self.assertEqual(len(st[("IU", "PET")]["chans"]), 2)
        self.assertEqual(st[("IU", "PET")]["chans"][1][2], "9999")

    def test_choose_loc(self):
        st = g.parse_station_text(self.TEXT)
        ch = st[("IU", "PET")]["chans"]
        self.assertEqual(g.choose_loc(ch, "2006-05-01T00:00:00"), "00")
        self.assertEqual(g.choose_loc(ch, "2012-05-01T00:00:00Z".rstrip("Z")), "10")
        self.assertIsNone(g.choose_loc(ch, "1990-01-01T00:00:00"))

    def test_select(self):
        st = g.parse_station_text(self.TEXT)
        times = ["2001-0%d-01T00:00:00Z" % k for k in range(1, 10)] * 3   # 27 events
        found, chosen = g.select_candidates(st, times)
        names = ["%s.%s" % (f["net"], f["sta"]) for f in found]
        self.assertIn("IU.PET", names)       # about 173 deg from the SSI centroid
        self.assertNotIn("II.FAR", names)    # (0, 0) is far outside 155-179.5 deg
        self.assertEqual(chosen, ["IU.PET"])
        found2, chosen2 = g.select_candidates(st, times[:5])
        self.assertEqual(chosen2, [])        # cover 5 < 20

    def test_members(self):
        self.assertEqual(g.parse_members("a(2001);b(2002); "), ["a", "b"])


if __name__ == "__main__":
    unittest.main()

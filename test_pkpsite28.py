# Unit tests for pkpsite28.py (pure helpers only; stdlib + numpy). Apache-2.0.
# Made with AI agents working in MV2's name (MV2, Munim, Inc.).
import math
import unittest

import pkpsite28 as P


class TestGeometry(unittest.TestCase):
    def test_quarter_circle_on_equator(self):
        d, az, baz = P.gc_dist_az(0.0, 0.0, 0.0, 90.0)
        self.assertAlmostEqual(d, 90.0, places=6)
        self.assertAlmostEqual(az, 90.0, places=6)
        self.assertAlmostEqual(baz, 270.0, places=6)

    def test_due_north(self):
        d, az, _ = P.gc_dist_az(10.0, 20.0, 30.0, 20.0)
        self.assertAlmostEqual(d, 20.0, places=6)
        self.assertAlmostEqual(az % 360.0, 0.0, places=6)

    def test_mururoa_brvk_rough(self):
        d, az, _ = P.gc_dist_az(P.MURUROA[0], P.MURUROA[1], P.BRVK[0], P.BRVK[1])
        self.assertTrue(140.5 < d < 143.0, d)
        self.assertTrue(320.0 < az < 340.0, az)

    def test_offset_point_roundtrip(self):
        lat, lon = P.offset_point(-21.8, -138.9, 331.0, 10.0)
        d, az, _ = P.gc_dist_az(-21.8, -138.9, lat, lon)
        self.assertAlmostEqual(d * 6371.0 * math.pi / 180.0, 10.0, places=3)
        self.assertAlmostEqual(az, 331.0, places=2)

    def test_along_across(self):
        lat, lon = P.offset_point(-21.8, -138.9, 331.0, 5.0)
        al, ac = P.along_across_km(-21.8, -138.9, lat, lon, 331.0)
        self.assertAlmostEqual(al, 5.0, delta=0.05)
        self.assertAlmostEqual(ac, 0.0, delta=0.05)
        lat2, lon2 = P.offset_point(-21.8, -138.9, 61.0, 3.0)   # 90 deg to the right
        al2, ac2 = P.along_across_km(-21.8, -138.9, lat2, lon2, 331.0)
        self.assertAlmostEqual(al2, 0.0, delta=0.05)
        self.assertAlmostEqual(ac2, 3.0, delta=0.05)


class TestHelpers(unittest.TestCase):
    def test_slope(self):
        self.assertAlmostEqual(P.slope([0, 1, 2, 3], [1, 3, 5, 7]), 2.0)
        self.assertAlmostEqual(P.slope([0, 1, 2], [1, float("nan"), 5]), 2.0)
        self.assertTrue(math.isnan(P.slope([1], [2])))

    def test_parse_iso(self):
        self.assertAlmostEqual(P.parse_iso("1970-01-01T00:00:01.5Z"), 1.5)
        self.assertAlmostEqual(P.parse_iso("1970-01-02 00:00:00"), 86400.0)

    def test_parse_zenodo_line_tagged(self):
        line = ("1984-12-02T03:17:00.005  0.02400    60729    470805420.00500   6  90   5 "
                "BRVK_s07Z 19841202.0317.BRV.SKM.s07Z0.024.txt")
        r = P.parse_zenodo_line(line)
        self.assertEqual(r["tag"], "SKM")
        self.assertEqual(r["channel"], "BRVK_s07Z")
        self.assertAlmostEqual(r["dt"], 0.024)
        self.assertEqual(r["nsamp"], 60729)
        self.assertAlmostEqual(r["ok"], 90.0)

    def test_parse_zenodo_line_untagged_and_kod(self):
        r = P.parse_zenodo_line("1977-08-17T04:26:59.999  0.02600  1000  -1.0  0 100 0 BRVK_sZ07 "
                                "19770817.042659.BRV.sZ070.026.txt")
        self.assertEqual(r["tag"], "")
        r2 = P.parse_zenodo_line("1966-12-18T04:57:30.000  0.03000    44992    -95886150.00000   7  88"
                                 "   6 BRVK_SHE 19661218.0457.brvk.KODB.SHE0.030.txt")
        self.assertEqual(r2["tag"], "KODB")
        self.assertIsNone(P.parse_zenodo_line("garbage line"))

    def test_match_by_time(self):
        self.assertEqual(P.match_by_time([0.0, 100.0, 700.0], 50.0, 60.0), [0, 1])

    def test_close_pairs(self):
        ev = [dict(lat=-21.80, lon=-138.90), dict(lat=-21.83, lon=-138.92),
              dict(lat=-22.24, lon=-138.74)]
        pairs = P.close_pairs(ev, 0.05)
        self.assertEqual([(i, j) for i, j, _ in pairs], [(0, 1)])

    def test_clean(self):
        self.assertEqual(P._clean({"a": float("nan"), "b": [1.0, float("inf")]}),
                         {"a": None, "b": [1.0, None]})


if __name__ == "__main__":
    unittest.main()

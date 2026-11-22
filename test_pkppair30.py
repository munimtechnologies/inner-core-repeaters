# Tests for pkppair30.py (Apache-2.0). Synthetic only; stdlib + numpy.
# Made with AI agents working in MV2's name (MV2, Munim, Inc.).
import math
import unittest

import pkppair30 as P


class TestParsing(unittest.TestCase):
    def test_coord_forms(self):
        self.assertAlmostEqual(P.parse_coord("21.83S", "lat"), -21.83)
        self.assertAlmostEqual(P.parse_coord("-21.83", "lat"), -21.83)
        self.assertAlmostEqual(P.parse_coord("21.83", "lat"), -21.83)
        self.assertAlmostEqual(P.parse_coord("138.92W", "lon"), -138.92)
        self.assertAlmostEqual(P.parse_coord("138.92", "lon"), -138.92)
        self.assertAlmostEqual(P.parse_coord("21°49.8'S", "lat"), -(21 + 49.8 / 60), places=6)
        self.assertIsNone(P.parse_coord("", "lat"))
        self.assertIsNone(P.parse_coord("abc", "lat"))

    def test_dates(self):
        self.assertEqual(P.parse_date("1977-11-24"), "1977-11-24")
        self.assertEqual(P.parse_date("1977/11/24"), "1977-11-24")
        self.assertEqual(P.parse_date("24/11/1977"), "1977-11-24")
        self.assertEqual(P.parse_date("19771124"), "1977-11-24")
        self.assertIsNone(P.parse_date("x"))

    def test_parse_s1_loose_headers(self):
        txt = ("Event date,Origin time (UT),Latitude,Longitude,mb,Name\n"
               "1977-11-24,16:59:59.9,21.879S,138.939W,5.9,Enee\n"
               "1982-07-25,18:02:00,-21.855,-138.922,6.0,Laios\n"
               "bad,,,,,\n")
        rows = P.parse_s1(txt)
        self.assertEqual(len(rows), 2)
        self.assertEqual(rows[0]["date"], "1977-11-24")
        self.assertAlmostEqual(rows[0]["lat"], -21.879)
        self.assertAlmostEqual(rows[0]["lon"], -138.939)
        self.assertAlmostEqual(rows[1]["mag"], 6.0)
        self.assertEqual(rows[1]["name"], "Laios")
        self.assertEqual(rows[0]["time"], "16:59:59.9")

    def test_tod(self):
        self.assertAlmostEqual(P.tod_seconds("17:31:58.48"), 17 * 3600 + 31 * 60 + 58.48)
        self.assertIsNone(P.tod_seconds(""))


class TestMatchPairs(unittest.TestCase):
    def setUp(self):
        self.shafts = [
            dict(name="A", date="1977-11-24", time_ut="16:59:59.9", lat=-21.879, lon=-138.939, approx=0, atoll="Mururoa"),
            dict(name="B", date="1982-07-25", time_ut="18:02:00.0", lat=-21.855, lon=-138.922, approx=1, atoll="Mururoa"),
            dict(name="C1", date="1990-06-02", time_ut="17:29:58.5", lat=-21.844, lon=-138.894, approx=1, atoll="Mururoa"),
            dict(name="C2", date="1990-06-02", time_ut="20:00:00", lat=-22.24, lon=-138.73, approx=0, atoll="Fangataufa"),
        ]

    def test_match(self):
        s1 = [dict(date="1977-11-24", time=""), dict(date="1982-07-26", time=""),
              dict(date="1990-06-02", time="17:30"), dict(date="1990-06-02", time=""),
              dict(date="1985-01-01", time="")]
        self.assertEqual(P.match_events(s1, self.shafts), [0, 1, 2, None, None])

    def test_pairs(self):
        ev = [dict(lat=-21.85, lon=-138.90), dict(lat=-21.87, lon=-138.93), dict(lat=-21.85, lon=-138.97),
              dict(lat=-22.24, lon=-138.73)]
        pairs = P.form_pairs(ev, 0.05)
        self.assertIn((0, 1), pairs)      # ~0.034 deg arc
        self.assertNotIn((0, 2), pairs)   # 0.07 deg of longitude apart, ~0.065 deg arc
        self.assertNotIn((0, 3), pairs)
        self.assertIn((1, 2), pairs)


class TestInversion(unittest.TestCase):
    def test_year_terms_recover(self):
        years = [1977, 1977, 1982, 1982, 1988, 1995]
        truth = {1977: 0.0, 1982: -0.1, 1988: -0.15, 1995: 0.2}
        off = 3.0   # a common constant cancels
        x = [off + truth[y] for y in years]
        pairs = [(i, j) for i in range(6) for j in range(i + 1, 6)]
        t = P.year_terms(x, years, pairs)
        for y, v in truth.items():
            self.assertAlmostEqual(t[y], v, places=9)
        ep = P.epoch_terms(t, years)
        self.assertAlmostEqual(ep["E1982_83"], -0.1)
        self.assertAlmostEqual(ep["E1988_90"], -0.15)
        self.assertAlmostEqual(ep["E1995"], 0.2)

    def test_year_terms_mean_and_missing(self):
        years = [1977, 1977, 1982, 1982, 1990]
        x = [0.0, 0.0, 1.0, 1.2, None]
        pairs = [(i, j) for i in range(5) for j in range(i + 1, 5)]
        t = P.year_terms(x, years, pairs)
        self.assertAlmostEqual(t[1982], 1.1)
        self.assertNotIn(1990, t)
        ep = P.epoch_terms(t, years)
        self.assertIsNone(ep["E1988_90"])

    def test_disconnected_year_left_out(self):
        years = [1977, 1982, 1995]
        x = [0.0, 0.1, 0.3]
        t = P.year_terms(x, years, [(0, 1)])
        self.assertIn(1982, t)
        self.assertNotIn(1995, t)

    def test_no_ref(self):
        self.assertEqual(P.year_terms([0.0, 1.0], [1982, 1995], [(0, 1)]), {})


class TestGrades(unittest.TestCase):
    def test_explains(self):
        g = P.grade({"E1982_83": -0.12, "E1988_90": -0.14, "E1995": 0.16})
        self.assertEqual(g["grade"], "explains")
        self.assertEqual(g["signs_match"], 3)

    def test_partly(self):
        # signs right everywhere, residual rms ~0.076 s: above 0.05, below 0.75 x rms(Zhou) ~0.114
        g = P.grade({"E1982_83": -0.10, "E1988_90": -0.12, "E1995": 0.05})
        self.assertEqual(g["grade"], "partly")
        # two signs right but too far off in size: not partly
        g = P.grade({"E1982_83": -0.08, "E1988_90": -0.06, "E1995": 0.0})
        self.assertEqual(g["grade"], "does-not-explain")

    def test_does_not(self):
        g = P.grade({"E1982_83": -0.036, "E1988_90": 0.018, "E1995": -0.09})
        self.assertEqual(g["grade"], "does-not-explain")
        g = P.grade({"E1982_83": 0.0, "E1988_90": 0.0, "E1995": 0.0})
        self.assertEqual(g["grade"], "does-not-explain")
        self.assertAlmostEqual(g["rms_res"], g["rms_zhou"])

    def test_small_values_have_no_sign(self):
        g = P.grade({"E1982_83": -0.01, "E1988_90": -0.15, "E1995": 0.175})
        self.assertEqual(g["signs_match"], 2)

    def test_untestable(self):
        self.assertEqual(P.grade({"E1982_83": -0.1, "E1988_90": None, "E1995": 0.1})["grade"], "untestable")

    def test_headline(self):
        gr = {"R": {"grade": "does-not-explain"}, "C": {"grade": "partly"}}
        self.assertEqual(P.headline({"distance_correction": "none"}, gr)["reading"], "R")
        self.assertEqual(P.headline({"distance_correction": "catalogue-1D"}, gr)["grade"], "partly")
        h = P.headline({}, gr)
        self.assertEqual(h["reading"], "R+C")
        self.assertIn("R:does-not-explain", h["grade"])


class TestGeometry(unittest.TestCase):
    def test_az_and_along(self):
        d, az = P.gc_dist_az(P.REF[0], P.REF[1], P.BRVK[0], P.BRVK[1])
        self.assertTrue(140.0 < d < 143.0)
        self.assertTrue(325.0 < az < 340.0)
        # a point 0.01 deg north of REF is ~1.1 km along an azimuth of 0
        self.assertAlmostEqual(P.along_km(P.REF[0], P.REF[1], P.REF[0] + 0.01, P.REF[1], 0.0), 1.11195, places=4)
        self.assertTrue(math.isfinite(az))


if __name__ == "__main__":
    unittest.main()

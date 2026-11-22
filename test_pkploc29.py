# Unit tests for pkploc29.py (synthetic inputs only; stdlib + numpy). Apache-2.0.
import csv
import math
import os
import tempfile
import unittest

import pkploc29 as P


class Geometry(unittest.TestCase):
    def test_path_to_brvk(self):
        d, az = P.gc_dist_az(P.REF[0], P.REF[1], P.BRVK[0], P.BRVK[1])
        self.assertTrue(140.5 < d < 143.0, d)
        self.assertTrue(325.0 < az < 340.0, az)

    def test_along_across(self):
        al, ac = P.along_across_km(0.0, 0.0, 0.01, 0.0, 0.0)
        self.assertAlmostEqual(al, 1.11195, places=4)
        self.assertAlmostEqual(ac, 0.0, places=6)
        al, ac = P.along_across_km(0.0, 0.0, 0.0, 0.01, 90.0)
        self.assertAlmostEqual(al, 1.11195, places=4)
        al, ac = P.along_across_km(0.0, 0.0, 0.0, 0.01, 0.0)
        self.assertAlmostEqual(ac, 1.11195, places=4)

    def test_north_rim_is_closer_to_brvk_than_southeast_rim(self):
        _, az = P.gc_dist_az(P.REF[0], P.REF[1], P.BRVK[0], P.BRVK[1])
        north, _ = P.along_across_km(P.REF[0], P.REF[1], -21.785, -138.866, az)
        southeast, _ = P.along_across_km(P.REF[0], P.REF[1], -21.875, -138.86, az)
        self.assertGreater(north, 2.0)
        self.assertLess(southeast, -5.0)


class Wikitext(unittest.TestCase):
    def test_coord_templates(self):
        w = ("x {{coord|21.87442|S|138.86477|W|name=Ulysse}} y "
             "{{Coord|21|51|10|S|138|53|39|W|type:event}} z {{coord|-22.233|-138.74}}")
        c = P.parse_coord_templates(w)
        self.assertEqual(len(c), 3)
        self.assertAlmostEqual(c[0][0], -21.87442, places=5)
        self.assertAlmostEqual(c[0][1], -138.86477, places=5)
        self.assertAlmostEqual(c[1][0], -(21 + 51 / 60 + 10 / 3600), places=5)
        self.assertAlmostEqual(c[1][1], -(138 + 53 / 60 + 39 / 3600), places=5)
        self.assertAlmostEqual(c[2][1], -138.74, places=5)

    def test_verify(self):
        rows = [dict(lat=-21.8, lon=-138.9), dict(lat=-21.9, lon=-138.95)]
        self.assertEqual(P.verify_rows(rows, [(-21.8003, -138.9002)]), [True, False])


class TextParsing(unittest.TestCase):
    def test_dates(self):
        self.assertEqual(P.find_date("1983/04/19 18:53:00.2 -21.85 -138.89"), "1983-04-19")
        self.assertEqual(P.find_date("19/04/1983 Eurytos"), "1983-04-19")
        self.assertEqual(P.find_date("19830419 185300"), "1983-04-19")
        self.assertEqual(P.find_date("1983 4 19 18 53"), "1983-04-19")
        self.assertIsNone(P.find_date("no date here 2015-01-01"))

    def test_latlon(self):
        self.assertEqual(P.find_latlon("1983/04/19 -21.852 -138.894 0.0 5.6"), (-21.852, -138.894))
        self.assertEqual(P.find_latlon("21.85S 138.89W mb 5.6"), (-21.85, -138.89))
        self.assertIsNone(P.find_latlon("37.1 -116.0 Nevada"))


class Epochs(unittest.TestCase):
    def test_epoch_of(self):
        self.assertEqual(P.epoch_of(1977), "E1977")
        self.assertEqual(P.epoch_of(1983), "E1982_83")
        self.assertEqual(P.epoch_of(1989), "E1988_90")
        self.assertEqual(P.epoch_of(1995), "E1995")
        self.assertIsNone(P.epoch_of(1986))

    def test_epoch_means(self):
        rows = [dict(date="1977-03-19", atoll="Mururoa", yield_kt=47.0, v=1.0),
                dict(date="1977-11-24", atoll="Mururoa", yield_kt=5.0, v=3.0),
                dict(date="1989-06-03", atoll="Mururoa", yield_kt=20.0, v=10.0),
                dict(date="1989-06-10", atoll="Fangataufa", yield_kt=80.0, v=99.0),
                dict(date="1995-10-27", atoll="Mururoa", yield_kt=39.0, v=float("nan"))]
        m = P.epoch_means(rows, "v")
        self.assertEqual(m["E1977"]["n"], 2)
        self.assertAlmostEqual(m["E1977"]["mean"], 2.0)
        self.assertAlmostEqual(m["E1988_90"]["minus_E1977"], 8.0)
        self.assertEqual(m["E1995"]["n"], 0)
        g = P.epoch_means(rows, "v", min_yield=20.0)
        self.assertAlmostEqual(g["E1977"]["mean"], 1.0)
        self.assertAlmostEqual(g["E1988_90"]["minus_E1977"], 9.0)


class Table(unittest.TestCase):
    def test_load_and_origin(self):
        with tempfile.TemporaryDirectory() as d:
            p = os.path.join(d, "s.csv")
            with open(p, "w", newline="", encoding="utf-8") as f:
                w = csv.writer(f)
                w.writerow(["name", "date", "time_ut", "atoll", "zone", "shaft", "lat", "lon", "approx",
                            "yield_kt", "wiki_page"])
                w.writerow(["Nestor", "1977-03-19", "23:00:59.9", "Mururoa", "z", "Zoe1", "-21.88009",
                            "-138.94102", "0", "47", "p"])
                w.writerow(["Didon", "1978-03-08", "", "Mururoa", "z", "", "-21.78", "-138.89", "0", "", "p"])
            rows = P.load_shafts(p)
        self.assertEqual(len(rows), 2)
        self.assertIsNone(rows[1]["yield_kt"])
        t = P.origin_seconds(rows[0]["date"], rows[0]["time_ut"])
        self.assertAlmostEqual(t % 60.0, 59.9, places=3)
        self.assertIsNone(P.origin_seconds("1978-03-08", ""))

    def test_shipped_table_is_sane(self):
        here = os.path.dirname(os.path.abspath(__file__))
        p = os.path.join(here, "mururoa_shafts_1977_1996.csv")
        if not os.path.exists(p):
            p = os.path.join(here, "data", "mururoa_shafts_1977_1996.csv")
        if not os.path.exists(p):
            self.skipTest("table not next to the tests")
        rows = P.load_shafts(p)
        self.assertGreater(len(rows), 140)
        for r in rows:
            self.assertTrue(-22.3 < r["lat"] < -21.7, r)
            self.assertTrue(-139.05 < r["lon"] < -138.7, r)
            self.assertTrue(r["atoll"] in ("Mururoa", "Fangataufa"), r)
            d, _ = P.gc_dist_az(r["lat"], r["lon"], -22.24, -138.74)
            if r["atoll"] == "Fangataufa":
                self.assertLess(d, 0.05, r)
            else:
                self.assertGreater(d, 0.3, r)
        self.assertEqual(sorted(r["date"] for r in rows), [r["date"] for r in rows])


if __name__ == "__main__":
    unittest.main()

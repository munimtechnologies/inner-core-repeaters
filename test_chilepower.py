# Copyright 2026 MV2, Munim, Inc.
# Licensed under the Apache License, Version 2.0 (Apache-2.0). See LICENSE.
# Made with AI agents working in MV2's name.
import os
import tempfile
import unittest

import chilepower as C


def ev(i, t, lat, lon, depth, mag):
    return dict(id="e%d" % i, time="%d-01-01T00:00:00Z" % int(t), t=float(t), lat=lat, lon=lon,
                depth=depth, mag=mag, magType="mww")


class TestGeometry(unittest.TestCase):
    def test_mkar_distances_match_pilot(self):
        # distances chilepilot-v2 reported for MKAR (chile_rows.csv)
        self.assertAlmostEqual(C.mkar_dist(-22.817, -66.684), 145.447, delta=0.05)
        self.assertAlmostEqual(C.mkar_dist(-37.1085, -74.0511), 160.041, delta=0.05)
        self.assertAlmostEqual(C.mkar_dist(-19.198, -68.902), 143.748, delta=0.05)

    def test_gate_examples(self):
        self.assertGreater(C.mkar_dist(-30.0, -71.0), 150.0)
        self.assertLess(C.mkar_dist(-10.0, -75.0), 143.0)

    def test_era(self):
        self.assertAlmostEqual(C.ERA0, 2000.6557, delta=0.001)
        self.assertAlmostEqual(C.ERA1, 2025.6192, delta=0.001)


class TestClock(unittest.TestCase):
    def test_reoccupied(self):
        self.assertAlmostEqual(C.reoccupied(2016.0), 2005.5, places=9)
        self.assertAlmostEqual(C.reoccupied(2008.5), 2008.5, places=9)

    def test_residual(self):
        self.assertAlmostEqual(C.clock_residual(2005.5, 2016.0), 0.0, places=9)
        self.assertAlmostEqual(C.clock_residual(2003.5, 2016.0), -2.0, places=9)
        self.assertAlmostEqual(C.clock_residual(2007.5, 2024.0), 5.2, places=9)


class TestPairs(unittest.TestCase):
    def setUp(self):
        self.rows = [
            ev(1, 2005.5, -30.0, -71.0, 30.0, 5.2),     # before
            ev(2, 2016.0, -30.05, -71.0, 35.0, 5.4),    # after, 5.56 km from e1
            ev(3, 2024.0, -30.0, -71.06, 30.0, 5.0),    # after, 5.78 km from e1/e4/e5
            ev(4, 2002.0, -30.0, -71.0, 30.0, 5.1),     # before
            ev(5, 2007.5, -30.0, -71.0, 30.0, 5.0),     # before
            ev(6, 2015.0, -30.0, -71.0, 30.0, 4.9),     # after, below M5.0
            ev(7, 1999.0, -30.0, -71.0, 30.0, 5.2),     # before MKAR era
            ev(8, 2026.0, -30.0, -71.0, 30.0, 5.2),     # after MKAR data end
            ev(9, 2004.0, -10.0, -75.0, 30.0, 5.2),     # outside distance gate
            ev(10, 2012.0, -10.0, -75.0, 30.0, 5.2),    # outside distance gate
            ev(11, 2012.0, -30.2, -71.0, 30.0, 5.2),    # 22 km from the cluster
        ]

    def test_eligible(self):
        ids = sorted(r["id"] for r in C.eligible(self.rows))
        self.assertEqual(ids, sorted(["e1", "e2", "e3", "e4", "e5", "e6", "e11"]))

    def test_main_definition(self):
        p = C.straddling_pairs(C.eligible(self.rows), 10.0, 0.3, 20.0, 5.0)
        got = sorted((x["id_before"], x["id_after"]) for x in p)
        self.assertEqual(got, sorted([("e1", "e2"), ("e1", "e3"), ("e4", "e2"), ("e4", "e3"), ("e5", "e3")]))
        s = C.summarise(p)
        self.assertEqual(s["pairs"], 5)
        self.assertEqual(s["near"], 2)    # e1-e2 (r 0), e4-e3 (r -0.3)
        self.assertEqual(s["far"], 1)     # e5-e3 (r 5.2)
        self.assertEqual(s["events"], 5)
        self.assertAlmostEqual(s["median_abs_r"], 3.2, places=6)   # |r|: 0, 0.3, 3.2, 3.5, 5.2

    def test_proxy_and_loose(self):
        ev_ = C.eligible(self.rows)
        self.assertEqual(len(C.straddling_pairs(ev_, 5.0, 0.2, 10.0, 5.0)), 0)
        self.assertEqual(len(C.straddling_pairs(ev_, 20.0, 0.3, 20.0, 5.0)), 5)
        self.assertEqual(len(C.straddling_pairs(ev_, 10.0, 0.3, 20.0, 5.5)), 0)

    def test_decide(self):
        self.assertEqual(C.decide(dict(pairs=60, near=15, far=15))[0], "mkar-multiplet-pilot")
        self.assertEqual(C.decide(dict(pairs=59, near=30, far=30))[0], "close-chile-path")
        self.assertEqual(C.decide(dict(pairs=100, near=14, far=30))[0], "close-chile-path")
        self.assertEqual(C.decide(dict(pairs=100, near=30, far=14))[0], "close-chile-path")

    def test_analyse(self):
        res = C.analyse(self.rows)
        self.assertEqual(res["decision"], "close-chile-path")
        self.assertEqual(res["main"]["pairs"], 5)
        self.assertEqual(res["eligible_m5.0"], dict(n=6, before=3, after=3))
        self.assertFalse(res["reproduce_pilot"]["ok"])
        self.assertIn("r_bins", res["main"])

    def test_read_catalogue_and_csv(self):
        txt = ("time,latitude,longitude,depth,mag,magType,id,type\n"
               "2003-08-05T09:04:38.790Z,-22.817,-66.684,217.3,4.8,mb,a1,earthquake\n"
               "2014-04-03T02:43:13.110Z,-20.5709,-70.4931,22.4,7.7,mww,a2,earthquake\n")
        with tempfile.TemporaryDirectory() as d:
            p = os.path.join(d, "c.csv")
            with open(p, "w", encoding="utf-8") as f:
                f.write(txt)
            rows = C.read_catalogue(p)
            res = C.analyse(rows)
            C.write_pairs_csv(res, os.path.join(d, "pairs.csv"))
            self.assertTrue(os.path.exists(os.path.join(d, "pairs.csv")))
        self.assertEqual(len(rows), 2)
        self.assertAlmostEqual(rows[1]["mag"], 7.7, places=9)
        self.assertLess(rows[0]["t"], 2008.5)


if __name__ == "__main__":
    unittest.main()

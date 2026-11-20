# Copyright 2026 MV2, Munim, Inc.
# Licensed under the Apache License, Version 2.0 (Apache-2.0). See LICENSE.
# Made with AI agents working in MV2's name.
import os
import tempfile
import unittest

import ssipairs as S


def ev(i, t, lat, lon, depth, mag):
    return dict(id="e%d" % i, time="%d-01-01T00:00:00Z" % int(t), t=float(t), lat=lat, lon=lon,
                depth=depth, mag=mag, magType="mww")


class TestGeometry(unittest.TestCase):
    def test_km_zero_and_degree(self):
        self.assertAlmostEqual(S.km_between(-58.0, -25.0, -58.0, -25.0), 0.0, places=6)
        self.assertAlmostEqual(S.km_between(-58.0, -25.0, -59.0, -25.0), 111.19, delta=0.05)

    def test_decimal_year(self):
        self.assertAlmostEqual(S.decimal_year("2008-01-01T00:00:00.000Z"), 2008.0, places=9)
        self.assertLess(S.decimal_year("2008-06-30T06:17:43.020Z"), 2008.5)
        self.assertGreater(S.decimal_year("2008-07-02T12:00:00Z"), 2008.5)


class TestPairs(unittest.TestCase):
    def setUp(self):
        # a: before, b: after 2 km away (dM 0.1, dz 5) -> proxy pair; c: after, 2.65 km, dM 0.4;
        # d: after, 8 km away; e: before, far away
        self.rows = [ev(1, 2000, -58.0, -25.0, 20.0, 6.2),
                     ev(2, 2015, -58.018, -25.0, 25.0, 6.1),
                     ev(3, 2016, -58.0, -25.045, 20.0, 6.6),
                     ev(4, 2017, -58.072, -25.0, 20.0, 6.2),
                     ev(5, 2001, -60.0, -22.0, 20.0, 6.2)]

    def test_proxy(self):
        r = S.close_pairs(self.rows, min_mag=6.0)
        self.assertEqual(r["n_events"], 5)
        self.assertEqual(r["straddling"], 1)
        self.assertEqual(r["events_in_straddling"], 2)

    def test_proxy_min_mag(self):
        r = S.close_pairs(self.rows, min_mag=6.15)
        self.assertEqual(r["straddling"], 0)

    def test_loose_list(self):
        lst = S.straddling_list(self.rows, 10.0, 6.0)
        self.assertEqual(len(lst), 3)
        self.assertEqual([p["id_after"] for p in lst], ["e2", "e3", "e4"])
        self.assertEqual([p["meets_dm_dz"] for p in lst], [True, False, True])
        self.assertLessEqual(lst[0]["km"], lst[1]["km"])

    def test_analyse_decision(self):
        res = S.analyse(self.rows)
        self.assertEqual(res["decision"], "close-near-antipodal-line")
        d = res["by_min_mag"]["6.0"]
        self.assertEqual(d["before_turn"], 2)
        self.assertEqual(d["after_turn"], 3)
        self.assertEqual(d["straddling_within_10km_any"], 3)
        self.assertEqual(d["straddling_within_10km_dm_dz"], 2)

    def test_read_catalogue(self):
        txt = ("time,latitude,longitude,depth,mag,magType,id,type\n"
               "2001-04-13T15:33:53.550Z,-59.723,-25.586,26,6.2,mwc,a1,earthquake\n"
               "2021-08-12T18:35:17.231Z,-58.3753,-25.2637,22.79,8.1,mww,a2,earthquake\n"
               "2021-08-12T18:35:17.231Z,-58.3753,-25.2637,0,3.0,ml,a3,ice quake\n")
        with tempfile.TemporaryDirectory() as d:
            p = os.path.join(d, "c.csv")
            with open(p, "w", encoding="utf-8") as f:
                f.write(txt)
            rows = S.read_catalogue(p)
        self.assertEqual(len(rows), 2)
        self.assertAlmostEqual(rows[1]["mag"], 8.1, places=9)
        self.assertLess(rows[0]["t"], 2008.5)


def find_catalogue():
    here = os.path.dirname(os.path.abspath(__file__))
    for p in (os.environ.get("SSIPAIRS_CATALOGUE", ""),
              os.path.join(here, "results", "snrpilot-v2", "snr_catalogue.csv"),
              os.path.join(here, "snr_catalogue.csv"),
              "/out/snr_catalogue.csv"):
        if p and os.path.exists(p):
            return p
    return None


class TestCommittedCatalogue(unittest.TestCase):
    """The snrpilot-v2 catalogue (ComCat, 300 km of the SSI centroid, M >= 4.5, 1990 to
    2026-11-19). Robust facts only: the nearest straddling M >= 6 pair (4.3 km) fails |dM|."""

    def test_m6(self):
        p = find_catalogue()
        if not p:
            self.skipTest("catalogue not present")
        res = S.analyse(S.read_catalogue(p))
        d = res["by_min_mag"]["6.0"]
        self.assertEqual(d["n_events"], 69)
        self.assertEqual(d["before_turn"], 28)
        self.assertEqual(d["after_turn"], 41)
        self.assertEqual(d["proxy"]["straddling"], 0)
        self.assertEqual(d["straddling_within_5km_any"], 1)
        self.assertEqual(res["decision"], "close-near-antipodal-line")


if __name__ == "__main__":
    unittest.main()

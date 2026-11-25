# Copyright 2026 MV2, Munim, Inc.
# Licensed under the Apache License, Version 2.0 (Apache-2.0). See LICENSE.
# Made with AI agents working in MV2's name.
"""Stdlib-only tests for dprk37probe.py (pure helpers)."""

import os
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import dprk37probe as P  # noqa: E402


class TestGeometry(unittest.TestCase):
    def test_antipode(self):
        self.assertAlmostEqual(P.gc_deg(41.3, 129.05, -41.3, 129.05 - 180.0), 180.0, places=6)

    def test_zero(self):
        self.assertAlmostEqual(P.gc_deg(10.0, 20.0, 10.0, 20.0), 0.0, places=9)
        self.assertAlmostEqual(P.gc_deg(P.SRC_LAT, P.SRC_LON, P.SRC_LAT, P.SRC_LON), 0.0, places=9)

    def test_quarter(self):
        self.assertAlmostEqual(P.gc_deg(0.0, 0.0, 0.0, 90.0), 90.0, places=9)
        self.assertAlmostEqual(P.gc_deg(0.0, 0.0, 90.0, 0.0), 90.0, places=9)

    def test_chile_in_range(self):
        # C1.AC02 (Maricunga, Chile) from the test site
        d = P.gc_deg(P.SRC_LAT, P.SRC_LON, -26.835517, -69.129067)
        self.assertTrue(P.DIST_RANGE[0] <= d <= P.DIST_RANGE[1], d)


class TestBranches(unittest.TestCase):
    def test_two_pkp(self):
        b = P.branches([("PKIKP", 1180.0), ("PKP", 1190.0), ("PKP", 1200.0)], 152.0)
        self.assertEqual(b, dict(df=1180.0, bc=1190.0, ab=1200.0))

    def test_one_pkp_far(self):
        b = P.branches([("PKP", 1230.0), ("PKIKP", 1200.0)], 165.0)
        self.assertEqual(b, dict(df=1200.0, bc=None, ab=1230.0))

    def test_one_pkp_near(self):
        b = P.branches([("PKP", 1150.0), ("PKIKP", 1140.0)], 148.0)
        self.assertEqual(b["bc"], 1150.0)
        self.assertIsNone(b["ab"])

    def test_no_df(self):
        self.assertIsNone(P.branches([("PKP", 1200.0)], 160.0)["df"])

    def test_windows_clear(self):
        self.assertTrue(P.windows_clear(dict(df=1200.0, bc=None, ab=1230.0)))
        self.assertFalse(P.windows_clear(dict(df=1200.0, bc=None, ab=1204.0)))
        self.assertFalse(P.windows_clear(dict(df=1200.0, bc=1202.0, ab=1230.0)))
        self.assertFalse(P.windows_clear(dict(df=1200.0, bc=None, ab=None)))


class TestChannel(unittest.TestCase):
    def test_prefers_bhz_00(self):
        self.assertEqual(P.choose_channel([("10", "HHZ"), ("00", "HHZ"), ("10", "BHZ"),
                                           ("00", "BHZ")]), ("00", "BHZ"))

    def test_empty_loc_before_10(self):
        self.assertEqual(P.choose_channel([("10", "BHZ"), ("", "BHZ")]), ("", "BHZ"))

    def test_hhz_only(self):
        self.assertEqual(P.choose_channel([("", "HHZ")]), ("", "HHZ"))

    def test_none(self):
        self.assertIsNone(P.choose_channel([]))
        self.assertIsNone(P.choose_channel([("00", "LHZ")]))


def rec(ev, st, df, ab, chan="00.BHZ", epoch="2010-01-01", clear=True, status="ok"):
    return {"event": ev, "station": st, "status": status, "channel": chan, "epoch_start": epoch,
            "windows_clear": clear, "snr_df_0.7-2.0": df, "snr_ab_0.7-2.0": ab, "dist": "165.0"}


class TestCounts(unittest.TestCase):
    def test_pair_counts(self):
        e = [x[0] for x in P.EVENTS]
        rows = [
            rec(e[2], "C1.AA", 8.0, 9.0), rec(e[3], "C1.AA", 7.0, 6.0),      # both, same epoch
            rec(e[2], "C1.BB", 8.0, 2.0), rec(e[3], "C1.BB", 7.0, 6.0),      # AB fails in e2
            rec(e[2], "C1.CC", 8.0, 9.0), rec(e[3], "C1.CC", 7.0, 6.0, epoch="2015-01-01"),
            rec(e[2], "C1.DD", 4.0, 9.0), rec(e[3], "C1.DD", 7.0, 6.0),      # DF fails
            rec(e[2], "C1.EE", 8.0, 9.0, clear=False), rec(e[3], "C1.EE", 7.0, 6.0, clear=False),
            rec(e[2], "C1.FF", 8.0, 9.0), rec(e[3], "C1.FF", "", "", status="no data"),
        ]
        out = P.pair_counts(rows, "0.7-2.0")
        self.assertEqual(len(out), 15)
        r = [x for x in out if x["a"] == e[2] and x["b"] == e[3]][0]
        self.assertEqual(r["n_df_both"], 4)          # AA, BB, CC, EE
        self.assertEqual(r["n_df_ab_both"], 3)       # AA, CC, EE
        self.assertEqual(r["n_df_ab_same_epoch"], 2)  # AA, EE
        self.assertEqual(r["n_df_ab_same_epoch_clear"], 1)  # AA
        r0 = [x for x in out if x["a"] == e[0] and x["b"] == e[1]][0]
        self.assertEqual(r0["n_df_both"], 0)

    def test_event_summary(self):
        e = [x[0] for x in P.EVENTS]
        rows = [rec(e[5], "C1.AA", 8.0, 9.0), rec(e[5], "GT.LPAZ", 6.0, 1.0),
                rec(e[5], "IU.TRQA", 3.0, 9.0), rec(e[5], "C1.ZZ", "", "", status="no data")]
        s = P.event_summary(rows, "0.7-2.0")["2017-09-03"]
        self.assertEqual((s["n_listed"], s["n_ok"], s["n_df"], s["n_df_ab"]), (4, 3, 2, 1))
        self.assertEqual(s["networks"], ["C1", "GT"])
        self.assertEqual(s["df_by_distance"], {"160-170": 2})

    def test_record_ok(self):
        self.assertTrue(P.record_ok({"snr_df_1.0-2.0": 5.0}, "1.0-2.0", "df"))
        self.assertFalse(P.record_ok({"snr_df_1.0-2.0": ""}, "1.0-2.0", "df"))
        self.assertFalse(P.record_ok({}, "1.0-2.0", "ab"))


class TestBlind(unittest.TestCase):
    def test_snr_bin(self):
        self.assertEqual(P.snr_bin(2.0), "<3")
        self.assertEqual(P.snr_bin(3.0), "3-5")
        self.assertEqual(P.snr_bin(4.99), "3-5")
        self.assertEqual(P.snr_bin(5.0), "5-10")
        self.assertEqual(P.snr_bin(19.9), "10-20")
        self.assertEqual(P.snr_bin(25.0), ">=20")
        self.assertEqual(P.snr_bin(""), "")
        self.assertEqual(P.snr_bin(None), "")
        self.assertEqual(P.snr_bin(float("nan")), "")

    def test_blind_rows(self):
        e = P.EVENTS[3][0]
        rows = [rec(e, "C1.AA", 8.123, 21.5)]
        out = P.blind_rows(rows)
        self.assertEqual(out[0]["snr_df_0.7-2.0"], "5-10")
        self.assertEqual(out[0]["snr_ab_0.7-2.0"], ">=20")
        self.assertEqual(out[0]["station"], "C1.AA")
        self.assertEqual(rows[0]["snr_df_0.7-2.0"], 8.123)   # original untouched


class TestEventsAndCsv(unittest.TestCase):
    def test_events(self):
        self.assertEqual(len(P.EVENTS), 6)
        dates = [e[1][:10] for e in P.EVENTS]
        self.assertEqual(dates, sorted(dates))
        for e in P.EVENTS:
            self.assertLess(P.gc_deg(P.SRC_LAT, P.SRC_LON, e[2], e[3]), 0.1)

    def test_write_csv(self):
        d = tempfile.mkdtemp()
        p = os.path.join(d, "x.csv")
        P.write_csv(p, [dict(a=1), dict(a=2, b=3)])
        with open(p, encoding="utf-8") as f:
            self.assertEqual(f.read().splitlines(), ["a,b", "1,", "2,3"])
        P.write_csv(p, [])
        with open(p, encoding="utf-8") as f:
            self.assertEqual(f.read().splitlines(), ["empty"])


if __name__ == "__main__":
    unittest.main()

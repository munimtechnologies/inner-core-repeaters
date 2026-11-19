# Copyright 2026 MV2, Munim, Inc.
# Licensed under the Apache License, Version 2.0 (Apache-2.0). See LICENSE.
# Made with AI agents working in MV2's name.
"""Unit tests for snrpilot.py (stdlib + numpy only)."""

import math
import unittest

import numpy as np

import snrpilot as sp


class TestHelpers(unittest.TestCase):
    def test_gc_deg(self):
        self.assertAlmostEqual(sp.gc_deg(0, 0, 0, 90), 90.0, places=6)
        self.assertAlmostEqual(sp.gc_deg(10, 20, -10, -160), 180.0, places=4)
        self.assertAlmostEqual(sp.gc_deg(5, 5, 5, 5), 0.0, places=4)

    def test_median(self):
        self.assertEqual(sp.median([3, 1, 2]), 2)
        self.assertAlmostEqual(sp.median([4, 1, 2, 3]), 2.5)
        self.assertTrue(math.isnan(sp.median([None, float("nan")])))

    def test_rms_and_bounds(self):
        x = np.ones(100) * 2.0
        self.assertAlmostEqual(sp.rms(x, 10, 50), 2.0, places=9)
        self.assertTrue(math.isnan(sp.rms(x, -1, 50)))
        self.assertTrue(math.isnan(sp.rms(x, 10, 101)))

    def test_decimal_year(self):
        self.assertAlmostEqual(sp.decimal_year("2001-01-01T00:00:00.000Z"), 2001.0, places=6)
        self.assertAlmostEqual(sp.decimal_year("2001-07-02T12:00:00Z"), 2001.5, places=3)


class TestWindows(unittest.TestCase):
    def test_window_snrs(self):
        x = np.ones(sp.N)
        i0, i1 = int(sp.idx(sp.SIG_DF[0])), int(sp.idx(sp.SIG_DF[1]))
        x[i0:i1] = 5.0
        df, ab, nl = sp.window_snrs(x, 30.0)
        self.assertAlmostEqual(df, 5.0, places=6)
        self.assertAlmostEqual(nl, 1.0, places=6)
        self.assertAlmostEqual(ab, 1.0, places=6)
        df2, ab2, nl2 = sp.window_snrs(x, float("nan"))
        self.assertTrue(math.isnan(ab2))

    def test_windows_inside_trace(self):
        self.assertGreaterEqual(sp.idx(sp.NULL[0]), 0)
        self.assertLessEqual(sp.NULL[1], sp.NOISE[0])
        self.assertLessEqual(sp.NOISE[1], sp.SIG_DF[0])
        self.assertLess(sp.idx(80.0 + sp.SIG_AB[1]), sp.N)


class TestSelection(unittest.TestCase):
    def _cat(self):
        rows = []
        for i in range(40):
            rows.append(dict(id="a%02d" % i, time="20%02d-01-01T00:00:00Z" % (i % 30), mag=4.6 + 0.05 * (i % 30),
                             t=2000.0 + (i % 30)))
        return rows

    def test_pick_events_caps_and_spreads(self):
        cat = self._cat()
        ev = sp.pick_events(cat, bins=((4.5, 5.0), (5.0, 9.9)), per_bin=5)
        b1 = [e for e in ev if e["bin"] == "4.5-5.0"]
        b2 = [e for e in ev if e["bin"] == "5.0-9.9"]
        self.assertEqual(len(b1), 5)
        self.assertEqual(len(b2), 5)
        times = [e["time"] for e in b1]
        self.assertEqual(times, sorted(times))
        self.assertEqual(sp.pick_events(cat, bins=((4.5, 5.0),), per_bin=5), b1)

    def test_pick_events_keeps_small_bins(self):
        cat = [dict(id="x", time="2001-01-01T00:00:00Z", mag=6.2, t=2001.0)]
        ev = sp.pick_events(cat)
        self.assertEqual(len(ev), 1)
        self.assertEqual(ev[0]["bin"], "6.0-9.9")

    def test_bin_counts(self):
        cat = [dict(mag=5.6, t=2000.0), dict(mag=5.7, t=2010.0), dict(mag=4.7, t=2009.0)]
        out = {b["bin"]: b for b in sp.bin_counts(cat)}
        self.assertEqual(out["5.5-6.0"]["n"], 2)
        self.assertEqual(out["5.5-6.0"]["before_turn"], 1)
        self.assertEqual(out["4.5-5.0"]["after_turn"], 1)

    def test_choose_loc(self):
        chans = [("10", "2000-01-01T00:00:00", "9999", 1.0, 2.0),
                 ("00", "2005-01-01T00:00:00", "2010-01-01T00:00:00", 3.0, 4.0)]
        self.assertEqual(sp.choose_loc(chans, "2006-05-05T00:00:00.000Z"), ("00", (3.0, 4.0)))
        self.assertEqual(sp.choose_loc(chans, "2012-05-05T00:00:00Z"), ("10", (1.0, 2.0)))
        self.assertEqual(sp.choose_loc(chans, "1999-05-05T00:00:00Z"), (None, None))


class TestFits(unittest.TestCase):
    def _rows(self, slope=0.8):
        rows = []
        for sta, off in (("A", 0.0), ("B", 0.5)):
            for m in (4.6, 5.1, 5.6, 6.1):
                rows.append(dict(station=sta, mag=m, bin=None, snr=10 ** (off - 4.0 * slope + slope * m)))
        return rows

    def test_fe_slope(self):
        s, n = sp.fe_slope(self._rows(0.8), "snr")
        self.assertAlmostEqual(s, 0.8, places=6)
        self.assertEqual(n, 8)

    def test_station_fit_mreq(self):
        fit = sp.station_fit(self._rows(1.0), "snr", req=3.0)
        # log10 snr = -4 + m at A -> m_req = 4 + log10 3
        self.assertAlmostEqual(fit["A"]["m_req"], 4.0 + math.log10(3.0), places=6)
        self.assertAlmostEqual(fit["B"]["m_req"], 3.5 + math.log10(3.0), places=6)

    def test_bin_table_and_verdict(self):
        rows = []
        for sta in ("A", "B", "C", "D"):
            for k in range(3):
                rows.append(dict(station=sta, bin="4.5-5.0", snr=1.2))
                rows.append(dict(station=sta, bin="5.5-6.0", snr=4.0 if sta != "D" else 2.0))
        tab = sp.bin_table(rows, "snr")
        v = sp.verdict(tab)
        self.assertEqual(v["first_bin"], "5.5-6.0")
        self.assertEqual(v["verdict"], "testable-from-M5.5")
        self.assertEqual(v["stations_passing"]["5.5-6.0"], ["A", "B", "C"])
        self.assertEqual(v["stations_passing"]["4.5-5.0"], [])
        none = sp.verdict(sp.bin_table([r for r in rows if r["bin"] == "4.5-5.0"], "snr"))
        self.assertEqual(none["verdict"], "none")

    def test_jsonable(self):
        out = sp.jsonable({"a": float("nan"), "b": np.float64(2.0), "c": (1, 2)})
        self.assertEqual(out, {"a": None, "b": 2.0, "c": [1, 2]})


if __name__ == "__main__":
    unittest.main()

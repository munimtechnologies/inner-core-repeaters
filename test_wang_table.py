"""Offline tests for wang_table.py (no network, no openpyxl needed).

Part of inner-core-repeaters by MV2, Munim, Inc. (made with AI agents). Apache-2.0.
"""
import datetime as dt
import unittest

import wang_table as W

UTC = dt.timezone.utc


class TestTimeText(unittest.TestCase):
    def test_iso(self):
        t, h = W.parse_time_text("1993-05-03T12:34:56.5")
        self.assertTrue(h)
        self.assertEqual(t, dt.datetime(1993, 5, 3, 12, 34, 56, 500000, tzinfo=UTC))

    def test_slash_space(self):
        t, h = W.parse_time_text("2003/11/21 07:05:12")
        self.assertTrue(h)
        self.assertEqual((t.year, t.month, t.day, t.hour, t.minute, t.second), (2003, 11, 21, 7, 5, 12))

    def test_date_only(self):
        t, h = W.parse_time_text("2009.06.30")
        self.assertFalse(h)
        self.assertEqual(t.date(), dt.date(2009, 6, 30))

    def test_compact(self):
        t, h = W.parse_time_text("ev_20210812_183544")
        self.assertTrue(h)
        self.assertEqual((t.year, t.month, t.day, t.hour, t.minute, t.second), (2021, 8, 12, 18, 35, 44))

    def test_julian(self):
        t, h = W.parse_time_text("1997.032 05:06:07")
        self.assertTrue(h)
        self.assertEqual((t.month, t.day, t.hour), (2, 1, 5))

    def test_garbage(self):
        self.assertEqual(W.parse_time_text("Multiplet 3"), (None, False))
        self.assertEqual(W.parse_time_text(""), (None, False))


class TestRowEvent(unittest.TestCase):
    def test_datetime_cell_and_coords(self):
        row = ["M01", dt.datetime(2001, 2, 3, 4, 5, 6), -56.1, -27.3, 33.0, 5.6]
        hdr = ["multiplet", "origin time", "lat", "lon", "depth", "mb"]
        t, h, lat, lon, mag, lab = W.row_event(row, hdr)
        self.assertTrue(h)
        self.assertEqual(t.hour, 4)
        self.assertEqual((lat, lon, mag, lab), (-56.1, -27.3, 5.6, "M01"))

    def test_date_plus_time_cells(self):
        row = [dt.datetime(2005, 6, 1), dt.time(10, 20, 30), -58.0, 333.5]
        t, h, lat, lon, mag, lab = W.row_event(row, None)
        self.assertTrue(h)
        self.assertEqual((t.hour, t.minute, t.second), (10, 20, 30))
        self.assertAlmostEqual(lon, -26.5)

    def test_split_columns(self):
        row = [1999, 7, 14, 1, 2, 3.5, -60.0, -26.0]
        hdr = ["year", "month", "day", "hour", "min", "sec", "lat", "lon"]
        t, h, lat, lon, mag, lab = W.row_event(row, hdr)
        self.assertTrue(h)
        self.assertEqual((t.year, t.day, t.second), (1999, 14, 3))

    def test_no_event(self):
        t, h, *_ = W.row_event(["Event", "Lat", "Lon"], None)
        self.assertIsNone(t)


class TestMatch(unittest.TestCase):
    def setUp(self):
        base = dt.datetime(2003, 3, 3, 3, 3, 3, tzinfo=UTC)
        self.cat = [dict(id="a", time=base, lat=-56.0, lon=-27.0, mag=5.5),
                    dict(id="b", time=base + dt.timedelta(seconds=40), lat=-59.0, lon=-25.0, mag=4.8),
                    dict(id="c", time=base + dt.timedelta(hours=5), lat=-56.1, lon=-27.1, mag=5.9)]
        self.base = base

    def test_time_match(self):
        r, how = W.match_event(dict(time=self.base + dt.timedelta(seconds=2), has_time=True, lat=None, lon=None),
                               self.cat)
        self.assertEqual(r["id"], "a")

    def test_time_match_prefers_close_location(self):
        ev = dict(time=self.base + dt.timedelta(seconds=30), has_time=True, lat=-59.0, lon=-25.0)
        r, _ = W.match_event(ev, self.cat)
        self.assertEqual(r["id"], "b")

    def test_no_match(self):
        r, how = W.match_event(dict(time=self.base + dt.timedelta(hours=2), has_time=True, lat=None, lon=None),
                               self.cat)
        self.assertIsNone(r)

    def test_date_only(self):
        r, _ = W.match_event(dict(time=self.base.replace(hour=0, minute=0, second=0), has_time=False,
                                  lat=None, lon=None), self.cat)
        self.assertEqual(r["id"], "c")  # largest that day
        r, _ = W.match_event(dict(time=self.base.replace(hour=0), has_time=False, lat=-59.0, lon=-25.0), self.cat)
        self.assertEqual(r["id"], "b")


class TestMisc(unittest.TestCase):
    def test_single_link(self):
        pts = [(-56.0, -27.0), (-56.05, -27.0), (-56.1, -27.0), (-59.0, -25.0)]
        cl = W.single_link(pts, 8.0)
        self.assertEqual(cl[0], cl[1])
        self.assertEqual(cl[1], cl[2])
        self.assertNotEqual(cl[0], cl[3])

    def test_gc_km(self):
        self.assertAlmostEqual(W.gc_km(0, 0, 1, 0), 111.19, places=1)

    def test_v1_cut(self):
        self.assertTrue(W.in_v1_catalog(2001.0, 5.0))
        self.assertFalse(W.in_v1_catalog(2001.0, 4.9))
        self.assertTrue(W.in_v1_catalog(2025.3, 4.7))
        self.assertFalse(W.in_v1_catalog(2025.3, None))


if __name__ == "__main__":
    unittest.main()

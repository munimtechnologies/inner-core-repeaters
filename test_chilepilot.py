# Copyright 2026 MV2, Munim, Inc.
# Licensed under the Apache License, Version 2.0 (Apache-2.0). See LICENSE.
# Made with AI agents working in MV2's name.
"""Unit tests for chilepilot.py (stdlib + numpy only)."""

import math
import unittest

import numpy as np

import chilepilot as cp


class TestHelpers(unittest.TestCase):
    def test_gc_deg_and_km(self):
        self.assertAlmostEqual(cp.gc_deg(0, 0, 0, 90), 90.0, places=6)
        self.assertAlmostEqual(cp.gc_deg(10, 20, -10, -160), 180.0, places=4)
        self.assertAlmostEqual(cp.km_between(0, 0, 0, 1), 6371.0 * math.pi / 180.0, places=6)

    def test_median_rms(self):
        self.assertEqual(cp.median([3, 1, 2]), 2)
        self.assertTrue(math.isnan(cp.median([None, float("nan")])))
        x = np.ones(100) * 2.0
        self.assertAlmostEqual(cp.rms(x, 10, 50), 2.0, places=9)
        self.assertTrue(math.isnan(cp.rms(x, -1, 50)))

    def test_decimal_year(self):
        self.assertAlmostEqual(cp.decimal_year("2001-07-02T12:00:00Z"), 2001.5, places=3)

    def test_windows_inside_trace(self):
        self.assertGreaterEqual(cp.idx(cp.NULL[0]), 0)
        self.assertLessEqual(cp.NULL[1], cp.NOISE[0])
        self.assertLessEqual(cp.NOISE[1], cp.SIG[0])
        self.assertLess(cp.idx(cp.SIG[1]), cp.N)
        self.assertLess(cp.idx(cp.CODA[1]), cp.N)

    def test_geometry_of_box(self):
        # the box centre is inside the accepted distance range for Makanchi
        d = cp.gc_deg(-30.0, -71.0, 46.79, 82.29)
        self.assertTrue(cp.DIST_RANGE[0] <= d <= cp.DIST_RANGE[1])


class TestBeam(unittest.TestCase):
    def _trace(self, rng, sig=0.0):
        x = rng.normal(size=cp.N)
        i0, i1 = int(cp.idx(cp.SIG[0])), int(cp.idx(cp.SIG[1]))
        x[i0:i1] += sig * np.sin(np.arange(i1 - i0) * 0.9)
        return cp.normalise(x)

    def test_window_snrs_constant(self):
        x = np.ones(cp.N)
        i0, i1 = int(cp.idx(cp.SIG[0])), int(cp.idx(cp.SIG[1]))
        x[i0:i1] = 5.0
        s, c, nl = cp.window_snrs(x)
        self.assertAlmostEqual(s, 5.0, places=6)
        self.assertAlmostEqual(nl, 1.0, places=6)

    def test_beam_gain(self):
        rng = np.random.default_rng(1)
        xs = [self._trace(rng, sig=1.0) for _ in range(9)]
        b = cp.beam(xs)
        el = cp.median([cp.window_snrs(x)[0] for x in xs])
        self.assertGreater(cp.window_snrs(b)[0], 1.5 * el)
        self.assertLess(abs(cp.window_snrs(b)[2] - 1.0), 0.5)

    def test_beam_needs_min_el(self):
        rng = np.random.default_rng(2)
        self.assertIsNone(cp.beam([self._trace(rng), None, self._trace(rng)]))
        self.assertIsNotNone(cp.beam([self._trace(rng) for _ in range(3)]))

    def test_normalise_dead(self):
        self.assertIsNone(cp.normalise(np.zeros(cp.N)))


class TestSelection(unittest.TestCase):
    def test_pick_events_sides(self):
        cat = [dict(id="a%02d" % i, time="%04d-01-01T00:00:00Z" % (1995 + i), t=1995.0 + i, mag=5.2)
               for i in range(30)]
        ev = cp.pick_events(cat, bins=((5.0, 5.5),), per_side=4)
        before = [e for e in ev if e["side"] == "before"]
        after = [e for e in ev if e["side"] == "after"]
        self.assertEqual(len(before), 4)
        self.assertEqual(len(after), 4)
        self.assertTrue(all(e["t"] < cp.TURN for e in before))
        self.assertTrue(all(e["t"] >= cp.TURN for e in after))
        self.assertEqual(before[0]["id"], "a00")
        self.assertEqual(after[-1]["id"], "a29")

    def test_bin_counts(self):
        cat = [dict(mag=5.6, t=2000.0), dict(mag=5.7, t=2010.0), dict(mag=4.7, t=2009.0)]
        out = {b["bin"]: b for b in cp.bin_counts(cat)}
        self.assertEqual(out["5.5-6.0"]["before_turn"], 1)
        self.assertEqual(out["4.5-5.0"]["after_turn"], 1)

    def test_close_pairs(self):
        cat = [dict(lat=-30.0, lon=-71.0, depth=30.0, mag=5.6, t=2000.0),
               dict(lat=-30.01, lon=-71.01, depth=32.0, mag=5.7, t=2015.0),   # straddling, M5.5
               dict(lat=-30.02, lon=-71.0, depth=31.0, mag=5.6, t=2003.0),    # same side as first
               dict(lat=-30.0, lon=-71.0, depth=90.0, mag=5.6, t=2016.0),     # depth too different
               dict(lat=-33.0, lon=-71.0, depth=30.0, mag=5.6, t=2016.0),     # too far
               dict(lat=-30.0, lon=-71.0, depth=30.0, mag=4.8, t=2016.0)]     # too small
        r = cp.close_pairs(cat)
        self.assertEqual(r["n_events"], 5)
        self.assertEqual(r["pairs"], 3)
        self.assertEqual(r["straddling"], 2)
        self.assertEqual(r["straddling_m55"], 2)
        self.assertEqual(r["events_in_straddling"], 3)

    def test_choose_chan(self):
        chans = [("SHZ", "", "2000-01-01T00:00:00", "9999", 1.0, 2.0),
                 ("BHZ", "00", "2005-01-01T00:00:00", "2010-01-01T00:00:00", 3.0, 4.0),
                 ("BHZ", "", "2005-01-01T00:00:00", "2010-01-01T00:00:00", 5.0, 6.0)]
        self.assertEqual(cp.choose_chan(chans, "2006-05-05T00:00:00.000Z"), ("BHZ", "", 5.0, 6.0))
        self.assertEqual(cp.choose_chan(chans, "2012-05-05T00:00:00Z"), ("SHZ", "", 1.0, 2.0))
        self.assertIsNone(cp.choose_chan(chans, "1999-05-05T00:00:00Z"))


class TestDecision(unittest.TestCase):
    def _rows(self, med_by_bin, side_ok=True):
        rows = []
        for b, v in med_by_bin.items():
            for k in range(4):
                for side in ("before", "after"):
                    rows.append(dict(bin=b, side=side, status="ok", beam_sig_primary=v if side_ok else None))
        return rows

    def test_array_verdict(self):
        rows = self._rows({"4.5-5.0": 1.5, "5.0-5.5": 3.5, "5.5-6.0": 8.0})
        v = cp.array_verdict(rows, "beam_sig_primary")
        self.assertEqual(v["verdict"], "testable-from-M5.0")
        few = [r for r in rows if r["bin"] == "5.0-5.5"][:3]
        self.assertEqual(cp.array_verdict(few, "beam_sig_primary")["verdict"], "none")

    def test_availability(self):
        rows = [dict(side="before", status="ok", beam_sig_primary=4.0),
                dict(side="before", status="nodata"),
                dict(side="after", status="out-of-range"),
                dict(side="after", status="ok", beam_sig_primary=2.0)]
        a = cp.availability(rows)
        self.assertAlmostEqual(a["before"]["frac"], 0.5)
        self.assertEqual(a["after"]["n"], 1)
        self.assertAlmostEqual(a["after"]["frac"], 1.0)

    def test_decide(self):
        good = {"before": {"frac": 0.8}, "after": {"frac": 0.9}}
        bad = {"before": {"frac": 0.2}, "after": {"frac": 0.9}}
        v = {"A": {"verdict": "testable-from-M5.0"}, "B": {"verdict": "testable-from-M4.5"}}
        d = cp.decide(v, {"A": good, "B": bad}, {"straddling_m55": 0})
        self.assertEqual(d["best_array"], "A")
        self.assertEqual(d["decision"], "build-catalogue")
        v2 = {"A": {"verdict": "testable-from-M5.5"}}
        self.assertEqual(cp.decide(v2, {"A": good}, {"straddling_m55": 40})["decision"], "build-catalogue-m5.5")
        self.assertEqual(cp.decide(v2, {"A": good}, {"straddling_m55": 5})["decision"], "close-chile-path")
        self.assertEqual(cp.decide(v, {"A": bad, "B": bad}, {})["decision"], "close-chile-path")

    def test_beam_gain_frac(self):
        rows = [dict(beam_sig_primary=5.0, el_sig_primary_median=2.0, n_el_ok=6),
                dict(beam_sig_primary=1.0, el_sig_primary_median=2.0, n_el_ok=6),
                dict(beam_sig_primary=9.0, el_sig_primary_median=2.0, n_el_ok=3)]
        g, n = cp.beam_gain_frac(rows)
        self.assertEqual(n, 2)
        self.assertAlmostEqual(g, 0.5)

    def test_jsonable(self):
        out = cp.jsonable({"a": float("nan"), "b": np.float64(2.0), "c": (1, 2)})
        self.assertEqual(out, {"a": None, "b": 2.0, "c": [1, 2]})


if __name__ == "__main__":
    unittest.main()

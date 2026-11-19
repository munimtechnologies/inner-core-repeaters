"""Tests for nea.py (pre-registration 16). Part of inner-core-repeaters by MV2, Munim, Inc.
(made with AI agents). Licence: Apache-2.0. Needs numpy only (no network, no obspy)."""
import math
import unittest

import numpy as np

import nea


def wavelet(kind, fs=nea.FS, dur=4.0, amp=15.0):
    t = np.arange(int(dur * fs)) / fs
    h = np.hanning(len(t))
    if kind == "a":
        w = np.sin(2 * np.pi * 0.9 * t) + 0.6 * np.sin(2 * np.pi * 1.4 * t + 0.5)
    else:
        w = np.sin(2 * np.pi * 1.4 * t) - 0.8 * np.sin(2 * np.pi * 0.7 * t + 1.0)
    w = h * w
    return amp * w / np.max(np.abs(w))


def trace(rng, df_kind, sep_s=20.0, shift_s=0.0):
    x = rng.normal(0.0, 1.0, nea.N)
    i = int(nea.PRE * nea.FS)
    w = wavelet(df_kind, amp=20.0)
    x[i:i + len(w)] += w
    j = i + int(sep_s * nea.FS)
    ab = wavelet("a", amp=25.0)
    x[j:j + len(ab)] += ab
    if shift_s:
        x = nea.fshift(x, shift_s * nea.FS)
    return x


def synth_rows(truth, stations=("S1", "S2"), n_mult=6, seed=0):
    rng = np.random.default_rng(seed)
    base = [1999.0, 2001.0, 2003.0, 2007.0, 2012.0, 2014.0, 2018.0, 2025.0]
    rows = []
    for st in stations:
        for m in range(n_mult):
            ts = [t + 0.3 * m for t in base]
            for a in range(len(ts)):
                for b in range(a + 1, len(ts)):
                    t1, t2 = ts[a], ts[b]
                    if truth == "B":
                        D = nea.d_b(t1, t2) + 0.01 * rng.random()
                    elif truth == "M":
                        D = nea.d_m(t1, t2) + 0.01 * rng.random()
                    else:
                        D = 0.5
                    rows.append(dict(station=st, mult="Z%d" % m, ev1="e%d_%d" % (m, a), ev2="e%d_%d" % (m, b),
                                     t1=t1, t2=t2, D=D, D_ab=0.5, ok=True,
                                     cross_turn=t1 < nea.TURN < t2))
    return rows


class TestPure(unittest.TestCase):
    def test_equiv_and_distances(self):
        self.assertEqual(nea.equiv(2000.0), 2000.0)
        self.assertAlmostEqual(nea.equiv(2026.0), 2008.5 - 17.5 / 2.5)
        self.assertAlmostEqual(nea.d_b(2001.5, 2026.0), 0.0)
        self.assertAlmostEqual(nea.d_m(2001.5, 2026.0), 24.5)

    def test_gc_and_antipode(self):
        self.assertAlmostEqual(nea.gc_deg(0, 0, 0, 90), 90.0, places=6)
        self.assertAlmostEqual(nea.gc_deg(nea.CENTROID[0], nea.CENTROID[1], nea.ANTIPODE[0], nea.ANTIPODE[1]),
                               180.0, places=4)
        self.assertAlmostEqual(nea.ANTIPODE[1], 154.3061, places=4)

    def test_decimal_year(self):
        self.assertAlmostEqual(nea.decimal_year("2001-01-01T00:00:00.000Z"), 2001.0)
        self.assertAlmostEqual(nea.decimal_year("2001-07-02T12:00:00Z"), 2001.5, places=3)

    def test_parse_members(self):
        self.assertEqual(nea.parse_members("usp000a1fh(2000-10-08);us7000ty3q(2026-11-13)"),
                         ["usp000a1fh", "us7000ty3q"])

    def test_choose_loc(self):
        chans = [("00", "2000-01-01T00:00:00", "2010-01-01T00:00:00"), ("10", "2005-01-01T00:00:00", "9999"),
                 ("", "1995-01-01T00:00:00", "2000-01-01T00:00:00")]
        self.assertEqual(nea.choose_loc(chans, "1998-05-05T00:00:00"), "")
        self.assertEqual(nea.choose_loc(chans, "2007-05-05T00:00:00"), "00")
        self.assertEqual(nea.choose_loc(chans, "2015-05-05T00:00:00.5"), "10")
        self.assertIsNone(nea.choose_loc(chans, "1990-01-01T00:00:00"))

    def test_spearman(self):
        self.assertAlmostEqual(nea.spearman([1, 2, 3, 4], [10, 20, 30, 40]), 1.0)
        self.assertAlmostEqual(nea.spearman([1, 2, 3, 4], [4, 3, 2, 1]), -1.0)
        self.assertTrue(math.isnan(nea.spearman([1, 2], [1, 2])))

    def test_block_normalise(self):
        rows = [dict(station="S", mult="Z", D=d) for d in (0.1, 0.4, 0.2, 0.3)]
        rows += [dict(station="S", mult="Y", D=d) for d in (0.1, 0.2, 0.3)]
        out = nea.block_normalise(rows, "D")
        self.assertEqual(len(out), 4)
        self.assertEqual([round(r["u"], 6) for r in out], [0.0, 1.0, round(1 / 3, 6), round(2 / 3, 6)])


class TestDecisions(unittest.TestCase):
    def test_clock_B(self):
        r = nea.clock_test(synth_rows("B"), "D", nboot=200)
        self.assertEqual(r["label"], "B")
        self.assertGreater(r["rho_B"], r["rho_M"])

    def test_clock_M(self):
        r = nea.clock_test(synth_rows("M"), "D", nboot=200)
        self.assertEqual(r["label"], "M")

    def test_clock_untestable(self):
        r = nea.clock_test(synth_rows("M", stations=("S1",), n_mult=1), "D", nboot=50)
        self.assertEqual(r["label"], "untestable")

    def test_calibration(self):
        c = nea.calibration(synth_rows("M"), "D", nboot=200)
        self.assertTrue(c["testable"])
        self.assertTrue(c["holds"])
        c0 = nea.calibration(synth_rows("flat"), "D", nboot=200)
        self.assertTrue(c0["testable"])
        self.assertFalse(c0["holds"])

    def test_headline(self):
        cal_ok = dict(testable=True, holds=True)
        self.assertEqual(nea.headline(dict(testable=False), {}, {}), "untestable")
        self.assertEqual(nea.headline(dict(testable=True, holds=False), dict(label="B"), {}), "no-time-signal")
        self.assertEqual(nea.headline(cal_ok, dict(label="B"), dict(label="inconclusive")), "B-supported")
        self.assertEqual(nea.headline(cal_ok, dict(label="M"), dict(label="M")), "source-confounded:M")
        self.assertEqual(nea.headline(cal_ok, dict(label="inconclusive"), dict(label="B")), "inconclusive")

    def test_new_event_matches(self):
        rows = [dict(station="S", ev1=p, ev2="n", t1=t, t2=2026.0, D=d)
                for p, t, d in (("p1", 2001.5, 0.05), ("p2", 2012.0, 0.3), ("p3", 2024.0, 0.2))]
        r = nea.new_event_matches(rows, "D")
        self.assertEqual(r["n_discriminating"], 1)
        self.assertEqual(r["n_match_B"], 1)
        self.assertEqual(r["cases"][0]["b_pick"], "p1")
        self.assertEqual(r["cases"][0]["m_pick"], "p3")


class TestSignal(unittest.TestCase):
    def test_family(self):
        rng = np.random.default_rng(3)
        traces = {"e1": trace(rng, "a"), "e2": trace(rng, "a", shift_s=0.6), "e3": trace(rng, "a"),
                  "e4": trace(rng, "b")}
        ab = {k: 20.0 for k in traces}
        info, per, al, onsets = nea.measure_family(traces, ab)
        self.assertEqual(info["status"], "ok", info)
        self.assertEqual(info["n_aligned"], 4)
        self.assertAlmostEqual(per["e2"]["lag_s"] - per["e1"]["lag_s"], 0.6, delta=0.06)
        self.assertTrue(-0.5 <= info["df_pick_minus_pred_s"] <= 1.5, info)
        evs = {k: dict(time="%d-01-01T00:00:00Z" % y, t=float(y), mag=5.0)
               for k, y in (("e1", 2000), ("e2", 2005), ("e3", 2015), ("e4", 2020))}
        rows = nea.family_pairs("S.X", "Z0", evs, al, per, onsets)
        self.assertEqual(len(rows), 6)
        d = {(r["ev1"], r["ev2"]): r for r in rows}
        self.assertLess(d[("e1", "e2")]["D"], 0.15)
        self.assertGreater(d[("e1", "e4")]["D"], d[("e1", "e2")]["D"] + 0.15)
        self.assertTrue(all(r["cc_ab"] >= 0.8 for r in rows))
        self.assertTrue(d[("e1", "e3")]["ok"])
        self.assertTrue(d[("e1", "e2")]["cross_turn"] is False and d[("e1", "e3")]["cross_turn"] is True)

    def test_ab_too_close(self):
        rng = np.random.default_rng(4)
        traces = {"e1": trace(rng, "a", sep_s=5.0), "e2": trace(rng, "a", sep_s=5.0)}
        info, per, al, onsets = nea.measure_family(traces, {"e1": 5.0, "e2": 5.0})
        self.assertEqual(info["status"], "ab-too-close")
        self.assertIsNone(onsets)

    def test_measure_pair_identical(self):
        # Identical traces. lag_scan is a windowed (non-circular) correlation: the windows at lags
        # -1 and +1 take in different edge samples, so the parabolic peak can sit ~1e-3 sample off
        # zero, in the PKPab reference lag and then in the DF lag. In v0.1.0's first test run this
        # gave dt_df_s = -4.7e-5 s (0.001 sample). The test bound is 0.005 s = 0.1 sample, far
        # below the 0.05 s sample interval and the 0.4 s DF search range.
        rng = np.random.default_rng(5)
        x = trace(rng, "a")
        i = int(nea.PRE * nea.FS)
        r = nea.measure_pair(x, x.copy(), i, i + int(20 * nea.FS))
        self.assertAlmostEqual(r["cc_ab"], 1.0, places=6)
        self.assertAlmostEqual(r["cc_df_best"], 1.0, places=6)
        self.assertLess(abs(r["lag_ref_s"]), 0.005)
        self.assertLess(abs(r["dt_df_s"]), 0.005)

    def test_measure_pair_shift(self):
        # A pure 0.25 s delay of the whole trace is absorbed by the PKPab alignment: the DF
        # residual delay is ~0 and the DF correlation stays ~1.
        rng = np.random.default_rng(6)
        x = trace(rng, "a")
        i = int(nea.PRE * nea.FS)
        r = nea.measure_pair(x, nea.fshift(x, 0.25 * nea.FS), i, i + int(20 * nea.FS))
        self.assertAlmostEqual(r["lag_ref_s"], 0.25, delta=0.02)
        self.assertLess(abs(r["dt_df_s"]), 0.02)
        self.assertGreater(r["cc_df_best"], 0.98)


if __name__ == "__main__":
    unittest.main()

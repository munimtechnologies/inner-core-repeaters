"""Synthetic tests for cola.py (numpy; the trace-processing test also needs ObsPy).

Part of inner-core-repeaters by MV2, Munim, Inc. (made with AI agents). Apache-2.0.
"""
import math
import unittest

try:
    import numpy as np
    import bcdf
    import cola
    import model
    from phases import Phases
except ImportError:  # pragma: no cover
    np = None

try:
    import obspy  # noqa: F401
    HAVE_OBSPY = True
except ImportError:  # pragma: no cover
    HAVE_OBSPY = False


def ricker(t, f=1.5):
    a = (np.pi * f * t) ** 2
    return (1 - 2 * a) * np.exp(-a)


@unittest.skipIf(np is None, "numpy missing")
class TestChannels(unittest.TestCase):
    def test_epochs(self):
        self.assertEqual(cola.channel_for("1985-06-01T12:00:00.000Z"), ("DW", "COL", "", "SHZ"))
        self.assertEqual(cola.channel_for("1991-03-04T00:00:01.000Z"), ("IU", "COL", "", "BHZ"))
        self.assertIsNone(cola.channel_for("1996-05-20T00:00:00.000Z"))  # gap between COL and COLA
        self.assertEqual(cola.channel_for("1997-01-01T00:00:00Z"), ("IU", "COLA", "", "BHZ"))
        self.assertEqual(cola.channel_for("2025-09-01T03:04:05.6Z"), ("IU", "COLA", "00", "BHZ"))
        self.assertEqual(cola.channel_for("2025-09-01T03:04:05.6Z", cola.SECONDARY), ("IU", "COLA", "10", "BHZ"))
        self.assertIsNone(cola.channel_for("1981-12-31T00:00:00Z"))
        self.assertEqual(cola.chan_id(("IU", "COLA", "00", "BHZ")), "IU.COLA.00.BHZ")

    def test_parse_catalog(self):
        rows = [dict(id="a", time="2025-01-02T00:00:00.000Z", latitude="-57", longitude="-26", depth="", mag="5.1"),
                dict(id="b", time="1990-01-02T00:00:00.000Z", latitude="-58", longitude="-25", depth="33", mag=""),
                dict(id="c", time="bad", latitude="x", longitude="-25", depth="33", mag="5")]
        evs = cola.parse_catalog(rows)
        self.assertEqual([e["id"] for e in evs], ["b", "a"])
        self.assertEqual(evs[1]["depth"], 10.0)
        self.assertTrue(math.isnan(evs[0]["mag"]))


@unittest.skipIf(np is None, "numpy missing")
class TestStats(unittest.TestCase):
    def rows(self, truth, b=0.02, seed=0, noise=0.01):
        rng = np.random.default_rng(seed)
        out = []
        ts = list(rng.uniform(1982.0, 2024.95, 160)) + list(rng.uniform(2025.0, 2026.85, 14))
        for t in ts:
            lat, lon, z = rng.uniform(-60, -55), rng.uniform(-29, -24), rng.uniform(10, 80)
            g = float(bcdf.time_term([t], truth)[0])
            res = 0.5 + 0.05 * (lat + 57.5) - 0.03 * (lon + 26.5) + b * g + noise * rng.standard_normal()
            out.append(dict(t=t, lat=lat, lon=lon, depth=z, res=res, new=t >= model.NEW_CUT, ok=True,
                            res_sec=res + 0.01 * rng.standard_normal(), cc_sec=0.95))
        return out

    def split(self, rr):
        return [r for r in rr if not r["new"]], [r for r in rr if r["new"]]

    def test_bk_truth(self):
        ins, outs = self.split(self.rows("Bk"))
        r = cola.compare(ins, outs, "L", "Bk", nboot=300)
        self.assertTrue(r["decisive"])
        self.assertEqual(r["verdict"], "Bk")
        r2 = cola.compare(ins, outs, "P", "Bk", nboot=300)
        self.assertEqual(r2["verdict"], "Bk")

    def test_pause_truth(self):
        ins, outs = self.split(self.rows("P", seed=4))
        self.assertEqual(cola.compare(ins, outs, "L", "P", nboot=300)["verdict"], "P")
        self.assertEqual(cola.compare(ins, outs, "P", "Bk", nboot=300)["verdict"], "P")

    def test_l_truth(self):
        ins, outs = self.split(self.rows("L", seed=2))
        self.assertEqual(cola.compare(ins, outs, "L", "Bk", nboot=300)["verdict"], "L")

    def test_not_decisive_when_flat(self):
        ins, outs = self.split(self.rows("L", b=0.0, seed=3, noise=0.15))
        r = cola.compare(ins, outs, "L", "Bk", nboot=300)
        self.assertTrue(r["testable"])
        self.assertLess(abs(r["pred_gap"]), 0.1)  # no time signal: the two forecasts barely differ

    def test_not_testable(self):
        rr = self.rows("L")[:10]
        self.assertFalse(cola.compare(rr, rr[:2], "L", "Bk")["testable"])

    def test_calibration(self):
        ins, _ = self.split(self.rows("L", seed=5))
        c = cola.calibration_c1(ins, nboot=200)
        self.assertTrue(c["testable"])
        self.assertTrue(c["holds"])
        self.assertAlmostEqual(c["b_s_per_yr"], 0.02, delta=0.004)
        neg = [dict(r, res=-r["res"]) for r in ins]
        self.assertFalse(cola.calibration_c1(neg, nboot=200)["holds"])

    def test_screen_and_run_tests(self):
        rr = self.rows("Bk", seed=6)
        rr[0] = dict(rr[0], res=rr[0]["res"] + 2.0)  # a cycle skip
        rr[1] = dict(rr[1], ok=False)
        use = cola.screen(rr)
        self.assertEqual(len(use), len(rr) - 2)
        out = cola.run_tests(use)
        self.assertEqual(out["Q1_in_sample"]["best"], "Bk")
        self.assertEqual(out["Q2_L_vs_Bk"]["verdict"], "Bk")
        self.assertIn("C1_calibration", out)

    def test_secondary(self):
        s = cola.secondary_consistency(self.rows("L", seed=7))
        self.assertGreater(s["n"], 100)
        self.assertLess(abs(s["median_diff"]), 0.01)

    def test_pause_term(self):
        g = bcdf.time_term([2000.0, 2008.5, 2020.0], "P")
        self.assertEqual(list(g), [-8.5, 0.0, 0.0])


def fake_grid():
    dists = np.round(np.arange(140.0, 160.0 + 1e-9, 0.1), 2)
    depths = np.array([0.0, 100.0])
    df = np.full((2, len(dists)), 1200.0)
    return dict(dists=dists, depths=depths, df=df, bc=df + 2.0, ab=df + 5.0, p_df=np.full_like(df, 1.9))


@unittest.skipIf(np is None or not HAVE_OBSPY, "numpy or obspy missing")
class TestTrace(unittest.TestCase):
    def test_process_and_measure(self):
        from obspy import Stream, Trace, UTCDateTime

        ph = Phases(fake_grid())
        ev = dict(id="x", time="2025-09-01T03:04:05.000Z", t=2025.67, lat=-57.0, lon=-26.0, depth=20.0, mag=5.5)
        o = UTCDateTime(ev["time"])
        sr = 20.0
        start = o + 1200.0 - cola.PRE - 30
        n = int(140 * sr)
        tt = np.arange(n) / sr + (start - o) - 1200.0  # time relative to ak135 DF
        rng = np.random.default_rng(0)
        x = 0.02 * rng.standard_normal(n)
        df_late, bc_df = 0.3, 2.3
        x += ricker(tt - df_late) + 1.4 * ricker(tt - df_late - bc_df) + 0.6 * ricker(tt - df_late - 5.0, 1.2)
        tr = Trace(data=(x * 1e4).astype(np.int32))
        tr.stats.sampling_rate = sr
        tr.stats.starttime = start
        tr.stats.network, tr.stats.station, tr.stats.location, tr.stats.channel = "IU", "COLA", "00", "BHZ"
        st = Stream([tr])
        d, why = cola.process_trace(st, ev, ph, cola.COORDS[("IU", "COLA")], cola.BAND)
        self.assertEqual(why, "ok")
        self.assertGreater(d["dist"], 145)
        self.assertLess(d["dist"], 155)
        self.assertAlmostEqual(d["dt_bc"], 2.0, places=6)
        da, _ = cola.process_trace(st, ev, ph, cola.COORDS[("IU", "COLA")], cola.BAND_ALT)
        row = cola.measure_row(ev, d, ("IU", "COLA", "00", "BHZ"), da["x"])
        self.assertIsNotNone(row)
        self.assertTrue(row["ok"])
        self.assertAlmostEqual(row["res"], 0.3, delta=0.05)
        self.assertAlmostEqual(row["res_alt"], 0.3, delta=0.08)
        self.assertTrue(row["new"])

    def test_gappy_rejected(self):
        from obspy import Stream, Trace, UTCDateTime

        ph = Phases(fake_grid())
        ev = dict(id="y", time="2010-01-01T00:00:00Z", t=2010.0, lat=-57.0, lon=-26.0, depth=20.0, mag=5.5)
        o = UTCDateTime(ev["time"])
        a = Trace(data=np.ones(1000, dtype=np.int32))
        a.stats.sampling_rate, a.stats.starttime = 20.0, o + 1130
        b = Trace(data=np.ones(1000, dtype=np.int32))
        b.stats.sampling_rate, b.stats.starttime = 20.0, o + 1200
        d, why = cola.process_trace(Stream([a, b]), ev, ph, cola.COORDS[("IU", "COLA")])
        self.assertIsNone(d)


if __name__ == "__main__":
    unittest.main()

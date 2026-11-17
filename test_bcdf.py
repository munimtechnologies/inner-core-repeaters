"""Synthetic tests for bcdf.py (numpy only; skipped if numpy is missing).

Part of inner-core-repeaters by MV2, Munim, Inc. (made with AI agents). Apache-2.0.
"""
import unittest

try:
    import numpy as np
    import bcdf
    import model
except ImportError:  # pragma: no cover
    np = None


def ricker(t, f=1.5):
    a = (np.pi * f * t) ** 2
    return (1 - 2 * a) * np.exp(-a)


@unittest.skipIf(np is None, "numpy missing")
class TestMeasure(unittest.TestCase):
    def make(self, df_early=0.5, bcdf_true=1.8, seed=1):
        fs, pre = bcdf.FS, bcdf.PRE
        n = int(80 * fs)
        t = np.arange(n) / fs - pre  # 0 = predicted DF
        rng = np.random.default_rng(seed)
        x = 0.03 * rng.standard_normal(n)
        tdf = -df_early
        x += ricker(t - tdf)
        x += 1.3 * ricker(t - tdf - bcdf_true)
        # AB: Hilbert-like (90-degree) pulse later
        x += 0.8 * np.gradient(ricker(t - tdf - 4.6)) / 4.0
        return x

    def test_bcdf_recovered(self):
        x = self.make()
        r = bcdf.measure_bcdf(x, dt_bc_pred=1.6, dt_ab_pred=4.5)
        self.assertIsNotNone(r)
        self.assertAlmostEqual(r["dt_obs_s"], 1.8, delta=0.03)
        self.assertGreater(r["cc"], 0.9)
        self.assertLess(abs(r["onset_s"] + 0.5 + 0.4), 0.4)  # onset a little before the pulse centre

    def test_too_close_not_measured(self):
        x = self.make(bcdf_true=0.9)
        self.assertIsNone(bcdf.measure_bcdf(x, dt_bc_pred=0.8, dt_ab_pred=3.0))

    def test_elements(self):
        x = self.make()
        X = np.vstack([self.make(seed=s) for s in range(2, 9)])
        r = bcdf.measure_bcdf(x, 1.6, 4.5)
        v, n = bcdf.measure_elements(X, r, 1.6, 4.5)
        self.assertEqual(n, 7)
        self.assertAlmostEqual(v, 1.8, delta=0.03)


@unittest.skipIf(np is None, "numpy missing")
class TestFit(unittest.TestCase):
    def rows(self, truth, seed=0):
        rng = np.random.default_rng(seed)
        out = []
        ts = list(rng.uniform(1995.0, 2024.95, 120)) + list(rng.uniform(2025.0, 2026.85, 14))
        for t in ts:
            lat, lon, z = rng.uniform(-60, -55), rng.uniform(-29, -24), rng.uniform(10, 80)
            g = (model.equiv_date(t) if truth == "Bk" else t) - model.T_TURN
            res = 0.4 + 0.05 * (lat + 57.5) - 0.03 * (lon + 26.5) + 0.02 * g + 0.01 * rng.standard_normal()
            out.append(dict(t=t, lat=lat, lon=lon, depth=z, res=res, new=t >= model.NEW_CUT))
        return out

    def test_backtracking_truth(self):
        rr = self.rows("Bk")
        res = bcdf.out_of_sample([r for r in rr if not r["new"]], [r for r in rr if r["new"]], nboot=300)
        self.assertTrue(res["testable"])
        self.assertTrue(res["decisive"])
        self.assertEqual(res["verdict"], "Bk")
        self.assertAlmostEqual(res["coef_Bk"]["b"], 0.02, delta=0.003)

    def test_monotone_truth(self):
        rr = self.rows("L", seed=3)
        res = bcdf.out_of_sample([r for r in rr if not r["new"]], [r for r in rr if r["new"]], nboot=300)
        self.assertEqual(res["verdict"], "L")

    def test_scan_recovers_turn(self):
        rr = [r for r in self.rows("Bk", seed=5) if not r["new"]]
        s = bcdf.scan_turn(rr)
        self.assertLessEqual(abs(s["best_T"] - 2008.5), 1.0)

    def test_not_testable(self):
        rr = self.rows("L")[:10]
        self.assertFalse(bcdf.out_of_sample(rr, rr[:2])["testable"])

    def test_spatial_correction_keeps_time_signal(self):
        rr = self.rows("L", seed=7)
        f = bcdf.fit(rr, "L")
        c = bcdf.spatially_corrected(f, rr)
        t = np.array([r["t"] for r in rr])
        slope = np.polyfit(t, c, 1)[0]
        self.assertAlmostEqual(slope, 0.02, delta=0.003)


if __name__ == "__main__":
    unittest.main()

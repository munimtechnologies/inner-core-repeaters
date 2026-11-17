"""Tests for xcorr.py (synthetic signals only, no network). Apache-2.0."""
import unittest

try:
    import numpy as np

    from xcorr import best_lag, cc, fshift, lag_scan, measure_pair, pick_onset, sta_lta

    HAVE = True
except ImportError:  # the publish lane may lack numpy; the compute job runs these for real
    HAVE = False

FS = 40.0


def wavelet(n, t0, f=1.5, fs=FS):
    t = np.arange(n) / fs
    return np.exp(-((t - t0) * f * 1.5) ** 2) * np.sin(2 * np.pi * f * (t - t0))


def causal_wavelet(n, t0, f=1.5, fs=FS):
    t = np.arange(n) / fs - t0
    w = np.where(t >= 0, np.sin(2 * np.pi * f * t) * np.exp(-t / 0.8), 0.0)
    return w


@unittest.skipUnless(HAVE, "numpy not installed")
class TestXcorr(unittest.TestCase):
    def test_cc_identical_and_negated(self):
        x = wavelet(400, 5.0)
        self.assertAlmostEqual(cc(x, x), 1.0, places=9)
        self.assertAlmostEqual(cc(x, -x), -1.0, places=9)

    def test_fshift_integer_exact(self):
        x = wavelet(400, 5.0)
        y = fshift(x, 7)
        self.assertTrue(np.allclose(y[7:], x[:-7], atol=1e-9))

    def test_fshift_roundtrip_fractional(self):
        x = wavelet(800, 10.0)
        y = fshift(fshift(x, 3.37), -3.37)
        self.assertTrue(np.allclose(y[100:700], x[100:700], atol=1e-6))

    def test_lag_scan_subsample(self):
        x = wavelet(800, 8.0)
        y = fshift(x, 3.3)
        lags, c = lag_scan(x, y, 200, 500, 20)
        l, cpk = best_lag(lags, c)
        self.assertLess(abs(l - 3.3), 0.1)
        self.assertGreater(cpk, 0.99)

    def test_lag_scan_negative_and_edges(self):
        x = wavelet(800, 8.0)
        y = fshift(x, -12.0)
        lags, c = lag_scan(x, y, 200, 500, 30)
        l, _ = best_lag(lags, c)
        self.assertLess(abs(l + 12.0), 0.1)
        # window near the start: impossible lags are NaN, not errors
        lags, c = lag_scan(x, y, 5, 100, 30)
        self.assertTrue(np.isnan(c[0]))

    def test_measure_pair_recovers_df_delay(self):
        n = 1200
        df = wavelet(n, 8.0, f=1.5)
        ref = 1.5 * wavelet(n, 12.0, f=1.2)
        x = df + ref
        # later event: DF 0.05 s late relative to its reference phases, whole record 1.0 s late
        y = fshift(wavelet(n, 8.05, f=1.5) + ref, 40.0)
        r = measure_pair(x, y, (int(10.5 * FS), int(15 * FS)), (int(7 * FS), int(9.5 * FS)), FS,
                         max_lag_ref_s=0.5, max_lag_df_s=0.25, lag0=38.0)
        self.assertIsNotNone(r)
        self.assertLess(abs(r["lag_ref_s"] - 1.0), 0.01)
        self.assertGreater(r["cc_ref"], 0.99)
        self.assertLess(abs(r["dt_df_s"] - 0.05), 0.01)
        self.assertGreater(r["cc_df_best"], 0.99)
        self.assertLess(r["cc_df0"], r["cc_df_best"])

    def test_changed_df_lowers_similarity_only_in_df(self):
        n = 1200
        ref = 1.5 * wavelet(n, 12.0, f=1.2)
        x = wavelet(n, 8.0, f=1.5) + ref
        y = wavelet(n, 8.3, f=1.5) * 0.6 + ref  # DF moved by 0.3 s (about half a cycle)
        r = measure_pair(x, y, (int(10.5 * FS), int(15 * FS)), (int(7 * FS), int(9.5 * FS)), FS)
        self.assertGreater(r["cc_ref"], 0.99)
        self.assertLess(r["cc_df0"], 0.8)

    def test_sta_lta_and_pick(self):
        rng = np.random.default_rng(3)
        n = 3200
        x = 0.05 * rng.standard_normal(n) + causal_wavelet(n, 40.0)
        r = sta_lta(x, 12, 240)
        self.assertEqual(len(r), n)
        trig = pick_onset(x, FS, center=int(40 * FS), search_s=6.0)
        self.assertIsNotNone(trig)
        self.assertGreaterEqual(trig, int(40.0 * FS) - 12)
        self.assertLessEqual(trig, int(40.25 * FS))
        raw = pick_onset(x, FS, center=int(40 * FS), search_s=6.0, refine=False)
        self.assertGreaterEqual(raw, trig)

    def test_pick_none_on_noise(self):
        rng = np.random.default_rng(4)
        x = rng.standard_normal(3200)
        self.assertIsNone(pick_onset(x, FS, center=1600, search_s=6.0, thr=8.0))


if __name__ == "__main__":
    unittest.main()

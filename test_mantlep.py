"""Synthetic tests for mantlep.py (inner-core-repeaters, MV2, Munim, Inc.; made with AI agents).
Licence: Apache-2.0. numpy only; no network.
"""
import os
import unittest

import numpy as np

import mantlep as mp


def band_noise(rng, n, w=3):
    x = rng.standard_normal(n)
    return np.convolve(x, np.ones(w) / w, mode="same")


class TestHelpers(unittest.TestCase):
    def test_ncc_and_shift(self):
        rng = np.random.default_rng(1)
        x = band_noise(rng, 500)
        self.assertAlmostEqual(mp.ncc(x, x), 1.0, places=9)
        self.assertAlmostEqual(mp.ncc(x, -x), -1.0, places=9)
        y = mp.shift(x, 5)
        self.assertTrue(np.allclose(y[:495], x[5:]))
        self.assertTrue(np.all(y[495:] == 0))
        z = mp.shift(x, -5)
        self.assertTrue(np.allclose(z[5:], x[:495]))

    def test_best_lag(self):
        rng = np.random.default_rng(2)
        ref = band_noise(rng, 800)
        x = mp.shift(ref, -7)  # x[i] = ref[i-7]
        lag, c = mp.best_lag(ref, x, 300, 400, 20)
        self.assertEqual(lag, 7)
        self.assertGreater(c, 0.999)

    def test_align_recovers_delay(self):
        rng = np.random.default_rng(3)
        fs, n, center = mp.FS, 1400, 600
        t = np.arange(n)
        env = np.exp(-0.5 * ((t - (center + 20)) / 12.0) ** 2)
        ref = band_noise(rng, n) * env + 0.01 * rng.standard_normal(n)
        for true in (17, -23, 4):
            x = mp.shift(ref, -true) + 0.01 * rng.standard_normal(n)
            lag, c = mp.align(ref, x, fs, center)
            self.assertEqual(lag, true)
            self.assertGreater(c, 0.95)

    def test_pair_D_windows(self):
        rng = np.random.default_rng(4)
        fs, n, onset = mp.FS, 1400, 600
        xa = band_noise(rng, n)
        xb = mp.shift(xa, -2)
        r = mp.pair_D(xa, xb, onset, fs)
        self.assertAlmostEqual(r["lag_s"], 2 / fs)
        self.assertLess(r["D_early"], 1e-6)
        self.assertLess(r["D_late"], 1e-6)
        xc = xb.copy()
        a, b = onset + int(3 * fs) + 2, onset + int(10 * fs) + 2
        xc[a:b] = band_noise(rng, b - a)
        r2 = mp.pair_D(xa, xc, onset, fs)
        self.assertLess(r2["D_early"], 1e-6)
        self.assertGreater(r2["D_late"], 0.5)

    def test_envelope_onset_and_snr(self):
        fs, n, center = mp.FS, 1400, 600
        x = np.zeros(n)
        x[center + 10:center + 30] = np.sin(np.arange(20))
        rng = np.random.default_rng(5)
        x += 0.01 * rng.standard_normal(n)
        on = mp.envelope_onset(x, fs, center)
        self.assertTrue(center + 5 <= on <= center + 12)
        self.assertGreater(mp.snr(x, fs, center), 20)


class TestPartition(unittest.TestCase):
    def block(self, rng, y=(1, 2, 3), n=8, within=0.01, between=0.1):
        D = np.zeros((n, n))
        for i in range(n):
            for j in range(i + 1, n):
                same = (i in y) == (j in y)
                D[i, j] = D[j, i] = (within if same else between) + 0.003 * rng.random()
        return D

    def test_block_split_ranks_first(self):
        rng = np.random.default_rng(6)
        D = self.block(rng)
        r, m, t, vals, combos = mp.partition_rank(D, (1, 2, 3))
        self.assertEqual(m, 56)
        self.assertEqual(r, 1)
        self.assertGreater(t, 0.05)
        self.assertEqual(combos.index((1, 2, 3)), vals.index(max(vals)))

    def test_random_matrix_rank_not_forced(self):
        rng = np.random.default_rng(7)
        D = self.block(rng, y=(0, 4, 5))
        r, m, _, _, _ = mp.partition_rank(D, (1, 2, 3))
        self.assertGreater(r, 1)
        r0, _, _, _, _ = mp.partition_rank(D, (0, 4, 5))
        self.assertEqual(r0, 1)

    def test_decide(self):
        self.assertEqual(mp.decide([1, 2])[0], "untestable")
        self.assertEqual(mp.decide([1, 1, 2, 40, 50])[0], "S")
        self.assertEqual(mp.decide([1, 30, 40, 50, 20])[0], "N")
        self.assertEqual(mp.decide([1, 2, 30, 40, 50, 20, 10, 9, 8])[0], "inconclusive")
        self.assertEqual(mp.decide([10, 20, 30])[0], "N")

    def test_network_combined(self):
        rng = np.random.default_rng(8)
        per = []
        for _ in range(4):
            D = self.block(rng)
            r, m, t, vals, combos = mp.partition_rank(D, (1, 2, 3))
            per.append(vals)
        out = mp.network_combined(per, combos.index((1, 2, 3)))
        self.assertEqual(out["rank"], 1)
        self.assertEqual(out["n_stations"], 4)

    def test_core_controls_file(self):
        path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "z004_core_D.csv")
        out = mp.core_controls(path)
        for k in ("yka_c310", "yka_c13", "ilar_long", "ilar_early"):
            self.assertIn(k, out)
            self.assertTrue(1 <= out[k]["rank"] <= 56)
            self.assertEqual(out[k]["of"], 56)

    def test_groups(self):
        self.assertEqual(mp.ygroup_idx(), (1, 2, 3))
        self.assertEqual(len(mp.IDS), 8)
        self.assertNotIn("us6000ittv", mp.IDS)


if __name__ == "__main__":
    unittest.main()

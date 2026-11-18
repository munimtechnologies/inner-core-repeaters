# Synthetic tests for mantlep2.py (inner-core-repeaters, MV2, Munim, Inc.; made with AI agents).
# Licence: Apache-2.0. No network.
import math
import os
import unittest

import numpy as np

import mantlep2 as m

FS = m.FS
N = int(round((m.PRE + m.POST) * FS))
C = int(round(m.PRE * FS))
HERE = os.path.dirname(os.path.abspath(__file__))


def source(onset, coda_seed=1, n=N, f=1.5):
    """A ringing P pulse starting at sample `onset` followed by a decaying 1-2 Hz-ish coda."""
    t = (np.arange(n) - onset) / FS
    x = np.zeros(n)
    k = t >= 0
    x[k] = 10.0 * np.sin(2 * np.pi * f * t[k]) * np.exp(-((t[k] - 1.2) / 0.9) ** 2)
    rng = np.random.default_rng(coda_seed)
    w = rng.standard_normal(n + 400)
    coda = (np.convolve(w, np.ones(5) / 5, mode="same") - np.convolve(w, np.ones(20) / 20, mode="same"))
    coda = coda[200 - onset + C: 200 - onset + C + n] if 0 <= 200 - onset + C <= 400 else coda[200:200 + n]
    decay = np.where(t >= 1.5, 3.0 * np.exp(-(t - 1.5) / 6.0), 0.0)
    return x + coda * decay


def noisy(x, seed, level=0.1):
    return x + level * np.random.default_rng(1000 + seed).standard_normal(len(x))


def raw_snrs(X):
    return {e: m.window_ratio(x, FS, C, m.RAW_SIG, m.RAW_NOISE, "peak") for e, x in X.items()}


class TestHelpers(unittest.TestCase):
    def test_align_recovers_late_onset(self):
        # Reference P 2.5 s after ak135 (as at QSPA); the other event 1.3 s later still.
        ref = noisy(source(C + 50), 1)
        x = noisy(source(C + 76), 2)
        anchor = m.onset_in(ref, FS, m.ONSET_SEARCH[0], m.ONSET_SEARCH[1], C)
        self.assertTrue(C + 45 <= anchor <= C + 65, anchor)
        lag, c = m.align2(ref, x, FS, anchor)
        self.assertLessEqual(abs(lag - 26), 1)
        self.assertGreater(c, 0.9)

    def test_align_large_shift(self):
        ref = noisy(source(C + 100), 3)     # 5 s late (CASY-like)
        x = noisy(source(C - 10), 4)
        anchor = m.onset_in(ref, FS, m.ONSET_SEARCH[0], m.ONSET_SEARCH[1], C)
        lag, c = m.align2(ref, x, FS, anchor)
        self.assertLessEqual(abs(lag + 110), 1)

    def test_window_ratio_off_grid_is_nan(self):
        x = noisy(source(C), 5)
        self.assertTrue(math.isnan(m.window_ratio(x, FS, 100, m.SIG, m.NOISE)))
        self.assertGreater(m.window_ratio(x, FS, C, m.SIG, m.NOISE), 20)

    def test_window_ratio_ignores_zero_padding(self):
        x = noisy(source(C), 6)
        y = x.copy()
        y[:300] = 0.0
        a = m.window_ratio(x, FS, C, m.SIG, m.NOISE)
        b = m.window_ratio(y, FS, C, m.SIG, m.NOISE)
        self.assertLess(abs(a - b) / a, 0.3)

    def test_spearman(self):
        self.assertAlmostEqual(m.spearman([1, 2, 3, 4, 5], [2, 4, 6, 8, 10]), 1.0)
        self.assertAlmostEqual(m.spearman([1, 2, 3, 4, 5], [5, 4, 3, 2, 1]), -1.0)
        self.assertTrue(math.isnan(m.spearman([1, 2, 3], [1, 2, 3])))
        r = m.rankdata([3.0, 1.0, 3.0, 2.0])
        self.assertEqual(list(r), [2.5, 0.0, 2.5, 1.0])

    def test_parse_station_text(self):
        txt = ("#Network | Station | Latitude | Longitude | Elevation | SiteName | StartTime | EndTime\n"
               "IU|QSPA|-89.9289|144.4382|2847.0|South Pole|2002-11-25T00:00:00|\n"
               "GT|BOSA|-28.614|25.256|1180.0|Boshof|1992-10-01T00:00:00|\n"
               "IU|QSPA|-89.9289|144.4382|2847.0|South Pole|2010-01-01T00:00:00|\n"
               "garbage line\n")
        out = m.parse_station_text(txt)
        self.assertEqual([(a, b) for a, b, _, _ in out], [("GT", "BOSA"), ("IU", "QSPA")])

    def test_gc_deg(self):
        self.assertAlmostEqual(m.gc_deg(m.SRC_LAT, m.SRC_LON, -89.9289, 144.4382), 31.0, delta=0.3)
        self.assertAlmostEqual(m.gc_deg(0, 0, 0, 90), 90.0, places=6)


class TestPartition(unittest.TestCase):
    def test_usable_subset(self):
        g = {e: True for e in m.IDS}
        sub, ok = m.usable_subset(g)
        self.assertTrue(ok)
        self.assertEqual(len(sub), 8)
        g2 = dict(g, usp000bt5x=False, us6000h2wz=False, us7000am0e=False)   # 5 left
        self.assertFalse(m.usable_subset(g2)[1])
        g3 = dict(g, usp000fz4d=False, usp000hute=False)                    # 6 left, 1 Y
        self.assertFalse(m.usable_subset(g3)[1])
        g4 = dict(g, usp000bt5x=False, usp000fz4d=False)                    # 6 left, 2 Y
        self.assertTrue(m.usable_subset(g4)[1])

    def test_rank_subset_threshold(self):
        sub = [e for e in m.IDS if e not in ("usp000bt5x", "usp000fz4d")]   # 6 events, 2 Y
        n = len(m.IDS)
        D = np.full((n, n), 0.1)
        for a in sub:
            for b in sub:
                if (a in m.GROUP_Y) != (b in m.GROUP_Y):
                    D[m.IDS.index(a), m.IDS.index(b)] = 0.5
        r = m.rank_subset(D, sub)
        self.assertEqual((r["rank"], r["of"]), (1, 15))
        self.assertTrue(r["hit"])                       # 1/15 = 0.067 <= 0.07

    def test_rank_full_56(self):
        n = len(m.IDS)
        rng = np.random.default_rng(7)
        D = rng.uniform(0.1, 0.2, (n, n))
        D = (D + D.T) / 2
        r = m.rank_subset(D, list(m.IDS))
        self.assertEqual(r["of"], 56)
        self.assertEqual(r["hit"], r["rank"] <= 3)

    def test_decide(self):
        self.assertEqual(m.decide([True, False])[0], "untestable")
        self.assertEqual(m.decide([True, True, False])[0], "S")
        self.assertEqual(m.decide([True, False, False])[0], "N")
        self.assertEqual(m.decide([True, True, False, False, False])[0], "inconclusive")
        self.assertEqual(m.decide([False] * 6)[0], "N")

    def test_core_file(self):
        core = m.load_core(os.path.join(HERE, "z004_core_D.csv"))
        self.assertEqual(len(core), 28)
        self.assertAlmostEqual(core[("us6000h2wz", "us6000mk65")]["yka_c310"], 0.0017842)


class TestStation(unittest.TestCase):
    lags = [0, 37, -44, 12, -21, 58, -9, 25]       # origin-time scatter, samples

    def build(self, split, drop=()):
        X = {}
        for k, e in enumerate(m.IDS):
            if e in drop:
                continue
            seed = 2 if (split and e in m.GROUP_Y) else 1
            X[e] = noisy(source(C + 50 + self.lags[k], coda_seed=seed), k)
        return X

    def test_identical_sources_aligned_and_alike(self):
        X = self.build(split=False)
        res, A, info, D = m.measure_station(X, raw_snrs(X))
        self.assertEqual(res["status"], "ok")
        self.assertTrue(res["usable"], res.get("gates"))
        # the stack sits in the reference event's time frame
        expected = 2.5 + self.lags[m.IDS.index(res["ref"])] / FS
        self.assertAlmostEqual(res["onset_minus_pred_s"], expected, delta=0.5)
        iu = np.triu_indices(len(m.IDS), 1)
        self.assertLess(np.nanmax(D["D_late"][iu]), 0.15)
        self.assertGreater(np.nanmin(D["cc_direct"][iu]), 0.95)

    def test_split_is_found(self):
        X = self.build(split=True)
        res, A, info, D = m.measure_station(X, raw_snrs(X))
        self.assertTrue(res["usable"])
        self.assertEqual(res["D_late"]["rank"], 1)
        self.assertTrue(res["D_late"]["hit"])
        pc = m.pair_contrast(D["D_late"], res["gated_events"])
        self.assertGreater(pc["ratio"], 3)

    def test_split_with_missing_events(self):
        X = self.build(split=True, drop=("usp000bt5x", "usp000fz4d"))     # 4 X + 2 Y
        res, A, info, D = m.measure_station(X, raw_snrs(X))
        self.assertTrue(res["usable"])
        self.assertEqual(res["D_late"]["of"], 15)
        self.assertEqual(res["D_late"]["rank"], 1)
        self.assertTrue(res["D_late"]["hit"])

    def test_three_plus_three_not_usable(self):
        X = self.build(split=True, drop=("usp000bt5x", "us6000h2wz"))     # 3 X + 3 Y
        res, A, info, D = m.measure_station(X, raw_snrs(X))
        self.assertFalse(res["usable"])
        self.assertEqual(res["D_late"]["of"], 10)                          # distinct splits
        self.assertEqual(res["D_late"]["rank"], 1)
        self.assertFalse(res["D_late"]["hit"])                             # 1/10 > 0.07

    def test_noise_only_event_gated_out(self):
        X = self.build(split=False)
        X["usp000bt5x"] = noisy(np.zeros(N), 99)
        res, A, info, D = m.measure_station(X, raw_snrs(X))
        self.assertFalse(res["gates"]["usp000bt5x"])
        self.assertTrue(res["usable"])
        self.assertEqual(len(res["gated_events"]), 7)


if __name__ == "__main__":
    unittest.main()

"""Synthetic tests for mantlep3.py (pre-registration 13).

Part of inner-core-repeaters by MV2, Munim, Inc. (made with AI agents). Licence: Apache-2.0.
"""
import random
import unittest

import numpy as np

import mantlep2 as m2
import mantlep3 as m3


def _coda(n, t, s):
    lr = np.random.default_rng(s)
    y = np.zeros(n)
    for _ in range(8):
        y += np.sin(2 * np.pi * lr.uniform(1.1, 1.9) * t + lr.uniform(0, 2 * np.pi))
    return y / 2.0


def synth(seed, late_seed, lag=0, noise=0.01):
    """A mantle-P-like record on mantlep2's grid: direct pulse 1 s after ak135 P, then a
    band-limited coda. The coda is shared by all records for its first 2 s and set by late_seed
    after 3 s (linear blend in between): sources that are alike, with a late coda that differs.
    The direct pulse dominates the first second, as in real mantle P at 1-2 Hz.

    (Second version of this test. In the first, the coda was 0.3 and differed from the first
    sample, and the stacked onset came out 1.85 s late; the likely cause, inferred from the code,
    is that the odd-one-out record became the reference and the envelope alignment of the other
    five locked onto a coda beat. On the real records of mantlep2-v1 the stack onset was within
    0.05 s of the reference onset at every usable station.)"""
    n = int(round((m2.PRE + m2.POST) * m2.FS))
    c = int(round(m2.PRE * m2.FS))
    t = (np.arange(n) - c) / m2.FS - 1.0
    on = (t >= 0).astype(float)
    direct = np.sin(2 * np.pi * 1.5 * t) * np.exp(-np.clip(t, 0, None) / 1.0) * on
    w = np.clip(3.0 - t, 0.0, 1.0)            # 1 before t = 2 s, 0 after t = 3 s
    late = (w * _coda(n, t, 1) + (1.0 - w) * _coda(n, t, late_seed)) * np.exp(-np.clip(t, 0, None) / 4.0) * on
    x = direct + 0.2 * late + np.random.default_rng(seed).normal(0, noise, n)
    return m2.shift(x, lag)


def make_rows(fn, n_fam=12, n_ev=6, seed=5):
    rng = random.Random(seed)
    rows = []
    for f in range(n_fam):
        ts = sorted(rng.uniform(2000.0, 2025.0) for _ in range(n_ev))
        ids = ["f%02de%d" % (f, k) for k in range(n_ev)]
        for i in range(n_ev):
            for j in range(i + 1, n_ev):
                t1, t2 = ts[i], ts[j]
                r = dict(array="yka", mult="Y%03d" % f, ev1=ids[i], ev2=ids[j], key=(ids[i], ids[j]),
                         t1=t1, t2=t2, dB=m3.d_B(t1, t2), dM=m3.d_M(t1, t2), cross=m3.is_cross(t1, t2),
                         m0324=False)
                r["D"] = fn(r, rng)
                rows.append(r)
    return rows


class TestClocks(unittest.TestCase):
    def test_equiv(self):
        self.assertEqual(m3.equiv(2003.0), 2003.0)
        self.assertAlmostEqual(m3.equiv(2021.0), 2008.5 - 12.5 / 2.5)
        # 2003.5 and 2008.5 + 2.5 * 5 = 2021.0 are at the same backtracking position
        self.assertAlmostEqual(m3.d_B(2003.5, 2021.0), 0.0)
        self.assertTrue(m3.is_cross(2003.0, 2010.0))
        self.assertFalse(m3.is_cross(2010.0, 2020.0))

    def test_decide13(self):
        self.assertEqual(m3.decide13(0.3, 0.1, 0.95, 0.05, 0.4, 0.2), "B")
        self.assertEqual(m3.decide13(0.1, 0.3, 0.05, 0.95, 0.2, 0.4), "M")
        self.assertEqual(m3.decide13(0.02, 0.01, 0.5, 0.5, 0.12, 0.11), "none")
        self.assertEqual(m3.decide13(0.3, 0.25, 0.7, 0.3, 0.4, 0.4), "inconclusive")
        self.assertEqual(m3.decide13(0.05, 0.05, 0.5, 0.5, 0.25, 0.1), "inconclusive")
        self.assertEqual(m3.decide13(float("nan"), 0.1, 0.5, 0.5, 0.2, 0.2), "untestable")


class TestGrade(unittest.TestCase):
    def scores(self, rows, S=0.1):
        return {r["key"]: dict(S=S, n_st=3, D_late_median=0.01) for r in rows}

    def test_backtracking_truth(self):
        rows = make_rows(lambda r, g: 0.02 * r["dB"] + g.gauss(0, 0.005))
        c = m3.clock(rows, "D", nboot=200)
        self.assertTrue(c["testable"])
        self.assertEqual(c["decision"], "B")

    def test_monotone_truth(self):
        rows = make_rows(lambda r, g: 0.02 * r["dM"] + g.gauss(0, 0.005))
        c = m3.clock(rows, "D", nboot=200)
        self.assertEqual(c["decision"], "M")

    def test_no_clock(self):
        rows = make_rows(lambda r, g: g.random(), n_fam=20)
        c = m3.clock(rows, "D", nboot=200)
        self.assertIn(c["decision"], ("none", "inconclusive"))

    def test_calibration_and_selection(self):
        # core D = source term only; S tracks it -> calibration passes; verified pairs show no clock
        rng = random.Random(9)
        rows = make_rows(lambda r, g: 0.0, n_fam=20)
        sc = {}
        for r in rows:
            s = rng.random()
            r["D"] = s + rng.gauss(0, 0.05)
            sc[r["key"]] = dict(S=s, n_st=3, D_late_median=s)
        out = m3.clock_test(rows, sc, "D", nboot=200)
        self.assertTrue(out["calibration"]["passed"])
        self.assertLess(out["clock"]["n_pairs"], out["n_scored"])
        self.assertNotIn(out["decision"], ("B", "M"))

    def test_calibration_fail(self):
        rows = make_rows(lambda r, g: g.random(), n_fam=20)
        rng = random.Random(3)
        sc = {r["key"]: dict(S=rng.random(), n_st=3, D_late_median=0.1) for r in rows}
        out = m3.clock_test(rows, sc, "D", nboot=200)
        self.assertFalse(out["calibration"]["passed"])
        self.assertEqual(out["decision"], "untestable-calibration")

    def test_m0324_excluded(self):
        rows = make_rows(lambda r, g: g.random())
        rows[0]["m0324"] = True
        sc = self.scores(rows)
        a = m3.clock_test(rows, sc, "D", nboot=50)
        b = m3.clock_test(rows, sc, "D", include_m0324=True, nboot=50)
        self.assertEqual(b["n_pairs"] - a["n_pairs"], 1)

    def test_min_stations(self):
        rows = make_rows(lambda r, g: g.random())
        sc = {r["key"]: dict(S=0.1, n_st=1, D_late_median=0.1) for r in rows}
        out = m3.clock_test(rows, sc, "D", nboot=50)
        self.assertEqual(out["n_scored"], 0)


class TestScores(unittest.TestCase):
    def test_percentiles(self):
        rows = []
        for st in ("A", "B", "C"):
            for k in range(25):
                rows.append(dict(station=st, ev1="e%02d" % k, ev2="z", D_late=k * (1.0 if st != "B" else 10.0)))
        rows.append(dict(station="D", ev1="e00", ev2="z", D_late=0.5))   # too few pairs at D
        rows.append(dict(station="A", ev1="z", ev2="e00", D_late=99.0))  # duplicate pair, ignored
        sc, used = m3.source_scores(rows)
        self.assertEqual(used, ["A", "B", "C"])
        k0 = tuple(sorted(("e00", "z")))
        self.assertEqual(sc[k0]["n_st"], 3)
        self.assertAlmostEqual(sc[k0]["S"], 1.0 / 25.0)
        k24 = tuple(sorted(("e24", "z")))
        self.assertAlmostEqual(sc[k24]["S"], 1.0)

    def test_read_pairs(self):
        text = ("mult,ev1,ev2,time1,time2,t1,t2,mag1,mag2,D,D_c13,D_c310,ok\r\n"
                "Y001,a,b,2001-01-01T00:00:00.000Z,2021-01-01T00:00:00Z,2001.0,2021.0,5,5,0.1,0.05,0.2,True\r\n"
                "Y001,a,c,2001-01-01T00:00:00.000Z,2002-01-01T00:00:00Z,2001.0,2002.0,5,5,0.1,0.05,0.2,False\r\n"
                "Y007,usp000bt5x,us7000am0e,x,y,2003.2,2020.5,5,5,0.1,0.05,0.2,True\r\n")
        rows = m3.read_pairs(text, "yka")
        self.assertEqual(len(rows), 2)
        self.assertTrue(rows[0]["cross"])
        self.assertAlmostEqual(rows[0]["dM"], 20.0)
        self.assertAlmostEqual(rows[0]["dB"], abs(m3.equiv(2021.0) - 2001.0))
        self.assertTrue(rows[1]["m0324"])
        il = m3.read_pairs(text, "ilar")
        self.assertEqual(il[0]["mult"], "I:Y001")

    def test_locate(self):
        cat = m3.read_catalog("time,latitude,longitude,depth,mag,id\n"
                              "2001-01-01T00:00:00.000Z,-58,-26,35,5,a\n"
                              "2002-01-01T00:00:01.000Z,-59,-25,,5,q\n")
        e, how = m3.locate("a", "", cat)
        self.assertEqual(how, "id")
        e, how = m3.locate("b", "2002-01-01T00:00:03.000Z", cat)
        self.assertEqual(how, "time")
        self.assertEqual(e["lat"], -59.0)
        ev, miss = m3.build_events([dict(ev1="a", ev2="b", time1="2001-01-01T00:00:00.000Z",
                                         time2="2002-01-01T00:00:03.000Z")], cat)
        self.assertEqual(ev["b"]["depth"], 33.0)
        self.assertEqual(miss, [])


class TestMeasureFamily(unittest.TestCase):
    def test_alike_and_different(self):
        ids = ["e0", "e1", "e2", "e3", "e4", "e5"]
        lags = [0, 7, -9, 4, -3, 6]
        X = {}
        for k, (e, lag) in enumerate(zip(ids, lags)):
            X[e] = synth(100 + k, late_seed=2 if e == "e5" else 1, lag=lag)
        c = int(round(m2.PRE * m2.FS))
        raw = {e: m2.window_ratio(x, m2.FS, c, m2.RAW_SIG, m2.RAW_NOISE, "peak") for e, x in X.items()}
        need = {tuple(sorted(p)) for p in [("e0", "e1"), ("e2", "e3"), ("e0", "e5"), ("e1", "e5"), ("e4", "zz")]}
        res, gates, pairs = m3.measure_family(X, raw, ids, need)
        self.assertEqual(res["status"], "ok", res)
        self.assertTrue(all(gates[e] for e in ids), gates)
        self.assertLess(pairs[("e0", "e1")]["D_late"], 0.1)
        self.assertLess(pairs[("e2", "e3")]["D_late"], 0.1)
        self.assertGreater(pairs[("e0", "e5")]["D_late"], 0.3)
        self.assertGreater(pairs[("e1", "e5")]["D_late"], 0.3)
        self.assertNotIn(("e4", "zz"), pairs)
        # the direct pulse starts 1 s after ak135 P, moved by the reference's lag (at most 0.45 s)
        self.assertLess(abs(res["onset_minus_pred_s"] - 1.0), 0.8, res)
        self.assertLess(abs(res["onset_minus_pred_s"] - res["ref_onset_minus_pred_s"]), 0.3, res)

    def test_noise_station_is_edge_or_gated(self):
        rng = np.random.default_rng(1)
        n = int(round((m2.PRE + m2.POST) * m2.FS))
        X = {e: rng.normal(0, 1, n) for e in ("a", "b", "c")}
        c = int(round(m2.PRE * m2.FS))
        raw = {e: m2.window_ratio(x, m2.FS, c, m2.RAW_SIG, m2.RAW_NOISE, "peak") for e, x in X.items()}
        res, gates, pairs = m3.measure_family(X, raw, ["a", "b", "c"], {("a", "b"), ("a", "c")})
        self.assertEqual(pairs, {})


if __name__ == "__main__":
    unittest.main()

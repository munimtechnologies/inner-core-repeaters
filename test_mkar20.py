# Copyright 2026 MV2, Munim, Inc.
# Licensed under the Apache License, Version 2.0 (Apache-2.0). See LICENSE.
# Made with AI agents working in MV2's name.
"""Unit tests for mkar20.py (standard library and numpy only)."""

import math
import os
import tempfile
import unittest

import numpy as np

import mkar20 as m


def ev(i, t, lat, lon, depth, mag):
    return dict(id=i, t=t, lat=lat, lon=lon, depth=depth, mag=mag, time="x")


class TestClocks(unittest.TestCase):
    def test_equiv_and_distances(self):
        self.assertAlmostEqual(m.equiv(2000.0), 2000.0)
        self.assertAlmostEqual(m.equiv(2013.5), 2006.5)
        self.assertAlmostEqual(m.d_B(2006.5, 2013.5), 0.0)
        self.assertAlmostEqual(m.d_M(2006.5, 2013.5), 7.0)
        self.assertAlmostEqual(m.d_B(2001.0, 2005.0), 4.0)

    def test_kind_of(self):
        self.assertEqual(m.kind_of(2001.0, 2010.0), "straddle")
        self.assertEqual(m.kind_of(2001.0, 2005.0), "before")
        self.assertEqual(m.kind_of(2009.0, 2015.0), "after")
        self.assertEqual(m.kind_of(2001.0, 2008.5), "straddle")


class TestStats(unittest.TestCase):
    def test_ranks_ties(self):
        self.assertEqual(m.ranks([1.0, 1.0, 2.0]), [1.5, 1.5, 3.0])
        self.assertEqual(m.ranks([3.0, 1.0, 2.0]), [3.0, 1.0, 2.0])

    def test_spearman(self):
        self.assertAlmostEqual(m.spearman([1, 2, 3, 4], [10, 20, 30, 40]), 1.0)
        self.assertAlmostEqual(m.spearman([1, 2, 3, 4], [4, 3, 2, 1]), -1.0)
        self.assertTrue(math.isnan(m.spearman([1, 2], [1, 2])))
        self.assertTrue(math.isnan(m.spearman([1, 1, 1], [1, 2, 3])))

    def test_median_pct(self):
        self.assertEqual(m.median([3, 1, float("nan"), 2]), 2.0)
        self.assertEqual(m.median([4, 1, 2, 3]), 2.5)
        self.assertTrue(math.isnan(m.median([])))
        self.assertEqual(m.pct([1, 2, 3, 4, 5], 0.0), 1)
        self.assertEqual(m.pct([1, 2, 3, 4, 5], 1.0), 5)


class TestPairs(unittest.TestCase):
    def events(self):
        return [
            ev("F", 2001.0, -30.02, -71.0, 25.0, 5.0),
            ev("A", 2005.0, -30.00, -71.0, 30.0, 5.2),
            ev("B", 2012.0, -30.05, -71.0, 35.0, 5.4),     # 5.6 km from A
            ev("C", 2015.0, -30.20, -71.0, 30.0, 5.8),     # 22 km from A; magnitude far from all
            ev("D", 2003.0, -30.00, -71.0, 60.0, 5.2),     # depth too different
            ev("E", 2006.0, -30.00, -71.0, 30.0, 4.8),     # below M5.0
            ev("G", 2016.0, -30.13, -71.0, 36.0, 5.4),     # 8.9 km from B, 14.5 km from A: after pair only
            ev("H", 2002.0, -20.00, -70.0, float("nan"), 5.5),
            ev("I", 2004.0, -20.03, -70.0, 100.0, 5.6),    # missing depth of H passes
        ]

    def test_candidate_pairs(self):
        p = m.candidate_pairs(self.events())
        got = sorted((x["ev1"], x["ev2"], x["kind"]) for x in p)
        self.assertEqual(got, [("A", "B", "straddle"), ("B", "G", "after"), ("F", "A", "before"),
                               ("H", "I", "before")])
        ab = [x for x in p if x["ev1"] == "A"][0]
        self.assertAlmostEqual(ab["sep_km"], 5.56, delta=0.05)
        self.assertAlmostEqual(ab["dB"], abs(m.equiv(2012.0) - 2005.0))
        self.assertAlmostEqual(ab["dM"], 7.0)
        self.assertEqual(ab["key"], ("A", "B"))

    def test_families(self):
        evs = self.events()
        p = m.candidate_pairs(evs)
        fam, rows = m.families(p, {e["id"]: e["t"] for e in evs})
        self.assertEqual(fam, {"F": "F001", "A": "F001", "B": "F001", "H": "F002", "I": "F002"})
        self.assertEqual(sorted((r["ev1"], r["ev2"]) for r in rows), [("A", "B"), ("F", "A"), ("H", "I")])
        ab = [r for r in rows if r["ev1"] == "A"][0]
        self.assertTrue(ab["cross"])
        self.assertFalse(ab["m0324"])
        self.assertEqual(ab["mult"], "F001")

    def test_after_pair_inside_family_is_kept(self):
        evs = self.events() + [ev("J", 2001.5, -30.08, -71.0, 34.0, 5.4)]   # 5.6 km from G: links it in
        p = m.candidate_pairs(evs)
        fam, rows = m.families(p, {e["id"]: e["t"] for e in evs})
        self.assertEqual(fam["G"], fam["A"])
        self.assertIn(("B", "G"), [(r["ev1"], r["ev2"]) for r in rows])


def wav(n, fs, c0, t0, f=1.5, w=0.4):
    t = (np.arange(n) - c0) / fs
    return np.exp(-((t - t0) / w) ** 2) * np.cos(2 * np.pi * f * (t - t0))


def synth(seed=1, n=3000, fs=20.0, c0=1200):
    rng = np.random.default_rng(seed)
    noise = np.convolve(rng.normal(size=n), np.ones(5) / 5.0, mode="same")
    return (0.2 * noise + 5 * wav(n, fs, c0, 0.0) + 8 * wav(n, fs, c0, 1.6) + 6 * wav(n, fs, c0, 4.0)
            + 3 * wav(n, fs, c0, 7.5))


class TestSignal(unittest.TestCase):
    def test_fshift(self):
        n = 400
        x = np.exp(-((np.arange(n) - 100.0) / 5.0) ** 2)
        self.assertEqual(int(np.argmax(m.fshift(x, 10))), 90)
        y = m.fshift(x, 2.5)
        c = float(np.sum(np.arange(n) * y) / np.sum(y))
        self.assertAlmostEqual(c, 97.5, delta=0.01)

    def test_parabolic(self):
        p = 0.3
        f = [-(k - p) ** 2 + 1.0 for k in (-1, 0, 1)]
        self.assertAlmostEqual(m.parabolic(*f), 0.3, places=9)
        self.assertEqual(m.parabolic(0.0, 0.0, 0.0), 0.0)
        self.assertEqual(m.parabolic(1.0, 0.0, 1.0), 0.0)

    def test_lag_scan(self):
        x = synth()
        b = np.roll(x, 7)                       # b[n] = x[n - 7]: b is late by 7 samples
        L, cc, _, _ = m.lag_scan(x, b, 1100, 1400, 20)
        self.assertEqual(L, 7)
        self.assertGreater(cc, 0.999)

    def test_ncc(self):
        a = np.array([1.0, 2.0, 3.0, 4.0])
        self.assertAlmostEqual(m.ncc(a, 2 * a + 1), 1.0)
        self.assertAlmostEqual(m.ncc(a, -a), -1.0)
        self.assertTrue(math.isnan(m.ncc(a, np.ones(4))))

    def test_mkar_pair_same_wave_shifted(self):
        xa = synth()
        xb = m.fshift(xa, -27.4)                # xb[n] = xa[n - 27.4]: late by 1.37 s
        r = m.mkar_pair(xa, xb, 20.0, 1200)
        self.assertIsNotNone(r)
        self.assertAlmostEqual(r["lag_s"], 1.37, delta=0.02)
        self.assertLess(r["D"], 0.01)
        self.assertLess(r["D_c310"], 0.01)
        self.assertGreater(r["cc_direct"], 0.99)
        self.assertLess(abs(r["dt_coda_s"]), 0.03)
        self.assertLess(abs(r["onset_minus_pred_s"]), 1.0)

    def test_mkar_pair_changed_coda(self):
        xa = synth()
        same = m.mkar_pair(xa, m.fshift(xa, -27.4), 20.0, 1200)
        xc = xa + 6 * wav(3000, 20.0, 1200, 3.0, f=1.2)
        diff = m.mkar_pair(xa, m.fshift(xc, -27.4), 20.0, 1200)
        self.assertIsNotNone(diff)
        self.assertGreater(diff["D"], same["D"] + 0.05)
        self.assertAlmostEqual(diff["lag_s"], 1.37, delta=0.05)

    def test_mkar_pair_none(self):
        self.assertIsNone(m.mkar_pair(None, synth(), 20.0, 1200))

    def test_gate_pair(self):
        good = dict(cc_direct=0.9, D=0.05)
        self.assertTrue(m.gate_pair(3.0, 5.0, good))
        self.assertFalse(m.gate_pair(2.9, 5.0, good))
        self.assertFalse(m.gate_pair(3.0, float("nan"), good))
        self.assertFalse(m.gate_pair(5.0, 5.0, dict(cc_direct=0.4, D=0.05)))
        self.assertFalse(m.gate_pair(5.0, 5.0, None))


class TestSecondary(unittest.TestCase):
    def test_positive_control_and_verified(self):
        rows, scores = [], {}
        for f in range(4):
            for j in range(1, 5):
                dM = j + 0.1 * f
                key = ("a%d%d" % (f, j), "b%d%d" % (f, j))
                rows.append(dict(key=key, mult="F%d" % f, kind="before", dM=dM, dB=dM, D=0.01 * dM))
                scores[key] = dict(S=0.2, n_st=3)
        rows.append(dict(key=("x", "y"), mult="F9", kind="before", dM=1.0, dB=1.0, D=0.5))
        scores[("x", "y")] = dict(S=0.9, n_st=3)                       # not verified
        V = m.verified(rows, scores)
        self.assertEqual(len(V), 16)
        pc = m.positive_control(V)
        self.assertEqual(pc["n_pairs"], 16)
        self.assertAlmostEqual(pc["rho"], 1.0)
        self.assertTrue(pc["passed"])
        self.assertFalse(m.positive_control(V[:10])["passed"])

    def test_near_far(self):
        V = [dict(kind="straddle", dB=1.0, D=0.02), dict(kind="straddle", dB=2.0, D=0.04),
             dict(kind="straddle", dB=5.0, D=0.08), dict(kind="before", dB=6.0, D=0.5)]
        nf = m.near_far(V)
        self.assertEqual((nf["n_straddle"], nf["n_near"], nf["n_far"]), (3, 2, 1))
        self.assertAlmostEqual(nf["median_D_near"], 0.03)
        self.assertAlmostEqual(nf["median_D_far"], 0.08)

    def test_stability_flags_step(self):
        rows = []
        for y in range(2001, 2005):
            for c in (0.78, 0.80, 0.82):
                rows.append(dict(t=y + 0.5, beam_snr=5.0, el_beam_cc=c))
        for c in (0.20, 0.22, 0.18):
            rows.append(dict(t=2005.5, beam_snr=5.0, el_beam_cc=c))
        rows.append(dict(t=2006.5, status="nodata"))
        st = m.stability(rows)
        self.assertEqual(st["step_flags"], [2005])
        self.assertEqual(len(st["years"]), 5)

    def test_jsonable_and_csv(self):
        self.assertEqual(m.jsonable({"a": float("nan"), "b": (1, 2), "c": {3}}), {"a": None, "b": [1, 2], "c": [3]})
        with tempfile.TemporaryDirectory() as d:
            p = os.path.join(d, "x.csv")
            m.write_csv(p, [dict(key=("a", "b"), D=0.123456789, S=float("nan")), dict(extra=1)])
            with open(p, encoding="utf-8") as f:
                txt = f.read()
            self.assertIn("a;b,0.12346,", txt)
            self.assertIn("extra", txt.splitlines()[0])


if __name__ == "__main__":
    unittest.main()

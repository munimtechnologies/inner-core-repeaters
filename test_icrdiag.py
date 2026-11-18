"""Tests for icrdiag.py (stdlib only, offline). MV2, Munim, Inc. (made with AI agents). Apache-2.0."""
import csv
import json
import math
import os
import random
import tempfile
import unittest

import icrdiag as d


class TestBasics(unittest.TestCase):
    def test_equiv_and_deltas(self):
        self.assertAlmostEqual(d.equiv(2026.0), 2008.5 - 17.5 / 2.5)
        self.assertAlmostEqual(d.delta_B(2001.5, 2026.0), 0.0)
        self.assertAlmostEqual(d.delta_M(2001.5, 2026.0), 24.5)

    def test_sep(self):
        self.assertAlmostEqual(d.sep_km(-56.0, -27.0, -57.0, -27.0), 111.19, delta=0.1)

    def test_ranks_ties(self):
        self.assertEqual(d.ranks([3.0, 1.0, 3.0, 2.0]), [3.5, 1.0, 3.5, 2.0])

    def test_mw(self):
        a, b = [1, 2, 3, 4, 5, 6], [7, 8, 9, 10, 11, 12]
        self.assertLess(d.mw_less(a, b), 0.01)
        self.assertGreater(d.mw_less(b, a), 0.99)

    def test_ols(self):
        X = [[1.0, x] for x in range(10)]
        y = [1.0 + 2.0 * x for x in range(10)]
        c = d.ols(X, y)
        self.assertAlmostEqual(c[0], 1.0)
        self.assertAlmostEqual(c[1], 2.0)

    def test_binom(self):
        self.assertAlmostEqual(d.binom_tail(5, 5), 1 / 32)
        self.assertEqual(d.binom_tail(0, 4), 1.0)

    def test_epoch_and_locate(self):
        cat = d.read_catalog("time,latitude,longitude,depth,mag,id\n"
                             "2001-05-29T06:22:55.100Z,-56.0,-27.0,30,5.0,aaa\n"
                             "2002-01-01T00:00:00.000Z,-57.0,-26.0,40,5.1,bbb\n")
        e, how = d.locate("aaa", "", cat)
        self.assertEqual(how, "id")
        e, how = d.locate("zzz", "2002-01-01T00:00:02.000Z", cat)
        self.assertEqual(how, "time")
        self.assertEqual(e["lat"], -57.0)
        e, how = d.locate("zzz", "2003-01-01T00:00:00Z", cat)
        self.assertIsNone(e)


def prow(mult, t1, t2, D, snr=20.0, ok=True):
    return dict(mult=mult, ev1="e%.3f" % t1, ev2="e%.3f" % t2, t1=t1, t2=t2, D=D, D_c13=D, D_c310=D,
                ok=ok, dt=abs(t2 - t1), dB=d.delta_B(t1, t2), minsnr=snr, sep_km=10.0, ddep_km=5.0,
                dmag=0.1, seq21=False, same_gen=True)


class TestPartY(unittest.TestCase):
    def simpson_rows(self):
        rows = []
        # family A: high baseline, many near pairs (D 0.9), few far (1.0)
        for k in range(10):
            rows.append(prow("A", 2021.7 + 0.01 * k, 2022.0 + 0.01 * k, 0.9 + 0.001 * k))
        for k in range(2):
            rows.append(prow("A", 2008.0 + 0.3 * k, 2026.0, 1.0))  # dt 18, dB > 5: 'far'
        # family B: low baseline, few near (0.05), many far (0.10-0.19)
        for k in range(2):
            rows.append(prow("B", 2003.0 + 0.1 * k, 2003.5, 0.05))
        for k in range(10):
            rows.append(prow("B", 1995.0 + 0.1 * k, 2008.0 + 0.1 * k, 0.1 + 0.01 * k))
        return rows

    def test_simpson(self):
        rows = self.simpson_rows()
        pooled = d.p4x(rows, "D", nboot=50)
        self.assertFalse(pooled["holds"])
        wf = d.within_family(rows, "D")
        self.assertEqual(wf["n_families"], 2)
        self.assertEqual(wf["n_holds"], 2)
        dm = d.p4x(d.demeaned(rows, "D"), "D_dm", nboot=50)
        self.assertTrue(dm["holds"])

    def test_far_definition(self):
        self.assertTrue(d.is_far(prow("A", 1995.0, 2008.0, 0.1)))
        # 2001 vs 2015: 14 yr apart, but only 4.9 yr apart on the backtracking clock: not 'far'
        self.assertFalse(d.is_far(prow("A", 2001.0, 2015.0, 0.1)))
        # 2001.5 vs 2026: 24.5 yr apart but the same backtracking-equivalent date: not 'far'
        self.assertFalse(d.is_far(prow("A", 2001.5, 2026.0, 0.1)))

    def test_rank_regression_finds_snr(self):
        rng = random.Random(0)
        rows = []
        for k in range(200):
            snr = rng.uniform(5, 100)
            r = prow("F%d" % (k % 20), 2000 + rng.uniform(0, 25), 2000 + rng.uniform(0, 25), 1.0 / snr)
            r.update(minsnr=snr, sep_km=rng.uniform(0, 50), ddep_km=rng.uniform(0, 20), dmag=rng.uniform(0, 1),
                     seq21=rng.random() < 0.3)
            rows.append(r)
        res = d.rank_regression(rows, "D", nboot=30)
        self.assertTrue(res["testable"])
        self.assertLess(res["coef"]["minsnr"]["b"], -0.8)


class TestPartB(unittest.TestCase):
    def test_em2(self):
        x = [3.4 + 0.01 * (k % 11 - 5) for k in range(100)] + [4.0 + 0.01 * (k % 11 - 5) for k in range(60)]
        m = d.em2(x)
        self.assertAlmostEqual(m["mu"][0], 3.4, delta=0.03)
        self.assertAlmostEqual(m["mu"][1], 4.0, delta=0.03)
        self.assertLess(m["bic2"] - m["bic1"], 0)

    def test_nn_agreement(self):
        rng = random.Random(1)
        pts = [(-55.0 - rng.random() * 5, -27.0 + rng.random()) for _ in range(80)]
        lab = [1 if p[0] > -57.5 else 0 for p in pts]
        res = d.nn_agreement(pts, lab, 1, nperm=200)
        self.assertGreater(res["observed"], 0.85)
        self.assertLess(res["p"], 0.02)
        rnd = [rng.randint(0, 1) for _ in pts]
        res2 = d.nn_agreement(pts, rnd, 1, nperm=200)
        self.assertLess(res2["observed"], 0.75)


def write_csv(path, rows, cols):
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(cols)
        for r in rows:
            w.writerow([r.get(c, "") for c in cols])


class TestRunOffline(unittest.TestCase):
    def test_main(self):
        rng = random.Random(2)
        with tempfile.TemporaryDirectory() as tmp:
            cat, pairs, bc = [], [], []
            evs = []
            for k in range(40):
                t = 2000.0 + k * 0.65
                iso = "%04d-06-01T00:00:%02d.000Z" % (int(t), k % 60)
                evs.append(dict(id="ev%02d" % k, t=t, iso=iso, lat=-56 - rng.random(), lon=-27 + rng.random()))
                cat.append(dict(time=iso, latitude=evs[-1]["lat"], longitude=evs[-1]["lon"], depth=30, mag=5, id="ev%02d" % k))
            for i in range(40):
                for j in range(i + 1, min(40, i + 25)):
                    a, b = evs[i], evs[j]
                    pairs.append(dict(mult="Y%03d" % (i % 4), ev1=a["id"], ev2=b["id"], time1=a["iso"], time2=b["iso"],
                                      t1=a["t"], t2=b["t"], mag1=5.0, mag2=5.1, gen1="legacy", gen2="new",
                                      snr_dir1=rng.uniform(5, 50), snr_dir2=rng.uniform(5, 50), ok="True",
                                      D=rng.uniform(0, 0.2), D_c13=rng.uniform(0, 0.2), D_c310=rng.uniform(0, 1),
                                      D_full=rng.uniform(0, 0.2)))
            for k in range(60):
                br = k % 3 == 0
                bc.append(dict(id="b%d" % k, t=1999 + k * 0.45, lat=-56 - rng.random() * 3, lon=-27 + rng.random() * 2,
                               depth=30 + rng.random() * 50, mag=5.0, res=3.4 + (0.6 if br else 0.0) + rng.gauss(0, 0.05),
                               res_el=3.4 + (0.6 if br else 0.0) + rng.gauss(0, 0.05), onset_s=-4.0, cc=0.9,
                               snr_df=6.0, ok="True"))
            pp, cp, bp = (os.path.join(tmp, x) for x in ("p.csv", "c.csv", "b.csv"))
            write_csv(pp, pairs, list(pairs[0].keys()))
            write_csv(cp, cat, list(cat[0].keys()))
            write_csv(bp, bc, list(bc[0].keys()))
            out = os.path.join(tmp, "out")
            d.main(["--out-dir", out, "--pairs", pp, "--catalog", cp, "--bcdf", bp, "--nboot", "20", "--nperm", "50"])
            with open(os.path.join(out, "icrdiag_summary.json")) as f:
                s = json.load(f)
            self.assertNotIn("part_y_error", s)
            self.assertNotIn("part_b_error", s)
            self.assertEqual(s["pairs_located"], len(pairs))
            self.assertIn("within_family", s["part_y"])
            self.assertAlmostEqual(s["part_b"]["mixture"]["gap_s"], 0.6, delta=0.1)
            self.assertTrue(os.path.exists(os.path.join(out, "bcdf_branches.csv")))


if __name__ == "__main__":
    unittest.main()

"""Synthetic tests for ilar_diag.py (skipped if numpy or scipy is missing).

Part of inner-core-repeaters by MV2, Munim, Inc. (made with AI agents). Apache-2.0.
"""
import csv
import os
import tempfile
import unittest

try:
    import numpy as np
    import scipy  # noqa: F401
    import ilar_diag
    import model
except ImportError:  # pragma: no cover
    np = None


def synth(truth="B", n_mult=14, seed=0):
    rng = np.random.default_rng(seed)
    rows = []
    for m in range(n_mult):
        ts = sorted(rng.uniform(2000.0, 2024.9, 5))
        for i in range(len(ts)):
            for j in range(i + 1, len(ts)):
                t1, t2 = ts[i], ts[j]
                g = model.delta_B(t1, t2) if truth == "B" else model.delta_M(t1, t2)
                D = 0.02 + 0.02 * g + 0.01 * rng.standard_normal()
                rows.append(dict(mult="M%03d" % m, ev1="a%d_%d" % (m, i), ev2="a%d_%d" % (m, j), t1=t1, t2=t2,
                                 D=D, Dx=D, cc_ref=0.97, snr=float(rng.uniform(3, 8)), ok=True))
    return rows


@unittest.skipIf(np is None, "numpy/scipy missing")
class TestDiag(unittest.TestCase):
    def test_backtracking_truth(self):
        res = ilar_diag.run(synth("B"), nboot_T=20)
        self.assertGreater(res["J_jackknife"]["full"]["diff"], 0)
        self.assertGreater(res["W_within_multiplet"]["n_B_better"], res["W_within_multiplet"]["n_M_better"])
        self.assertGreater(res["S_snr"]["partial_rho_B"], res["S_snr"]["partial_rho_M"])
        self.assertLessEqual(abs(res["T_turn_bootstrap"]["best"]["best_T"] - 2008.5), 2.0)

    def test_monotone_truth(self):
        res = ilar_diag.run(synth("M", seed=3), nboot_T=20)
        self.assertLess(res["J_jackknife"]["full"]["diff"], 0)
        self.assertGreater(res["W_within_multiplet"]["n_M_better"], res["W_within_multiplet"]["n_B_better"])

    def test_epoch_classes(self):
        self.assertEqual(ilar_diag.epoch_class(dict(t1=2005.0, t2=2010.0)), "pre")
        self.assertEqual(ilar_diag.epoch_class(dict(t1=2005.0, t2=2020.0)), "cross")
        self.assertEqual(ilar_diag.epoch_class(dict(t1=2019.0, t2=2020.0)), "post")

    def test_load_local_csv(self):
        keys = ["mult", "ev1", "ev2", "t1", "t2", "D", "Dx", "cc_ref", "snr_df1", "snr_df2", "ok"]
        with tempfile.TemporaryDirectory() as d:
            p = os.path.join(d, "pairs.csv")
            with open(p, "w", newline="") as f:
                w = csv.DictWriter(f, fieldnames=keys)
                w.writeheader()
                w.writerow(dict(mult="M1", ev1="a", ev2="b", t1=2001.0, t2=2010.0, D=0.1, Dx=0.05, cc_ref=0.95,
                                snr_df1=4.0, snr_df2=3.5, ok="True"))
                w.writerow(dict(mult="M1", ev1="a", ev2="c", t1=2001.0, t2=2012.0, D="nan", Dx=0.05, cc_ref=0.95,
                                snr_df1=4.0, snr_df2=3.5, ok="False"))
            rows = ilar_diag.load_pairs(p, d)
            self.assertEqual(len(rows), 1)
            self.assertTrue(rows[0]["ok"])
            self.assertEqual(rows[0]["snr"], 3.5)


if __name__ == "__main__":
    unittest.main()

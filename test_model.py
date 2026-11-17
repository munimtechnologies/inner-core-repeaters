"""Tests for model.py: the rotation models and the pre-registered grading. Apache-2.0."""
import itertools
import unittest

try:
    import numpy as np
    import scipy  # noqa: F401

    HAVE = True
except ImportError:  # the publish lane may lack numpy/scipy; the compute job runs these for real
    HAVE = False

import model


class TestDates(unittest.TestCase):
    def test_equiv_date(self):
        self.assertAlmostEqual(model.equiv_date(2000.0), 2000.0)
        self.assertAlmostEqual(model.equiv_date(2008.5), 2008.5)
        self.assertAlmostEqual(model.equiv_date(2025.0), 2008.5 - 16.5 / 2.5)  # 2001.9
        self.assertAlmostEqual(model.delta_B(2001.9, 2025.0), 0.0, places=9)
        self.assertAlmostEqual(model.delta_B(2025.0, 2001.9), 0.0, places=9)
        self.assertAlmostEqual(model.delta_B(2010.0, 2015.0), 2.0)  # after the turn: 5 yr / 2.5
        self.assertAlmostEqual(model.delta_M(2010.0, 2015.0), 5.0)

    def test_decimal_year(self):
        self.assertAlmostEqual(model.decimal_year("2025-01-01T00:00:00Z"), 2025.0)
        self.assertAlmostEqual(model.decimal_year("2024-07-02T00:00:00.000Z"), 2024.5)  # leap year
        self.assertAlmostEqual(model.decimal_year("2023-07-02T12:00:00"), 2023.5)


def synth(kind, seed=1, n_mult=30):
    rng = np.random.default_rng(seed)
    rows = []
    for k in range(n_mult):
        n = int(rng.integers(3, 7))
        ts = sorted(rng.uniform(1999.0, 2026.8, size=n))
        ids = ["m%d_e%d" % (k, i) for i in range(n)]
        for (i, t1), (j, t2) in itertools.combinations(list(enumerate(ts)), 2):
            d = model.delta_B(t1, t2) if kind == "B" else model.delta_M(t1, t2)
            D = 0.03 + 0.04 * d + rng.normal(0, 0.03)
            rows.append(dict(mult="M%d" % k, ev1=ids[i], ev2=ids[j], t1=t1, t2=t2, D=D, Dx=D, Ds=D, ok=True))
    return rows


@unittest.skipUnless(HAVE, "numpy/scipy not installed")
class TestGrading(unittest.TestCase):
    def test_backtracking_data_pass_P1(self):
        res = model.evaluate(synth("B"), nboot=300)
        p1 = res["P1_in_sample"]
        self.assertTrue(p1["testable"])
        self.assertTrue(p1["holds"], p1)
        self.assertGreater(p1["rho_B"], p1["rho_M"])
        self.assertTrue(res["P4_calibration"]["holds"])
        fit = res["P5_fit_all"]
        self.assertLess(abs(fit["best_T"] - 2008.5), 2.0, fit)

    def test_monotone_data_fail_P1(self):
        res = model.evaluate(synth("M"), nboot=300)
        p1 = res["P1_in_sample"]
        self.assertTrue(p1["testable"])
        self.assertFalse(p1["holds"], p1)
        self.assertGreater(p1["rho_M"], p1["rho_B"])

    def test_p3_logic(self):
        # new event in 2025.5 (equivalent date ~2001.7) with partners 2001.5 and 2007.5
        rows = [dict(mult="A", ev1="a", ev2="n", t1=2001.5, t2=2025.5, D=0.05, ok=True),
                dict(mult="A", ev1="b", ev2="n", t1=2007.5, t2=2025.5, D=0.30, ok=True),
                dict(mult="B", ev1="c", ev2="m", t1=2002.0, t2=2026.0, D=0.40, ok=True),
                dict(mult="B", ev1="d", ev2="m", t1=2006.0, t2=2026.0, D=0.10, ok=True)]
        r = model.p3(rows)
        self.assertEqual(r["n_discriminating"], 2)
        self.assertEqual(r["n_match_B"], 1)
        self.assertEqual(r["n_match_M"], 1)
        self.assertFalse(r["holds"])

    def test_untestable_when_empty(self):
        res = model.evaluate([], nboot=10)
        self.assertFalse(res["P1_in_sample"]["testable"])
        self.assertIsNone(res["P2_out_of_sample"]["holds"])


if __name__ == "__main__":
    unittest.main()

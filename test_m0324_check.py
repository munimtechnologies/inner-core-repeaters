"""Synthetic tests for m0324_check.py (numpy + scipy; no network, no ObsPy needed).

Part of inner-core-repeaters by MV2, Munim, Inc. (made with AI agents). Apache-2.0.
"""
import csv
import math
import os
import tempfile
import unittest

import numpy as np

import m0324_check as mc

HERE = os.path.dirname(os.path.abspath(__file__))


def ricker(t, f=1.5):
    a = (np.pi * f * t) ** 2
    return (1 - 2 * a) * np.exp(-a)


def synth_event(df_sign, rng, n_el=6):
    """ILAR-like processed event: big 'reference' phases at +8..+19 s, a DF pulse at 0 s."""
    t = np.arange(mc.N) / mc.FS - mc.PRE
    ref = ricker(t - 8.0, 1.0) + 0.8 * ricker(t - 12.0, 1.0) - 0.6 * ricker(t - 19.0, 0.8)
    X = np.array([4.0 * ref + df_sign * ricker(t, 1.5) + 0.05 * rng.standard_normal(mc.N) for _ in range(n_el)])
    X = X / np.sqrt(np.mean(X[:, mc.N0:mc.N1] ** 2, axis=1))[:, None]
    beam = X.mean(axis=0)
    beam = beam / np.sqrt(np.mean(beam[mc.N0:mc.N1] ** 2))
    return dict(names=np.array(["E%d" % k for k in range(n_el)]), X=X, beam=beam, dist=150.0, dt_bc=7.2,
                dt_ab=18.3, snr_group=10.0)


def synth_set():
    rng = np.random.default_rng(1)
    data = {"A": synth_event(1.0, rng), "B": synth_event(-1.0, rng), "C": synth_event(1.0, rng)}
    evs = {"A": dict(id="A", t=2003.0, time="2003-01-01T00:00:00Z"),
           "B": dict(id="B", t=2011.0, time="2011-01-01T00:00:00Z"),
           "C": dict(id="C", t=2025.5, time="2025-07-01T00:00:00Z")}
    edges = {("A", "B"): (0.99, 0.0), ("A", "C"): (0.99, 0.0), ("B", "C"): (0.99, 0.0)}
    return data, evs, edges


class TestHelpers(unittest.TestCase):
    def test_beam_from_normalised(self):
        rng = np.random.default_rng(0)
        X = rng.standard_normal((5, mc.N))
        b = mc.beam_from(X, [0, 1, 2, 2])
        self.assertAlmostEqual(float(np.sqrt(np.mean(b[mc.N0:mc.N1] ** 2))), 1.0, places=9)

    def test_common_elements(self):
        data = {"a": dict(names=np.array(["IL01", "IL02", "IL03"])), "b": dict(names=np.array(["IL02", "IL03"]))}
        self.assertEqual(mc.common_elements(data, ["a", "b"]), ["IL02", "IL03"])

    def test_window(self):
        self.assertEqual(mc.window(1600, (-0.3, 1.2)), (1588, 1648))

    def test_frac_greater(self):
        self.assertEqual(mc.frac_greater([1, 2, float("nan")], [0, 3, 1]), 0.5)

    def test_read_edges_converts_seconds(self):
        with tempfile.TemporaryDirectory() as d:
            p = os.path.join(d, "e.csv")
            with open(p, "w", newline="") as f:
                w = csv.writer(f)
                w.writerow(["ev1", "ev2", "cc", "lag_s"])
                w.writerow(["a", "b", "0.95", "-0.25"])
            e = mc.read_edges(p)
            self.assertEqual(e[("a", "b")], (0.95, -10.0))

    def test_jsonable(self):
        out = mc.jsonable({"x": float("nan"), "y": np.float64(1.5), "z": (np.int64(2), np.bool_(True))})
        self.assertEqual(out, {"x": None, "y": 1.5, "z": [2, True]})


class TestShippedFiles(unittest.TestCase):
    def test_z004_files(self):
        evs = mc.read_events(os.path.join(HERE, "z004_events.csv"))
        edges = mc.read_edges(os.path.join(HERE, "z004_edges.csv"))
        pub = mc.read_published(os.path.join(HERE, "z004_published_D.csv"))
        self.assertEqual(len(evs), 9)
        self.assertEqual(len(edges), 36)
        self.assertEqual(len(pub), 36)
        pairs = set(mc.ordered_pairs(list(evs), evs))
        self.assertEqual(set(pub), pairs)
        self.assertEqual(set(edges), pairs)
        self.assertIn(mc.KEY_A, pub)
        self.assertIn(mc.KEY_B, pub)
        self.assertGreater(pub[mc.KEY_A]["D"], pub[mc.KEY_B]["D"])  # the anomaly being checked


class TestPipeline(unittest.TestCase):
    def test_frame_and_D(self):
        data, evs, edges = synth_set()
        fr = mc.ilar_frame(["A", "B", "C"], data, edges)
        self.assertIsNotNone(fr)
        self.assertLess(abs(fr["pick_minus_pred_s"]), 1.0)
        D = mc.all_D({m: data[m]["beam"] for m in data}, fr, evs)
        self.assertLess(D[("A", "C")], 0.05)
        self.assertGreater(D[("A", "B")], 1.5)
        Dl = mc.all_D({m: data[m]["beam"] for m in data}, fr, evs, dfw=mc.window(fr["onset"], (-0.3, 6.0)))
        self.assertLess(Dl[("A", "C")], Dl[("A", "B")])

    def test_bootstrap_and_grade(self):
        data, evs, edges = synth_set()
        fr = mc.ilar_frame(["A", "B", "C"], data, edges)
        common, boot = mc.element_bootstrap(data, fr, evs, nboot=12, seed=3)
        self.assertEqual(common, ["E%d" % k for k in range(6)])
        loo = mc.leave_one_out(data, fr, evs, common)
        g = mc.grade_c1(boot, loo, a=("A", "B"), b=("A", "C"))
        self.assertTrue(g["testable"])
        self.assertTrue(g["robust"])
        self.assertEqual(g["frac_boot"], 1.0)
        g2 = mc.grade_c1(boot, loo, a=("A", "C"), b=("A", "B"))
        self.assertFalse(g2["robust"])


class TestYkaGrade(unittest.TestCase):
    def test_y1_y0(self):
        rows = [dict(ev1=mc.E2003, ev2=mc.NEW, D=0.30, ok=True), dict(ev1=mc.E2011, ev2=mc.NEW, D=0.10, ok=True),
                dict(ev1=mc.E2003, ev2=mc.E2011, D=0.4, ok=True), dict(ev1=mc.E2003, ev2=mc.E2020, D=0.1, ok=True),
                dict(ev1=mc.E2003, ev2=mc.E2024, D=0.2, ok=False)]
        g = mc.grade_y(rows)
        self.assertTrue(g["Y1"]["testable"])
        self.assertTrue(g["Y1"]["replicated"])
        self.assertFalse(g["Y0"]["testable"])
        self.assertIsNone(g["Y0"]["holds"])
        g = mc.grade_y(rows[:1])
        self.assertFalse(g["Y1"]["testable"])
        self.assertIsNone(g["Y1"]["replicated"])


if __name__ == "__main__":
    unittest.main()

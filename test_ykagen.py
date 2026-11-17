"""Tests for ykagen.py. MV2, Munim, Inc. (made with AI agents). Apache-2.0.

The grading tests are stdlib only. The synthetic comparison test needs yka.py v0.2.0 and
numpy: in the repository yka.py sits next to this file; in the compete job it is the pinned
download in /data/mods_<commit>.
"""
import math
import unittest

import ykagen

ykagen.add_mod_path("/data")
try:
    import numpy as np
    import yka

    HAVE_YKA = getattr(yka, "VERSION", "") == "0.2.0"
except Exception:  # pragma: no cover - environment without the pinned modules
    HAVE_YKA = False


def okrow(D, cc_ref=0.95, flip=False, off_s=0.0, ssi=True, ok=True):
    return dict(status="ok", ok=ok, D=D, cc_ref=cc_ref, flip=flip, off_s=off_s, ssi=ssi, D_c13=D, D_c310=D,
                D_full=D, Ds_median=D)


class TestGrade(unittest.TestCase):
    def test_median_quantile(self):
        self.assertEqual(ykagen.median([3, 1, 2]), 2)
        self.assertEqual(ykagen.median([4, 1, 2, 3]), 2.5)
        self.assertTrue(math.isnan(ykagen.median([])))
        self.assertAlmostEqual(ykagen.quantile([0, 10], 0.25), 2.5)

    def test_in_ssi(self):
        self.assertTrue(ykagen.in_ssi(-58.0, -26.0))
        self.assertFalse(ykagen.in_ssi(-60.3, -46.4))  # 2013-11-17 Scotia Sea M7.7 is outside the SSI box

    def test_all_hold(self):
        g = ykagen.grade([okrow(0.02 + 0.001 * k) for k in range(8)])
        self.assertTrue(g["testable"])
        self.assertEqual(g["H1_floor_small"], "holds")
        self.assertEqual(g["H2_gate_survives"], "holds")
        self.assertEqual(g["H3_no_flip"], "holds")

    def test_failures(self):
        rows = [okrow(0.2, cc_ref=0.8, flip=True) for _ in range(6)]
        g = ykagen.grade(rows)
        self.assertEqual(g["H1_floor_small"], "fails")
        self.assertEqual(g["H2_gate_survives"], "fails")
        self.assertEqual(g["H3_no_flip"], "fails")

    def test_untestable_and_gates(self):
        rows = [okrow(0.01) for _ in range(4)] + [okrow(0.01, ok=False) for _ in range(5)]
        rows.append(dict(status="one-generation-only"))
        g = ykagen.grade(rows)
        self.assertEqual(g["n_ok"], 4)
        self.assertFalse(g["testable"])
        self.assertEqual(g["H1_floor_small"], "untestable")


def synth(shift_s=0.0, flip=False, noise_seed=0, coda_change=0.0):
    """A 1-2 Hz-like PKIKP pulse plus coda at PRE, in unit-rms noise (yka.py's beam scaling)."""
    FS, PRE, N = yka.FS, yka.PRE, yka.N
    rng = np.random.default_rng(noise_seed)
    t = np.arange(N) / FS - PRE - shift_s
    sig = np.zeros(N)
    on = t >= 0
    sig[on] = 60.0 * np.sin(2 * np.pi * 1.5 * t[on]) * np.exp(-t[on] / 0.6)
    coda = (t >= 1.0) & (t <= 8.0)
    tc = t[coda]
    # broadband coda: 12 tones in 0.9-2.1 Hz with fixed random phases (the "path"); a changed
    # coda draws new phases for every tone
    crng = np.random.default_rng(99)
    freqs = np.linspace(0.9, 2.1, 12)
    phases = crng.uniform(0, 2 * np.pi, 12)
    if coda_change:  # half of the tones get new phases
        phases[::2] = np.random.default_rng(1234).uniform(0, 2 * np.pi, 6)
    amp = 15.0 / math.sqrt(6.0)
    c = sum(np.sin(2 * np.pi * f * tc + p) for f, p in zip(freqs, phases))
    sig[coda] += amp * c * np.exp(-(tc - 1.0) / 4.0)
    b = sig + rng.normal(0.0, 1.0, N)
    return dict(beam=-b if flip else b, names=np.array([]))


@unittest.skipUnless(HAVE_YKA, "yka.py v0.2.0 not importable here")
class TestCompare(unittest.TestCase):
    def test_same_event_recovers_offset(self):
        r = ykagen.compare(yka, synth(0.0, noise_seed=1), synth(0.1, noise_seed=2))
        self.assertEqual(r["status"], "ok")
        self.assertFalse(r["flip"])
        self.assertAlmostEqual(r["off_s"], 0.1, delta=0.03)
        self.assertLess(r["D"], 0.05)
        self.assertGreater(r["cc_ref"], 0.9)
        self.assertTrue(r["ok"])

    def test_flip_detected(self):
        r = ykagen.compare(yka, synth(0.0, noise_seed=3), synth(0.0, flip=True, noise_seed=4))
        self.assertEqual(r["status"], "ok")
        self.assertTrue(r["flip"])
        self.assertLess(r["D"], 0.05)

    def test_changed_coda_raises_D(self):
        same = ykagen.compare(yka, synth(0.0, noise_seed=5), synth(0.0, noise_seed=6))
        diff = ykagen.compare(yka, synth(0.0, noise_seed=5), synth(0.0, noise_seed=6, coda_change=1.0))
        self.assertGreater(diff["D"], same["D"] + 0.1)


if __name__ == "__main__":
    unittest.main()

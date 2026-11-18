"""Tests for yka_repick.py (MV2, Munim, Inc.; made with AI agents; Apache-2.0)."""
import os
import sys
import unittest

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
for p in (HERE, os.path.join(HERE, "v080")):
    if p not in sys.path:
        sys.path.insert(0, p)

import yka_repick as yr  # noqa: E402

try:
    from xcorr import pick_onset  # repo root, or v080/ inside a compete job
except Exception:  # pragma: no cover
    pick_onset = None

FS = 40.0
PRE = 40.0
CENTER = int(PRE * FS)


def burst(n, i0, dur_s, amp, f=1.5, seed=1, decay_s=None):
    """Band-limited random-phase burst starting at sample i0 (sum of 1-2 Hz sines)."""
    rng = np.random.default_rng(seed)
    t = np.arange(n) / FS
    x = np.zeros(n)
    m = int(dur_s * FS)
    tt = t[:m]
    sig = np.zeros(m)
    for fr in np.linspace(1.0, 2.0, 9):
        sig += np.sin(2 * np.pi * fr * tt + rng.uniform(0, 2 * np.pi))
    sig /= np.max(np.abs(sig))
    env = np.exp(-tt / decay_s) if decay_s else np.ones(m)
    ramp = np.minimum(1.0, tt / 0.15)
    x[i0:i0 + m] = amp * sig * env * ramp
    return x


def synthetic_yka_stack(seed=0):
    """Quiet stack: weak coherent precursor at ak135 - 4.8 s, strong PKIKP at ak135 - 0.8 s."""
    n = int(80 * FS)
    rng = np.random.default_rng(seed)
    x = 0.01 * rng.standard_normal(n)
    x += burst(n, CENTER - int(4.8 * FS), 1.2, 0.2, seed=2)
    x += burst(n, CENTER - int(0.8 * FS), 8.0, 4.0, seed=3, decay_s=1.5)
    return x


class TestRules(unittest.TestCase):
    def test_pick_valid(self):
        self.assertTrue(yr.pick_valid(0.0))
        self.assertTrue(yr.pick_valid(-2.0))
        self.assertTrue(yr.pick_valid("1.5"))
        self.assertFalse(yr.pick_valid(-4.8))
        self.assertFalse(yr.pick_valid(2.01))
        self.assertFalse(yr.pick_valid(None))
        self.assertFalse(yr.pick_valid(float("nan")))

    def test_wrapper_overrides_search_and_forwards(self):
        calls = []

        def fake(x, fs, center, search_s=6.0, sta_s=0.5, lta_s=8.0, thr=6.0):
            calls.append(dict(center=center, search_s=search_s, sta_s=sta_s, thr=thr))
            return center

        p = yr.constrained_picker(fake, 2.0)
        out = p(np.zeros(10), 40.0, center=5, search_s=6.0, sta_s=0.4, thr=7.0)
        self.assertEqual(out, 5)
        self.assertEqual(calls[0], dict(center=5, search_s=2.0, sta_s=0.4, thr=7.0))

    def test_envelope_onset(self):
        x = synthetic_yka_stack()
        on, mx = yr.envelope_onset(x, FS, CENTER)
        self.assertIsNotNone(on)
        self.assertTrue(-1.0 <= (on - CENTER) / FS <= 0.0, (on - CENTER) / FS)
        self.assertTrue(-1.0 <= (mx - CENTER) / FS <= 3.0, (mx - CENTER) / FS)

    def test_envelope_onset_empty(self):
        self.assertEqual(yr.envelope_onset(np.zeros(100), FS, 50), (None, None))


@unittest.skipIf(pick_onset is None, "xcorr.py not importable here")
class TestPicker(unittest.TestCase):
    def test_original_picker_takes_precursor(self):
        """Reproduces the icr-m0324chk-v1 failure: +/-6 s picks the weak early energy."""
        i = pick_onset(synthetic_yka_stack(), FS, CENTER, search_s=6.0)
        self.assertIsNotNone(i)
        self.assertLess((i - CENTER) / FS, -4.0)

    def test_constrained_picker_takes_pkikp(self):
        rec = []
        p = yr.constrained_picker(pick_onset, 2.0, record=rec)
        i = p(synthetic_yka_stack(), FS, center=CENTER, search_s=6.0, sta_s=0.5, lta_s=8.0, thr=6.0)
        self.assertIsNotNone(i)
        rel = (i - CENTER) / FS
        self.assertTrue(-1.3 <= rel <= -0.2, rel)
        self.assertTrue(yr.pick_valid(rel))
        self.assertEqual(len(rec), 1)
        self.assertLess(rec[0]["pick_original_s"], -4.0)
        self.assertAlmostEqual(rec[0]["pick_constrained_s"], round(rel, 3), places=3)


if __name__ == "__main__":
    unittest.main()

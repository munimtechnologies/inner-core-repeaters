# Regression tests for brvpkp31.py 0.1.1 (the despike fix). Synthetic data only; stdlib + numpy.
# Part of munimtechnologies/inner-core-repeaters. Licence: Apache-2.0.
# Made with AI agents working in MV2's name (MV2, Munim, Inc.).
import unittest

import numpy as np

import brvpkp31 as B
import brvpkp31_core as C


class TestDespikeFix(unittest.TestCase):
    def test_wrapper_is_used_by_pipeline(self):
        self.assertIs(C.despike, B.despike)
        self.assertEqual(B.VERSION, "0.1.1")
        self.assertEqual(C.VERSION, "0.1.1")
        self.assertIs(B.PARAMS, C.PARAMS)

    def test_quantised_trace(self):
        # Integer counts of a smooth wave: the running-median residual is mostly exact zeros.
        t = np.arange(4000)
        clean = np.round(100.0 * np.sin(2 * np.pi * t / 40.0))
        x = clean.copy()
        x[1234] += 5000.0
        y, cnt = B.despike(x)
        self.assertEqual(cnt, 1)
        self.assertEqual(int(np.count_nonzero(y != x)), 1)
        self.assertLess(abs(y[1234] - clean[1234]), 10.0)

    def test_no_false_alarms_on_strong_signal(self):
        # 0.8 Hz at 20 samples/s (25 samples per period), amplitude 1000: nothing is a glitch.
        t = np.arange(5000)
        x = 1000.0 * np.sin(2 * np.pi * t / 25.0)
        y, cnt = B.despike(x)
        self.assertEqual(cnt, 0)
        np.testing.assert_array_equal(y, x)

    def test_long_run_kept(self):
        # Runs longer than spike_max_run (3) are not glitches by the rule and stay untouched.
        t = np.arange(3000)
        x = np.sin(2 * np.pi * t / 100.0)
        x[1000:1005] += 30.0
        y, cnt = B.despike(x)
        self.assertEqual(cnt, 0)
        np.testing.assert_array_equal(y, x)

    def test_flat_trace(self):
        y, cnt = B.despike(np.zeros(100))
        self.assertEqual(cnt, 0)
        y, cnt = B.despike(np.ones(3))
        self.assertEqual(cnt, 0)


if __name__ == "__main__":
    unittest.main()

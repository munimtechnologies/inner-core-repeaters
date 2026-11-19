# test_ilresp.py - unit tests for ilresp.py (stdlib + numpy; scipy tests skip if absent).
# Part of github.com/munimtechnologies/inner-core-repeaters. Apache-2.0.
# By MV2 of Munim, Inc.; made with AI agents.
#
# 2026-11-19: test_max_cc_identity_and_shift asserted an exactly-zero lag for identical
# inputs; the parabolic peak refinement over a non-circular window gives ~1.5e-4 samples.
# The tolerance is now 0.01 samples (0.5 ms at 20 Hz). ilresp.py is unchanged.
import importlib.util
import math
import unittest

import numpy as np

import ilresp

HAVE_SCIPY = importlib.util.find_spec("scipy") is not None


class FIRResponseStage:
    def __init__(self, coeffs, symmetry, fs, corr=0.0, gain=1.0):
        self.coefficients = coeffs
        self.symmetry = symmetry
        self.decimation_input_sample_rate = fs
        self.decimation_correction = corr
        self.stage_gain = gain


class CoefficientsTypeResponseStage:
    def __init__(self, num, den, fs=None, gain=1.0):
        self.numerator = num
        self.denominator = den
        self.decimation_input_sample_rate = fs
        self.decimation_correction = 0.0
        self.stage_gain = gain


class PolesZerosResponseStage:
    def __init__(self, zeros, poles, a0=1.0, gain=1.0):
        self.zeros = zeros
        self.poles = poles
        self.normalization_factor = a0
        self.pz_transfer_function_type = "LAPLACE (RADIANS/SECOND)"
        self.stage_gain = gain
        self.decimation_input_sample_rate = None
        self.decimation_correction = None


class TestResponses(unittest.TestCase):
    def test_expand_fir(self):
        self.assertEqual(list(ilresp.expand_fir([1, 2], "EVEN")), [1, 2, 2, 1])
        self.assertEqual(list(ilresp.expand_fir([1, 2, 3], "ODD")), [1, 2, 3, 2, 1])
        self.assertEqual(list(ilresp.expand_fir([1, 2, 3], "NONE")), [1, 2, 3])

    def test_unit_delay_phase(self):
        f = np.array([0.5, 1.0, 2.0])
        h = ilresp.digital_response([0.0, 1.0], [1.0], 20.0, f)
        self.assertTrue(np.allclose(np.abs(h), 1.0))
        self.assertTrue(np.allclose(np.angle(h), -2 * np.pi * f / 20.0))

    def test_fir_with_full_correction_is_zero_phase(self):
        fs = 100.0
        st = FIRResponseStage([0.25, 0.25], "EVEN", fs, corr=1.5 / fs)
        h = ilresp.stage_response(st, np.array([0.5, 1.0, 2.0]))
        self.assertLess(float(np.max(np.abs(np.angle(h)))), 1e-9)
        st2 = FIRResponseStage([0.25, 0.25], "EVEN", fs, corr=0.0)
        h2 = ilresp.stage_response(st2, np.array([1.0]))
        self.assertAlmostEqual(float(np.angle(h2[0])), -2 * np.pi * 1.0 * 1.5 / fs, places=9)

    def test_constant_coefficient_stage(self):
        st = CoefficientsTypeResponseStage([1.0], [], fs=None, gain=4.0)
        h = ilresp.stage_response(st, np.array([1.0, 3.0]))
        self.assertTrue(np.allclose(h, 4.0))

    def test_pz_single_pole(self):
        st = PolesZerosResponseStage([], [complex(-2 * np.pi, 0)])
        h = ilresp.stage_response(st, np.array([1.0]))
        self.assertAlmostEqual(abs(h[0]), 1.0 / (2 * np.pi * math.sqrt(2)), places=9)


class TestAnalysis(unittest.TestCase):
    def test_ratio_metrics_recovers_delay(self):
        f = np.fft.rfftfreq(ilresp.NFFT, 1.0 / ilresp.FS)[1:]
        r = np.exp(-2j * np.pi * f * 0.03)
        m = ilresp.ratio_metrics(f, r, (1.0, 2.0))
        self.assertAlmostEqual(m["delay_s"], 0.03, places=6)
        self.assertLess(m["resid_phase_rms_deg"], 1e-6)
        self.assertLess(m["amp_dev_max_pct"], 1e-9)
        self.assertLess(abs(m["phase_intercept_deg"]), 1e-4)

    def test_remove_delay(self):
        f = np.fft.rfftfreq(ilresp.NFFT, 1.0 / ilresp.FS)
        r = np.exp(-2j * np.pi * f * 0.0237)
        self.assertTrue(np.allclose(ilresp.remove_delay(f, r, 0.0237), 1.0))

    def test_apply_ratio_identity(self):
        x = np.random.default_rng(1).standard_normal(ilresp.NFFT)
        r = np.ones(ilresp.NFFT // 2 + 1, dtype=complex)
        self.assertTrue(np.allclose(ilresp.apply_ratio(x, r), x))

    def test_max_cc_identity_and_shift(self):
        rng = np.random.default_rng(2)
        a = np.convolve(rng.standard_normal(500), np.hanning(9), mode="same")
        cc, lag = ilresp.max_cc(a, a, 20)
        self.assertGreater(cc, 1 - 1e-9)
        self.assertLess(abs(lag), 0.01)
        b = np.roll(a, 3)
        cc, lag = ilresp.max_cc(a, b, 20)
        self.assertGreater(cc, 0.999)
        self.assertLess(abs(lag - 3.0), 0.2)

    def test_ricker_peak(self):
        x = ilresp.ricker(101, 20.0, 1.5, 2.5)
        self.assertEqual(int(np.argmax(x)), 50)

    @unittest.skipUnless(HAVE_SCIPY, "scipy not installed")
    def test_decorrelation_identity(self):
        r = np.ones(ilresp.NFFT // 2 + 1, dtype=complex)
        d = ilresp.decorrelation(r, (1.0, 2.0), np.random.default_rng(3), n_noise=5)
        self.assertLess(max(d["noise"]), 1e-9)
        self.assertLess(d["ricker"], 1e-9)

    @unittest.skipUnless(HAVE_SCIPY, "scipy not installed")
    def test_decorrelation_sees_shape_change(self):
        f = np.fft.rfftfreq(ilresp.NFFT, 1.0 / ilresp.FS)
        # a strong 90-degree phase rotation changes the shape: D must be clearly > 0
        r = np.where(f > 0, np.exp(0.5j * np.pi), 1.0).astype(complex)
        d = ilresp.decorrelation(r, (1.0, 2.0), np.random.default_rng(4), n_noise=5)
        self.assertGreater(float(np.median(d["noise"])), 0.01)


if __name__ == "__main__":
    unittest.main()

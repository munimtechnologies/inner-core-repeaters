"""Tests for yka.py (synthetic signals only; no obspy, no network). Apache-2.0.

By MV2, Munim, Inc. (made with AI agents).
"""
import unittest

try:
    import numpy as np

    import yka
    from xcorr import fshift

    HAVE = True
except ImportError:  # the publish lane may lack numpy; the compute job runs these for real
    HAVE = False

FS = 40.0
PRE = 40.0
N = 3200


def event(seed_coda, amp=1.0, noise=0.02, seed_noise=0, coda_scale=0.5):
    """Direct PKIKP pulse at PRE s, followed by a 6-s scattered coda, plus weak noise."""
    t = np.arange(N) / FS
    tau = t - PRE
    direct = np.where(tau >= 0, np.sin(2 * np.pi * 1.5 * tau) * np.exp(-tau / 0.35), 0.0)
    rc = np.random.default_rng(seed_coda)
    coda = np.zeros(N)
    for _ in range(25):
        t0 = PRE + rc.uniform(1.6, 7.0)  # keep the coda clear of the direct-pulse window
        f = rc.uniform(1.1, 1.9)
        a = rc.normal(0, 1) * np.exp(-(t0 - PRE) / 4.0)
        coda += a * np.exp(-((t - t0) * f * 1.8) ** 2) * np.sin(2 * np.pi * f * (t - t0))
    rn = np.random.default_rng(seed_noise)
    return amp * (direct + coda_scale * coda) + noise * rn.standard_normal(N)


@unittest.skipUnless(HAVE, "numpy not installed")
class TestYka(unittest.TestCase):
    def test_canonical_names(self):
        self.assertEqual(yka.canonical("YKAB3"), "B3")
        self.assertEqual(yka.canonical("YKB3"), "B3")
        self.assertEqual(yka.canonical("YKAR7"), "R7")
        self.assertEqual(yka.canonical("YKR7"), "R7")
        self.assertIsNone(yka.canonical("YKAW1"))
        self.assertIsNone(yka.canonical("YKB5"))  # no B5 element
        self.assertIsNone(yka.canonical("IL01"))
        self.assertEqual(yka.generation("YKAR1"), "new")
        self.assertEqual(yka.generation("YKR1"), "legacy")

    def test_providers(self):
        self.assertEqual(yka.providers_for(2001.0)[0], yka.NRCAN)
        self.assertEqual(yka.providers_for(2025.5)[0], yka.EARTHSCOPE)

    def test_parse_members(self):
        s = "us7000abcd(2005-06-01);usp000xyz(2001-02-03); "
        self.assertEqual(yka.parse_members(s), ["us7000abcd", "usp000xyz"])
        self.assertEqual(yka.parse_members(""), [])

    def test_envelope_peaks_at_pulse(self):
        x = event(1, coda_scale=0.0, noise=0.0)
        e = yka.envelope(x, FS, smooth_s=0.1)
        self.assertLess(abs(int(np.argmax(e)) - int(PRE * FS)), int(0.5 * FS))

    def test_align_member_recovers_offset(self):
        x = event(1)
        y = fshift(event(1, amp=0.7, seed_noise=5), 57.3)  # 1.43 s later, smaller
        a = yka.align_member(x, y, FS, PRE)
        self.assertLess(abs(a["off"] - 57.3), 0.3, a)
        self.assertGreater(a["cc_fine"], 0.95)
        self.assertFalse(a["flip"])

    def test_align_member_flip(self):
        x = event(2)
        y = -fshift(event(2, seed_noise=7), -31.0)
        a = yka.align_member(x, y, FS, PRE, allow_flip=True)
        self.assertTrue(a["flip"], a)
        self.assertLess(abs(a["off"] + 31.0), 0.3, a)
        b = yka.align_member(x, y, FS, PRE, allow_flip=False)
        self.assertFalse(b["flip"])

    def test_pair_same_vs_changed_coda(self):
        x = event(3)
        same = fshift(event(3, seed_noise=9), 22.0)
        changed = fshift(event(4, seed_noise=9), 22.0)  # same direct pulse, different coda
        onset = int(PRE * FS)
        r1 = yka.yka_pair(x, same, onset, 0.0, 22.0, FS)
        r2 = yka.yka_pair(x, changed, onset, 0.0, 22.0, FS)
        self.assertGreater(r1["cc_ref"], 0.95)
        self.assertGreater(r2["cc_ref"], 0.95)
        self.assertGreater(r1["cc_df0"], 0.95)
        self.assertLess(r2["cc_df0"], 0.8)  # two independent random codas
        self.assertLess(r2["cc_df0"], r1["cc_df0"] - 0.15)
        self.assertLess(abs(r1["lag_total"] - 22.0), 0.3)

    def test_pair_offsets_in_member_frames(self):
        # root at 0, a delayed by 10 samples, b delayed by 30: the pair must still line up
        base = event(5)
        xa = fshift(base, 10.0)
        xb = fshift(event(5, seed_noise=3), 30.0)
        r = yka.yka_pair(xa, xb, int(PRE * FS), 10.0, 30.0, FS)
        self.assertIsNotNone(r)
        self.assertGreater(r["cc_df0"], 0.95)
        self.assertLess(abs(r["lag_total"] - 20.0), 0.3)

    def test_pair_out_of_range(self):
        x = event(6)
        self.assertIsNone(yka.yka_pair(x, x, 5, 0.0, 0.0, FS))


if __name__ == "__main__":
    unittest.main()

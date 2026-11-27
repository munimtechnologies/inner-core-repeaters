"""Tests for m2measure.py (pre-registration 38, arm M, measurement front end). Expected values
derived by hand (findings/2026-11-27-prereg38-armM-m2measure-frozen.md). Standard library only.
Part of inner-core-repeaters by MV2, Munim, Inc. (made with AI agents). Licence: Apache-2.0."""
import math
import unittest

import m2measure as m


class TestEvents(unittest.TestCase):
    def test_decimal_year_matches_z004(self):
        # z004_events.csv / yka_z004_pairs.csv: 2003-03-21T10:47:47.620Z -> 2003.2176708403094
        self.assertAlmostEqual(m.decimal_year("2003-03-21T10:47:47.620Z"), 2003.2176708403094, places=9)
        self.assertAlmostEqual(m.decimal_year("2003-03-21T10:47:47.620"), 2003.2176708403094, places=9)
        self.assertAlmostEqual(m.decimal_year("2020-01-01T00:00:00"), 2020.0)

    def test_measured_letters(self):
        self.assertEqual(len(m.MEASURED), 14)
        self.assertNotIn("I", m.MEASURED)
        self.assertNotIn("M", m.MEASURED)

    def test_multiplet_A_uses_734(self):
        a = m.multiplet_events("A")
        self.assertEqual(len(a), 7)
        self.assertIn("usp000fz4d", a)          # Wang 734, 2008-02-17
        self.assertEqual(a["usp000fz4d"]["wang"], 734)

    def test_unmatched_event(self):
        self.assertIsNone(m.event_dict(711))
        self.assertEqual(sorted(m.multiplet_events("K")), ["usp0007zxb", "usp000dr5u"])

    def test_no_event_in_two_multiplets(self):
        seen = {}
        for k in m.MULTIPLETS:
            for i in m.multiplet_events(k):
                self.assertNotIn(i, seen, "%s in %s and %s" % (i, seen.get(i), k))
                seen[i] = k

    def test_member_pairs_earlier_first(self):
        evs = {"x": dict(t=2010.0), "y": dict(t=2001.0), "z": dict(t=2005.0)}
        self.assertEqual(m.member_pairs(evs), [("y", "z"), ("y", "x"), ("z", "x")])


class TestIlar(unittest.TestCase):
    def test_gates(self):
        self.assertTrue(m.ilar_ok(0.1, 3.0, 3.0, 5.0, 5.0, 3.0, 3.0))
        self.assertFalse(m.ilar_ok(0.1, 2.9, 3.0, 5.0, 5.0, 3.0, 3.0))
        self.assertFalse(m.ilar_ok(0.1, 3.0, 3.0, 4.9, 5.0, 3.0, 3.0))
        self.assertFalse(m.ilar_ok(0.1, 3.0, 3.0, 5.0, 5.0, 3.0, 2.9))
        self.assertFalse(m.ilar_ok(float("nan"), 9, 9, 9, 9, 9, 9))

    def test_rows(self):
        evs = {"x": dict(t=2001.0), "y": dict(t=2005.0), "z": dict(t=2009.0)}
        good = (4.0, 6.0, 4.0)
        rows = m.ilar_rows("Q", evs, {("x", "y"): 0.08, ("z", "y"): 0.2}, {"x": good, "y": good})
        by = {(r["ev1"], r["ev2"]): r for r in rows}
        self.assertEqual(len(rows), 3)
        self.assertTrue(by[("x", "y")]["ok"])
        self.assertAlmostEqual(by[("x", "y")]["D"], 0.08)
        self.assertAlmostEqual(by[("y", "z")]["D"], 0.2)       # stored reversed, found anyway
        self.assertFalse(by[("y", "z")]["ok"])                  # z has no SNR values
        self.assertTrue(math.isnan(by[("x", "z")]["D"]))
        self.assertFalse(by[("x", "z")]["ok"])
        self.assertEqual(by[("x", "y")]["array"], "ILAR")
        self.assertEqual(by[("x", "y")]["t2"], 2005.0)

    def test_frame_snr(self):
        beam = [1.0, 1.0, 1.0, 1.0, 2.0, 2.0, 2.0, 2.0]
        self.assertAlmostEqual(m.frame_snr(beam, 0, (4, 8), 0, 4, 8), 2.0)
        self.assertTrue(math.isnan(m.frame_snr(beam, 1, (4, 8), 0, 4, 8)))
        self.assertAlmostEqual(m.frame_snr(beam, -1, (4, 8), 0, 4, 8), math.sqrt(3.25))
        self.assertTrue(math.isnan(m.frame_snr([0.0] * 8, 0, (4, 8), 0, 4, 8)))


class TestYka(unittest.TestCase):
    def test_rows(self):
        evs = {"x": dict(t=2001.0), "y": dict(t=2005.0), "z": dict(t=2009.0)}
        raw = [dict(ev1="y", ev2="x", D_c310="0.05", D="0.02", D_c13="0.03", ok="True"),
               dict(ev1="x", ev2="z", D_c310="0.2", D="0.04", D_c13="0.1", ok="False")]
        by = {(r["ev1"], r["ev2"]): r for r in m.yka_rows("Q", evs, raw)}
        self.assertTrue(by[("x", "y")]["ok"])
        self.assertAlmostEqual(by[("x", "y")]["D"], 0.05)
        self.assertAlmostEqual(by[("x", "y")]["D_direct"], 0.02)
        self.assertFalse(by[("x", "z")]["ok"])
        self.assertAlmostEqual(by[("x", "z")]["D"], 0.2)
        self.assertTrue(math.isnan(by[("y", "z")]["D"]))
        self.assertFalse(by[("y", "z")]["ok"])
        self.assertEqual(by[("y", "z")]["array"], "YKA")


class TestGate(unittest.TestCase):
    def test_compare(self):
        ours = [dict(ev1="a", ev2="b", D=0.1005), dict(ev1="a", ev2="d", D=0.3)]
        pub = [dict(ev1="b", ev2="a", D_reproduced="0.1"), dict(ev1="a", ev2="c", D_reproduced="")]
        c = m.compare(ours, pub, "D", "D_reproduced")
        self.assertEqual(c["n_published"], 1)
        self.assertEqual(c["n_common"], 1)
        self.assertAlmostEqual(c["max_abs_diff"], 0.0005)
        self.assertTrue(m.gate_ok(c))
        c2 = m.compare([dict(ev1="a", ev2="b", D=0.102)], pub, "D", "D_reproduced")
        self.assertFalse(m.gate_ok(c2))
        c3 = m.compare([], pub, "D", "D_reproduced")
        self.assertIsNone(c3["max_abs_diff"])
        self.assertFalse(m.gate_ok(c3))
        self.assertFalse(m.gate_ok({}))

    def test_coverage(self):
        self.assertFalse(m.gate_ok(dict(n_published=36, n_common=32, max_abs_diff=0.0)))
        self.assertTrue(m.gate_ok(dict(n_published=36, n_common=33, max_abs_diff=0.0)))


if __name__ == "__main__":
    unittest.main()

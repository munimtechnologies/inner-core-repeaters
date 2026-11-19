"""Tests for patch_nea.py (addendum T to pre-registration 16).
Part of inner-core-repeaters by MV2, Munim, Inc. (made with AI agents). Licence: Apache-2.0."""
import importlib.util
import unittest

import patch_nea

SRC = '''VERSION = "0.1.0"
NAN = float("nan")


def main():
    for x in y:
        if True:
            try:
                arr = model.get_travel_times(source_depth_in_km=max(ev["depth"], 0.0), distance_in_degree=dist,
                                             phase_list=["PKIKP", "PKPab"])
            except Exception as ex:  # noqa: BLE001
                continue
            tdf = next((a.time for a in arr if a.name == "PKIKP"), NAN)
            tab = next((a.time for a in arr if a.name == "PKPab"), NAN)
    ax.set_xlabel("time relative to ak135 PKIKP (s); red = DF window, blue = PKPab reference window")
'''


class Arr:
    def __init__(self, name, t):
        self.name, self.time = name, t


class TestPatch(unittest.TestCase):
    def test_applies_exactly_two_changes(self):
        out = patch_nea.apply_patch(SRC)
        self.assertIn(patch_nea.NEW_PHASES, out)
        self.assertIn("tab = " + patch_nea.NEW_AB, out)
        self.assertNotIn('"PKPab"', out)
        # everything else untouched
        self.assertIn('tdf = next((a.time for a in arr if a.name == "PKIKP"), NAN)', out)
        self.assertIn("blue = PKPab reference window", out)
        self.assertIn('VERSION = "0.1.0"', out)
        self.assertEqual(len(out.splitlines()), len(SRC.splitlines()))

    def test_single_quotes_also_matched(self):
        out = patch_nea.apply_patch(SRC.replace('"PKPab"', "'PKPab'").replace('"PKIKP", ', "'PKIKP', "))
        self.assertIn(patch_nea.NEW_PHASES, out)

    def test_missing_pattern_refused(self):
        with self.assertRaises(ValueError):
            patch_nea.apply_patch(SRC.replace('phase_list=["PKIKP", "PKPab"]', 'phase_list=["PKIKP"]'))

    def test_duplicate_pattern_refused(self):
        dup = SRC + '        tab = next((a.time for a in arr if a.name == "PKPab"), NAN)\n'
        with self.assertRaises(ValueError):
            patch_nea.apply_patch(dup)

    def test_already_patched_refused(self):
        with self.assertRaises(ValueError):
            patch_nea.apply_patch(patch_nea.apply_patch(SRC))

    def test_ab_time_is_latest_pkp(self):
        arr = [Arr("PKIKP", 1180.0), Arr("PKP", 1183.5), Arr("PKP", 1191.25), Arr("PKiKP", 1181.0)]
        self.assertEqual(patch_nea.ab_time(arr), 1191.25)
        nan = patch_nea.ab_time([Arr("PKIKP", 1180.0)])
        self.assertTrue(nan != nan)

    def test_inline_expression_matches_ab_time(self):
        arr = [Arr("PKIKP", 1180.0), Arr("PKP", 1183.5), Arr("PKP", 1191.25)]
        got = eval(patch_nea.NEW_AB, {"arr": arr, "NAN": float("nan")})  # noqa: S307 (fixed literal)
        self.assertEqual(got, patch_nea.ab_time(arr))


@unittest.skipUnless(importlib.util.find_spec("obspy") is not None, "obspy not installed")
class TestTauP(unittest.TestCase):
    def test_pkp_gives_ab_after_df(self):
        from obspy.taup import TauPyModel

        m = TauPyModel("ak135")
        for dist in (156.0, 165.0, 175.0, 179.5):
            arr = m.get_travel_times(source_depth_in_km=20.0, distance_in_degree=dist, phase_list=["PKIKP", "PKP"])
            tdf = next((a.time for a in arr if a.name == "PKIKP"), float("nan"))
            tab = patch_nea.ab_time(arr)
            self.assertTrue(tdf == tdf and tab == tab, "missing phase at %.1f deg" % dist)
            self.assertGreater(tab - tdf, 0.0)
            self.assertLess(tab - tdf, 60.0)


if __name__ == "__main__":
    unittest.main()

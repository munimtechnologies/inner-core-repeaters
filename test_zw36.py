# Copyright 2026 MV2, Munim, Inc.
# Licensed under the Apache License, Version 2.0 (Apache-2.0). See LICENSE.
# Made with AI agents working in MV2's name.
"""Unit tests for zw36.py (stdlib only)."""

import math
import os
import sys
import unittest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import zw36  # noqa: E402

A, B, C = zw36.F2005_DATE, zw36.F2012_DATE, zw36.F2022_DATE
N1, N2, N3 = "2024-04-19", "2025-09-23", "2026-05-02"


def trow(a, b, c, grade, rec="YKA", other="ILAR"):
    return dict(recovery_channel=rec, a=a, b=b, c=c, other_channel=other, grade=grade)


def all_n(*pairs):
    return {zw36.key(x, y): "N" for x, y in pairs}


class TestEvents(unittest.TestCase):
    def test_counts_and_dates(self):
        self.assertEqual(len(zw36.F_OLD), 3)
        self.assertEqual(len(zw36.F_NEW), 11)
        self.assertEqual(len(zw36.F_REPORT), 2)
        dates = [zw36.date_of(e) for e in zw36.F_OLD + zw36.F_NEW + zw36.F_REPORT]
        self.assertEqual(len(set(dates)), 16)
        self.assertEqual(len(zw36.graded_dates()), 14)
        self.assertNotIn("2025-12-17", zw36.new_dates())
        self.assertIn(N2, zw36.new_dates())

    def test_new_events_within_rule(self):
        # 15.0 km of the F centroid, depth <= 45 km, mag >= 4.7 (pre-registration 36)
        lat0, lon0 = -58.7887, -25.3637
        for e in zw36.F_NEW:
            dy = (e[2] - lat0) * 111.19
            dx = (e[3] - lon0) * 111.19 * math.cos(math.radians(lat0))
            self.assertLessEqual(math.hypot(dx, dy), 15.0, e[0])
            self.assertLessEqual(e[4], 45.0)
            self.assertGreaterEqual(e[5], 4.7)
        for e in zw36.F_REPORT:
            dy = (e[2] - lat0) * 111.19
            dx = (e[3] - lon0) * 111.19 * math.cos(math.radians(lat0))
            self.assertGreater(math.hypot(dx, dy), 15.0, e[0])

    def test_pins(self):
        for c in (zw36.ZW33_COMMIT, zw36.CLOSURE35_COMMIT, zw36.MANTLEP2_COMMIT):
            self.assertEqual(len(c), 40)
        self.assertEqual(zw36.VERSION, "0.1.0")


class TestSource(unittest.TestCase):
    def test_score_and_decide(self):
        refs = zw36.OLD8[:5]
        gated = {e: True for e in refs}
        gated.update(x=True, y=True)
        D = {}
        for i in range(5):
            for j in range(i + 1, 5):
                D[zw36.key(refs[i], refs[j])] = 0.1 + 0.01 * (i + j)
        D[zw36.key("x", "y")] = 0.9
        sc = zw36.station_pair_score(D, gated, refs, "x", "y")
        self.assertTrue(sc["counted"])
        self.assertEqual(sc["n_ref"], 10)
        self.assertEqual(sc["q"], 1.0)
        self.assertTrue(sc["hit"])
        D[zw36.key("x", "y")] = 0.0
        self.assertFalse(zw36.station_pair_score(D, gated, refs, "x", "y")["hit"])
        gated["y"] = False
        self.assertFalse(zw36.station_pair_score(D, gated, refs, "x", "y")["counted"])

    def test_decide(self):
        self.assertEqual(zw36.decide([True, False])[0], "untestable")
        self.assertEqual(zw36.decide([True, True, False, False])[0], "S")
        self.assertEqual(zw36.decide([True, False, False, False, False])[0], "N")
        self.assertEqual(zw36.decide([True, True, False, False, False, False, False])[0],
                         "inconclusive")

    def test_source_decisions(self):
        sc = {"S1": {(A, N1): dict(counted=True, hit=False, q=0.1)},
              "S2": {(A, N1): dict(counted=True, hit=False, q=0.2)},
              "S3": {(A, N1): dict(counted=True, hit=True, q=0.95)},
              "S4": {(A, N1): dict(counted=False)}}
        d = zw36.source_decisions(sc)
        self.assertEqual(d[(A, N1)]["decision"], "N")
        self.assertEqual(d[(A, N1)]["n_stations"], 3)

    def test_combine_external(self):
        self.assertEqual(zw36.combine_external("N", None), "N")
        self.assertEqual(zw36.combine_external("N", "N"), "N")
        self.assertEqual(zw36.combine_external("N", "S"), "inconclusive")
        self.assertEqual(zw36.combine_external("untestable", "S"), "untestable")


class TestP36a(unittest.TestCase):
    new = {N1, N2, N3}

    def test_old_only_triple_does_not_count(self):
        rows = [trow(A, B, C, "violated")]
        q, s = zw36.p36a_grade(rows, all_n((A, B), (A, C)), self.new)
        self.assertEqual(q, [])
        self.assertEqual(s["grade"], "untestable")

    def test_source_gate(self):
        rows = [trow(A, B, N1, "violated")]
        src = {zw36.key(A, B): "N", zw36.key(A, N1): "S"}
        self.assertEqual(zw36.p36a_grade(rows, src, self.new)[1]["grade"], "untestable")
        src[zw36.key(A, N1)] = {"decision": "N"}
        self.assertEqual(zw36.p36a_grade(rows, src, self.new)[1]["n_qualifying"], 1)

    def test_broken_forward(self):
        rows = [trow(A, B, N1, "violated"), trow(A, B, N2, "violated"),
                trow(A, B, N3, "holds")]
        src = all_n((A, B), (A, N1), (A, N2), (A, N3))
        s = zw36.p36a_grade(rows, src, self.new)[1]
        self.assertEqual((s["V"], s["H"]), (2, 1))
        self.assertEqual(s["grade"], "broken forward")

    def test_holds_forward(self):
        rows = [trow(A, B, N1, "holds"), trow(A, B, N2, "holds"),
                trow(A, B, N3, "untestable")]
        src = all_n((A, B), (A, N1), (A, N2), (A, N3))
        s = zw36.p36a_grade(rows, src, self.new)[1]
        self.assertEqual(s["grade"], "holds forward")
        self.assertEqual(s["n_qualifying"], 3)
        self.assertEqual(s["n_measured"], 2)

    def test_mixed_and_m_events(self):
        rows = [trow(A, B, N1, "violated"), trow(A, N1, N2, "holds")]
        src = all_n((A, B), (A, N1), (A, N2))
        s = zw36.p36a_grade(rows, src, self.new)[1]
        self.assertEqual(s["M_events"], [N1])
        self.assertEqual(s["grade"], "mixed")

    def test_implied_grades_not_measured(self):
        rows = [trow(A, B, N1, "violated-implied"), trow(A, B, N2, "holds-implied")]
        src = all_n((A, B), (A, N1), (A, N2))
        self.assertEqual(zw36.p36a_grade(rows, src, self.new)[1]["grade"], "untestable")


class TestHardwareAndP36b(unittest.TestCase):
    def test_hw_flag(self):
        lab = {(A, C): "changed", (B, N1): "changed", (C, N1): "same", (N1, N2): "same",
               (A, B): "nd"}
        self.assertTrue(zw36.hardware_flag(lab)[0])
        self.assertEqual(zw36.apply_hw("broken forward", True),
                         "broken forward, ILAR hardware not excluded")
        self.assertEqual(zw36.apply_hw("mixed", True), "mixed")
        lab[(A, B)] = "changed"   # within the old epoch: hardware no longer the only reading
        self.assertFalse(zw36.hardware_flag(lab)[0])
        lab = {(A, C): "changed", (B, C): "same"}
        self.assertFalse(zw36.hardware_flag(lab)[0])
        self.assertFalse(zw36.hardware_flag({(C, N1): "same"})[0])
        self.assertTrue(zw36.hardware_flag({(A, C): "changed"})[0])

    def test_p36b(self):
        new = {N1, N2, N3}
        src = all_n((C, N1), (C, N2), (C, N3))
        yka = {(C, N1): "same", (C, N2): "same", (C, N3): "nd"}
        self.assertEqual(zw36.p36b_grade(yka, src, new)[1]["grade"], "holds")
        yka = {(C, N1): "changed", (C, N2): "changed", (C, N3): "same"}
        self.assertEqual(zw36.p36b_grade(yka, src, new)[1]["grade"], "broken")
        yka = {(C, N1): "changed", (C, N2): "same", (C, N3): "same"}
        self.assertEqual(zw36.p36b_grade(yka, src, new)[1]["grade"], "mixed")
        yka = {(C, N1): "changed"}
        self.assertEqual(zw36.p36b_grade(yka, src, new)[1]["grade"], "untestable")
        src[zw36.key(C, N2)] = "S"
        yka = {(C, N1): "same", (C, N2): "same"}
        s = zw36.p36b_grade(yka, src, new)[1]
        self.assertEqual((s["n_source_N"], s["grade"]), (2, "untestable"))


class TestPlumbing(unittest.TestCase):
    def test_labels_and_closure_input(self):
        pr = [dict(band="1.0-2.0", array="ILAR", a=A, b=C, label_D_c310="changed"),
              dict(band="0.5-2.0", array="ILAR", a=A, b=C, label_D_c310="same"),
              dict(band="1.0-2.0", array="ILAR", a=A, b="2025-12-17", label_D_c310="same")]
        lab = zw36.labels_from_pairs(pr, "ILAR", "1.0-2.0", "D_c310", set(zw36.graded_dates()))
        self.assertEqual(lab, {(A, C): "changed"})
        ci = zw36.closure_input(pr, "1.0-2.0", "D_c310", set(zw36.graded_dates()))
        self.assertEqual(ci, [dict(channel="ILAR", event_i=A, event_j=C, label="changed")])

    def test_jsonable(self):
        o = zw36.jsonable({(A, B): {"x": float("nan"), "s": {2, 1}}})
        self.assertEqual(o, {A + "|" + B: {"x": None, "s": [1, 2]}})


if __name__ == "__main__":
    unittest.main()

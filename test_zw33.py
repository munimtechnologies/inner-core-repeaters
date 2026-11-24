# Copyright 2026 MV2, Munim, Inc.
# Licensed under the Apache License, Version 2.0 (Apache-2.0). See LICENSE.
# Made with AI agents working in MV2's name.
"""Unit tests for zw33.py (stdlib only)."""

import os
import sys
import unittest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import zw33  # noqa: E402

try:
    import closure35  # noqa: E402
except ImportError:  # the job folder and the repo both have it; skip only if absent
    closure35 = None

NEEDS_CLOSURE = unittest.skipIf(closure35 is None, "closure35.py not on the path")


def all_agree_labels():
    lab = {(a, x, y): e for a, x, y, e in zw33.E1_LABELS}
    dv = {}
    for arr, a, b, c in zw33.E1_RECOVERY:
        dv[(arr, a, b)] = 0.30
        dv[(arr, a, c)] = 0.10
    return lab, dv


class TestParsing(unittest.TestCase):
    def test_versions(self):
        self.assertEqual(zw33.VERSION, "0.1.0")

    @NEEDS_CLOSURE
    def test_rule_matches_closure35(self):
        self.assertEqual(closure35.VERSION, "0.1.0")
        self.assertEqual(zw33.SAME_MAX_RATIO, closure35.SAME_MAX_RATIO)
        self.assertEqual(zw33.CHANGED_MIN_RATIO, closure35.CHANGED_MIN_RATIO)
        self.assertEqual(zw33.FLOOR_MIN, closure35.FLOOR_MIN)
        for d, f in ((0.02, 0.01), (0.04, 0.01), (0.025, 0.01), (0.009, 0.001), (-1, 0.01)):
            self.assertEqual(zw33.label_pair(True, d, f), closure35.label_from_ratio(d, f))

    def test_event_folder(self):
        p = zw33.parse_event_folder("Event_2005_06_10_13_08_12")
        self.assertEqual(p["iso"], "2005-06-10T13:08:12")
        self.assertEqual(p["date"], "2005-06-10")
        self.assertIsNone(zw33.parse_event_folder("README.md"))
        self.assertIsNone(zw33.parse_event_folder("Event_2005_06_10"))

    def test_canonical_element(self):
        c = zw33.canonical_element
        self.assertEqual(c("YKAB0"), "YKB0")
        self.assertEqual(c("YKAR9"), "YKR9")
        self.assertEqual(c("YKB3"), "YKB3")
        self.assertEqual(c("YKR1"), "YKR1")
        self.assertIsNone(c("YKW1"))
        self.assertIsNone(c("YKAW3"))
        self.assertEqual(c("IL01"), "IL01")
        self.assertEqual(c("IL19"), "IL19")
        self.assertIsNone(c("IL31"))
        self.assertEqual(zw33.element_array("YKB0"), "YKA")
        self.assertEqual(zw33.element_array("IL07"), "ILAR")
        self.assertIsNone(zw33.element_array(None))

    def test_filename(self):
        p = zw33.parse_filename("CN.YKAB0..SHZ.sac")
        self.assertEqual((p["net"], p["sta"], p["loc"], p["cha"]), ("CN", "YKAB0", "", "SHZ"))
        self.assertIsNone(zw33.parse_filename("README.md"))

    def test_iso_epoch(self):
        self.assertEqual(zw33.iso_to_epoch("1970-01-01T00:01:00"), 60.0)
        self.assertAlmostEqual(zw33.iso_to_epoch("1970-01-01T00:00:01.500Z"), 1.5)
        self.assertEqual(zw33.epoch_to_iso(60.0), "1970-01-01T00:01:00")

    def test_pick_comcat(self):
        text = ("time,latitude,longitude,depth,mag,magType,id,place\n"
                "2005-06-10T13:08:40.000Z,-58.1,-26.2,40.0,5.0,mb,far1,\"South Sandwich Islands, region\"\n"
                "2005-06-10T13:08:12.590Z,-59.0,-25.7,35.0,5.6,mb,near1,\"South Sandwich Islands, region\"\n")
        c = zw33.pick_comcat(text, "2005-06-10T13:08:12")
        self.assertEqual(c["id"], "near1")
        self.assertEqual(c["lat"], -59.0)
        self.assertEqual(c["depth"], 35.0)
        self.assertAlmostEqual(c["dt"], 0.59, places=3)
        self.assertIsNone(zw33.pick_comcat(text, "2005-06-10T14:08:12"))
        self.assertIsNone(zw33.pick_comcat("", "2005-06-10T13:08:12"))


class TestLabels(unittest.TestCase):
    def test_label_pair(self):
        self.assertEqual(zw33.label_pair(True, 0.02, 0.01), "same")      # r = 2.0
        self.assertEqual(zw33.label_pair(True, 0.04, 0.01), "changed")   # r = 4.0
        self.assertEqual(zw33.label_pair(True, 0.025, 0.01), "nd")       # r = 2.5
        self.assertEqual(zw33.label_pair(False, 0.04, 0.01), "nd")       # gate failed
        self.assertEqual(zw33.label_pair(True, 0.009, 0.001), "same")    # floor -> 0.005, r = 1.8

    def test_recovery_result(self):
        self.assertEqual(zw33.recovery_result(0.30, 0.10)[0], "agree")     # 0.333
        self.assertEqual(zw33.recovery_result(0.30, 0.30)[0], "disagree")  # 1.0
        self.assertEqual(zw33.recovery_result(0.30, 0.24)[0], "nd")        # 0.8
        self.assertEqual(zw33.recovery_result(None, 0.10)[0], "nd")
        self.assertEqual(zw33.recovery_result(0.0, 0.10)[0], "nd")


class TestE1(unittest.TestCase):
    def test_all_agree(self):
        lab, dv = all_agree_labels()
        rows, s = zw33.e1_grade(lab, dv)
        self.assertEqual(s["n_graded"], 8)
        self.assertEqual(s["n_agree"], 8)
        self.assertEqual(s["headline"], "reproduced")
        self.assertEqual(len(rows), 10)   # 6 labels + 2 recoveries + 2 report-only

    def test_two_disagree(self):
        lab, dv = all_agree_labels()
        lab[("YKA", "2008-02-17", "2011-02-20")] = "same"
        lab[("ILAR", "2008-02-17", "2011-02-20")] = "changed"
        s = zw33.e1_grade(lab, dv)[1]
        self.assertEqual((s["n_agree"], s["n_disagree"], s["n_nd"]), (6, 2, 0))
        self.assertEqual(s["headline"], "not-reproduced")

    def test_one_nd_is_mostly(self):
        lab, dv = all_agree_labels()
        lab[("YKA", "2008-02-17", "2011-02-20")] = "nd"
        s = zw33.e1_grade(lab, dv)[1]
        self.assertEqual((s["n_agree"], s["n_disagree"], s["n_nd"]), (7, 0, 1))
        self.assertEqual(s["headline"], "mostly")

    def test_three_nd_is_inconclusive(self):
        lab, dv = all_agree_labels()
        del lab[("YKA", "2003-03-21", "2008-02-17")]
        lab[("ILAR", "2003-03-21", "2008-02-17")] = "nd"
        dv.pop(("ILAR", "2003-03-21", "2022-03-07"))
        s = zw33.e1_grade(lab, dv)[1]
        self.assertEqual((s["n_agree"], s["n_disagree"], s["n_nd"]), (5, 0, 3))
        self.assertEqual(s["headline"], "inconclusive")


class TestP35(unittest.TestCase):
    a, b, c = zw33.F_EVENTS

    def present(self):
        return {(arr, d): True for arr in ("YKA", "ILAR") for d in zw33.F_EVENTS}

    def labels(self, ilar_ac):
        a, b, c = self.a, self.b, self.c
        return {("YKA", a, b): "changed", ("YKA", a, c): "same", ("YKA", b, c): "changed",
                ("ILAR", a, b): "same", ("ILAR", a, c): ilar_ac, ("ILAR", b, c): "changed"}

    def test_violated(self):
        g = zw33.p35_grade(self.labels("changed"), self.present())
        self.assertTrue(g["G1"])
        self.assertEqual(g["grade"], "violated")

    def test_holds_names_failed_paper_label(self):
        g = zw33.p35_grade(self.labels("same"), self.present())
        self.assertEqual(g["grade"], "holds")
        self.assertEqual(len(g["paper_ilar_labels_not_reproduced"]), 0)
        lab = self.labels("same")
        lab[("ILAR", self.b, self.c)] = "same"
        g = zw33.p35_grade(lab, self.present())
        self.assertEqual(len(g["paper_ilar_labels_not_reproduced"]), 1)

    def test_nd_untestable(self):
        g = zw33.p35_grade(self.labels("nd"), self.present())
        self.assertEqual(g["grade"], "untestable")

    def test_g1_fails(self):
        lab = self.labels("changed")
        lab[("YKA", self.a, self.c)] = "changed"
        g = zw33.p35_grade(lab, self.present())
        self.assertFalse(g["G1"])
        self.assertEqual(g["grade"], "untestable")

    def test_missing_data(self):
        pr = self.present()
        pr[("ILAR", self.c)] = False
        g = zw33.p35_grade(self.labels("changed"), pr)
        self.assertEqual(g["grade"], "untestable")
        self.assertTrue(g["reason"].startswith("data"))

    @NEEDS_CLOSURE
    def test_closure_agrees_with_p35(self):
        rows = [dict(channel=k[0], event_i=k[1], event_j=k[2], label=v)
                for k, v in self.labels("changed").items()]
        s = closure35.grade(closure35.build_labels(rows))[1]
        self.assertEqual(s["headline"], "violated")
        rows = [dict(channel=k[0], event_i=k[1], event_j=k[2], label=v)
                for k, v in self.labels("same").items()]
        s = closure35.grade(closure35.build_labels(rows))[1]
        self.assertEqual(s["headline"], "holds")

    def test_closure_rows(self):
        pr = [dict(set="F", band="1.0-2.0", array="YKA", a=self.a, b=self.b, label_D_c310="same"),
              dict(set="F", band="0.5-2.0", array="YKA", a=self.a, b=self.b, label_D_c310="nd"),
              dict(set="M0324", band="1.0-2.0", array="YKA", a="x", b="y", label_D_c310="same")]
        out = zw33.closure_rows(pr, "F", "1.0-2.0", "D_c310")
        self.assertEqual(out, [dict(channel="YKA", event_i=self.a, event_j=self.b, label="same")])


if __name__ == "__main__":
    unittest.main()

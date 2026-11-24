# Copyright 2026 MV2, Munim, Inc.
# Licensed under the Apache License, Version 2.0 (Apache-2.0). See LICENSE.
# Made with AI agents working in MV2's name. Stdlib only; every expected value is hand-counted.
import csv
import json
import os
import tempfile
import unittest

import closure35 as C

HERE = os.path.dirname(os.path.abspath(__file__))

# Pre-registration 11 split of multiplet M03-24 (results/mantlep2-v1), seen at both arrays.
X = ["2003-03-21", "2020-07-11", "2022-03-07", "2024-03-21", "2025-12-13"]
Y = ["2008-02-17", "2011-02-20", "2014-12-09"]


def xy_rows(channels, flip=None):
    """All pairs at each channel: same within X or within Y, changed across.
    flip: optional (channel, a, b) whose label is inverted."""
    evs = sorted(X + Y)
    rows = []
    for ch in channels:
        for i, a in enumerate(evs):
            for b in evs[i + 1:]:
                lab = "same" if ((a in X) == (b in X)) else "changed"
                if flip and flip[0] == ch and {flip[1], flip[2]} == {a, b}:
                    lab = "changed" if lab == "same" else "same"
                rows.append(dict(channel=ch, event_i=a, event_j=b, label=lab))
    return rows


class TestRatioRule(unittest.TestCase):
    def test_thresholds(self):
        self.assertEqual(C.label_from_ratio(0.02, 0.01), "same")      # r = 2.0
        self.assertEqual(C.label_from_ratio(0.03, 0.01), "changed")   # r = 3.0
        self.assertEqual(C.label_from_ratio(0.025, 0.01), "nd")       # r = 2.5
        self.assertEqual(C.label_from_ratio(0.01, 0.001), "same")     # floor 0.005 -> r = 2
        self.assertEqual(C.label_from_ratio(0.016, 0.0), "changed")   # floor 0.005 -> r = 3.2
        self.assertEqual(C.label_from_ratio(float("nan"), 0.01), "nd")
        self.assertEqual(C.label_from_ratio(None, 0.01), "nd")
        self.assertEqual(C.label_from_ratio(-0.1, 0.01), "nd")

    def test_build_from_D(self):
        rows = [dict(channel="A", event_i="2005-01-01", event_j="2001-01-01", D="0.04",
                     D_floor="0.01")]
        lab = C.build_labels(rows)
        self.assertEqual(lab, {"A": {("2001-01-01", "2005-01-01"): "changed"}})

    def test_conflict_becomes_nd(self):
        rows = [dict(channel="A", event_i="2001-01-01", event_j="2002-01-01", label="same"),
                dict(channel="A", event_i="2002-01-01", event_j="2001-01-01", label="changed")]
        self.assertEqual(C.build_labels(rows)["A"][("2001-01-01", "2002-01-01")], "nd")


class TestClosure(unittest.TestCase):
    def test_xy_split_holds(self):
        # 12 recovery triples per channel (a = 2003-03-21, b in Y, c in X after b: 3 x 4),
        # 24 in all, every one confirmed by the other channel.
        rows, s = C.grade(C.build_labels(xy_rows(["ILAR", "YKA"])))
        self.assertEqual(s["n_triples"], 24)
        self.assertEqual(s["triples_holds"], 24)
        self.assertEqual(s["triples_violated"], 0)
        self.assertEqual(s["headline"], "holds")
        self.assertEqual(len(rows), 24)

    def test_one_channel_breaks_closure(self):
        # ILAR sees 2003-03-21 vs 2020-07-11 as changed. YKA triples with a = 2003-03-21 and
        # c = 2020-07-11 (b = 2008, 2011, 2014: 3 triples) are violated at ILAR; the other 9
        # YKA triples hold. ILAR has 12 triples of its own: b in Y with c in 2022, 2024, 2025
        # (9), plus b = 2020-07-11 (now changed) with c in 2022, 2024, 2025 (3); YKA
        # confirms all 12. So 24 triples: 3 violated, 21 hold -> 'mixed'.
        rows, s = C.grade(C.build_labels(xy_rows(["ILAR", "YKA"],
                                                 flip=("ILAR", "2003-03-21", "2020-07-11"))))
        self.assertEqual(s["n_triples"], 24)
        self.assertEqual(s["triples_violated"], 3)
        self.assertEqual(s["triples_holds"], 21)
        self.assertEqual(s["headline"], "mixed")

    def test_only_violations(self):
        rows = [dict(channel="A", event_i="2001-01-01", event_j="2002-01-01", label="changed"),
                dict(channel="A", event_i="2001-01-01", event_j="2003-01-01", label="same"),
                dict(channel="B", event_i="2001-01-01", event_j="2003-01-01", label="changed")]
        out, s = C.grade(C.build_labels(rows))
        self.assertEqual(s["headline"], "violated")
        self.assertEqual(s["triples_violated"], 1)
        self.assertEqual(out[0]["grade"], "violated")

    def test_paper_labels_multiF(self):
        # Zhang & Wen (2026) labels for multiplet F as printed: one recovery triple at YKAR
        # (2005-06-10, 2012-12-03, 2022-05-19); ILAR has no measured 2005-2022 label, but
        # same(2005, 2012) and changed(2012, 2022) imply changed(2005, 2022).
        with open(os.path.join(HERE, "labels", "closure35", "zw26_multiF_paper_labels.csv"), newline="",
                  encoding="utf-8") as f:
            rows = list(csv.DictReader(f))
        out, s = C.grade(C.build_labels(rows))
        self.assertEqual(s["n_triples"], 1)
        self.assertEqual(s["triples_violated_implied"], 1)
        self.assertEqual(s["headline"], "untestable")   # no measured ILAR 2005-2022 label yet
        self.assertEqual(out[0]["recovery_channel"], "YKAR")
        self.assertEqual(out[0]["other_channel"], "ILAR")
        self.assertEqual(out[0]["grade"], "violated-implied")

    def test_run_writes_files(self):
        with tempfile.TemporaryDirectory() as d:
            p = os.path.join(d, "pairs.csv")
            with open(p, "w", newline="", encoding="utf-8") as f:
                w = csv.DictWriter(f, fieldnames=["channel", "event_i", "event_j", "label"])
                w.writeheader()
                for r in xy_rows(["ILAR", "YKA"]):
                    w.writerow(r)
            o = os.path.join(d, "out")
            self.assertEqual(C.main(["--pairs", p, "--out-dir", o]), 0)
            with open(os.path.join(o, "closure35_summary.json"), encoding="utf-8") as f:
                s = json.load(f)
            self.assertEqual(s["headline"], "holds")
            with open(os.path.join(o, "closure35_triples.csv"), newline="",
                      encoding="utf-8") as f:
                self.assertEqual(len(list(csv.DictReader(f))), 24)


if __name__ == "__main__":
    unittest.main()

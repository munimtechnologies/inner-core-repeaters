# Copyright 2026 MV2, Munim, Inc.
# Licensed under the Apache License, Version 2.0 (Apache-2.0). See LICENSE.
# Made with AI agents working in MV2's name.
"""Stdlib-only tests for wsc39.py (pre-registration 39, amendment 1). Expected values were derived by
hand before any code ran; see findings/2026-11-27-prereg39-wang-source-census.md and
findings/2026-11-27-prereg39-amendment1-multiplet-A.md."""

import unittest

import wsc39 as w


class TestTables(unittest.TestCase):
    def test_version(self):
        self.assertEqual(w.VERSION, "0.2.0")

    def test_members_match_comcat_table(self):
        members = {m for ms in w.MULTIPLETS.values() for m in ms}
        self.assertEqual(len(members), 65)
        located = {k for k, v in w.WANG_EVENTS.items() if v}
        self.assertEqual(len(located), 64)
        self.assertEqual(members - {711}, located)
        self.assertIsNone(w.WANG_EVENTS[711])

    def test_amendment_multiplet_A(self):
        # Table 3's own pair list for A (P40, P52, P102, P104, P124, P129) uses 734, never 724.
        self.assertIn(734, w.MULTIPLETS["A"])
        self.assertNotIn(724, w.MULTIPLETS["A"])
        self.assertNotIn(724, w.WANG_EVENTS)
        self.assertEqual(w.cid(734), "usp000fz4d")
        self.assertEqual(w.located("A"), [712, 718, 734, 742, 848, 854, 859])

    def test_ids_unique(self):
        ids = [v[0] for v in w.WANG_EVENTS.values() if v]
        self.assertEqual(len(ids), len(set(ids)))

    def test_graded_set(self):
        self.assertEqual(w.GRADED, ("A", "B", "C", "D", "E", "F", "G", "H", "J", "K", "L", "N", "O", "P"))

    def test_located_order(self):
        self.assertEqual(w.located("K"), [705, 729])
        self.assertEqual(w.located("F"), [709, 730, 744, 856])
        self.assertEqual(w.located("G"), [802, 715, 830, 758, 860])
        self.assertEqual(w.located("N"), [707, 809, 756])

    def test_pairs(self):
        # A has five M03-24 events (718, 734, 742, 848, 854): 21 pairs minus 10 yardstick pairs.
        self.assertEqual(len(w.member_pairs("A")), 11)
        self.assertNotIn((718, 742), w.member_pairs("A"))
        self.assertNotIn((718, 734), w.member_pairs("A"))
        self.assertIn((712, 718), w.member_pairs("A"))
        self.assertIn((734, 859), w.member_pairs("A"))
        self.assertEqual(len(w.member_pairs("K")), 1)
        self.assertEqual(sum(len(w.member_pairs(k)) for k in w.MULTIPLETS), 99)

    def test_old8_overlap(self):
        a = [w.cid(m) for m in w.MULTIPLETS["A"]]
        self.assertEqual(sum(1 for x in a if x in w.OLD8), 5)
        f = [w.cid(m) for m in (730, 744, 856)]
        self.assertEqual(f, ["usp000dsnt", "usp000jwcf", "us6000hm6g"])

    def test_catalogue_sep(self):
        km, dd = w.catalogue_sep(730, 744)
        self.assertTrue(26.5 < km < 27.7, km)
        self.assertAlmostEqual(dd, 31.4, places=6)
        km, dd = w.catalogue_sep(718, 734)
        self.assertTrue(15.5 < km < 16.6, km)
        self.assertAlmostEqual(dd, 33.5, places=6)


class TestScores(unittest.TestCase):
    def setUp(self):
        self.ref10 = [0.1 * i for i in range(1, 11)]
        self.g = {"a": True, "c": True}

    def test_q_and_hit(self):
        r = w.station_pair_score(self.ref10, {("a", "c"): 0.95}, self.g, "a", "c")
        self.assertTrue(r["counted"])
        self.assertAlmostEqual(r["q"], 0.9)
        self.assertFalse(r["hit"])
        r = w.station_pair_score(self.ref10, {("a", "c"): 1.05}, self.g, "c", "a")
        self.assertAlmostEqual(r["q"], 1.0)
        self.assertTrue(r["hit"])
        ref15 = [float(i) for i in range(1, 16)]
        self.assertTrue(w.station_pair_score(ref15, {("a", "c"): 14.5}, self.g, "a", "c")["hit"])
        self.assertFalse(w.station_pair_score(ref15, {("a", "c"): 13.5}, self.g, "a", "c")["hit"])

    def test_not_counted(self):
        r = w.station_pair_score(self.ref10[:9], {("a", "c"): 2.0}, self.g, "a", "c")
        self.assertFalse(r["counted"])
        self.assertEqual(r["reason"], "too_few_reference_pairs")
        r = w.station_pair_score(self.ref10, {("a", "c"): 2.0}, {"a": True, "c": False}, "a", "c")
        self.assertEqual(r["reason"], "member_not_gated")
        r = w.station_pair_score(self.ref10, {}, self.g, "a", "c")
        self.assertEqual(r["reason"], "no_pair_D")

    def test_reference_values(self):
        D = {w.key(w.OLD8[0], w.OLD8[1]): 0.1, w.key(w.OLD8[0], w.OLD8[2]): 0.2,
             w.key(w.OLD8[1], w.OLD8[2]): float("nan")}
        g = {w.OLD8[0]: True, w.OLD8[1]: True, w.OLD8[2]: True}
        self.assertEqual(sorted(w.reference_values(D, g)), [0.1, 0.2])
        g[w.OLD8[2]] = False
        self.assertEqual(w.reference_values(D, g), [0.1])


class TestDecisions(unittest.TestCase):
    def test_decide(self):
        T, F = True, False
        self.assertEqual(w.decide([T, T, T, F, F, F, F]), ("inconclusive", 3))
        self.assertEqual(w.decide([T, T, T, T, F, F, F]), ("S", 4))
        self.assertEqual(w.decide([T, F, F, F, F, F, F]), ("N", 1))
        self.assertEqual(w.decide([T, T]), ("untestable", 2))
        self.assertEqual(w.decide([T, T, F]), ("S", 2))
        self.assertEqual(w.decide([T, F, F]), ("N", 1))

    def test_verdict(self):
        self.assertEqual(w.multiplet_verdict({1: "S", 2: "N", 3: "S"}), "odd")
        self.assertEqual(w.multiplet_verdict({1: "N", 2: "N", 3: "inconclusive"}), "clean")
        self.assertEqual(w.multiplet_verdict({1: "N", 2: "inconclusive", 3: "inconclusive"}), "unresolved")
        self.assertEqual(w.multiplet_verdict({1: "untestable", 2: "untestable"}), "untestable")
        self.assertEqual(w.multiplet_verdict({1: "N", 2: "untestable"}), "clean")
        self.assertEqual(w.multiplet_verdict({}), "untestable")

    def test_odd_members(self):
        self.assertEqual(w.odd_members({(730, 744): "S", (730, 856): "N", (744, 856): "S"}), [744])
        self.assertEqual(w.odd_members({(1, 2): "S", (1, 3): "S", (2, 3): "S"}), [])
        d = {(1, 2): "N", (1, 3): "N", (2, 3): "N", (1, 4): "S", (2, 4): "S", (3, 4): "inconclusive"}
        self.assertEqual(w.odd_members(d), [4])

    def test_census(self):
        v = {k: "untestable" for k in w.GRADED}
        v.update(A="odd", F="odd", B="clean", C="clean", D="clean", E="clean")
        self.assertEqual(w.census(v)[0], "some source problems")
        self.assertEqual(w.census(v)[2], 6)
        v["A"] = "clean"
        self.assertEqual(w.census(v)[0], "source-clean")
        v = {k: "untestable" for k in w.GRADED}
        v.update(A="odd", F="odd", B="odd", C="clean", D="clean", E="clean", G="clean", H="clean")
        self.assertEqual(w.census(v)[0], "source problems common")
        for k in ("J", "K", "L", "N", "O"):
            v[k] = "clean"
        self.assertEqual(w.census(v)[2], 13)
        self.assertEqual(w.census(v)[0], "some source problems")
        v = {k: "untestable" for k in w.GRADED}
        v.update(A="clean", B="clean", C="odd", D="unresolved")
        self.assertEqual(w.census(v)[0], "untestable")

    def test_control(self):
        self.assertEqual(w.control_status("odd"), "passed")
        self.assertEqual(w.control_status("clean"), "failed")
        self.assertEqual(w.control_status("unresolved"), "inconclusive")
        self.assertEqual(w.final_headline("source-clean", "failed"), "source-clean [positive control failed]")
        self.assertEqual(w.final_headline("source-clean", "passed"), "source-clean")

    def test_replication(self):
        f = dict(pairs={(730, 744): {"decision": "S"}, (730, 856): {"decision": "N"},
                        (744, 856): {"decision": "inconclusive"}})
        r = w.replication(f)
        self.assertTrue(r["730-744"]["agree"])
        self.assertTrue(r["730-856"]["agree"])
        self.assertFalse(r["744-856"]["agree"])

    def test_grade_multiplet_synthetic(self):
        T, F = True, False
        hits = {(716, 737): [T, T, F], (716, 852): [F, F, F], (737, 852): [T, T, T]}
        scores = {}
        for i, st in enumerate(("X.A", "X.B", "X.C")):
            scores[st] = {p: dict(counted=True, q=(1.0 if h[i] else 0.5), hit=h[i]) for p, h in hits.items()}
        g = w.grade_multiplet("J", scores)
        self.assertEqual(g["pairs"][(716, 737)]["decision"], "S")
        self.assertEqual(g["pairs"][(716, 852)]["decision"], "N")
        self.assertEqual(g["pairs"][(737, 852)]["decision"], "S")
        self.assertEqual(g["verdict"], "odd")
        self.assertEqual(g["odd_members"], [737])

    def test_summarise_missing_control(self):
        graded = {k: dict(verdict="clean") for k in w.GRADED}
        s = w.summarise(graded)
        self.assertEqual(s["headline"], "source-clean")
        self.assertEqual(s["positive_control"]["status"], "inconclusive")
        self.assertEqual(s["n_testable"], 14)


if __name__ == "__main__":
    unittest.main()

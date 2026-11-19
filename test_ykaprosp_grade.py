"""Tests for ykaprosp_grade.py on synthetic families (no network, stdlib only).

Part of inner-core-repeaters by MV2, Munim, Inc. (made with AI agents). Licence: Apache-2.0.
"""
import csv
import json
import os
import shutil
import tempfile
import unittest

import ykaprosp_grade as g

COLS = ["array", "mult", "ev1", "ev2", "time1", "time2", "t1", "t2", "snr_dir1", "snr_dir2", "ok", "D", "D_c13",
        "D_c310"]


def iso(y, mo, d, hh=12):
    return "%04d-%02d-%02dT%02d:00:00.000Z" % (y, mo, d, hh)


class World:
    """Synthetic South Sandwich families. Per family f: pre-turn o1 (1996), o2 (2002), o3 (2007);
    mid (2012); n25 (2025); new events n1 (2027-06) and n2 (2028-03). For n1 and n2 the backtracking
    clock gives ~2000.9 and ~2000.6, so the B pick is o2 and the M pick is o3."""

    def __init__(self, n_fam=4, d_b=0.05, d_m=0.30, d_o1=0.20, d_mid=0.30, d_near=0.04, extras=True):
        self.events = {}  # id -> (time, lat, lon)
        self.rows = []
        for f in range(n_fam):
            day = 1 + f
            lat0, lon0 = -56.0 - 0.5 * f, -27.0
            ev = dict(o1=iso(1996, 6, day), o2=iso(2002, 6, day), o3=iso(2007, 6, day), mid=iso(2012, 6, day),
                      n25=iso(2025, 6, day), n1=iso(2027, 6, day), n2=iso(2028, 3, day))
            for k, (name, t) in enumerate(sorted(ev.items())):
                self.events["F%d_%s" % (f, name)] = (t, lat0 + 0.01 * k, lon0)
            dval = dict(o1=d_o1, o2=d_b, o3=d_m, mid=d_mid, n25=d_near)
            for new in ("n1", "n2"):
                for part in ("o1", "o2", "o3", "mid", "n25"):
                    self.pair(f, part, new, dval[part])
            self.pair(f, "n1", "n2", d_near)
            self.pair(f, "o1", "o2", 0.01)  # not prospective: ignored
        if extras:
            lat0 = -56.0
            self.events["F0_far"] = (iso(2004, 6, 1), lat0 + 0.3, -27.0)  # ~33 km away
            self.events["F0_lowsnr"] = (iso(2003, 6, 1), lat0, -27.0)
            self.events["F1_x"] = (iso(2004, 6, 2), -56.5, -27.0)
            self.pair(0, "far", "n1", 0.0)
            self.pair(0, "lowsnr", "n1", 0.0, snr2=11.0)
            self.pair(1, "x", "n2", 0.0, ok=False)
            # partner not in the catalogue, at a time no catalogue event has
            self.rows.append(self.row("F2", "ghost", "F2_n1", "2005-03-03T03:03:03.000Z", self.events["F2_n1"][0],
                                      0.0, 20.0, 20.0, True))

    def row(self, mult, e1, e2, t1, t2, d, s1, s2, ok):
        return dict(array="YKA", mult=mult, ev1=e1, ev2=e2, time1=t1, time2=t2,
                    t1=g.decyear(g.epoch(t1)), t2=g.decyear(g.epoch(t2)), snr_dir1=s1, snr_dir2=s2,
                    ok=str(bool(ok)), D=d, D_c13=d, D_c310=d)

    def pair(self, f, a, b, d, snr2=20.0, ok=True):
        ea, eb = "F%d_%s" % (f, a), "F%d_%s" % (f, b)
        self.rows.append(self.row("F%d" % f, ea, eb, self.events[ea][0], self.events[eb][0], d, 20.0, snr2, ok))

    def write(self, d):
        pairs = os.path.join(d, "yka_pairs.csv")
        with open(pairs, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=COLS)
            w.writeheader()
            for r in self.rows:
                w.writerow(r)
        cat = os.path.join(d, "comcat.csv")
        with open(cat, "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["time", "latitude", "longitude", "depth", "mag", "id"])
            for k, (t, la, lo) in sorted(self.events.items()):
                w.writerow([t, la, lo, 30.0, 5.0, k])
        return pairs, cat


class Base(unittest.TestCase):
    def setUp(self):
        self.d = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.d, ignore_errors=True)

    def grade(self, world, look="final", nboot=500):
        pairs, cat = world.write(self.d)
        with open(cat) as fh:
            c = g.read_catalog(fh.read())
        return g.grade_all(g.read_pairs(pairs), c, look, nboot=nboot)


class TestHelpers(unittest.TestCase):
    def test_decyear_matches_yka_b(self):
        # results/ykab-v2/yka_pairs.csv row 1: 2001-05-29T06:22:55.100Z -> 2001.4062079876967
        self.assertAlmostEqual(g.decyear(g.epoch("2001-05-29T06:22:55.100Z")), 2001.4062079876967, places=7)
        self.assertAlmostEqual(g.decyear(g.epoch("2002-02-02T07:29:27.270Z")), 2002.0885263593989, places=7)

    def test_epoch_forms(self):
        a = g.epoch("2026-11-19T00:00:00Z")
        self.assertEqual(a, g.epoch("2026-11-19T00:00:00+00:00"))
        self.assertEqual(a, g.epoch("2026-11-19"))
        self.assertAlmostEqual(g.epoch("2026-11-19T00:00:00.25Z") - a, 0.25)

    def test_equiv(self):
        self.assertAlmostEqual(g.equiv(2027.5), 2000.9)
        self.assertEqual(g.equiv(2001.0), 2001.0)

    def test_binom(self):
        self.assertAlmostEqual(g.binom_tail(8, 8), 1 / 256.0)
        self.assertEqual(g.binom_tail(0, 0), 1.0)

    def test_look_gate(self):
        self.assertFalse(g.look_open("interim", g.epoch("2027-11-18T23:59:59Z")))
        self.assertTrue(g.look_open("interim", g.epoch("2027-11-19T00:00:00Z")))
        self.assertFalse(g.look_open("final", g.epoch("2028-11-18T00:00:00Z")))

    def test_sep(self):
        self.assertAlmostEqual(g.sep_km(-56.0, -27.0, -55.7, -27.0), 33.36, places=1)


class TestGrade(Base):
    def test_supported(self):
        res = self.grade(World())
        pr = res["primary"]
        self.assertEqual(pr["n_prospective"], 48)
        self.assertEqual(pr["n_selected"], 44)
        self.assertEqual(pr["excluded"], {"separation": 1, "snr": 1, "not_ok": 1, "no_location": 1})
        self.assertEqual(pr["P4"]["n_near"], 12)
        self.assertEqual(pr["P4"]["n_far"], 16)
        self.assertTrue(pr["P4"]["holds"])
        f1 = pr["F1"]
        self.assertEqual((f1["n_discriminating"], f1["n_families"], f1["n_match_B"], f1["n_match_M"]), (8, 4, 8, 0))
        self.assertLess(f1["p_binom_B"], 0.025)
        self.assertEqual(f1["boot_frac_B_ahead"], 1.0)
        self.assertEqual(pr["decision"], "supported")
        for c in res["cases"]:
            self.assertEqual(c["B_pick"][-2:], "o2")
            self.assertEqual(c["M_pick"][-2:], "o3")

    def test_contradicted(self):
        res = self.grade(World(d_b=0.20, d_m=0.02))
        pr = res["primary"]
        self.assertTrue(pr["P4"]["holds"])
        self.assertEqual(pr["F1"]["n_match_M"], 8)
        self.assertEqual(pr["decision"], "contradicted")

    def test_p4_fails_blocks_vote(self):
        res = self.grade(World(d_near=0.5))
        pr = res["primary"]
        self.assertFalse(pr["P4"]["holds"])
        self.assertEqual(pr["F1"]["verdict"], "supported")
        self.assertEqual(pr["decision"], "not interpreted (P4-prosp fails)")

    def test_untestable_interim_and_final(self):
        res = self.grade(World(n_fam=2, extras=False), look="interim")
        self.assertEqual(res["primary"]["F1"]["n_discriminating"], 4)
        self.assertEqual(res["primary"]["decision"], "not graded at interim (vote untestable)")
        res = self.grade(World(n_fam=2, extras=False), look="final")
        self.assertEqual(res["primary"]["decision"], "untestable")

    def test_inconclusive_when_mixed(self):
        # Two families favour B, two favour M: 4 vs 4.
        w = World()
        for r in w.rows:
            if r["mult"] in ("F2", "F3") and r["ev2"].endswith(("n1", "n2")):
                if r["ev1"].endswith("o2"):
                    r["D"] = r["D_c13"] = 0.20
                elif r["ev1"].endswith("o3"):
                    r["D"] = r["D_c13"] = 0.02
        pr = self.grade(w)["primary"]
        self.assertEqual((pr["F1"]["n_match_B"], pr["F1"]["n_match_M"]), (4, 4))
        self.assertEqual(pr["decision"], "inconclusive")

    def test_extras(self):
        res = self.grade(World())
        self.assertFalse(res["extras"]["D"]["graded"])
        self.assertEqual(res["extras"]["D"]["decision"], "supported")
        pool = res["extras"]["pool_drop"]
        self.assertEqual(pool["P4"]["n_near"], 4)  # only the n1-n2 pairs remain near
        self.assertEqual(pool["P4"]["verdict"], "untestable")


class TestMain(Base):
    def test_too_early(self):
        pairs, cat = World().write(self.d)
        out = os.path.join(self.d, "out")
        code = g.main(["--pairs", pairs, "--catalog", cat, "--out-dir", out, "--look", "interim",
                       "--now", "2027-06-01T00:00:00Z"])
        self.assertEqual(code, 2)
        with open(os.path.join(out, "ykaprosp_grade.json")) as fh:
            self.assertTrue(json.load(fh)["status"].startswith("too-early"))
        self.assertFalse(os.path.exists(os.path.join(out, "ykaprosp_cases.csv")))

    def test_final_run_writes_outputs(self):
        pairs, cat = World().write(self.d)
        out = os.path.join(self.d, "out")
        code = g.main(["--pairs", pairs, "--catalog", cat, "--out-dir", out, "--look", "final",
                       "--now", "2028-11-20T00:00:00Z", "--nboot", "300"])
        self.assertEqual(code, 0)
        with open(os.path.join(out, "ykaprosp_grade.json")) as fh:
            res = json.load(fh)
        self.assertEqual(res["status"], "graded")
        self.assertEqual(res["primary"]["decision"], "supported")
        with open(os.path.join(out, "ykaprosp_cases.csv"), newline="") as fh:
            self.assertEqual(len(list(csv.DictReader(fh))), 8)
        with open(os.path.join(out, "ykaprosp_pairs.csv"), newline="") as fh:
            self.assertEqual(len(list(csv.DictReader(fh))), 44)


if __name__ == "__main__":
    unittest.main()

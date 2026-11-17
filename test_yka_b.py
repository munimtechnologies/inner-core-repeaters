"""Tests for yka_b.py (stdlib only). MV2, Munim, Inc. (made with AI agents). Apache-2.0."""
import csv
import os
import tempfile
import unittest

import yka_b


def row(mult, ev1, t1, ev2, t2, D, ok=True, snr1=10.0, Dx=None):
    nan = float("nan")
    return dict(mult=mult, ev1=ev1, ev2=ev2, t1=t1, t2=t2, D=D, Dx=D if Dx is None else Dx, D_c13=D, D_c310=D,
                D_full=D, cc_ref=0.95, snr_dir1=snr1, snr_dir2=10.0, snr_coda1=3.0, snr_coda2=3.0, ok=ok,
                ignored=nan)


class TestBasics(unittest.TestCase):
    def test_equiv(self):
        self.assertAlmostEqual(yka_b.equiv(2026.0), 2008.5 - 17.5 / 2.5)
        self.assertEqual(yka_b.equiv(2003.0), 2003.0)

    def test_parse(self):
        s = "a1(2001-01-01);b2(2025-02-03)"
        self.assertEqual(yka_b.parse_members(s), ["a1", "b2"])
        self.assertEqual(yka_b.member_dates(s)["b2"], "(2025-02-03)")

    def test_binom(self):
        self.assertAlmostEqual(yka_b.binom_tail(6, 6), 1 / 64)
        self.assertAlmostEqual(yka_b.binom_tail(0, 5), 1.0)
        self.assertEqual(yka_b.binom_tail(0, 0), 1.0)


class TestMerge(unittest.TestCase):
    def test_union(self):
        fams = [dict(name="I1", ids=["a", "b"], dates={"a": "(2001)", "b": "(2025)"}),
                dict(name="Z1", ids=["b", "c"], dates={"b": "(2025)", "c": "(2011)"}),
                dict(name="F1", ids=["d", "e"], dates={})]
        m = yka_b.merge_families(fams)
        self.assertEqual(len(m), 2)
        self.assertEqual(sorted(m[0]["ids"]), ["a", "b", "c"])
        self.assertEqual(m[0]["sources"], ["I1", "Z1"])
        self.assertEqual(m[0]["mult"], "Y000")

    def test_write_and_read_back(self):
        fams = [dict(name="I1", ids=["a", "b"], dates={"a": "(2001-01-01)", "b": "(2025-01-01)"})]
        m = yka_b.merge_families(fams)
        with tempfile.TemporaryDirectory() as d:
            p = os.path.join(d, "m.csv")
            yka_b.write_merged(p, m)
            back = yka_b.read_families(p, "X")
            self.assertEqual(back[0]["ids"], ["a", "b"])
            self.assertEqual(back[0]["name"], "XY000")

    def test_ok_pairs_first_source_wins(self):
        with tempfile.TemporaryDirectory() as d:
            p1, p2, out = (os.path.join(d, x) for x in ("1.csv", "2.csv", "o.csv"))
            for p, D in ((p1, "0.1"), (p2, "0.9")):
                with open(p, "w", newline="") as f:
                    w = csv.writer(f)
                    w.writerow(["mult", "ev1", "ev2", "D", "ok"])
                    w.writerow(["M", "a", "b", D, "True"])
                    w.writerow(["M", "a", "c", D, "False"])
            self.assertEqual(yka_b.merge_ok_pairs([p1, p2], out), 1)
            with open(out, newline="") as f:
                r = list(csv.DictReader(f))
            self.assertEqual(float(r[0]["D"]), 0.1)


class TestVote(unittest.TestCase):
    def test_case_logic(self):
        # new event 2026.0 -> e = 2001.5; partners 2001.6 (B pick) and 2007.0 (M pick)
        rows = [row("Y1", "p01", 2001.6, "n1", 2026.0, 0.05), row("Y1", "p07", 2007.0, "n1", 2026.0, 0.20),
                row("Y1", "p01", 2001.6, "p07", 2007.0, 0.3)]
        c = yka_b.cases_from_pairs(rows)
        self.assertEqual(len(c), 1)
        self.assertTrue(c[0]["discriminating"])
        self.assertTrue(c[0]["matches_B"])
        self.assertFalse(c[0]["matches_M"])
        self.assertAlmostEqual(c[0]["chance"], 0.5)

    def test_not_ok_and_single_partner_ignored(self):
        rows = [row("Y1", "p01", 2001.6, "n1", 2026.0, 0.05, ok=False), row("Y1", "p07", 2007.0, "n1", 2026.0, 0.2)]
        self.assertEqual(yka_b.cases_from_pairs(rows), [])

    def _cases(self, pattern):
        rows = []
        for k, b_wins in enumerate(pattern):
            m = "Y%d" % k
            n = "n%d" % k
            rows.append(row(m, "a%d" % k, 2001.5, n, 2026.0, 0.05 if b_wins else 0.3))
            rows.append(row(m, "b%d" % k, 2007.5, n, 2026.0, 0.3 if b_wins else 0.05))
        return yka_b.cases_from_pairs(rows)

    def test_supported(self):
        g = yka_b.grade_vote(self._cases([True] * 7), nboot=200)
        self.assertTrue(g["testable"])
        self.assertEqual(g["verdict"], "supported")
        self.assertEqual(g["n_match_B"], 7)

    def test_contradicted(self):
        g = yka_b.grade_vote(self._cases([False] * 7), nboot=200)
        self.assertEqual(g["verdict"], "contradicted")

    def test_inconclusive_and_untestable(self):
        self.assertEqual(yka_b.grade_vote(self._cases([True, False] * 4), nboot=200)["verdict"], "inconclusive")
        self.assertEqual(yka_b.grade_vote(self._cases([True] * 3), nboot=200)["verdict"], "untestable")

    def test_grade_all_from_csv(self):
        rows = []
        for k in range(7):
            rows.append(row("Y%d" % k, "a%d" % k, 2001.5, "n%d" % k, 2026.0, 0.05))
            rows.append(row("Y%d" % k, "b%d" % k, 2007.5, "n%d" % k, 2026.0, 0.3))
        rows.append(row("Y0", "a0", 2001.5, "us6000runq", 2025.95, 0.1))
        with tempfile.TemporaryDirectory() as d:
            p = os.path.join(d, "yka_pairs.csv")
            with open(p, "w", newline="") as f:
                w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
                w.writeheader()
                for r in rows:
                    w.writerow(r)
            res, cases = yka_b.grade_all(p, nboot=100)
        self.assertEqual(res["F1_Y"]["verdict"], "supported")
        self.assertTrue(res["F1_Y_robust"])
        self.assertEqual(len(res["anomaly_event_pairs"]), 1)

    def test_empty_file(self):
        with tempfile.TemporaryDirectory() as d:
            p = os.path.join(d, "yka_pairs.csv")
            open(p, "w").close()
            res, cases = yka_b.grade_all(p, nboot=10)
        self.assertEqual(res["F1_Y"]["verdict"], "untestable")


if __name__ == "__main__":
    unittest.main()

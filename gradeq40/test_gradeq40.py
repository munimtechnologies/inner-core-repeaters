# Copyright 2026 MV2, Munim, Inc.
# Licensed under the Apache License, Version 2.0 (Apache-2.0). See LICENSE.
# Made with AI agents working in MV2's name. Expected values below were derived by hand.
"""Tests for gradeq40.py 0.1.1 (pre-registration 40). Standard library only, no network."""
import csv
import io
import math
import unittest

import gradeq40 as g


def row(**kw):
    base = dict(array="YKA", mult="Z1", ev1="a", ev2="b", t1="2021.0", t2="2023.0",
                m0324="False", S="0.2", n_st="3", D="0.01", D_c310="0.05")
    base.update({k: str(v) for k, v in kw.items()})
    return base


class TestBasics(unittest.TestCase):
    def test_version(self):
        self.assertEqual(g.VERSION, "0.1.1")

    def test_blob_sha_known(self):
        # git hash-object of an empty file and of "hello\n"
        self.assertEqual(g.blob_sha(b""), "e69de29bb2d1d6434b8b29ae775ad8c2e48c5391")
        self.assertEqual(g.blob_sha(b"hello\n"), "ce013625030ba8dba906f756967f9e9ca394464a")

    def test_decimal_year(self):
        # 2003-03-21 is day 80: (79 + 10.79656/24) / 365 = 0.2176708...
        self.assertAlmostEqual(g.decimal_year("2003-03-21T10:47:47.620Z"), 2003.2176708403, places=7)
        # 2020-04-08 is day 99 of a leap year: 98/366
        self.assertAlmostEqual(g.POST_START, 2020 + 98 / 366.0, places=9)
        # 2021-12-01 is day 335: 334/365 ; 2022-04-01 is day 91: 90/365
        self.assertAlmostEqual(g.ERA2[0], 2021 + 334 / 365.0, places=9)
        self.assertAlmostEqual(g.ERA2[1], 2022 + 90 / 365.0, places=9)

    def test_ranks_ties(self):
        self.assertEqual(g.ranks([10, 20, 20, 30]), [1.0, 2.5, 2.5, 4.0])
        self.assertEqual(g.ranks([3, 1, 2]), [3.0, 1.0, 2.0])

    def test_spearman(self):
        self.assertAlmostEqual(g.spearman([1, 2, 3, 4], [10, 20, 30, 40]), 1.0)
        self.assertAlmostEqual(g.spearman([1, 2, 3], [3, 2, 1]), -1.0)
        self.assertTrue(math.isnan(g.spearman([1, 2, 3], [5, 5, 5])))
        self.assertTrue(math.isnan(g.spearman([1, 2], [1, 2])))
        # x = 1..4, y = 1,3,2,4: d = 0,-1,1,0 -> rho = 1 - 6*2/(4*15) = 0.8
        self.assertAlmostEqual(g.spearman([1, 2, 3, 4], [1, 3, 2, 4]), 0.8)

    def test_nearest_rank_and_median(self):
        v = list(range(1, 21))  # 1..20
        self.assertEqual(g.nearest_rank(v, 0.95), 19)  # ceil(19.0) = 19th value
        self.assertEqual(g.nearest_rank(v, 0.05), 1)
        self.assertEqual(g.median([3, 1, 2]), 2)
        self.assertEqual(g.median([4, 1, 3, 2]), 2.5)


class TestSelect(unittest.TestCase):
    def test_post_filters(self):
        rows = [
            row(),                                   # kept: dt 2.0
            row(t1="2022.0"),                        # in era 2 -> out
            row(t1="2020.1"),                        # before 2020-04-08 (2020.2678) -> out
            row(S="0.5"),                            # S > 1/3 -> out
            row(S=""),                               # no S -> out
            row(n_st="1"),                           # one station -> out
            row(m0324="True"),                       # M03-24 -> out by default
            row(array="ILAR"),                       # other array -> out
            row(D_c310="nan"),                       # no metric -> out
            row(t1="2020.5", t2="2026.5", mult="Z2"),  # kept: dt 6.0
        ]
        p = g.select(rows, "post")
        self.assertEqual([(x["mult"], x["dt"]) for x in p], [("Z1", 2.0), ("Z2", 6.0)])
        p2 = g.select(rows, "post", include_m0324=True)
        self.assertEqual(len(p2), 3)
        p3 = g.select(rows, "post", metric="D")
        self.assertEqual([x["v"] for x in p3], [0.01, 0.01, 0.01])  # nan D_c310 row has D = 0.01

    def test_array_case(self):
        # 0.1.1: the pinned table writes "yka"; any case of YKA is kept, ILAR in any case is not.
        rows = [row(array="yka"), row(array=" YKA "), row(array="Yka"),
                row(array="ilar"), row(array="ILAR"), row(array="")]
        self.assertEqual(len(g.select(rows, "post")), 3)
        self.assertEqual(len(g.select([row(array="yka", t1="2005.0", t2="2007.0")], "cal")), 1)

    def test_real_table_format(self):
        # Exact header of results/h13-v1/h13_pairs_scored.csv, lower-case array, CRLF line ends,
        # synthetic values. Row 2 is M03-24. Row 3 is a calibration pair (2005-2007, dt 2.0).
        text = ("array,mult,ev1,ev2,t1,t2,dB,dM,cross,m0324,D,D_c310,sep_km,S,n_st\r\n"
                "yka,Y900,e1,e2,2021.0,2023.0,0.8,2.0,False,False,0.01,0.05,5.0,0.2,3\r\n"
                "yka,Y900,e1,e3,2021.0,2024.0,1.2,3.0,False,True,0.01,0.05,5.0,0.2,3\r\n"
                "yka,Y901,e4,e5,2005.0,2007.0,0.8,2.0,False,False,0.01,0.07,5.0,0.1,4\r\n")
        rows = list(csv.DictReader(io.StringIO(text)))
        self.assertEqual(len(g.select(rows, "post")), 1)
        self.assertEqual(len(g.select(rows, "post", include_m0324=True)), 2)
        c = g.select(rows, "cal")
        self.assertEqual(len(c), 1)
        self.assertAlmostEqual(c[0]["v"], 0.07)
        self.assertAlmostEqual(c[0]["dt"], 2.0)
        self.assertEqual(sum(1 for r in rows if g.is_yka(r)), 3)

    def test_cal_filters(self):
        rows = [
            row(t1="2005.0", t2="2007.0"),   # kept, dt 2.0
            row(t1="2004.0", t2="2007.5"),   # dt 3.5 > 2.65 -> out
            row(t1="2007.0", t2="2009.0"),   # after the turn -> out
            row(t1="1989.5", t2="1990.5"),   # before 1990 -> out
            row(t1="2006.0", t2="2008.4"),   # kept, dt 2.4
        ]
        c = g.select(rows, "cal")
        self.assertEqual(len(c), 2)
        self.assertAlmostEqual(c[0]["dt"], 2.0)
        self.assertAlmostEqual(c[1]["dt"], 2.4)


def monotone_pairs(n_fam, per_fam=3, base=0.0):
    """Family k has pairs with dt = k + 0.1j and v = dt: every resample has rho = 1."""
    out = []
    for k in range(n_fam):
        for j in range(1, per_fam + 1):
            dt = k + 0.1 * j + base
            out.append(dict(mult="F%d" % k, ev1="e%d_0" % k, ev2="e%d_%d" % (k, j),
                            t1=2021.0, t2=2021.0 + dt, dt=dt, v=dt))
    return out


class TestGrade(unittest.TestCase):
    def test_monotone_boot(self):
        s = g.summarize(monotone_pairs(4), seed=0, nboot=200)
        self.assertAlmostEqual(s["rho"], 1.0)
        self.assertEqual(s["n"], 12)
        self.assertEqual(s["n_fam"], 4)
        self.assertEqual(s["n_long"], 6)  # families 2 and 3 have dt >= 2
        self.assertAlmostEqual(s["f_pos"], 1.0)
        self.assertAlmostEqual(s["boot_p95"], 1.0)

    def _P(self, **kw):
        d = dict(n=30, n_fam=6, n_long=8, rho=0.05, f_pos=0.6, boot_p95=0.2)
        d.update(kw)
        return d

    def _C(self, **kw):
        d = dict(n=30, n_fam=6, rho=0.4, f_pos=0.99)
        d.update(kw)
        return d

    def test_headlines(self):
        self.assertEqual(g.headline(self._P(n=10), self._C()), "untestable")
        self.assertEqual(g.headline(self._P(n_fam=3), self._C()), "untestable")
        self.assertEqual(g.headline(self._P(n_long=4), self._C()), "untestable")
        self.assertEqual(g.headline(self._P(rho=0.3, f_pos=0.97), self._C(f_pos=0.1)), "still-changing")
        self.assertEqual(g.headline(self._P(), self._C()), "quiet-since-2020")
        self.assertEqual(g.headline(self._P(), self._C(rho=0.15)), "unresolved")
        self.assertEqual(g.headline(self._P(), self._C(f_pos=0.5)), "insensitive")
        self.assertEqual(g.headline(self._P(), self._C(n=19)), "insensitive")
        self.assertEqual(g.headline(self._P(rho=0.3, f_pos=0.94, boot_p95=0.5), self._C()), "unresolved")

    def test_grade_all_end_to_end(self):
        for arr in ("YKA", "yka"):
            rows = []
            for p in monotone_pairs(5, per_fam=5):  # 25 post pairs, dt 0.1 .. 4.5, rho 1
                # t1 = 2022.5 is after era 2 ends (2022.2466), so no pair touches it
                rows.append(row(array=arr, mult=p["mult"], ev1=p["ev1"], ev2=p["ev2"], t1=2022.5,
                                t2=2022.5 + p["dt"], D_c310=p["v"]))
            r = g.grade_all(rows, nboot=100)
            self.assertEqual(r["post"]["n"], 25)
            self.assertTrue(r["testable_post"])
            self.assertEqual(r["headline"], "still-changing")
            self.assertEqual(r["cal"]["n"], 0)
            self.assertFalse(r["calibrated"])

    def test_perm_p_range(self):
        p = g.perm_p(monotone_pairs(4), nperm=50, seed=0)
        self.assertTrue(0 < p <= 1)


if __name__ == "__main__":
    unittest.main()

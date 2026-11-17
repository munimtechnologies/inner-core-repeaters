"""Synthetic tests for fdscan.py (numpy; grading also needs scipy).

Part of inner-core-repeaters by MV2, Munim, Inc. (made with AI agents). Apache-2.0.
"""
import math
import unittest

try:
    import numpy as np
    import fdscan
    import model
except ImportError:  # pragma: no cover
    np = None

try:
    import scipy  # noqa: F401
    HAVE_SCIPY = True
except ImportError:  # pragma: no cover
    HAVE_SCIPY = False


def ev(i, t, lat=-59.0, lon=-25.7, mag=4.8):
    return dict(id=i, time="x", lat=lat, lon=lon, depth=30.0, mag=mag, magType="mb", t=t)


@unittest.skipIf(np is None, "numpy missing")
class TestSelect(unittest.TestCase):
    def test_parse_rows(self):
        rows = [dict(id="a", time="2025-03-01T00:00:00Z", latitude="-59", longitude="-25", depth="", mag="4.6",
                     magType="mb", type="earthquake"),
                dict(id="a", time="2025-03-01T00:00:00Z", latitude="-59", longitude="-25", depth="10", mag="4.6",
                     magType="mb", type="earthquake"),
                dict(id="q", time="2001-03-01T00:00:00Z", latitude="-59", longitude="-25", depth="10", mag="4.6",
                     magType="mb", type="quarry blast")]
        evs = fdscan.parse_rows(rows)
        self.assertEqual([e["id"] for e in evs], ["a"])
        self.assertEqual(evs[0]["depth"], 10.0)

    def test_cases(self):
        evs = [ev("old_b", 2001.6), ev("old_m", 2006.0), ev("old_early", 1999.2),
               ev("far_b", 2001.6, lat=-56.0), ev("far_m", 2006.0, lat=-56.0),
               ev("new1", 2025.4, mag=5.0), ev("new_far", 2025.5, lat=-56.0, mag=5.2),
               ev("new_lonely", 2025.6, lat=-60.5), ev("mid", 2015.0)]
        evs.sort(key=lambda e: e["t"])
        cases = fdscan.select_cases(evs, r_km=25.0, b_half=1.0)
        self.assertEqual([c["new"]["id"] for c in cases], ["new_far", "new1"])  # by magnitude
        c = cases[1]
        self.assertAlmostEqual(c["e2"], model.equiv_date(2025.4))
        self.assertEqual([o["id"] for o in c["b"]], ["old_b"])
        self.assertEqual([o["id"] for o in c["m"]], ["old_m"])  # 1999.2 is not "later", 2015 is post-turn
        order = [e["id"] for e in fdscan.event_order(cases)]
        self.assertEqual(order, ["new_far", "far_b", "far_m", "new1", "old_b", "old_m"])
        self.assertEqual(len(fdscan.event_order(cases, max_events=4)), 4)

    def test_binom(self):
        self.assertAlmostEqual(fdscan.binom_sf(0, 5), 1.0)
        self.assertAlmostEqual(fdscan.binom_sf(5, 5), 1 / 32.0)
        self.assertAlmostEqual(fdscan.binom_sf(9, 10), 11 / 1024.0)
        self.assertEqual(fdscan.binom_sf(0, 0), 1.0)


def synth_rows(prefer, n_mult=5, seed=1):
    """n_mult multiplets, each with partners at 2001.6 (B pick) and 2006.5 (M pick) and two new
    events. prefer = 'B' or 'M' sets which partner is the better DF match."""
    rng = np.random.default_rng(seed)
    rows = []
    for k in range(n_mult):
        mid = "F%03d" % k
        ts = {"b%d" % k: 2001.6, "m%d" % k: 2006.5, "n%da" % k: 2025.3, "n%db" % k: 2026.1}
        names = sorted(ts, key=ts.get)
        for i, p in enumerate(names):
            for q in names[i + 1:]:
                if prefer == "B":
                    D = 0.03 + 0.05 * model.delta_B(ts[p], ts[q])
                else:
                    D = 0.03 + 0.05 * model.delta_M(ts[p], ts[q])
                D += rng.normal(0, 0.003)
                rows.append(dict(mult=mid, ev1=p, ev2=q, t1=ts[p], t2=ts[q], D=D, Dx=D, ok=True,
                                 snr_df1=5.0 + rng.random()))
    return rows


@unittest.skipUnless(np is not None and HAVE_SCIPY, "numpy/scipy missing")
class TestGrade(unittest.TestCase):
    def evs(self, rows):
        out = {}
        for r in rows:
            for i, t in ((r["ev1"], r["t1"]), (r["ev2"], r["t2"])):
                out[i] = ev(i, t, lat=-57.0)  # far from M03-24
        return out

    def test_backtracking(self):
        rows = synth_rows("B")
        res = fdscan.grade_all(rows, self.evs(rows), nboot=300)
        F1 = res["F1"]
        self.assertEqual(F1["n_discriminating"], 10)
        self.assertEqual(F1["n_match_B"], 10)
        self.assertEqual(F1["verdict"], "backtracking-supported", F1)
        self.assertTrue(res["F1_robust"])
        self.assertEqual(res["F1e_without_M0324"]["n_discriminating"], 10)
        self.assertAlmostEqual(F1["chance_expected_B"], 10 * 0.5)

    def test_monotone(self):
        rows = synth_rows("M")
        res = fdscan.grade_all(rows, self.evs(rows), nboot=300)
        self.assertEqual(res["F1"]["verdict"], "backtracking-contradicted", res["F1"])
        self.assertEqual(res["F1"]["n_match_M"], 10)

    def test_untestable_and_exclusion(self):
        rows = synth_rows("B", n_mult=2)
        evs = self.evs(rows)
        for i in evs:
            evs[i]["lat"], evs[i]["lon"] = -59.03, -25.70  # all inside the M03-24 patch
        res = fdscan.grade_all(rows, evs, nboot=100)
        self.assertEqual(res["F1"]["verdict"], "untestable")  # 4 cases, 2 multiplets
        self.assertEqual(res["F1e_without_M0324"]["n_discriminating"], 0)

    def test_not_ok_rows_ignored(self):
        rows = synth_rows("B")
        for r in rows:
            if r["t1"] < 2002:
                r["ok"] = False  # remove every B partner -> one partner per new event -> no cases
        res = fdscan.grade_all(rows, self.evs(rows), nboot=100)
        self.assertEqual(res["F1"]["n_cases"], 0)
        self.assertEqual(res["F1"]["verdict"], "untestable")
        self.assertTrue(math.isnan(res["F1"]["frac_best_is_top_snr"]))


if __name__ == "__main__":
    unittest.main()

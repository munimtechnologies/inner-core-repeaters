"""Synthetic tests for m0324.py (numpy; grading tests also need scipy).

Part of inner-core-repeaters by MV2, Munim, Inc. (made with AI agents). Apache-2.0.
"""
import unittest

try:
    import numpy as np
    import m0324
    import model
except ImportError:  # pragma: no cover
    np = None

try:
    import scipy  # noqa: F401
    HAVE_SCIPY = True
except ImportError:  # pragma: no cover
    HAVE_SCIPY = False


@unittest.skipIf(np is None, "numpy missing")
class TestGeometry(unittest.TestCase):
    def test_centroid_inside_members(self):
        la, lo = m0324.centroid()
        self.assertAlmostEqual(la, -59.036, places=2)
        self.assertAlmostEqual(lo, -25.694, places=2)
        for _, _, mla, mlo, _, _ in m0324.M0324:
            self.assertLess(m0324.km_between(la, lo, mla, mlo), 20.0)

    def test_km_between(self):
        self.assertAlmostEqual(m0324.km_between(0, 0, 1, 0), 111.195, places=2)
        self.assertAlmostEqual(m0324.km_between(-59, -25, -59, -25), 0.0, places=6)

    def test_forecast_dates(self):
        # the backtracking forecast this patch tests: 2025-2026 map to 2001.5-2001.9
        self.assertAlmostEqual(model.equiv_date(2025.0), 2001.9, places=6)
        self.assertAlmostEqual(model.equiv_date(2026.0), 2001.5, places=6)


@unittest.skipIf(np is None, "numpy missing")
class TestCatalog(unittest.TestCase):
    def test_merge_and_patch(self):
        cat = [dict(id="usp000bt5x", time="2003-03-21T10:47:47.620Z", latitude="-59.07", longitude="-25.47",
                    depth="33", mag="4.9", magType="mb", type="earthquake"),
               dict(id="near", time="2025-10-13T02:40:33.881Z", latitude="-59.079", longitude="-25.4283",
                    depth="35", mag="4.9", magType="mb", type="earthquake"),
               dict(id="far", time="2025-03-14T23:42:32Z", latitude="-55.8", longitude="-27.07",
                    depth="13", mag="5.9", magType="mww", type="earthquake"),
               dict(id="blast", time="2025-03-15T00:00:00Z", latitude="-59.0", longitude="-25.6",
                    depth="0", mag="4.5", magType="mb", type="explosion")]
        evs = m0324.merge_catalog(cat, m0324.member_rows())
        ids = [e["id"] for e in evs]
        self.assertEqual(len(ids), len(set(ids)))
        self.assertNotIn("blast", ids)
        self.assertEqual(len(evs), 6 + 2)  # 6 members (one from the catalogue) + near + far
        self.assertEqual(next(e for e in evs if e["id"] == "usp000bt5x")["mag"], 4.9)  # catalogue row wins
        self.assertTrue(all(evs[i]["t"] <= evs[i + 1]["t"] for i in range(len(evs) - 1)))
        la, lo = m0324.centroid()
        keep = {e["id"] for e in m0324.patch_events(evs, la, lo, 35.0)}
        self.assertIn("near", keep)
        self.assertNotIn("far", keep)
        self.assertTrue({m[0] for m in m0324.M0324} <= keep)


@unittest.skipIf(np is None, "numpy missing")
class TestZ0(unittest.TestCase):
    def test_reproduced(self):
        mem = [m[0] for m in m0324.M0324]
        comps = [mem[:5] + ["x1", "x2"], ["y1", "y2"]]
        r = m0324.z0_membership(comps, mem, mem)
        self.assertTrue(r["testable"])
        self.assertTrue(r["holds"])
        self.assertEqual(r["n_in_largest_shared_component"], 5)

    def test_split(self):
        mem = [m[0] for m in m0324.M0324]
        comps = [mem[:3], mem[3:]]
        r = m0324.z0_membership(comps, mem, mem)
        self.assertFalse(r["holds"])

    def test_untestable(self):
        mem = [m[0] for m in m0324.M0324]
        r = m0324.z0_membership([mem[:3]], mem, mem[:3])
        self.assertFalse(r["testable"])
        self.assertIsNone(r["holds"])


@unittest.skipUnless(np is not None and HAVE_SCIPY, "numpy/scipy missing")
class TestGrade(unittest.TestCase):
    def test_backtracking_patch(self):
        # a patch whose new events best match their ~2001-2002 partners (model B)
        rng = np.random.default_rng(3)
        # equivalent dates of n1..n3: 2001.82, 2001.58, 2001.38 -> B picks "b", M picks "d"
        ts = {"a": 2000.0, "b": 2001.7, "c": 2003.2, "d": 2008.1, "e": 2011.1, "f": 2020.5,
              "n1": 2025.2, "n2": 2025.8, "n3": 2026.3}
        rows = []
        names = sorted(ts, key=ts.get)
        for i, p in enumerate(names):
            for q in names[i + 1:]:
                D = 0.03 + 0.05 * model.delta_B(ts[p], ts[q]) + rng.normal(0, 0.005)
                rows.append(dict(mult="Z000", ev1=p, ev2=q, t1=ts[p], t2=ts[q], D=D, Dx=D, Ds=D, ok=True))
        res = m0324.grade(rows, [names], names, nboot=100)
        z1 = res["Z1_new_event_best_match"]
        self.assertTrue(z1["testable"], z1)
        self.assertTrue(z1["holds"], z1)
        self.assertEqual(z1["n_match_M"], 0)
        self.assertIsNone(res["Z0_reproduce_M03_24"]["holds"])  # no published ids in this synthetic set
        self.assertTrue(res["Z3_calibration"]["holds"] in (True, None))


if __name__ == "__main__":
    unittest.main()

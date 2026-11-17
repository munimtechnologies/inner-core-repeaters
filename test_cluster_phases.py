"""Tests for cluster.py and the interpolation in phases.py (no obspy, no network). Apache-2.0."""
import unittest

from cluster import components, max_spanning_offsets

try:
    import numpy as np

    from phases import Phases, gc_deg

    HAVE = True
except ImportError:
    HAVE = False


class TestCluster(unittest.TestCase):
    def test_components(self):
        comps = components([0, 1, 2, 3, 4], {(0, 1): (0.95, 1.0), (1, 2): (0.92, 0.0), (3, 4): (0.9, 0.0)})
        self.assertEqual(sorted(comps), [[0, 1, 2], [3, 4]])

    def test_offsets_follow_strongest_edges(self):
        edges = {(0, 1): (0.95, 5.0), (1, 2): (0.93, -2.0), (0, 2): (0.50, 99.0)}
        off = max_spanning_offsets([0, 1, 2], edges, root=0)
        self.assertAlmostEqual(off[0], 0.0)
        self.assertAlmostEqual(off[1], 5.0)
        self.assertAlmostEqual(off[2], 3.0)
        off1 = max_spanning_offsets([0, 1, 2], edges, root=1)
        self.assertAlmostEqual(off1[0], -5.0)
        self.assertAlmostEqual(off1[2], -2.0)


@unittest.skipUnless(HAVE, "numpy not installed")
class TestPhases(unittest.TestCase):
    def grid(self):
        dists = np.array([140.0, 150.0, 160.0])
        depths = np.array([0.0, 100.0])
        df = np.array([[1000.0, 1100.0, 1200.0], [990.0, 1090.0, 1190.0]])
        bc = np.array([[np.nan, 1103.0, np.nan], [np.nan, 1093.0, np.nan]])
        return dict(dists=dists, depths=depths, df=df, bc=bc, ab=df + 6.0, p_df=df * 0 + 1.9)

    def test_bilinear(self):
        ph = Phases(self.grid())
        self.assertAlmostEqual(ph.t("df", 145.0, 0.0), 1050.0)
        self.assertAlmostEqual(ph.t("df", 150.0, 50.0), 1095.0)
        self.assertAlmostEqual(ph.t("df", 140.0, -5.0), 1000.0)  # depth clipped
        self.assertTrue(np.isnan(ph.t("df", 139.0, 10.0)))
        self.assertTrue(np.isnan(ph.t("bc", 155.0, 10.0)))  # branch missing -> NaN
        self.assertAlmostEqual(ph.t("bc", 150.0, 0.0), 1103.0)

    def test_gc(self):
        self.assertAlmostEqual(float(gc_deg(0, 0, 0, 90)), 90.0, places=6)
        self.assertAlmostEqual(float(gc_deg(-90, 0, 90, 0)), 180.0, places=6)
        d = float(gc_deg(-58.0, -26.5, 64.7714, -146.8866))  # SSI to ILAR
        self.assertTrue(148.0 < d < 155.0, d)


if __name__ == "__main__":
    unittest.main()

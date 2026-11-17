"""ak135 travel times of PKIKP (DF), PKPbc and PKPab on a distance x depth grid.

Part of inner-core-repeaters by MV2, Munim, Inc. (made with AI agents). Apache-2.0.

The two PKP arrivals returned by TauP under the name "PKP" are told apart by time:
beyond the B caustic (~145 deg) the BC branch arrives before the AB branch; beyond the
C cusp (~153-155 deg) only AB remains.
"""
import math
import os

import numpy as np

DISTS = np.round(np.arange(140.0, 160.0 + 1e-9, 0.1), 2)
DEPTHS = np.array([0.0, 10.0, 20.0, 35.0, 60.0, 100.0])
KEYS = ("df", "bc", "ab", "p_df")


def build_grid(model="ak135", dists=DISTS, depths=DEPTHS):
    from obspy.taup import TauPyModel

    m = TauPyModel(model=model)
    g = {k: np.full((len(depths), len(dists)), np.nan) for k in KEYS}
    for a, z in enumerate(depths):
        for b, d in enumerate(dists):
            try:
                arr = m.get_travel_times(source_depth_in_km=float(z), distance_in_degree=float(d),
                                         phase_list=["PKIKP", "PKP"])
            except Exception:
                continue
            df = [x for x in arr if x.name == "PKIKP"]
            pkp = sorted(x.time for x in arr if x.name == "PKP")
            if df:
                g["df"][a, b] = df[0].time
                g["p_df"][a, b] = df[0].ray_param_sec_degree
            if len(pkp) >= 2:
                g["bc"][a, b] = pkp[0]
                g["ab"][a, b] = pkp[-1]
            elif len(pkp) == 1:
                g["ab"][a, b] = pkp[0]
    return dict(dists=np.asarray(dists, dtype=float), depths=np.asarray(depths, dtype=float), **g)


def load_or_build(path):
    if path and os.path.exists(path):
        z = np.load(path)
        return {k: z[k] for k in z.files}
    g = build_grid()
    if path:
        np.savez(path, **g)
    return g


class Phases:
    """Bilinear interpolation of the grid; NaN outside it or where a branch does not exist."""

    def __init__(self, grid):
        self.g = grid
        self.d = np.asarray(grid["dists"], dtype=float)
        self.z = np.asarray(grid["depths"], dtype=float)

    def t(self, key, dist, depth):
        d, z, A = self.d, self.z, np.asarray(self.g[key], dtype=float)
        dist = float(dist)
        depth = float(depth) if depth is not None and math.isfinite(float(depth)) else 10.0
        if not (d[0] <= dist <= d[-1]):
            return float("nan")
        depth = min(max(depth, z[0]), z[-1])
        ib = int(np.clip(np.searchsorted(d, dist) - 1, 0, len(d) - 2))
        fb = (dist - d[ib]) / (d[ib + 1] - d[ib])
        ia = int(np.clip(np.searchsorted(z, depth) - 1, 0, len(z) - 2))
        fa = (depth - z[ia]) / (z[ia + 1] - z[ia])
        terms = (((1 - fa) * (1 - fb), A[ia, ib]), ((1 - fa) * fb, A[ia, ib + 1]),
                 (fa * (1 - fb), A[ia + 1, ib]), (fa * fb, A[ia + 1, ib + 1]))
        v = 0.0
        for w, a in terms:
            if w > 1e-12:
                if not math.isfinite(a):
                    return float("nan")
                v += w * a
        return float(v)


def gc_deg(lat1, lon1, lat2, lon2):
    """Great-circle distance in degrees (haversine); works on numpy arrays."""
    p1, p2 = np.radians(lat1), np.radians(lat2)
    dp = p2 - p1
    dl = np.radians(np.asarray(lon2) - np.asarray(lon1))
    h = np.sin(dp / 2) ** 2 + np.cos(p1) * np.cos(p2) * np.sin(dl / 2) ** 2
    return np.degrees(2 * np.arcsin(np.sqrt(np.clip(h, 0, 1))))

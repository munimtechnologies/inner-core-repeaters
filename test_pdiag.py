"""Tests for pdiag.py (stdlib + numpy only). Part of inner-core-repeaters by MV2, Munim, Inc.
(made with AI agents). Licence: Apache-2.0."""
import math
import random
import unittest

import numpy as np

import pdiag as P


def fP(strike, dip, rake, takeoff, az):
    """Aki & Richards (4.89) P radiation pattern, written independently of sdr_to_ned."""
    s, d, r = math.radians(strike), math.radians(dip), math.radians(rake)
    i, f = math.radians(takeoff), math.radians(az)
    return (math.cos(r) * math.sin(d) * math.sin(i) ** 2 * math.sin(2 * (f - s))
            - math.cos(r) * math.cos(d) * math.sin(2 * i) * math.cos(f - s)
            + math.sin(r) * math.sin(2 * d) * (math.cos(i) ** 2 - math.sin(i) ** 2 * math.sin(f - s) ** 2)
            + math.sin(r) * math.cos(2 * d) * math.sin(2 * i) * math.sin(f - s))


class TestRadiation(unittest.TestCase):
    def test_thrust_vertical_ray(self):
        M = P.sdr_to_ned(0, 45, 90)
        self.assertAlmostEqual(P.p_amp(M, 0.0, 123.0), 1.0, places=9)

    def test_strike_slip_lobes_and_nodes(self):
        M = P.sdr_to_ned(0, 90, 0)
        self.assertAlmostEqual(P.p_amp(M, 90.0, 45.0), 1.0, places=9)
        self.assertAlmostEqual(P.p_amp(M, 90.0, 0.0), 0.0, places=9)
        self.assertAlmostEqual(P.p_amp(M, 90.0, 90.0), 0.0, places=9)

    def test_matches_aki_richards_pattern(self):
        rng = random.Random(4)
        for _ in range(300):
            s, d, r = rng.uniform(0, 360), rng.uniform(1, 89), rng.uniform(-180, 180)
            i, a = rng.uniform(0, 90), rng.uniform(0, 360)
            A = P.p_amp(P.sdr_to_ned(s, d, r, m0=3.7e17), i, a)
            self.assertAlmostEqual(A, fP(s, d, r, i, a), places=9)
            self.assertLessEqual(abs(A), 1.0 + 1e-9)

    def test_trace_free_and_m0(self):
        M = P.sdr_to_ned(37, 61, -112, m0=2.5)
        self.assertAlmostEqual(float(np.trace(M)), 0.0, places=12)
        self.assertAlmostEqual(P.scalar_moment(M), 2.5, places=12)

    def test_rtp_roundtrip(self):
        M = P.sdr_to_ned(200, 30, 75, m0=1e24)
        c = P.ned_to_rtp(M)
        M2 = P.mt_ned(c["mrr"], c["mtt"], c["mpp"], c["mrt"], c["mrp"], c["mtp"])
        self.assertTrue(np.allclose(M, M2, rtol=0, atol=1e12))

    def test_mw(self):
        self.assertAlmostEqual(P.mw_from_dyncm(10 ** (1.5 * 5.0 + 16.1)), 5.0, places=9)


NDK = """PDE  2014/12/09 19:23:36.4 -58.99  -25.60  50.8 5.1 0.0 SOUTH SANDWICH ISLANDS REGION
C201412091923A   B: 30   43  40 S: 31   62  50 M:  0    0   0 CMT: 1 TRIHD:  0.8
CENTROID:      1.1 0.2 -59.10 0.02  -25.37 0.02  33.6  1.3 FREE O-20150302120000
24  0.838 0.201-0.005 0.231-0.833 0.270  1.050 0.121-0.369 0.161  0.044 0.240
V10   1.581 56  12  -0.537 23 140  -1.044 24 241   1.312   9 29  142 133 72   66
PDE  2005/01/01 01:20:05.4  13.78  -88.78 193.1 5.0 0.0 EL SALVADOR
C200501010120A   B:  4    4  40 S: 27   33  50 M:  0    0   0 CMT: 1 TRIHD:  0.6
CENTROID:     -0.3 0.9  13.76 0.06  -89.08 0.09 162.8 12.5 FREE S-20050322125201
23  0.838 0.201 -0.005 0.231 -0.833 0.270  1.050 0.121 -0.369 0.161  0.044 0.240
V10   1.581 56  12  -0.537 23 140  -1.044 24 241   1.312   9 29  142 133 72   66
"""


class TestNdk(unittest.TestCase):
    def test_parse_box_and_values(self):
        got = P.parse_ndk(NDK)
        self.assertEqual(len(got), 1)
        c = got[0]
        self.assertEqual(c["name"], "C201412091923A")
        self.assertAlmostEqual(c["lat"], -58.99)
        for k, v in (("mrr", 0.838), ("mtt", -0.005), ("mpp", -0.833), ("mrt", 1.050), ("mrp", -0.369),
                     ("mtp", 0.044)):
            self.assertAlmostEqual(c[k] / 1e24, v, places=9)
        self.assertAlmostEqual(c["epoch"], P.iso_epoch("2014-12-09T19:23:36.4"), places=6)

    def test_match_prefers_final_and_tolerances(self):
        t = P.iso_epoch("2014-12-09T19:23:36.400")
        cm = [dict(name="quick", epoch=t + 0.5, lat=-59.0, lon=-25.6, quick=True),
              dict(name="final", epoch=t + 3.0, lat=-59.0, lon=-25.6, quick=False),
              dict(name="far", epoch=t + 1.0, lat=-55.0, lon=-25.6, quick=False)]
        self.assertEqual(P.match_gcmt(t, -58.99, -25.60, cm)["name"], "final")
        self.assertEqual(P.match_gcmt(t, -58.99, -25.60, cm[:1])["name"], "quick")
        self.assertIsNone(P.match_gcmt(t + 100, -58.99, -25.60, cm))


class TestMeasure(unittest.TestCase):
    def grid(self, amp, t_on, seed=1):
        rng = np.random.RandomState(seed)
        n = int(round((P.GRID_PRE + P.GRID_POST) * P.FS))
        x = rng.normal(0, 1, n)
        c = int(round(P.GRID_PRE * P.FS))
        t = (np.arange(n) - c) / P.FS
        m = (t >= t_on) & (t < t_on + 4)
        x[m] += amp * np.sin(2 * np.pi * 1.5 * (t[m] - t_on)) * np.exp(-(t[m] - t_on) / 1.5)
        return x

    def test_clear_arrival(self):
        m = P.measure(self.grid(40.0, 2.0))
        self.assertGreater(m["snr"], 15)
        self.assertGreaterEqual(m["t_pick"], 1.5)
        self.assertLessEqual(m["t_pick"], 2.6)
        self.assertGreaterEqual(m["t_envmax"], 1.8)
        self.assertLessEqual(m["t_envmax"], 4.0)
        self.assertGreater(m["raw_snr_m2win"], 15)

    def test_noise_only(self):
        m = P.measure(self.grid(0.0, 2.0, seed=7))
        self.assertLess(m["snr"], 6)
        self.assertTrue(math.isnan(m["t_pick"]))

    def test_off_grid(self):
        self.assertIsNone(P.measure(np.ones(100)))


class TestStats(unittest.TestCase):
    def test_rank_ties(self):
        self.assertEqual(list(P.rankdata([3, 1, 3, 2])), [2.5, 0.0, 2.5, 1.0])

    def test_spearman(self):
        self.assertAlmostEqual(P.spearman([1, 2, 3, 4, 5], [2, 4, 9, 16, 30]), 1.0)
        self.assertAlmostEqual(P.spearman([1, 2, 3, 4, 5], [5, 4, 3, 2, 1]), -1.0)
        self.assertTrue(math.isnan(P.spearman([1, 2, 3], [1, 2, 3])))

    def test_median(self):
        self.assertEqual(P.median([3, "", 1, float("nan"), 2]), 2.0)
        self.assertTrue(math.isnan(P.median([])))


def _row(st, ev, snr, A, band="1-2", resp=True, first=True, loc="00", t_pick=0.3, pP=8.0):
    # pP may be "" (no ak135 depth phase), exactly as pdiag writes it; sP follows it.
    sP = (pP + 1.5) if isinstance(pP, (int, float)) else ""
    return dict(station=st, event=ev, dist_deg=50.0, az_deg=90.0, A_pred=A, loc=loc, first_loc=first,
                band=band, resp_requested=resp, resp_removed=resp, status="ok", snr=snr, t_pick=t_pick,
                noise_rms=1e-9, peak=snr * 1e-9, pP_minus_P_s=pP, sP_minus_P_s=sP)


class TestAnalyse(unittest.TestCase):
    def test_radiation_and_rescues(self):
        rows = []
        for k in range(6):
            A = 0.1 + 0.15 * k
            for e in range(4):
                rows.append(_row("S%d" % k, "e%d" % e, 2.0 + 20 * A + e * 0.1, A))
        # station S0 fails at 1-2 Hz but works at 0.5-1 Hz and at location 10
        for e in range(4):
            rows.append(_row("S0", "e%d" % e, 9.0, 0.1, band="0.5-1"))
            rows.append(_row("S0", "e%d" % e, 8.0, 0.1, loc="10", first=False))
            rows.append(_row("S0", "e%d" % e, 3.9, 0.1, resp=False))
        res, tab = P.analyse(rows)
        R = res["R_radiation"]
        self.assertGreater(R["rho_rows_absA_logsnr"], 0.8)
        self.assertAlmostEqual(R["rho_stations"], 1.0)
        self.assertLess(R["median_absA_bad"], R["median_absA_good"])
        self.assertIn("S0", res["failing_stations"])
        self.assertIn("S0", res["F_band_rescues"])
        self.assertIn("S0", res["L_location_rescues"])
        self.assertNotIn("S0", res["X_response_changes"])
        self.assertEqual(len(tab), 6)

    def test_depth_phase_lock(self):
        rows = [_row("A", "e%d" % e, 30.0, "", t_pick=4.0 + e, pP=4.0 + e) for e in range(6)]
        rows += [_row("B", "e%d" % e, 30.0, "", t_pick=0.2, pP=4.0 + e) for e in range(6)]
        rows.append(_row("C", "e0", 30.0, "", t_pick=5.0, pP=""))
        res, tab = P.analyse(rows)
        D = res["D_depth_phase"]
        self.assertEqual(D["n_clear"], 13)
        self.assertEqual(D["n_late"], 7)
        self.assertAlmostEqual(D["late_median_abs_minus_pP"], 0.0)
        a = [d for d in tab if d["station"] == "A"][0]
        self.assertAlmostEqual(a["rho_tpick_pP"], 1.0)
        self.assertEqual(res["R_radiation"]["n_rows"], 0)


if __name__ == "__main__":
    unittest.main()

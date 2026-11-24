# Copyright 2026 MV2, Munim, Inc.
# Licensed under the Apache License, Version 2.0 (Apache-2.0). See LICENSE.
# Made with AI agents working in MV2's name.
"""Synthetic tests for ks32f.py (pre-registration 32, amendment 1). numpy + stdlib only."""

import math
import unittest

import numpy as np

import ks32f as F

AZS = {"COLA": 15.0, "ULN": 310.0, "INK": 20.0, "YAK": 330.0}
# (dT/dDelta s/deg, dT/dh s/km), roughly iasp91-like values near 100-106 deg
PART = {"SKS": (6.0, -0.2155), "P1": (4.44, -0.1174), "PP": (7.6, -0.103),
        "S1": (8.3, -0.209), "SS": (13.7, -0.1836)}
PARTS = {s: dict(PART) for s in AZS}
TRUTH = dict(X=-0.8, uN=4.0, uE=-2.0, dh=3.0)


def dT(phase, az, t=TRUTH):
    p, h = PART[phase]
    a = math.radians(az)
    return (-p * (180.0 / math.pi) / F.R_EARTH * (t["uN"] * math.cos(a) + t["uE"] * math.sin(a))
            + h * t["dh"])


def synth_rows(t=TRUTH, edge_yak_ss=False, band="0.01-0.05"):
    rows = []
    for sta, az in AZS.items():
        o = 0.37
        for ph in ["SKS"] + F.REFS:
            tot = dT(ph, az, t) + o + (t["X"] if ph == "SKS" else 0.0)
            edge = "True" if (edge_yak_ss and sta == "YAK" and ph == "SS") else "False"
            rows.append(dict(A="E1997", B="E2018", net="XX", sta=sta, band=band, phase=ph,
                             d=repr(tot - 0.25), dp="0.25", cc="0.97", edge=edge,
                             snr_a="12", snr_b="12"))
    return rows


class FakeArr:
    def __init__(self, name, time):
        self.name = name
        self.time = time


class FakeModel:
    def get_travel_times(self, source_depth_in_km, distance_in_degree, phase_list):
        d, h = distance_in_degree, source_depth_in_km
        out = []
        if "SKS" in phase_list:
            out.append(FakeArr("SKS", 1400.0 + 6.0 * d - 0.2 * h))
        if "P" in phase_list:
            out.append(FakeArr("P", 800.0 + 4.5 * d - 0.12 * h))
        if "Pdiff" in phase_list:
            out.append(FakeArr("Pdiff", 820.0 + 4.4 * d - 0.12 * h))  # later than P here
        return out


class TestGeometry(unittest.TestCase):
    def test_dist_az(self):
        d, az = F.dist_az(0.0, 0.0, 0.0, 90.0)
        self.assertAlmostEqual(d, 90.0, places=6)
        self.assertAlmostEqual(az, 90.0, places=6)
        d, az = F.dist_az(0.0, 0.0, 10.0, 0.0)
        self.assertAlmostEqual(d, 10.0, places=6)
        self.assertAlmostEqual(az, 0.0, places=6)
        d, az = F.dist_az(10.0, 0.0, 0.0, 0.0)
        self.assertAlmostEqual(az, 180.0, places=6)

    def test_usgs_vector(self):
        uN, uE, dh = F.usgs_vector()
        self.assertTrue(abs(uN - 5.148) < 0.002)
        self.assertTrue(abs(uE - 0.8697) < 0.002)
        self.assertAlmostEqual(dh, 6.7, places=9)


class TestRows(unittest.TestCase):
    def rows(self):
        base = dict(A="E1997", B="E2018", sta="COLA", band="0.01-0.05")
        r = [
            dict(base, phase="SKS", d="-0.9", dp="-0.1", cc="0.95", edge="False", snr_a="10", snr_b="10"),
            dict(base, phase="P1", d="0.1", dp="-0.1", cc="0.95", edge="False", snr_a="10", snr_b="10"),
            dict(base, phase="PP", d="0.5", dp="0.0", cc="0.95", edge="True", snr_a="10", snr_b="10"),
            dict(base, phase="S1", d="0.2", dp="0", cc="0.85", edge="False", snr_a="10", snr_b="10"),
            dict(base, phase="SS", d="-0.3", dp="0.1", cc="0.95", edge="False", snr_a="10", snr_b="nan"),
            dict(base, phase="SKS", band="0.10-0.50", d="5", dp="0", cc="0.99", edge="False", snr_a="10", snr_b="10"),
            dict(base, sta="ULN", phase="SKS", d="-1", dp="0", cc="0.5", edge="False", snr_a="10", snr_b="10"),
            dict(base, sta="ULN", phase="P1", d="0", dp="0", cc="0.99", edge="False", snr_a="10", snr_b="10"),
        ]
        return r

    def test_parse_and_gate(self):
        p = F.parse_rows(self.rows(), "E1997", "E2018", "0.01-0.05")
        self.assertAlmostEqual(p["COLA"]["SKS"]["tot"], -1.0, places=9)
        self.assertTrue(p["COLA"]["PP"]["edge"])
        self.assertFalse(p["COLA"]["P1"]["edge"])
        d90 = F.raw_diffs(p, 0.90, 3.0)
        self.assertEqual(len(d90), 1)
        self.assertEqual(d90[0][:2], ("COLA", "P1"))
        self.assertAlmostEqual(d90[0][2], -1.0, places=9)
        d80 = F.raw_diffs(p, 0.80, 3.0)
        self.assertEqual([x[1] for x in d80], ["P1", "S1"])
        self.assertAlmostEqual(d80[1][2], -1.2, places=9)

    def test_ok_nan(self):
        m = dict(tot=float("nan"), cc=0.99, edge=False, snr_a=10.0, snr_b=10.0)
        self.assertFalse(F.ok(m, 0.9, 3.0))
        m = dict(tot=0.1, cc=0.99, edge=False, snr_a=float("inf"), snr_b=10.0)
        self.assertTrue(F.ok(m, 0.9, 3.0))
        self.assertFalse(F.ok({}, 0.9, 3.0))


class TestFit(unittest.TestCase):
    def test_design_row(self):
        row = F.design_row(PART["SKS"], PART["P1"], 0.0)
        dp = (6.0 - 4.44) * 180.0 / math.pi / F.R_EARTH
        self.assertEqual(row[0], 1.0)
        self.assertAlmostEqual(row[1], -dp, places=12)
        self.assertLess(row[1], 0.0)
        self.assertAlmostEqual(row[2], 0.0, places=12)
        self.assertAlmostEqual(row[3], -0.0981, places=12)

    def test_analyse_recovers(self):
        out = F.analyse(synth_rows(), PARTS, AZS, "0.01-0.05", 0.90)
        self.assertEqual(out["n"], 16)
        self.assertEqual(out["n_stations"], 4)
        self.assertAlmostEqual(out["X"], -0.8, places=6)
        self.assertAlmostEqual(out["uN_km"], 4.0, places=5)
        self.assertAlmostEqual(out["uE_km"], -2.0, places=5)
        self.assertAlmostEqual(out["dh_km"], 3.0, places=5)
        self.assertTrue(np.isfinite(out["x_fixed_usgs"]))
        self.assertEqual(out["grade"], F.grade_f(16, 4, out["cond"], out["X"], out["sigma_X"]))

    def test_analyse_edge_excluded(self):
        out = F.analyse(synth_rows(edge_yak_ss=True), PARTS, AZS, "0.01-0.05", 0.90)
        self.assertEqual(out["n"], 15)
        self.assertAlmostEqual(out["X"], -0.8, places=6)

    def test_analyse_other_band_empty(self):
        out = F.analyse(synth_rows(), PARTS, AZS, "0.10-0.50", 0.80)
        self.assertEqual(out["n"], 0)
        self.assertEqual(out["grade"], "untestable")

    def test_too_few(self):
        rows = [r for r in synth_rows() if r["sta"] == "COLA" and r["phase"] in ("SKS", "P1", "SS")]
        out = F.analyse(rows, PARTS, AZS, "0.01-0.05", 0.90)
        self.assertEqual(out["n"], 2)
        self.assertEqual(out["grade"], "untestable")


class TestGrade(unittest.TestCase):
    def test_cases(self):
        g = F.grade_f
        self.assertEqual(g(16, 4, 20.0, -0.8, 0.1), "core-term")
        self.assertEqual(g(16, 4, 20.0, -0.6, 0.35), "intermediate")
        self.assertEqual(g(16, 4, 20.0, -0.4, 0.1), "intermediate")
        self.assertEqual(g(16, 4, 20.0, -0.1, 0.1), "no-core-term")
        self.assertEqual(g(16, 4, 20.0, 0.7, 0.1), "no-core-term")
        self.assertEqual(g(6, 4, 20.0, -0.8, 0.1), "untestable")
        self.assertEqual(g(16, 2, 20.0, -0.8, 0.1), "untestable")
        self.assertEqual(g(16, 4, 150.0, -0.8, 0.1), "untestable")
        self.assertEqual(g(16, 4, 20.0, -0.8, 0.6), "untestable")
        self.assertEqual(g(16, 4, 20.0, float("nan"), 0.1), "untestable")
        self.assertEqual(g(16, 4, float("inf"), -0.8, 0.1), "untestable")


class TestPartials(unittest.TestCase):
    def test_numeric(self):
        m = FakeModel()
        p = F.partials(m, 108.3, 100.0, ["SKS"], 0.05, 2.0)
        self.assertAlmostEqual(p[0], 6.0, places=6)
        self.assertAlmostEqual(p[1], -0.2, places=6)
        p = F.partials(m, 108.3, 100.0, ["P", "Pdiff"], 0.05, 2.0)
        self.assertAlmostEqual(p[0], 4.5, places=6)
        self.assertAlmostEqual(p[1], -0.12, places=6)


if __name__ == "__main__":
    unittest.main()

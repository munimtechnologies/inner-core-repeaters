# Copyright 2026 MV2, Munim, Inc.
# Licensed under the Apache License, Version 2.0 (Apache-2.0). See LICENSE.
# Made with AI agents working in MV2's name.
"""Unit tests for mkar21.py (standard library only; no network)."""

import math
import os
import tempfile
import types
import unittest

import mkar21 as m


def ranks(v):
    order = sorted(range(len(v)), key=lambda i: v[i])
    r = [0.0] * len(v)
    for k, i in enumerate(order):
        r[i] = k + 1.0
    return r


def spearman(a, b):
    if len(a) < 3:
        return float("nan")
    ra, rb = ranks(a), ranks(b)
    n = len(a)
    ma, mb = sum(ra) / n, sum(rb) / n
    sab = sum((x - ma) * (y - mb) for x, y in zip(ra, rb))
    saa = sum((x - ma) ** 2 for x in ra)
    sbb = sum((y - mb) ** 2 for y in rb)
    return sab / math.sqrt(saa * sbb)


def boot_same(rows, stat):
    """Stand-in bootstrap: ten copies of the full-sample statistic."""
    return [stat(rows)] * 10 if len({r["mult"] for r in rows}) >= 2 else []


def pct(xs, q):
    v = sorted(x for x in xs if m.isf(x))
    if not v:
        return float("nan")
    return v[int(round(q * (len(v) - 1)))]


class TestPatch(unittest.TestCase):
    def fake(self):
        calls = []

        def candidate_pairs(ev, km=10.0, dm=0.3, dz=20.0, min_mag=5.0):
            calls.append(min_mag)
            return []

        written = []

        def write_csv(path, rows):
            written.append((path, len(rows)))

        mod = types.SimpleNamespace(candidate_pairs=candidate_pairs, PAIR=dict(km=10.0, dm=0.3, dz=20.0, min_mag=5.0),
                                    PLAN="old", PLAN_TOKEN=b"pre-registration 20", write_csv=write_csv)
        return mod, calls, written

    def test_patch_changes_only_floor_and_plan(self):
        mod, calls, written = self.fake()
        cap = m.patch(mod)
        mod.candidate_pairs([])
        self.assertEqual(calls, [4.5])
        mod.candidate_pairs([], min_mag=6.0)
        self.assertEqual(calls[-1], 6.0)
        self.assertEqual(mod.PAIR, dict(km=10.0, dm=0.3, dz=20.0, min_mag=4.5))
        self.assertEqual(mod.PLAN, "findings/2026-11-20-prereg21-mkar-m45-plan.md")
        self.assertEqual(mod.PLAN_TOKEN, b"pre-registration 21")
        rows = [dict(a=1), dict(a=2)]
        mod.write_csv("/out/mkar20_pairs.csv", rows)
        self.assertEqual(written, [("/out/mkar20_pairs.csv", 2)])
        rows[0]["a"] = 99  # the capture is a copy
        self.assertEqual(cap.rows["mkar20_pairs.csv"], [dict(a=1), dict(a=2)])

    def test_token_matches_plan_title(self):
        here = os.path.dirname(os.path.abspath(__file__))
        p = os.path.join(here, "findings", "2026-11-20-prereg21-mkar-m45-plan.md")
        if not os.path.exists(p):
            p = os.path.join(here, "2026-11-20-prereg21-mkar-m45-plan.md")
        if os.path.exists(p):
            with open(p, "rb") as f:
                self.assertIn(m.PLAN_TOKEN, f.read().lower())


class TestPure(unittest.TestCase):
    def test_is_new(self):
        mags = {"a": 5.0, "b": 4.9, "c": 5.2, "d": 4.5}
        self.assertFalse(m.is_new(dict(ev1="a", ev2="c"), mags))
        self.assertTrue(m.is_new(dict(ev1="a", ev2="b"), mags))
        self.assertTrue(m.is_new(dict(ev1="d", ev2="b"), mags))
        self.assertFalse(m.is_new(dict(ev1="a", ev2="zz"), mags))  # unknown magnitude: not new

    def test_scores_from_rows(self):
        rows = [dict(ev1="x", ev2="y", S="0.25", n_st="3", D_late_median="0.4"),
                dict(ev1="p", ev2="q", S="", n_st="0")]
        sc = m.scores_from_rows(rows)
        self.assertEqual(list(sc), [("x", "y")])
        self.assertAlmostEqual(sc[("x", "y")]["S"], 0.25)
        self.assertEqual(sc[("x", "y")]["n_st"], 3)

    def test_mags_and_kinds(self):
        self.assertEqual(m.mags_of([dict(id="a", mag=4.7), dict(id="b"), dict(id="c", mag="")]), {"a": 4.7})
        self.assertEqual(m.by_kind([dict(kind="before"), dict(kind="before"), dict(kind="after")]),
                         dict(straddle=0, before=2, after=1))

    def test_long_span_pc(self):
        V = [dict(kind="before", dM=0.5 + 0.1 * i, D=0.1 + 0.02 * i, mult="F%02d" % i) for i in range(16)]
        V += [dict(kind="before", dM=0.01, D=0.9, mult="F99"),  # short span: excluded
              dict(kind="straddle", dM=3.0, D=0.0, mult="F98")]  # not before: excluded
        out = m.long_span_pc(V, spearman, boot_same, pct)
        self.assertEqual(out["n_pairs"], 16)
        self.assertEqual(out["n_families"], 16)
        self.assertGreater(out["rho"], 0.99)
        self.assertTrue(out["passed"])
        out2 = m.long_span_pc(V[:14], spearman, boot_same, pct)
        self.assertEqual(out2["n_pairs"], 14)
        self.assertFalse(out2["passed"])
        V3 = [dict(r, D=1.0 - r["D"]) for r in V[:16]]
        self.assertFalse(m.long_span_pc(V3, spearman, boot_same, pct)["passed"])

    def test_forecasts_and_replication(self):
        f = m.forecasts("B", "untestable-calibration", False, 80, True)
        self.assertEqual(f, dict(H1_primary_testable=True, H2_primary_decision="B", H3_new_pairs_testable=False,
                                 H4_long_span_control_passes=False, H5_verified_ge_80=True,
                                 H6_calibration_passes=True))
        self.assertFalse(m.forecasts("untestable", "none", True, 79, False)["H1_primary_testable"])
        self.assertEqual(m.replication("B", "B"), "replicated")
        self.assertEqual(m.replication("M", "inconclusive"), "not replicated in new pairs")
        self.assertEqual(m.replication("none", "none"), "n/a")

    def test_ensure_m20_uses_existing_file(self):
        with tempfile.TemporaryDirectory() as d:
            with open(os.path.join(d, "mkar20.py"), "wb") as f:
                f.write(b'VERSION = "0.1.0"\n')
            info = m.ensure_m20(d)
            self.assertEqual(info["commit"], m.M20[1])
            self.assertEqual(len(info["sha256"]), 64)


if __name__ == "__main__":
    unittest.main()

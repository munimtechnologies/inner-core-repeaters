# Copyright 2026 MV2, Munim, Inc.
# Licensed under the Apache License, Version 2.0 (Apache-2.0). See LICENSE.
# Made with AI agents working in MV2's name.
"""Tests for p22parts.py (standard library only)."""

import hashlib
import json
import unittest

import p22parts as pp


def split(body, k):
    lines = body.splitlines(True)
    h, rows = lines[0], lines[1:]
    n = (len(rows) + k - 1) // k
    parts = [h + b"".join(rows[i * n:(i + 1) * n]) for i in range(k)]
    man = dict(sha256=hashlib.sha256(body).hexdigest(), bytes=len(body),
               parts=[dict(name="p%d.csv" % (i + 1), sha256=hashlib.sha256(p).hexdigest())
                      for i, p in enumerate(parts)])
    return man, {"p%d.csv" % (i + 1): p for i, p in enumerate(parts)}


BODY = b"ev1,ev2,t1\n" + b"".join(b"a%d,b%d,%d.5\n" % (i, i, 2000 + i) for i in range(10))
BODY_NO_EOL = BODY.rstrip(b"\n")


class TestJoin(unittest.TestCase):
    def test_roundtrip(self):
        for body in (BODY, BODY_NO_EOL):
            for k in (1, 2, 3, 4):
                man, bodies = split(body, k)
                self.assertEqual(pp.join_parts(man, bodies), body)

    def test_bad_part_sha(self):
        man, bodies = split(BODY, 3)
        bodies["p2.csv"] = bodies["p2.csv"] + b"x\n"
        with self.assertRaises(pp.PartsError):
            pp.join_parts(man, bodies)

    def test_missing_part(self):
        man, bodies = split(BODY, 3)
        del bodies["p3.csv"]
        with self.assertRaises(pp.PartsError):
            pp.join_parts(man, bodies)

    def test_wrong_total_sha(self):
        man, bodies = split(BODY, 2)
        man["sha256"] = "0" * 64
        with self.assertRaises(pp.PartsError):
            pp.join_parts(man, bodies)

    def test_empty_manifest(self):
        with self.assertRaises(pp.PartsError):
            pp.join_parts({"parts": []}, {})


class TestWaitFor(unittest.TestCase):
    def setUp(self):
        man, bodies = split(BODY, 3)
        self.store = {pp.MANIFEST: json.dumps(man).encode()}
        for name, b in bodies.items():
            self.store[pp.part_path(name)] = b

    def fake_clock(self):
        t = [0.0]

        def clock():
            return t[0]

        def sleep(s):
            t[0] += s
        return clock, sleep

    def test_serves_pairs_from_parts(self):
        clock, sleep = self.fake_clock()
        rec = {}
        wf = pp.make_wait_for(lambda *a: ("orig", 0), lambda p: self.store[p], "results/x.csv", rec,
                              sleep=sleep, clock=clock)
        body, mins = wf("results/x.csv", b"ev1,", 10)
        self.assertEqual(body, BODY)
        self.assertEqual(mins, 0.0)
        self.assertEqual(len(rec["parts"]), 3)

    def test_other_paths_go_to_orig(self):
        wf = pp.make_wait_for(lambda p, m, w: ("orig:" + p, 1.0), lambda p: b"", "results/x.csv", {})
        self.assertEqual(wf("findings/plan.md", b"x", 5), ("orig:findings/plan.md", 1.0))

    def test_times_out_when_absent(self):
        clock, sleep = self.fake_clock()
        rec = {}

        def fetch(p):
            raise IOError("404")
        wf = pp.make_wait_for(None, fetch, "results/x.csv", rec, sleep=sleep, clock=clock)
        body, mins = wf("results/x.csv", b"ev1,", 5)
        self.assertIsNone(body)
        self.assertTrue(mins >= 5)
        self.assertIn("404", rec["error"])

    def test_waits_then_succeeds(self):
        clock, sleep = self.fake_clock()
        calls = {"n": 0}

        def fetch(p):
            if clock() < 300:
                calls["n"] += 1
                raise IOError("not yet")
            return self.store[p]
        wf = pp.make_wait_for(None, fetch, "results/x.csv", {}, sleep=sleep, clock=clock)
        body, mins = wf("results/x.csv", b"ev1,", 30)
        self.assertEqual(body, BODY)
        self.assertTrue(calls["n"] >= 2)


if __name__ == "__main__":
    unittest.main()

# Copyright 2026 MV2, Munim, Inc.
# Licensed under the Apache License, Version 2.0 (Apache-2.0). See LICENSE.
# Made with AI agents working in MV2's name.
"""Unit tests for fdsnfix.py (stdlib only; no network, no obspy)."""

import unittest

import fdsnfix


class FakeClient:
    def __init__(self, base="X", timeout=120, _discover_services=True):
        self.base = base
        self.timeout = timeout
        self.discover = _discover_services


class TestArgs(unittest.TestCase):
    def test_split(self):
        own, rest = fdsnfix.split_args(["--sha", "abc1234", "--", "--data-dir", "/data"])
        self.assertEqual(own, ["--sha", "abc1234"])
        self.assertEqual(rest, ["--data-dir", "/data"])

    def test_split_no_dashes(self):
        own, rest = fdsnfix.split_args(["--sha", "abc1234"])
        self.assertEqual(own, ["--sha", "abc1234"])
        self.assertEqual(rest, [])

    def test_split_first_only(self):
        own, rest = fdsnfix.split_args(["a", "--", "b", "--", "c"])
        self.assertEqual(own, ["a"])
        self.assertEqual(rest, ["b", "--", "c"])


class TestUrl(unittest.TestCase):
    def test_raw_url(self):
        u = fdsnfix.raw_url("3df9ac8dec3a46d75b9be75265f7a9a7938acd7f", "chilepilot.py")
        self.assertEqual(u, "https://raw.githubusercontent.com/munimtechnologies/inner-core-repeaters/"
                            "3df9ac8dec3a46d75b9be75265f7a9a7938acd7f/chilepilot.py")

    def test_raw_url_rejects_branch_and_paths(self):
        for sha, path in (("main", "a.py"), ("", "a.py"), ("abc12", "a.py"),
                          ("abc1234", "/etc/passwd"), ("abc1234", "x/../a.py")):
            with self.assertRaises(ValueError):
                fdsnfix.raw_url(sha, path)

    def test_sha256(self):
        self.assertEqual(fdsnfix.sha256(b"abc"),
                         "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad")


class TestPatch(unittest.TestCase):
    def test_default_off_and_override_kept(self):
        class C(FakeClient):
            pass
        fdsnfix.default_no_discovery(C)
        self.assertFalse(C("B").discover)
        self.assertEqual(C("B", timeout=5).timeout, 5)
        self.assertTrue(C("B", _discover_services=True).discover)
        self.assertEqual(C("B").base, "B")

    def test_idempotent(self):
        class C(FakeClient):
            pass
        fdsnfix.default_no_discovery(C)
        first = C.__init__
        fdsnfix.default_no_discovery(C)
        self.assertIs(C.__init__, first)
        self.assertFalse(C().discover)

    def test_base_class_untouched(self):
        class C(FakeClient):
            pass
        fdsnfix.default_no_discovery(C)
        self.assertTrue(FakeClient().discover)


if __name__ == "__main__":
    unittest.main()

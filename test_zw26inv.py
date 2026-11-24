# Copyright 2026 MV2, Munim, Inc.
# Licensed under the Apache License, Version 2.0 (Apache-2.0). See LICENSE.
# Made with AI agents working in MV2's name.
"""Unit tests for zw26inv.py (pure helpers and safe extraction; no network)."""

import io
import os
import tempfile
import unittest
import zipfile

import zw26inv as z


class TestHelpers(unittest.TestCase):
    def test_safe_member(self):
        self.assertEqual(z.safe_member("a/b.SAC"), "a/b.SAC")
        self.assertEqual(z.safe_member("a\\b.SAC"), "a/b.SAC")
        self.assertIsNone(z.safe_member("../evil"))
        self.assertIsNone(z.safe_member("/abs/x"))
        self.assertIsNone(z.safe_member("C:/x"))
        self.assertIsNone(z.safe_member(""))

    def test_kind_of(self):
        self.assertEqual(z.kind_of("x/README"), "text")
        self.assertEqual(z.kind_of("x/list.txt"), "text")
        self.assertEqual(z.kind_of("x/IL01.BHZ.SAC"), "sac")
        self.assertEqual(z.kind_of("x/a.sac"), "sac")
        self.assertEqual(z.kind_of("x/a.mseed"), "mseed")
        self.assertEqual(z.kind_of("x/inner.zip"), "zip")
        self.assertEqual(z.kind_of("x/IL01.z"), "other")

    def test_origin_key(self):
        self.assertEqual(z.origin_key("2003-03-21T03:30:00.000000Z", 17.6),
                         "2003-03-21T03:30:18")
        self.assertEqual(z.origin_key("2003-03-21T03:30:00Z", 0.2), "2003-03-21T03:30:00")
        self.assertEqual(z.origin_key("2003-03-21T23:59:59.900000Z", 0.0),
                         "2003-03-22T00:00:00")
        self.assertEqual(z.origin_key("2003-03-21T03:30:00Z", None), "")

    def test_summarise_events(self):
        rows = [dict(status="ok", bundle="B", origin="2003-03-21T03:30:18", sta="IL01",
                     cha="SHZ", sr=20.0, member="e1/IL01.SHZ.SAC", sac_kevnm="M03"),
                dict(status="ok", bundle="B", origin="2003-03-21T03:30:18", sta="IL02",
                     cha="SHZ", sr=20.0, member="e1/IL02.SHZ.SAC", sac_kevnm=""),
                dict(status="ok", bundle="B", origin="", sta="IL01", cha="SHZ", sr=40.0,
                     member="e2/IL01.SHZ.SAC"),
                dict(status="unreadable: x", bundle="B", member="e3/bad")]
        ev = z.summarise_events(rows)
        self.assertEqual(len(ev), 2)
        e1 = [e for e in ev if e["origin"] == "2003-03-21T03:30:18"][0]
        self.assertEqual(e1["n_traces"], 2)
        self.assertEqual(e1["n_stations"], 2)
        self.assertEqual(e1["stations"], "IL01 IL02")
        self.assertEqual(e1["rates"], "20")
        self.assertEqual(e1["kevnm"], "M03")
        e2 = [e for e in ev if e["origin"].startswith("folder:")][0]
        self.assertEqual(e2["origin"], "folder:e2")
        self.assertEqual(e2["rates"], "40")

    def test_quantiles_and_rel(self):
        q = z._quantiles([3.0, None, 1.0, 2.0])
        self.assertEqual(q["n"], 3)
        self.assertEqual(q["min"], 1.0)
        self.assertEqual(q["median"], 2.0)
        self.assertEqual(q["max"], 3.0)
        self.assertEqual(z._quantiles([]), {})
        r = dict(origin="2003-03-21T03:30:18", start="2003-03-21T03:48:18.500000Z")
        self.assertAlmostEqual(z._rel(r), 1080.5, places=6)
        self.assertIsNone(z._rel(dict(origin="", start="x")))


class TestExtract(unittest.TestCase):
    def test_extract_nested_and_unsafe(self):
        with tempfile.TemporaryDirectory() as d:
            inner = io.BytesIO()
            with zipfile.ZipFile(inner, "w") as zi:
                zi.writestr("in/a.txt", "hello")
            outer = os.path.join(d, "o.zip")
            with zipfile.ZipFile(outer, "w") as zo:
                zo.writestr("top/readme.txt", "x" * 10)
                zo.writestr("top/inner.zip", inner.getvalue())
                zo.writestr("../evil.txt", "no")
            rows = []
            got = z.extract(outer, os.path.join(d, "x"), "o.zip", rows)
            names = sorted(m for m, _ in got)
            self.assertEqual(names, ["top/inner.zip/in/a.txt", "top/readme.txt"])
            self.assertFalse(os.path.exists(os.path.join(d, "evil.txt")))
            self.assertEqual(len(rows), 4)  # 3 outer members + 1 nested
            budget = [z.TEXT_TOTAL_MAX]
            out = os.path.join(d, "out")
            cp = z.copy_texts("o.zip", got, out, budget)
            self.assertTrue(all(c["copied"] for c in cp))
            self.assertTrue(os.path.exists(os.path.join(out, "text", "o", "top", "readme.txt")))
            self.assertEqual(budget[0], z.TEXT_TOTAL_MAX - 15)


if __name__ == "__main__":
    unittest.main()

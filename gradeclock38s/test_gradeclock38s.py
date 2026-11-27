"""Stdlib-only, offline tests for gradeclock38s.py (pre-registration 38, arm WS). Expected values were
derived by hand. Part of inner-core-repeaters by MV2, Munim, Inc.; Apache-2.0; made with AI agents."""
import unittest

import gradeclock38s as s

WSC = {
    "version": "0.2.0",
    "finished": "2026-12-01T00:00:00Z",
    "census": {"headline": "some source problems"},
    "multiplets": {
        "A": {"verdict": "clean", "odd_members": [],
              "pairs": {"712-718": {"decision": "N"}, "734-859": {"decision": "inconclusive"}}},
        "F": {"verdict": "odd", "odd_members": [744],
              "pairs": {"709-730": {"decision": "N"}, "730-744": {"decision": "S"},
                        "730-856": {"decision": "N"}, "744-856": {"decision": "S"}}},
        "B": {"verdict": "odd", "odd_members": [],
              "pairs": {"702-708": {"decision": "S"}, "702-723": {"decision": "untestable"}}},
    },
}


def row(e1, e2, mult="F", yka=0, ilar=0, pair="P00"):
    return dict(pair=pair, e1=e1, e2=e2, mult=mult, YKA=yka, ILAR=ilar, t1=2000.0, t2=2010.0)


class TestFlags(unittest.TestCase):
    def setUp(self):
        self.odd, self.dec = s.wsc_flags(WSC)

    def test_parse(self):
        self.assertEqual(self.odd, {"744"})
        self.assertEqual(self.dec[frozenset(("730", "744"))], "S")
        self.assertEqual(self.dec[frozenset(("723", "702"))], "untestable")
        self.assertEqual(len(self.dec), 8)

    def test_flag(self):
        f = lambda a, b: s.flag(row(a, b), self.odd, self.dec)  # noqa: E731
        self.assertEqual(f("709", "744"), "odd-member")      # 744 odd even without a pair decision
        self.assertEqual(f("744", "856"), "odd-member")
        self.assertEqual(f("730", "709"), "N")                 # order does not matter
        self.assertEqual(f("702", "708"), "S")                 # S pair, no odd member
        self.assertEqual(f("702", "723"), "untestable")
        self.assertEqual(f("734", "859"), "inconclusive")
        self.assertEqual(f("718", "742"), "yardstick")         # both M03-24, not graded by wsc39
        self.assertEqual(f("709", "855"), "untested")          # 855 not a Table 3 member
        self.assertEqual(f("848", "859"), "untested")          # only one M03-24 event, no decision

    def test_select(self):
        rows = [row("709", "744"), row("702", "708"), row("709", "730"), row("718", "742"),
                row("702", "723"), row("709", "855")]
        flags = [s.flag(r, self.odd, self.dec) for r in rows]
        ws1 = s.select(rows, flags, drop=s.DROP_WS1)
        self.assertEqual([(r["e1"], r["e2"]) for r in ws1],
                         [("709", "730"), ("718", "742"), ("702", "723"), ("709", "855")])
        ws2 = s.select(rows, flags, keep=s.KEEP_WS2)
        self.assertEqual([(r["e1"], r["e2"]) for r in ws2], [("709", "730"), ("718", "742")])


class TestGate(unittest.TestCase):
    def test_gate(self):
        self.assertTrue(s.wsc_gate(WSC))
        self.assertFalse(s.wsc_gate(dict(WSC, version="0.1.0")))
        self.assertFalse(s.wsc_gate(dict(WSC, finished=None)))
        bad = dict(WSC, multiplets={"A": {"pairs": {"712-724": {"decision": "S"}}}})
        self.assertFalse(s.wsc_gate(bad))


class TestReading(unittest.TestCase):
    def test_read_array(self):
        self.assertEqual(s.read_array(("structure", "turn-needed")), "survives")
        self.assertEqual(s.read_array(("structure", "weak")), "weakened")
        self.assertEqual(s.read_array(("weak", "turn-needed")), "weakened")
        self.assertEqual(s.read_array(("structure", "not-needed")), "lost")
        self.assertEqual(s.read_array(("no-structure", "turn-needed")), "lost")
        self.assertEqual(s.read_array(("untestable", "untestable")), "untestable")

    def test_headline(self):
        self.assertEqual(s.headline({"YKA": "survives", "ILAR": "survives"}), "survives")
        self.assertEqual(s.headline({"YKA": "lost", "ILAR": "lost"}), "lost")
        self.assertEqual(s.headline({"YKA": "survives", "ILAR": "lost"}), "mixed")
        self.assertEqual(s.headline({"YKA": "survives", "ILAR": "untestable"}), "mixed")

    def test_constants(self):
        self.assertEqual(s.REF_WANG, frozenset(("718", "734", "742", "848", "854")))
        self.assertEqual(s.ARM_W_EXPECTED["n_comp_YKA"], 193)


if __name__ == "__main__":
    unittest.main()

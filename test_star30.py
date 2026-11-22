# Tests for star30.py (synthetic only). Apache-2.0. Made with AI agents working in MV2's name.
import unittest

import star30

TXT = """Event Date,Event Time,Latitude,Longitude,Magnitude
1977-11-24*,16:59:59.92,-21.884,-138.886,5.86

1985-05-08,20:28:00.24,-21.831,-138.981,5.64
 1995-12-27* ,21:29:58.46,-21.792,-139.042,5.1
"""


class T(unittest.TestCase):
    def test_keeps_header_and_starred(self):
        out = star30.starred_only(TXT).strip().splitlines()
        self.assertEqual(len(out), 3)
        self.assertTrue(out[0].startswith("Event Date"))
        self.assertTrue(out[1].startswith("1977-11-24*"))
        self.assertTrue(out[2].strip().startswith("1995-12-27*"))

    def test_drops_unstarred(self):
        self.assertNotIn("1985-05-08", star30.starred_only(TXT))

    def test_empty(self):
        self.assertEqual(star30.starred_only(""), "")


if __name__ == "__main__":
    unittest.main()

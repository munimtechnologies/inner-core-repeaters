"""Tests for get_frozen.py. Expected blob SHAs are git's own well-known values. Standard library only.
Part of inner-core-repeaters by MV2, Munim, Inc. (made with AI agents). Licence: Apache-2.0."""
import unittest

import get_frozen as g


class TestBlobSha(unittest.TestCase):
    def test_known_values(self):
        # `printf '' | git hash-object --stdin` and `echo hello | git hash-object --stdin`
        self.assertEqual(g.blob_sha(b""), "e69de29bb2d1d6434b8b29ae775ad8c2e48c5391")
        self.assertEqual(g.blob_sha(b"hello\n"), "ce013625030ba8dba906f756967f9e9ca394464a")

    def test_pins(self):
        self.assertEqual(g.COMMIT, "a7b58bdb1d27cfe0fce09fd54ced519d557e5346")
        self.assertEqual(sorted(g.FILES), ["m2measure/m2measure.py", "m2measure/test_m2measure.py"])
        for v in g.FILES.values():
            self.assertEqual(len(v), 40)


if __name__ == "__main__":
    unittest.main()

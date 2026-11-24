# Tests for stage31h.py (data route of pre-registration 31). Synthetic tars only; stdlib only.
# Part of munimtechnologies/inner-core-repeaters. Licence: Apache-2.0.
# Made with AI agents working in MV2's name (MV2, Munim, Inc.).
import gzip
import io
import os
import subprocess
import tarfile
import tempfile
import unittest

import stage31h as S


def _tar_bytes(files):
    """{name: bytes} -> bytes of a plain tar."""
    buf = io.BytesIO()
    with tarfile.open(fileobj=buf, mode="w") as tf:
        for n, b in files.items():
            ti = tarfile.TarInfo(n)
            ti.size = len(b)
            tf.addfile(ti, io.BytesIO(b))
    return buf.getvalue()


class TestParse(unittest.TestCase):
    def test_file_member(self):
        u = S.parse_member("Exp/19771124-41-TS.tar.Z")
        self.assertEqual((u["kind"], u["ymd"], u["site"], u["seq"], u["sys"], u["ext"]),
                         ("file", "19771124", "4", 1, "TS", ".Z"))
        self.assertEqual(u["key"], "Exp/19771124-41-TS.tar.Z")

    def test_dir_member(self):
        u = S.parse_member("./FRENCH/19830525-42-SS/brv.wfdisc")
        self.assertEqual((u["kind"], u["ymd"], u["seq"], u["sys"]), ("dir", "19830525", 2, "SS"))
        self.assertEqual(u["key"], "FRENCH/19830525-42-SS")

    def test_other(self):
        self.assertIsNone(S.parse_member("README"))
        self.assertIsNone(S.parse_member("FRENCH/notes.txt"))

    def test_units_distinct(self):
        names = ["A/19830525-41-SS/a.wfdisc", "A/19830525-41-SS/a.w", "A/19830525-41-TS.tar.Z"]
        us = S.units_from_members(names)
        self.assertEqual(sorted(u["key"] for u in us), ["A/19830525-41-SS", "A/19830525-41-TS.tar.Z"])


class TestChoose(unittest.TestCase):
    def setUp(self):
        self.units = S.units_from_members([
            "E/19800720-41-SS.tar.Z", "E/19800720-41-TS.tar.Z",
            "E/19881123-61-TS.tar.Z",                       # Soviet test on the same day: never chosen
            "E/19881123-41-TS.tar.Z", "E/19881123-42-TS.tar.Z",
            "E/19851126-41-TS.tar.Z", "E/19851125-41-TS.tar.Z",
        ])

    def test_next_day_name(self):
        u, off, _ = S.choose_source(self.units, "19800719", "TS")
        self.assertEqual((u["key"], off), ("E/19800720-41-TS.tar.Z", 1))

    def test_site_and_seq(self):
        u, off, ties = S.choose_source(self.units, "19881123", "TS")
        self.assertEqual((u["key"], off, ties), ("E/19881123-41-TS.tar.Z", 0, 1))

    def test_exact_day_beats_neighbour(self):
        u, off, _ = S.choose_source(self.units, "19851126", "TS")
        self.assertEqual((u["key"], off), ("E/19851126-41-TS.tar.Z", 0))

    def test_none(self):
        u, off, ties = S.choose_source(self.units, "19900101", "TS")
        self.assertIsNone(u)
        self.assertIsNone(off)
        u, _, _ = S.choose_source(self.units, "19881123", "SS")
        self.assertIsNone(u)

    def test_day_offset(self):
        self.assertEqual(S.day_offset("19801231", "19810101"), 1)
        self.assertEqual(S.day_offset("19810101", "19801231"), -1)


class TestWrite(unittest.TestCase):
    def test_copy_Z_bytes_exactly(self):
        payload = b"\x1f\x9d\x90compressed-ish"
        master = _tar_bytes({"E/19771124-41-TS.tar.Z": payload})
        with tempfile.TemporaryDirectory() as d:
            mp = os.path.join(d, "m.tar")
            with open(mp, "wb") as f:
                f.write(master)
            with tarfile.open(mp) as tf:
                names = tf.getnames()
                u = S.units_from_members(names)[0]
                dest = os.path.join(d, "19771124-41-TS.tar.Z")
                S.write_unit(tf, u, names, dest)
            with open(dest, "rb") as f:
                self.assertEqual(f.read(), payload)

    def test_plain_tar_member_gzipped_under_Z_name(self):
        inner = _tar_bytes({"x.wfdisc": b"row\n", "x.w": b"\x00\x01"})
        master = _tar_bytes({"E/19771124-41-TS.tar": inner})
        with tempfile.TemporaryDirectory() as d:
            mp = os.path.join(d, "m.tar")
            with open(mp, "wb") as f:
                f.write(master)
            with tarfile.open(mp) as tf:
                names = tf.getnames()
                u = S.units_from_members(names)[0]
                dest = os.path.join(d, "19771124-41-TS.tar.Z")
                S.write_unit(tf, u, names, dest)
            with gzip.open(dest, "rb") as f:
                self.assertEqual(f.read(), inner)

    def test_directory_unit_packed(self):
        master = _tar_bytes({"F/19830525-41-SS/a.wfdisc": b"row\n", "F/19830525-41-SS/a.w": b"\x07" * 9,
                             "F/19830525-41-TS/b.w": b"no"})
        with tempfile.TemporaryDirectory() as d:
            mp = os.path.join(d, "m.tar")
            with open(mp, "wb") as f:
                f.write(master)
            with tarfile.open(mp) as tf:
                names = tf.getnames()
                u = [x for x in S.units_from_members(names) if x["sys"] == "SS"][0]
                dest = os.path.join(d, "19830525-41-SS.tar.Z")
                S.write_unit(tf, u, names, dest)
            with gzip.open(dest, "rb") as g, tarfile.open(fileobj=io.BytesIO(g.read())) as t2:
                self.assertEqual(sorted(t2.getnames()), ["19830525-41-SS/a.w", "19830525-41-SS/a.wfdisc"])
                self.assertEqual(t2.extractfile("19830525-41-SS/a.w").read(), b"\x07" * 9)

    def test_gzip_cli_reads_it(self):
        # The measurement code decompresses .tar.Z with `gzip -dc`; check that route when gzip exists.
        inner = _tar_bytes({"x.wfdisc": b"row\n"})
        master = _tar_bytes({"E/19771124-41-TS.tar": inner})
        with tempfile.TemporaryDirectory() as d:
            mp = os.path.join(d, "m.tar")
            with open(mp, "wb") as f:
                f.write(master)
            with tarfile.open(mp) as tf:
                names = tf.getnames()
                dest = os.path.join(d, "19771124-41-TS.tar.Z")
                S.write_unit(tf, S.units_from_members(names)[0], names, dest)
            try:
                out = subprocess.run(["gzip", "-dc", dest], capture_output=True, check=True).stdout
            except (OSError, subprocess.CalledProcessError):
                self.skipTest("gzip not available")
            self.assertEqual(out, inner)


if __name__ == "__main__":
    unittest.main()

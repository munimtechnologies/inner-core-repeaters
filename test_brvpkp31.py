# Unit tests for brvpkp31.py (pre-registration 31). Synthetic data only; stdlib + numpy.
# Part of munimtechnologies/inner-core-repeaters. Licence: Apache-2.0.
# Made with AI agents working in MV2's name (MV2, Munim, Inc.).
import calendar
import io
import math
import os
import tarfile
import tempfile
import unittest

import numpy as np

import brvpkp31 as B


class TestTable(unittest.TestCase):
    def test_parse_s1(self):
        text = ("Event Date,Event Time,Latitude,Longitude,Magnitude\n"
                "1977-11-24*,16:59:59.92,-21.884,-138.886,5.86\n"
                "1985-05-08,20:28:00.24,-21.831,-138.981,5.64\n")
        ev = B.parse_s1(text)
        self.assertEqual(len(ev), 2)
        self.assertTrue(ev[0]["starred"])
        self.assertFalse(ev[1]["starred"])
        self.assertEqual(ev[0]["ymd"], "19771124")
        self.assertEqual(ev[0]["year"], 1977)
        self.assertAlmostEqual(ev[0]["origin"], calendar.timegm((1977, 11, 24, 16, 59, 0, 0, 0, 0)) + 59.92, places=6)
        self.assertAlmostEqual(ev[1]["lon"], -138.981)

    def test_form_pairs(self):
        ev = [dict(lat=0.0, lon=0.0), dict(lat=0.0, lon=0.04), dict(lat=0.0, lon=0.2)]
        self.assertEqual(B.form_pairs(ev, 0.05), [(0, 1)])


class TestIndex(unittest.TestCase):
    def test_links(self):
        html = ('<a href="ftp://ftp.ldeo.columbia.edu/archive/wykim/Brv_archive/Exp/19771124-41-TS.tar.Z">TS</a> '
                '<a href=Exp/19771124-41-SS.tar.Z>SS</a> <a href="x.html">x</a>')
        links = B.parse_index_links(html, "https://www.example.edu/a/b.html")
        self.assertEqual(links, ["ftp://ftp.ldeo.columbia.edu/archive/wykim/Brv_archive/Exp/19771124-41-TS.tar.Z",
                                 "https://www.example.edu/a/Exp/19771124-41-SS.tar.Z"])
        self.assertEqual(B.links_for(links, "19771124", "TS"), [links[0]])
        self.assertEqual(B.links_for(links, "19771124", "SS"), [links[1]])
        self.assertEqual(B.links_for(links, "19820725", "TS"), [])
        self.assertEqual(B.system_of(links[0]), "TS")
        self.assertEqual(B.url_variants("ftp://h/p.tar.Z"), ["ftp://h/p.tar.Z", "https://h/p.tar.Z", "http://h/p.tar.Z"])


class TestWfdisc(unittest.TestCase):
    def test_fixed_roundtrip(self):
        rec = dict(sta="BRV", chan="sz", time=375382798.6, endtime=375384598.56875, nsamp=57600, samprate=32.0,
                   datatype="s4", dir=".", dfile="19771124.w", foff=1024, calib=0.5)
        line = B.format_wfdisc_line(rec)
        self.assertEqual(len(line), 283)
        p = B.parse_wfdisc_line(line)
        self.assertEqual(p["sta"], "BRV")
        self.assertEqual(p["chan"], "sz")
        self.assertAlmostEqual(p["time"], 375382798.6, places=4)
        self.assertEqual(p["nsamp"], 57600)
        self.assertAlmostEqual(p["samprate"], 32.0)
        self.assertEqual(p["datatype"], "s4")
        self.assertEqual(p["dir"], ".")
        self.assertEqual(p["dfile"], "19771124.w")
        self.assertEqual(p["foff"], 1024)
        self.assertAlmostEqual(p["calib"], 0.5)

    def test_whitespace_fallback(self):
        line = "BRV sz 375382798.6 1 1 1977328 375384598.6 54000 30.0 1.0 1.0 - - s4 - . f.w 0 -1 -"
        p = B.parse_wfdisc_line(line)
        self.assertEqual(p["chan"], "sz")
        self.assertEqual(p["nsamp"], 54000)
        self.assertAlmostEqual(p["samprate"], 30.0)
        self.assertEqual(p["dfile"], "f.w")
        self.assertEqual(p["foff"], 0)
        self.assertEqual(p["datatype"], "s4")
        self.assertIsNone(B.parse_wfdisc_line(""))
        self.assertIsNone(B.parse_wfdisc_line("garbage line"))

    def test_decode(self):
        np.testing.assert_array_equal(B.decode_samples(np.array([1, -2, 300000], ">i4").tobytes(), "s4"),
                                      [1, -2, 300000])
        np.testing.assert_array_equal(B.decode_samples(np.array([5, -7], ">i2").tobytes(), "s2"), [5, -7])
        np.testing.assert_array_equal(B.decode_samples(np.array([9, -3], "<i4").tobytes(), "i4"), [9, -3])
        np.testing.assert_array_equal(B.decode_samples(bytes([0, 0, 1, 255, 255, 255, 128, 0, 0]), "s3"),
                                      [1, -1, -8388608])
        np.testing.assert_array_equal(B.decode_samples(np.array([1.5], ">f4").tobytes(), "t4"), [1.5])
        with self.assertRaises(ValueError):
            B.decode_samples(b"\x00\x00", "e1")

    def test_read_wf(self):
        with tempfile.TemporaryDirectory() as tmp:
            os.makedirs(os.path.join(tmp, "sub"))
            with open(os.path.join(tmp, "sub", "d.w"), "wb") as f:
                f.write(b"\x00" * 16 + np.array([10, 20, 30], ">i4").tobytes())
            base = dict(sta="BRV", chan="sz", time=100.0, endtime=100.1, nsamp=3, samprate=20.0, datatype="s4",
                        dfile="d.w", foff=16)
            for d in ("sub", "/abs/elsewhere"):
                rec = B.parse_wfdisc_line(B.format_wfdisc_line(dict(base, dir=d)))
                t0, dt, x, note = B.read_wf(tmp, tmp, rec)
                self.assertAlmostEqual(t0, 100.0)
                self.assertAlmostEqual(dt, 0.05)
                np.testing.assert_array_equal(x, [10, 20, 30])
                self.assertEqual(note, "")

    def test_safe_extract(self):
        with tempfile.TemporaryDirectory() as tmp:
            tp = os.path.join(tmp, "t.tar")
            with tarfile.open(tp, "w") as tf:
                for name, data in (("db/x.wfdisc", b"abc"), ("../evil.txt", b"no"), ("/abs/y.w", b"yy")):
                    ti = tarfile.TarInfo(name)
                    ti.size = len(data)
                    tf.addfile(ti, io.BytesIO(data))
            dest = os.path.join(tmp, "out")
            os.makedirs(dest)
            n = B.safe_extract(tp, dest)
            self.assertEqual(n, 2)
            self.assertTrue(os.path.isfile(os.path.join(dest, "db", "x.wfdisc")))
            self.assertTrue(os.path.isfile(os.path.join(dest, "abs", "y.w")))
            self.assertFalse(os.path.exists(os.path.join(tmp, "evil.txt")))


class TestSignal(unittest.TestCase):
    def test_despike(self):
        n = 2000
        t = np.arange(n)
        rng = np.random.default_rng(0)
        clean = np.sin(2 * np.pi * t / 100.0) + 0.01 * rng.standard_normal(n)
        x = clean.copy()
        x[500] += 50.0
        x[1200] -= 40.0
        x[1201] -= 40.0
        y, cnt = B.despike(x)
        self.assertEqual(cnt, 3)
        self.assertEqual(int(np.count_nonzero(y != x)), 3)
        self.assertLess(float(np.max(np.abs(y - clean))), 0.1)

    def test_clipped_mask(self):
        t = np.arange(400)
        s = 1000.0 * np.sin(2 * np.pi * t / 100.0)
        self.assertEqual(int(B.clipped_mask(s).sum()), 0)
        c = np.clip(s, -800.0, 800.0)
        m = B.clipped_mask(c)
        self.assertTrue(m[25])
        self.assertTrue(m[75])
        self.assertFalse(m[0])
        self.assertFalse(m[50])

    def test_envelope(self):
        t = np.arange(1000) / 100.0
        env = B.envelope(np.sin(2 * np.pi * 5.0 * t))
        self.assertLess(abs(env[500] - 1.0), 1e-6)

    def test_seg_and_lag(self):
        dt = 0.025
        t = np.arange(0.0, 200.0, dt)

        def wav(tc):
            return np.exp(-((t - tc) / 1.5) ** 2) * np.sin(2 * np.pi * 0.8 * (t - tc))
        A, Bt = (0.0, dt, wav(100.0)), (0.0, dt, wav(100.23))
        r = B.lag_cc(A, Bt, 95.0, 95.0, 10.0, 0.6, 100.0)
        self.assertTrue(r["ok"])
        self.assertAlmostEqual(r["lag"], 0.23, delta=0.01)
        self.assertGreater(r["cc"], 0.99)
        self.assertFalse(r["edge"])
        self.assertIsNone(B.seg(A, -1.0, 10, 100.0))
        self.assertIsNone(B.seg(A, 199.95, 10, 100.0))

    def test_lag_edge(self):
        dt = 0.025
        t = np.arange(0.0, 200.0, dt)
        A = (0.0, dt, np.exp(-((t - 100.0) / 0.5) ** 2))
        Bt = (0.0, dt, np.exp(-((t - 101.5) / 0.5) ** 2))
        r = B.lag_cc(A, Bt, 97.0, 97.0, 6.0, 0.6, 100.0)
        self.assertTrue(r["edge"])

    def test_measure_pair(self):
        dt = 0.025
        t = np.arange(0.0, 400.0, dt)

        def wv(tc):
            return np.exp(-((t - tc) / 1.5) ** 2) * np.sin(2 * np.pi * 0.8 * (t - tc))
        xa = wv(100.0) + 0.8 * wv(220.0) + 0.5 * wv(330.0)
        xb = wv(103.1) + 0.8 * wv(223.0) + 0.5 * wv(332.95)
        A, Bt = (0.0, dt, xa), (0.0, dt, xb)
        ea = (0.0, dt, B.smooth(B.envelope(xa), 40))
        eb = (0.0, dt, B.smooth(B.envelope(xb), 40))
        times = {"PKIKP": 100.0, "PP": 220.0, "PPP": 330.0}
        m = B.measure_pair(A, Bt, ea, eb, times, times)
        self.assertFalse(m["edge_c"])
        self.assertAlmostEqual(m["c"], 3.1, delta=0.1)
        self.assertAlmostEqual(m["dK"], 0.10, delta=0.01)
        self.assertAlmostEqual(m["dP3"], -0.05, delta=0.01)
        self.assertGreater(m["cc_K"], 0.99)

    def test_choose_channel(self):
        sa = {"BRV.sz": {"eligible": True, "snr_K": 10}, "BRV.bz": {"eligible": True, "snr_K": 30},
              "BRV.lz": {"eligible": False, "snr_K": 99}}
        sb = {"BRV.sz": {"eligible": True, "snr_K": 20}, "BRV.bz": {"eligible": True, "snr_K": 5},
              "BRV.lz": {"eligible": True, "snr_K": 99}}
        self.assertEqual(B.choose_channel(sa, sb), "BRV.sz")
        self.assertIsNone(B.choose_channel({}, sb))


class TestInversion(unittest.TestCase):
    A = {1977: 0.0, 1982: -0.1, 1988: -0.15, 1995: 0.2}

    def test_pair_year_terms(self):
        years = [1977, 1982, 1988, 1995, 1982]
        pairs = [(i, j, self.A[years[j]] - self.A[years[i]]) for i in range(5) for j in range(i + 1, 5)]
        yt = B.pair_year_terms(pairs, years)
        self.assertEqual(sorted(yt), [1977, 1982, 1988, 1995])
        for y, v in self.A.items():
            self.assertAlmostEqual(yt[y], v, places=12)
        years2 = [1977, 1982, 1990, 1991]
        yt2 = B.pair_year_terms([(0, 1, -0.1), (2, 3, 0.05)], years2)
        self.assertEqual(sorted(yt2), [1977, 1982])
        self.assertAlmostEqual(yt2[1982], -0.1)
        self.assertEqual(B.pair_year_terms([(1, 2, 0.1)], years2), {})

    def test_epoch_terms(self):
        terms = {1977: 0.0, 1982: -0.1, 1983: -0.2, 1988: -0.1, 1989: -0.2, 1990: -0.3, 1995: 0.2}
        ep = B.epoch_terms(terms, [1977, 1982, 1983, 1983, 1988, 1989, 1990, 1995])
        self.assertAlmostEqual(ep["E1982_83"], -0.5 / 3)
        self.assertAlmostEqual(ep["E1988_90"], -0.2)
        self.assertAlmostEqual(ep["E1995"], 0.2)
        ep2 = B.epoch_terms({1977: 0.0, 1982: -0.1}, [1977, 1982])
        self.assertIsNone(ep2["E1995"])

    def test_bootstrap_exact(self):
        years = [1977, 1982, 1988, 1995]
        pairs = [(i, j, {"dK": self.A[years[j]] - self.A[years[i]]}) for i in range(4) for j in range(i + 1, 4)]
        boot = B.bootstrap_epochs(pairs, years, ["dK"], 30, 31)
        self.assertEqual(len(boot), 30)
        truth = {"E1982_83": -0.1, "E1988_90": -0.15, "E1995": 0.2}
        n_full = 0
        for s in boot:
            for e, v in s["dK"].items():
                if v is not None:
                    self.assertAlmostEqual(v, truth[e], places=9)
            n_full += all(v is not None for v in s["dK"].values())
        self.assertGreater(n_full, 0)


class TestGrades(unittest.TestCase):
    def test_grade_repro(self):
        self.assertEqual(B.grade_repro(dict(B.ZHOU))["grade"], "reproduced")
        g = B.grade_repro({"E1982_83": -0.1, "E1988_90": -0.1, "E1995": 0.05})
        self.assertEqual(g["grade"], "partly")
        self.assertEqual(g["signs_match"], 3)
        self.assertAlmostEqual(g["rms_res"], math.sqrt(0.00625), places=9)
        self.assertEqual(B.grade_repro({"E1982_83": 0.0, "E1988_90": 0.0, "E1995": 0.0})["grade"], "not-reproduced")
        self.assertEqual(B.grade_repro({"E1982_83": -0.1, "E1988_90": None, "E1995": 0.1})["grade"], "untestable")

    def test_hypotheses(self):
        ratios = {"core": 0.0, "position": -0.5, "time-base": -1.5}
        k = dict(B.ZHOU)
        p = {e: -0.5 * v for e, v in k.items()}
        rng = np.random.default_rng(1)
        bk = [{e: v + 0.01 * rng.standard_normal() for e, v in k.items()} for _ in range(200)]
        bp = [{e: v + 0.01 * rng.standard_normal() for e, v in p.items()} for _ in range(200)]
        chi = B.hypothesis_chi2(k, p, bk, bp, ratios)
        self.assertLess(chi["position"], 1e-9)
        self.assertGreater(chi["core"], 50.0)
        self.assertGreater(chi["time-base"], 50.0)
        self.assertEqual(B.grade_discriminant(chi), "position")
        self.assertIsNone(B.hypothesis_chi2(dict(k, E1995=None), p, bk, bp, ratios))
        self.assertIsNone(B.hypothesis_chi2(k, p, bk[:10], bp[:10], ratios))
        self.assertEqual(B.grade_discriminant(None), "untestable")
        self.assertEqual(B.grade_discriminant({"a": 1.0, "b": 2.0}), "unresolved")
        self.assertEqual(B.grade_discriminant({"a": 10.0, "b": 20.0}), "none-fits")

    def test_ratios_and_slopes(self):
        first = {"PKIKP": {"t": 1200.0, "p": 1.9}, "PP": {"t": 1320.0, "p": 5.6}, "PPP": {"t": 1500.0, "p": 7.5}}
        r = B.ratios_from(first)
        self.assertEqual(r["core"], 0.0)
        self.assertAlmostEqual(r["position"], 1.9 / -3.7)
        self.assertAlmostEqual(r["time-base"], -1.5)
        self.assertAlmostEqual(B.slope_through_zero([1, 2, 3], [2, 4, 6]), 2.0)
        self.assertAlmostEqual(B.tls_slope([1, 2, 3], [2, 4, 6]), 2.0)

    def test_clean(self):
        c = B._clean({"a": float("nan"), "_x": 1, "b": np.float64(2.0), "c": (1, 2), "d": np.bool_(True)})
        self.assertEqual(c, {"a": None, "b": 2.0, "c": [1, 2], "d": True})


if __name__ == "__main__":
    unittest.main()

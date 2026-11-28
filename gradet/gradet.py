#!/usr/bin/env python3
# Copyright 2026 MV2, Munim, Inc.
# Licensed under the Apache License, Version 2.0 (Apache-2.0). See LICENSE.
# Made with AI agents working in MV2's name.
"""gradet.py 0.1.0 - pre-registration 38, arm T: is the YKA disagreement between the two camps
bigger than the scatter of short-interval repeats?

Pre-registered in findings/2026-11-28-prereg38-armT-plan.md, written before any YKA record of a
Wang et al. (2024) multiplet other than M03-24 (Wang's A) was fetched or measured.

The three disputed YKA pairs of multiplet M03-24 (Zhang & Wen 2026 call them changed, Wang et al.
2024 grade them similar) have fixed, already published D_c310 values
(results/m0324chk-v2/yka_z004_pairs.csv). The yardstick is new: D_c310 of every usable YKA pair
less than 2 years apart in the other Wang multiplets, as measured by the frozen arm M job m2-v1
(m2measure 0.1.0), outside YKA's timing-error eras. Standard library only.
"""

import argparse
import calendar
import csv
import io
import json
import math
import os
import sys
import urllib.request

VERSION = "0.1.0"
NAN = float("nan")

REPO_RAW = "https://raw.githubusercontent.com/munimtechnologies/inner-core-repeaters/main/"
PAIRS_URL = REPO_RAW + "results/m2-v1/m2_pairs.csv"
M2_SUMMARY_URL = REPO_RAW + "results/m2-v1/m2measure_summary.json"
PUBLISHED_URL = REPO_RAW + "results/m0324chk-v2/yka_z004_pairs.csv"
UA = {"User-Agent": "inner-core-repeaters (MV2, Munim, Inc.)"}

MAX_DT_YR = 2.0          # yardstick pairs: |t2 - t1| < 2.0 yr
MIN_REF = 8              # fewer usable yardstick pairs -> untestable
Q_HI = 0.95              # nearest-rank percentile that a disputed pair must exceed
EXCLUDE_MULT = ("A",)    # multiplet A is M03-24 itself


def decimal_year(iso):
    """m2measure 0.1.0 convention: year + (day_of_year - 1 + fraction of day) / days."""
    s = iso.strip().rstrip("Z")
    date, hms = s.split("T")
    y, mo, d = (int(v) for v in date.split("-"))
    h, mi, sec = hms.split(":")
    t0 = calendar.timegm((y, 1, 1, 0, 0, 0, 0, 0, 0))
    t = calendar.timegm((y, mo, d, int(h), int(mi), 0, 0, 0, 0)) + float(sec)
    days = 366 if calendar.isleap(y) else 365
    return y + (t - t0) / 86400.0 / days


# YKA timing-error eras, exactly as in arm M (m2core v0.1.0; Ackerley 2026, SRL)
ERAS = ((decimal_year("2014-12-01T00:00:00Z"), decimal_year("2020-04-08T00:00:00Z")),
        (decimal_year("2021-12-01T00:00:00Z"), decimal_year("2022-04-01T00:00:00Z")))

# The disputed pairs, fixed: (ev1, ev2) -> (label, published D_c310 at YKA)
DISPUTED = {
    ("usp000fz4d", "usp000hute"): ("2008-02-17/2011-02-20", 0.017939),
    ("usp000bt5x", "us7000am0e"): ("2003-03-21/2020-07-11", 0.022813),
    ("us7000am0e", "us6000h2wz"): ("2020-07-11/2022-03-07", 0.020867),
}
DISPUTED_T = {"usp000bt5x": "2003-03-21T00:00:00Z", "usp000fz4d": "2008-02-17T00:00:00Z",
              "usp000hute": "2011-02-20T00:00:00Z", "us7000am0e": "2020-07-11T00:00:00Z",
              "us6000h2wz": "2022-03-07T00:00:00Z"}
# Report-only: pairs both camps call similar (published D_c310)
AGREED_SIMILAR = {
    ("usp000bt5x", "us6000h2wz"): ("2003-03-21/2022-03-07", 0.014711),
    ("usp000bt5x", "us6000mk65"): ("2003-03-21/2024-03-21", 0.012998),
}
# Report-only: the cut ranges each camp's YKA calls imply on the published D_c310
CUTS = {"zw26_low": 0.014711, "zw26_high": 0.017939, "wang_low": 0.022813, "wang_high": 0.12732,
        "armM_cut": 0.098}


def truthy(x):
    return x is True or str(x).strip().lower() in ("true", "1", "yes")


def fnum(x):
    try:
        v = float(x)
    except (TypeError, ValueError):
        return NAN
    return v


def isf(x):
    return isinstance(x, float) and math.isfinite(x)


def in_era(t):
    return any(lo <= t <= hi for lo, hi in ERAS)


def nearest_rank(vals, q):
    s = sorted(vals)
    if not s:
        return NAN
    k = max(0, int(math.ceil(q * len(s))) - 1)
    return s[min(k, len(s) - 1)]


def reference(rows):
    """Yardstick D values: YKA, not multiplet A, ok, finite D, |dt| < 2 yr, neither event in an era."""
    out = []
    for r in rows:
        if str(r.get("array")) != "YKA" or str(r.get("mult")) in EXCLUDE_MULT:
            continue
        if not truthy(r.get("ok")):
            continue
        d, t1, t2 = fnum(r.get("D")), fnum(r.get("t1")), fnum(r.get("t2"))
        if not (isf(d) and isf(t1) and isf(t2)):
            continue
        if abs(t2 - t1) >= MAX_DT_YR or in_era(t1) or in_era(t2):
            continue
        out.append(dict(mult=str(r.get("mult")), ev1=str(r.get("ev1")), ev2=str(r.get("ev2")),
                        dt=abs(t2 - t1), D=d))
    return out


def grade(ref_D, disputed=None):
    """ref_D: list of yardstick D values. Returns the arm T result dict."""
    disputed = DISPUTED if disputed is None else disputed
    n = len(ref_D)
    res = dict(n_ref=n, min_ref=MIN_REF)
    if n < MIN_REF:
        res.update(headline="untestable", reason="fewer than %d yardstick pairs" % MIN_REF, pairs=[])
        return res
    q50, q95 = nearest_rank(ref_D, 0.50), nearest_rank(ref_D, Q_HI)
    res.update(q50=q50, q95=q95)
    pairs, n_above = [], 0
    for (a, b), (label, d) in sorted(disputed.items(), key=lambda kv: kv[1][0]):
        above = d > q95
        n_above += int(above)
        era = any(in_era(decimal_year(DISPUTED_T[x])) for x in (a, b) if x in DISPUTED_T)
        frac_ge = sum(1 for v in ref_D if v >= d) / float(n)
        pairs.append(dict(pair=label, ev1=a, ev2=b, D=d, over_q50=d / q50 if q50 > 0 else NAN,
                          frac_ref_ge=frac_ge, verdict="above" if above else "within", era_pair=era))
    res["pairs"] = pairs
    res["n_above"] = n_above
    res["headline"] = "beyond-scatter" if n_above >= 2 else "within-scatter"
    # report-only
    res["report_only"] = dict(
        agreed_similar=[dict(pair=lab, D=d, verdict="above" if d > q95 else "within")
                        for (a, b), (lab, d) in sorted(AGREED_SIMILAR.items(), key=lambda kv: kv[1][0])],
        frac_ref_above_cut={k: sum(1 for v in ref_D if v > c) / float(n) for k, c in sorted(CUTS.items())},
        key_pair_2020_2022=[p["verdict"] for p in pairs if p["pair"] == "2020-07-11/2022-03-07"])
    return res


def read_csv_text(text):
    return list(csv.DictReader(io.StringIO(text))) if text else []


def fetch(url, tries=3):
    for _ in range(tries):
        try:
            req = urllib.request.Request(url, headers=UA)
            with urllib.request.urlopen(req, timeout=60) as r:
                return r.read().decode("utf-8")
        except Exception:  # noqa: BLE001
            pass
    return None


def check_published(rows):
    """The fixed disputed values must match the published rows to 1e-5."""
    by = {tuple(sorted((r.get("ev1", ""), r.get("ev2", "")))): fnum(r.get("D_c310")) for r in rows}
    out = {}
    for (a, b), (label, d) in list(DISPUTED.items()) + list(AGREED_SIMILAR.items()):
        v = by.get(tuple(sorted((a, b))), NAN)
        out[label] = dict(fixed=d, published=v, ok=isf(v) and abs(v - d) <= 1e-5)
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--pairs", default=PAIRS_URL, help="m2_pairs.csv path or URL")
    ap.add_argument("--m2-summary", default=M2_SUMMARY_URL, help="m2measure_summary.json path or URL")
    ap.add_argument("--out-dir", default="/out")
    a = ap.parse_args(argv)
    os.makedirs(a.out_dir, exist_ok=True)

    def load(p):
        if p.startswith("http"):
            return fetch(p)
        with open(p, encoding="utf-8") as f:
            return f.read()

    s = dict(version=VERSION, pairs_source=a.pairs)
    pub = check_published(read_csv_text(fetch(PUBLISHED_URL)))
    s["published_check"] = pub
    m2s = load(a.m2_summary)
    gate_passed = False
    if m2s:
        try:
            gate_passed = bool(json.loads(m2s).get("gate", {}).get("passed"))
        except ValueError:
            gate_passed = False
    s["m2_gate_passed"] = gate_passed
    rows = read_csv_text(load(a.pairs))
    ref = reference(rows)
    s["checks_ok"] = bool(rows) and all(v["ok"] for v in pub.values())
    if not gate_passed:
        s["result"] = dict(headline="untestable", reason="m2-v1 M03-24 reproduction gate did not pass")
    elif not s["checks_ok"]:
        s["result"] = dict(headline="not graded", reason="input checks failed")
    else:
        s["result"] = grade([r["D"] for r in ref])
    with open(os.path.join(a.out_dir, "gradet_reference.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["mult", "ev1", "ev2", "dt", "D"])
        w.writeheader()
        w.writerows(ref)
    with open(os.path.join(a.out_dir, "gradet_summary.json"), "w") as f:
        json.dump(s, f, indent=1)
    print(json.dumps(s.get("result"), indent=1)[:3000])
    return 0


if __name__ == "__main__":
    sys.exit(main())

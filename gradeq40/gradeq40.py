#!/usr/bin/env python3
# Copyright 2026 MV2, Munim, Inc.
# Licensed under the Apache License, Version 2.0 (Apache-2.0). See LICENSE.
# Made with AI agents working in MV2's name; no human expert has checked it.
"""gradeq40.py 0.1.0 - pre-registration 40: has the inner-core coda kept changing since 2020?

Part of inner-core-repeaters by MV2, Munim, Inc. (mv2@munimtech.com). Licence: Apache-2.0.
Pre-registered in findings/2026-11-28-prereg40-post2020-plan.md, written before any row of the
input table with both events after 2020-04-08 was looked at for this question.

Input: results/h13-v1/h13_pairs_scored.csv at commit 4c44a95 (pre-registration 13, 2026-11-19),
checked against its git blob SHA-1. It holds YKA late-coda differences (D_c310) and a mantle-P
source score S for repeat pairs of 83 South Sandwich Islands families found by our own clustering.

P (post-2020): YKA pairs with both events after the YKA firmware fix of 2020-04-08 and neither in
the 2021-12-01..2022-04-01 timing-error era. Every one-turn clock with its turn before 2020 maps
these times linearly, so the rank correlation of D_c310 with time apart is the same for all of them.
C (calibration): YKA pairs with both events in 1990.0..2008.5 (forward motion, agreed by both
rotation camps) and at most 2.65 years apart (6.6 yr / 2.5: the clock span of P under Wang et al.'s
2.5x slower return). Both sets: source-verified (S <= 1/3 at >= 2 stations), M03-24 excluded.
Standard library only.
"""
import argparse
import calendar
import csv
import hashlib
import io
import json
import math
import os
import random
import sys
import time
import urllib.request

VERSION = "0.1.0"
COMMIT = "4c44a952299eb29882930de2373c1aa16133e0ef"
IN_PATH = "results/h13-v1/h13_pairs_scored.csv"
IN_BLOB = "3bc1eacff2ce851426e9112083f98ae58afd9722"
RAW = "https://raw.githubusercontent.com/munimtechnologies/inner-core-repeaters/"
IN_URL = RAW + COMMIT + "/" + IN_PATH
PREREG_URL = RAW + "main/findings/2026-11-28-prereg40-post2020-plan.md"
UA = {"User-Agent": "inner-core-repeaters (MV2, Munim, Inc.)"}
NAN = float("nan")


def decimal_year(iso):
    """Project convention: year + (day_of_year - 1 + fraction of day) / days in year."""
    s = iso.strip().rstrip("Z")
    date, hms = s.split("T")
    y, mo, d = (int(v) for v in date.split("-"))
    h, mi, sec = hms.split(":")
    t0 = calendar.timegm((y, 1, 1, 0, 0, 0, 0, 0, 0))
    t = calendar.timegm((y, mo, d, int(h), int(mi), 0, 0, 0, 0)) + float(sec)
    days = 366 if calendar.isleap(y) else 365
    return y + (t - t0) / 86400.0 / days


POST_START = decimal_year("2020-04-08T00:00:00Z")      # YKA digitiser firmware fix (Ackerley 2026)
ERA2 = (decimal_year("2021-12-01T00:00:00Z"), decimal_year("2022-04-01T00:00:00Z"))
CAL_START = 1990.0
TURN = 2008.5
CAL_MAX_DT = 2.65
S_MAX = 1.0 / 3.0
MIN_ST = 2
MIN_N = 20
MIN_FAM = 4
LONG_DT = 2.0
MIN_LONG = 5
NBOOT = 2000
NPERM = 2000
F_POS = 0.95
BINS = ((0.0, 1.0), (1.0, 2.0), (2.0, 4.0), (4.0, 7.0))


def blob_sha(body):
    h = hashlib.sha1()
    h.update(b"blob %d\x00" % len(body))
    h.update(body)
    return h.hexdigest()


def fnum(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return NAN


def isf(x):
    return isinstance(x, float) and math.isfinite(x)


def truthy(x):
    return x is True or str(x).strip().lower() in ("true", "1", "yes")


def ranks(v):
    """Average ranks, 1-based."""
    idx = sorted(range(len(v)), key=lambda i: v[i])
    r = [0.0] * len(v)
    i = 0
    while i < len(idx):
        j = i
        while j + 1 < len(idx) and v[idx[j + 1]] == v[idx[i]]:
            j += 1
        avg = (i + j) / 2.0 + 1.0
        for k in range(i, j + 1):
            r[idx[k]] = avg
        i = j + 1
    return r


def pearson(x, y):
    n = len(x)
    if n < 2:
        return NAN
    mx, my = sum(x) / n, sum(y) / n
    sxx = sum((a - mx) ** 2 for a in x)
    syy = sum((b - my) ** 2 for b in y)
    sxy = sum((a - mx) * (b - my) for a, b in zip(x, y))
    if sxx <= 0 or syy <= 0:
        return NAN
    return sxy / math.sqrt(sxx * syy)


def spearman(x, y):
    if len(x) < 3:
        return NAN
    return pearson(ranks(x), ranks(y))


def nearest_rank(vals, q):
    s = sorted(vals)
    if not s:
        return NAN
    k = max(0, int(math.ceil(q * len(s))) - 1)
    return s[min(k, len(s) - 1)]


def median(vals):
    s = sorted(vals)
    n = len(s)
    if n == 0:
        return NAN
    return s[n // 2] if n % 2 else 0.5 * (s[n // 2 - 1] + s[n // 2])


def in_era2(t):
    return ERA2[0] <= t <= ERA2[1]


def select(rows, which, metric="D_c310", include_m0324=False):
    """Rows of the h13 table -> list of pair dicts for set 'post' or 'cal'."""
    out = []
    for r in rows:
        if str(r.get("array", "")).strip() != "YKA":
            continue
        if truthy(r.get("m0324")) and not include_m0324:
            continue
        s, nst = fnum(r.get("S")), fnum(r.get("n_st"))
        if not (isf(s) and s <= S_MAX and isf(nst) and nst >= MIN_ST):
            continue
        t1, t2, v = fnum(r.get("t1")), fnum(r.get("t2")), fnum(r.get(metric))
        if not (isf(t1) and isf(t2) and isf(v)):
            continue
        lo, hi, dt = min(t1, t2), max(t1, t2), abs(t2 - t1)
        if which == "post":
            if lo < POST_START or in_era2(t1) or in_era2(t2):
                continue
        elif which == "cal":
            if lo < CAL_START or hi >= TURN or dt > CAL_MAX_DT:
                continue
        else:
            raise ValueError(which)
        out.append(dict(mult=str(r.get("mult")), ev1=str(r.get("ev1")), ev2=str(r.get("ev2")),
                        t1=t1, t2=t2, dt=dt, v=v))
    return out


def boot(pairs, nboot=NBOOT, seed=0):
    """Whole-family bootstrap of Spearman(dt, v). Returns the list of rho (nan kept)."""
    fams = sorted(set(p["mult"] for p in pairs))
    by = {f: [p for p in pairs if p["mult"] == f] for f in fams}
    rng = random.Random(seed)
    out = []
    for _ in range(nboot):
        xs, ys = [], []
        for _k in range(len(fams)):
            for p in by[fams[rng.randrange(len(fams))]]:
                xs.append(p["dt"])
                ys.append(p["v"])
        out.append(spearman(xs, ys))
    return out


def summarize(pairs, seed=0, nboot=NBOOT):
    n = len(pairs)
    fams = set(p["mult"] for p in pairs)
    res = dict(n=n, n_fam=len(fams), n_long=sum(1 for p in pairs if p["dt"] >= LONG_DT))
    res["rho"] = spearman([p["dt"] for p in pairs], [p["v"] for p in pairs])
    if n >= 3 and len(fams) >= 1:
        b = boot(pairs, nboot, seed)
        fin = [x for x in b if isf(x)]
        res["f_pos"] = sum(1 for x in fin if x > 0) / float(nboot)
        res["boot_p05"] = nearest_rank(fin, 0.05)
        res["boot_p95"] = nearest_rank(fin, 0.95)
        res["n_boot_finite"] = len(fin)
    else:
        res.update(f_pos=NAN, boot_p05=NAN, boot_p95=NAN, n_boot_finite=0)
    res["bins"] = [dict(lo=lo, hi=hi, n=len([p for p in pairs if lo <= p["dt"] < hi]),
                        median=median([p["v"] for p in pairs if lo <= p["dt"] < hi]))
                   for lo, hi in BINS]
    return res


def testable_post(P):
    return P["n"] >= MIN_N and P["n_fam"] >= MIN_FAM and P["n_long"] >= MIN_LONG


def calibrated(C):
    return (C["n"] >= MIN_N and C["n_fam"] >= MIN_FAM and isf(C["rho"]) and C["rho"] > 0
            and isf(C["f_pos"]) and C["f_pos"] >= F_POS)


def headline(P, C):
    if not testable_post(P) or not isf(P["rho"]):
        return "untestable"
    if P["rho"] > 0 and isf(P["f_pos"]) and P["f_pos"] >= F_POS:
        return "still-changing"
    if not calibrated(C):
        return "insensitive"
    if isf(P["boot_p95"]) and P["boot_p95"] < C["rho"]:
        return "quiet-since-2020"
    return "unresolved"


def perm_p(pairs, nperm=NPERM, seed=0):
    """Report-only: permute event times within each family; one-sided p for rho >= observed."""
    vals = [p["v"] for p in pairs]
    obs = spearman([p["dt"] for p in pairs], vals)
    if not isf(obs):
        return NAN
    fam = {}
    for p in pairs:
        e = fam.setdefault(p["mult"], {})
        e.setdefault(p["ev1"], p["t1"])
        e.setdefault(p["ev2"], p["t2"])
    rng = random.Random(seed)
    ge = 0
    for _ in range(nperm):
        nt = {}
        for m in sorted(fam):
            ids = sorted(fam[m])
            ts = [fam[m][i] for i in ids]
            rng.shuffle(ts)
            for i, t in zip(ids, ts):
                nt[(m, i)] = t
        d = [abs(nt[(p["mult"], p["ev1"])] - nt[(p["mult"], p["ev2"])]) for p in pairs]
        r = spearman(d, vals)
        if isf(r) and r >= obs:
            ge += 1
    return (1 + ge) / (1.0 + nperm)


def grade_all(rows, metric="D_c310", include_m0324=False, nboot=NBOOT):
    P = summarize(select(rows, "post", metric, include_m0324), 0, nboot)
    C = summarize(select(rows, "cal", metric, include_m0324), 1, nboot)
    return dict(metric=metric, include_m0324=include_m0324, post=P, cal=C,
                testable_post=testable_post(P), calibrated=calibrated(C), headline=headline(P, C))


def clean(o):
    if isinstance(o, float):
        return o if math.isfinite(o) else None
    if isinstance(o, dict):
        return {k: clean(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [clean(v) for v in o]
    return o


def fetch(url, tries=4):
    for k in range(tries):
        try:
            req = urllib.request.Request(url, headers=UA)
            with urllib.request.urlopen(req, timeout=60) as r:
                if r.status == 200:
                    return r.read()
        except Exception as ex:  # noqa: BLE001
            print("fetch error", url, repr(ex)[:200], flush=True)
        time.sleep(5 * (k + 1))
    return None


def wait_prereg(max_min):
    """The plan must be public on main before grading; poll every 2 minutes."""
    t0 = time.time()
    while True:
        if fetch(PREREG_URL, tries=1) is not None:
            return True, round((time.time() - t0) / 60.0, 1)
        if time.time() - t0 > max_min * 60:
            return False, round((time.time() - t0) / 60.0, 1)
        time.sleep(120)


def write_pairs(path, pairs):
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["mult", "ev1", "ev2", "t1", "t2", "dt", "v"])
        w.writeheader()
        w.writerows(pairs)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--out-dir", default="/out")
    ap.add_argument("--wait-min", type=float, default=40.0)
    a = ap.parse_args(argv)
    os.makedirs(a.out_dir, exist_ok=True)
    s = dict(version=VERSION, input=dict(commit=COMMIT, path=IN_PATH, blob_expected=IN_BLOB),
             started=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()))
    pub, waited = wait_prereg(a.wait_min)
    s["prereg_public_at_start"] = pub
    s["prereg_wait_min"] = waited
    body = fetch(IN_URL)
    got = blob_sha(body) if body is not None else None
    s["input"]["blob_got"] = got
    s["checks_ok"] = bool(pub and got == IN_BLOB)
    if not s["checks_ok"]:
        s["result"] = dict(headline="not graded",
                           reason="plan not public" if not pub else "input missing or blob mismatch")
    else:
        rows = list(csv.DictReader(io.StringIO(body.decode("utf-8"))))
        s["n_rows"] = len(rows)
        s["n_rows_yka"] = sum(1 for r in rows if str(r.get("array", "")).strip() == "YKA")
        s["result"] = grade_all(rows)
        post = select(rows, "post")
        s["report_only"] = dict(
            metric_D=grade_all(rows, metric="D"),
            with_m0324=grade_all(rows, include_m0324=True),
            perm_p_post=perm_p(post),
            perm_p_cal=perm_p(select(rows, "cal")))
        write_pairs(os.path.join(a.out_dir, "gradeq40_post_pairs.csv"), post)
        write_pairs(os.path.join(a.out_dir, "gradeq40_cal_pairs.csv"), select(rows, "cal"))
    s["finished"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    with open(os.path.join(a.out_dir, "gradeq40_summary.json"), "w") as f:
        json.dump(clean(s), f, indent=1)
    res = s.get("result", {})
    print(json.dumps(clean(dict(headline=res.get("headline"), post=res.get("post"), cal=res.get("cal"))),
                     indent=1)[:4000], flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())

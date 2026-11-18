#!/usr/bin/env python3
"""icrdiag.py v0.1.0 -- exploratory diagnostics after icr-ykab-v2 and icr-bcdf-v2.

Part Y: why did YKA stage B's pre-registered calibration (P4-Y) fail in icr-ykab-v2?
        Pooled vs within-family calibration (Simpson's paradox), the 2021 M8.1 sequence,
        direct-pulse SNR, epicentral separation, magnitude difference.
Part B: are the two branches of the ILAR BC-DF residuals (about 0.6 s apart) cycle skips
        or source geography?

Everything here is EXPLORATORY (findings/2026-11-18-diag2-plan.md). It changes no
pre-registered verdict. Stdlib only; matplotlib is optional (plots).
MV2, Munim, Inc. (made with AI agents). Apache-2.0.
"""
import argparse
import bisect
import csv
import datetime
import io
import json
import math
import os
import random
import time
import urllib.request

VERSION = "0.1.0"
RAW = "https://raw.githubusercontent.com/munimtechnologies/inner-core-repeaters/main/"
USGS = ("https://earthquake.usgs.gov/fdsnws/event/1/query?format=csv&starttime=1990-01-01"
        "&endtime=2026-11-19&minlatitude=-61.5&maxlatitude=-54&minlongitude=-31&maxlongitude=-22"
        "&minmagnitude=4.3&orderby=time-asc&limit=20000")
TURN, RATIO = 2008.5, 2.5
SEQ21 = (2021.6, 2023.0)  # 2021-08-12 M8.1 South Sandwich Islands sequence and aftershocks
NAN = float("nan")


# ----------------------------------------------------------------------------- basics
def equiv(t):
    return t if t <= TURN else TURN - (t - TURN) / RATIO


def delta_M(t1, t2):
    return abs(t2 - t1)


def delta_B(t1, t2):
    return abs(equiv(t2) - equiv(t1))


def fnum(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return NAN


def isf(x):
    try:
        return math.isfinite(x)
    except TypeError:
        return False


def median(xs):
    xs = sorted(x for x in xs if isf(x))
    n = len(xs)
    if n == 0:
        return NAN
    m = n // 2
    return xs[m] if n % 2 else 0.5 * (xs[m - 1] + xs[m])


def mean(xs):
    xs = [x for x in xs if isf(x)]
    return sum(xs) / len(xs) if xs else NAN


def ranks(xs):
    idx = sorted(range(len(xs)), key=lambda i: xs[i])
    r = [0.0] * len(xs)
    i = 0
    while i < len(idx):
        j = i
        while j + 1 < len(idx) and xs[idx[j + 1]] == xs[idx[i]]:
            j += 1
        avg = 0.5 * (i + j) + 1.0
        for k in range(i, j + 1):
            r[idx[k]] = avg
        i = j + 1
    return r


def pearson(x, y):
    n = len(x)
    if n < 3:
        return NAN
    mx, my = sum(x) / n, sum(y) / n
    sxx = sum((a - mx) ** 2 for a in x)
    syy = sum((b - my) ** 2 for b in y)
    if sxx <= 0 or syy <= 0:
        return NAN
    return sum((a - mx) * (b - my) for a, b in zip(x, y)) / math.sqrt(sxx * syy)


def spearman(x, y):
    pr = [(a, b) for a, b in zip(x, y) if isf(a) and isf(b)]
    if len(pr) < 3:
        return NAN, len(pr)
    xs, ys = zip(*pr)
    return pearson(ranks(list(xs)), ranks(list(ys))), len(pr)


def mw_less(a, b):
    """One-sided Mann-Whitney p for 'a tends to be smaller than b' (normal approx.)."""
    n1, n2 = len(a), len(b)
    if n1 < 3 or n2 < 3:
        return NAN
    r = ranks(list(a) + list(b))
    u1 = sum(r[:n1]) - n1 * (n1 + 1) / 2.0
    mu = n1 * n2 / 2.0
    sd = math.sqrt(n1 * n2 * (n1 + n2 + 1) / 12.0)
    z = (u1 - mu + 0.5) / sd
    return 0.5 * math.erfc(-z / math.sqrt(2.0))


def binom_tail(k, n):
    """P(X >= k), X ~ Bin(n, 1/2)."""
    if n == 0:
        return 1.0
    return sum(math.comb(n, j) for j in range(k, n + 1)) / 2.0 ** n


def sep_km(lat1, lon1, lat2, lon2):
    p1, p2 = math.radians(lat1), math.radians(lat2)
    dp, dl = p2 - p1, math.radians(lon2 - lon1)
    h = math.sin(dp / 2) ** 2 + math.cos(p1) * math.cos(p2) * math.sin(dl / 2) ** 2
    return 2 * 6371.0 * math.asin(min(1.0, math.sqrt(h)))


def solve(A, b):
    n = len(b)
    M = [list(A[i]) + [b[i]] for i in range(n)]
    for c in range(n):
        p = max(range(c, n), key=lambda i: abs(M[i][c]))
        if abs(M[p][c]) < 1e-12:
            raise ValueError("singular")
        M[c], M[p] = M[p], M[c]
        for i in range(n):
            if i != c:
                f = M[i][c] / M[c][c]
                for j in range(c, n + 1):
                    M[i][j] -= f * M[c][j]
    return [M[i][n] / M[i][i] for i in range(n)]


def ols(X, y):
    p = len(X[0])
    A = [[sum(row[i] * row[j] for row in X) for j in range(p)] for i in range(p)]
    b = [sum(row[i] * v for row, v in zip(X, y)) for i in range(p)]
    return solve(A, b)


def iso_to_epoch(s):
    s = s.strip().rstrip("Z")
    fmt = "%Y-%m-%dT%H:%M:%S.%f" if "." in s else "%Y-%m-%dT%H:%M:%S"
    return datetime.datetime.strptime(s, fmt).replace(tzinfo=datetime.timezone.utc).timestamp()


def clean(o):
    if isinstance(o, float):
        return o if math.isfinite(o) else None
    if isinstance(o, dict):
        return {k: clean(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [clean(v) for v in o]
    return o


# ----------------------------------------------------------------------------- I/O
def fetch(url, tries=10, wait=60, timeout=90):
    last = None
    for k in range(tries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "mv2-inner-core-repeaters/" + VERSION})
            with urllib.request.urlopen(req, timeout=timeout) as r:
                return r.read().decode("utf-8")
        except Exception as ex:  # noqa: BLE001
            last = ex
            if k + 1 < tries:
                time.sleep(wait)
    raise RuntimeError("fetch failed: %s (%r)" % (url, last))


def get_text(local, url, **kw):
    if local and os.path.exists(local):
        with open(local, encoding="utf-8") as f:
            return f.read()
    return fetch(url, **kw)


def read_catalog(text):
    by_id, items = {}, []
    for r in csv.DictReader(io.StringIO(text)):
        try:
            e = dict(lat=float(r["latitude"]), lon=float(r["longitude"]), depth=fnum(r.get("depth")),
                     epoch=iso_to_epoch(r["time"]))
        except (KeyError, ValueError):
            continue
        by_id[r.get("id", "")] = e
        items.append(e)
    items.sort(key=lambda e: e["epoch"])
    return dict(by_id=by_id, items=items, keys=[e["epoch"] for e in items])


def locate(ev_id, iso, cat, tol=5.0):
    if ev_id in cat["by_id"]:
        return cat["by_id"][ev_id], "id"
    if not iso:
        return None, "none"
    try:
        t = iso_to_epoch(iso)
    except ValueError:
        return None, "none"
    i = bisect.bisect_left(cat["keys"], t)
    best = None
    for j in (i - 1, i):
        if 0 <= j < len(cat["items"]):
            d = abs(cat["keys"][j] - t)
            if d <= tol and (best is None or d < best[0]):
                best = (d, cat["items"][j])
    return (best[1], "time") if best else (None, "none")


def read_pairs(text):
    rows = []
    for r in csv.DictReader(io.StringIO(text)):
        x = dict(r)
        for k in ("t1", "t2", "mag1", "mag2", "snr_dir1", "snr_dir2", "D", "D_c13", "D_c310", "D_full"):
            x[k] = fnum(r.get(k))
        x["ok"] = str(r.get("ok", "")).strip() == "True"
        rows.append(x)
    return rows


def enrich(rows, cat):
    n_loc = 0
    for r in rows:
        e1, m1 = locate(r["ev1"], r.get("time1"), cat)
        e2, m2 = locate(r["ev2"], r.get("time2"), cat)
        if e1 and e2:
            n_loc += 1
            r["sep_km"] = sep_km(e1["lat"], e1["lon"], e2["lat"], e2["lon"])
            r["ddep_km"] = abs(e1["depth"] - e2["depth"]) if isf(e1["depth"]) and isf(e2["depth"]) else NAN
        else:
            r["sep_km"] = r["ddep_km"] = NAN
        r["loc"] = m1 + "/" + m2
        r["dt"] = delta_M(r["t1"], r["t2"])
        r["dB"] = delta_B(r["t1"], r["t2"])
        r["minsnr"] = min(r["snr_dir1"], r["snr_dir2"])
        r["dmag"] = abs(r["mag1"] - r["mag2"])
        r["seq21"] = any(SEQ21[0] <= t < SEQ21[1] for t in (r["t1"], r["t2"]))
        r["same_gen"] = r.get("gen1") == r.get("gen2")
    return n_loc


# ----------------------------------------------------------------------------- Part Y
def is_near(r):
    return r["dt"] < 2.0


def is_far(r):
    return r["dt"] > 10.0 and r["dB"] > 5.0


def p4x(rows, key="D", nboot=500, seed=1):
    rows = [r for r in rows if isf(r[key])]
    a = [r[key] for r in rows if is_near(r)]
    b = [r[key] for r in rows if is_far(r)]
    out = dict(n_near=len(a), n_far=len(b), median_near=median(a), median_far=median(b),
               p_mw_near_lower=mw_less(a, b),
               holds=(median(a) < median(b)) if len(a) >= 3 and len(b) >= 3 else None)
    by = {}
    for r in rows:
        by.setdefault(r["mult"], []).append(r)
    fams = sorted(by)
    rng = random.Random(seed)
    wins = tot = 0
    for _ in range(nboot if fams else 0):
        aa, bb = [], []
        for f in (rng.choice(fams) for _ in fams):
            for r in by[f]:
                if is_near(r):
                    aa.append(r[key])
                elif is_far(r):
                    bb.append(r[key])
        if len(aa) >= 3 and len(bb) >= 3:
            tot += 1
            wins += median(aa) < median(bb)
    out["boot_frac_holds"] = wins / tot if tot else NAN
    out["n_families_near"] = len({r["mult"] for r in rows if is_near(r)})
    out["n_families_far"] = len({r["mult"] for r in rows if is_far(r)})
    return out


def within_family(rows, key="D", min_each=2):
    by = {}
    for r in rows:
        if isf(r[key]):
            by.setdefault(r["mult"], []).append(r)
    fam = []
    for f, rr in sorted(by.items()):
        a = [r[key] for r in rr if is_near(r)]
        b = [r[key] for r in rr if is_far(r)]
        if len(a) >= min_each and len(b) >= min_each:
            fam.append(dict(mult=f, n_near=len(a), n_far=len(b), med_near=median(a), med_far=median(b),
                            holds=median(a) < median(b)))
    k = sum(1 for x in fam if x["holds"])
    return dict(n_families=len(fam), n_holds=k, p_sign=binom_tail(k, len(fam)), families=fam)


def demeaned(rows, key="D"):
    by = {}
    for r in rows:
        if isf(r[key]):
            by.setdefault(r["mult"], []).append(r[key])
    med = {f: median(v) for f, v in by.items()}
    out = []
    for r in rows:
        if isf(r[key]):
            x = dict(r)
            x[key + "_dm"] = r[key] - med[r["mult"]]
            out.append(x)
    return out


def composition(rows):
    out = {}
    for name, sel in (("near", is_near), ("far", is_far)):
        rr = [r for r in rows if sel(r)]
        out[name] = dict(n=len(rr), n_families=len({r["mult"] for r in rr}),
                         frac_seq21=mean([1.0 if r["seq21"] else 0.0 for r in rr]),
                         frac_same_gen=mean([1.0 if r["same_gen"] else 0.0 for r in rr]),
                         median_minsnr=median([r["minsnr"] for r in rr]),
                         median_sep_km=median([r["sep_km"] for r in rr]),
                         median_ddep_km=median([r["ddep_km"] for r in rr]),
                         median_dmag=median([r["dmag"] for r in rr]))
    return out


PRED = ("dt", "sep_km", "ddep_km", "minsnr", "dmag", "seq21")


def rank_regression(rows, key="D", nboot=300, seed=3):
    rr = [r for r in rows if isf(r[key]) and all(isf(float(r[p])) for p in PRED)]
    if len(rr) < 30:
        return dict(n=len(rr), testable=False)

    def design(sub):
        n = len(sub)
        cols = [ranks([float(r[p]) for r in sub]) for p in PRED]
        y = ranks([r[key] for r in sub])
        X = [[1.0] + [cols[j][i] / n for j in range(len(PRED))] for i in range(n)]
        return X, [v / n for v in y]

    X, y = design(rr)
    try:
        coef = ols(X, y)
    except ValueError:
        return dict(n=len(rr), testable=False)
    by = {}
    for r in rr:
        by.setdefault(r["mult"], []).append(r)
    fams = sorted(by)
    rng = random.Random(seed)
    boots = []
    for _ in range(nboot):
        sub = [r for f in (rng.choice(fams) for _ in fams) for r in by[f]]
        try:
            Xb, yb = design(sub)
            boots.append(ols(Xb, yb))
        except ValueError:
            continue
    res = dict(n=len(rr), testable=True, n_boot=len(boots), coef={})
    for j, p in enumerate(PRED):
        bs = sorted(b[j + 1] for b in boots)
        lo = bs[int(0.05 * (len(bs) - 1))] if bs else NAN
        hi = bs[int(0.95 * (len(bs) - 1))] if bs else NAN
        res["coef"][p] = dict(b=coef[j + 1], p05=lo, p95=hi)
    return res


def part_y(rows, nboot=500):
    ok = [r for r in rows if r["ok"]]
    out = dict(n_rows=len(rows), n_ok=len(ok))
    out["p4_pooled_as_graded"] = p4x(ok, "D", nboot)
    out["composition"] = composition(ok)
    out["spearman_D"] = {p: dict(zip(("rho", "n"), spearman([float(r[p]) for r in ok], [r["D"] for r in ok])))
                         for p in PRED}
    near = [r for r in ok if is_near(r)]
    out["spearman_D_near_only"] = {p: dict(zip(("rho", "n"), spearman([float(r[p]) for r in near],
                                                                      [r["D"] for r in near])))
                                   for p in ("sep_km", "ddep_km", "minsnr", "dmag", "seq21")}
    out["rank_regression_D"] = rank_regression(ok, "D", nboot=min(nboot, 300))
    snr_med = median([r["minsnr"] for r in ok])
    sep_med = median([r["sep_km"] for r in ok])
    strata = {
        "all": ok,
        "no_seq21": [r for r in ok if not r["seq21"]],
        "hi_snr": [r for r in ok if r["minsnr"] >= snr_med],
        "close": [r for r in ok if isf(r["sep_km"]) and r["sep_km"] <= sep_med],
        "hi_snr_close": [r for r in ok if r["minsnr"] >= snr_med and isf(r["sep_km"]) and r["sep_km"] <= sep_med],
        "same_gen": [r for r in ok if r["same_gen"]],
        "no_seq21_hi_snr": [r for r in ok if not r["seq21"] and r["minsnr"] >= snr_med],
    }
    out["strata_cut"] = dict(minsnr_median=snr_med, sep_km_median=sep_med)
    out["p4_strata"] = {}
    for name, rr in strata.items():
        out["p4_strata"][name] = {k: p4x(rr, k, nboot) for k in ("D", "D_c13", "D_c310")}
    out["within_family"] = {k: within_family(ok, k) for k in ("D", "D_c13", "D_c310")}
    out["within_family_no_seq21"] = {k: within_family([r for r in ok if not r["seq21"]], k)
                                     for k in ("D", "D_c13")}
    out["p4_family_demeaned"] = {}
    for k in ("D", "D_c13"):
        dm = demeaned(ok, k)
        out["p4_family_demeaned"][k] = p4x(dm, k + "_dm", nboot)
        out["p4_family_demeaned"][k + "_no_seq21"] = p4x([r for r in dm if not r["seq21"]], k + "_dm", nboot)
    return out


# ----------------------------------------------------------------------------- Part B
def em2(x, iters=300):
    n = len(x)
    xs = sorted(x)
    m = mean(x)
    s = math.sqrt(sum((v - m) ** 2 for v in x) / n) or 0.1
    mu = [xs[n // 4], xs[(3 * n) // 4]]
    if mu[0] == mu[1]:
        mu[1] = mu[0] + s
    sd = [max(0.02, s / 2)] * 2
    w = [0.5, 0.5]
    g = [0.5] * n
    for _ in range(iters):
        g = []
        for v in x:
            p = [w[k] * math.exp(-0.5 * ((v - mu[k]) / sd[k]) ** 2) / sd[k] for k in (0, 1)]
            tot = p[0] + p[1]
            g.append(p[1] / tot if tot > 0 else (1.0 if abs(v - mu[1]) < abs(v - mu[0]) else 0.0))
        n1 = sum(g)
        n0 = n - n1
        if n1 < 1 or n0 < 1:
            break
        mu = [sum((1 - gi) * v for gi, v in zip(g, x)) / n0, sum(gi * v for gi, v in zip(g, x)) / n1]
        sd = [max(0.02, math.sqrt(sum((1 - gi) * (v - mu[0]) ** 2 for gi, v in zip(g, x)) / n0)),
              max(0.02, math.sqrt(sum(gi * (v - mu[1]) ** 2 for gi, v in zip(g, x)) / n1))]
        w = [n0 / n, n1 / n]
    if mu[0] > mu[1]:
        mu, sd, w, g = mu[::-1], sd[::-1], w[::-1], [1 - gi for gi in g]

    def ll2():
        t = 0.0
        for v in x:
            p = sum(w[k] * math.exp(-0.5 * ((v - mu[k]) / sd[k]) ** 2) / (sd[k] * math.sqrt(2 * math.pi))
                    for k in (0, 1))
            t += math.log(max(p, 1e-300))
        return t

    ll1 = sum(-0.5 * ((v - m) / s) ** 2 - math.log(s * math.sqrt(2 * math.pi)) for v in x)
    return dict(mu=mu, sd=sd, w=w, post=g, gap=mu[1] - mu[0],
                bic1=-2 * ll1 + 2 * math.log(n), bic2=-2 * ll2() + 5 * math.log(n))


def nn_agreement(pts, labels, k=1, nperm=1000, seed=4):
    n = len(pts)
    nn = []
    for i in range(n):
        d = sorted((sep_km(pts[i][0], pts[i][1], pts[j][0], pts[j][1]), j) for j in range(n) if j != i)
        nn.append([j for _, j in d[:k]])

    def agree(lab):
        return sum(sum(1 for j in nn[i] if lab[j] == lab[i]) / float(len(nn[i])) for i in range(n)) / n

    obs = agree(labels)
    rng = random.Random(seed)
    lab = list(labels)
    null = []
    for _ in range(nperm):
        rng.shuffle(lab)
        null.append(agree(lab))
    null.sort()
    return dict(k=k, observed=obs, null_median=median(null), null_p95=null[int(0.95 * (len(null) - 1))],
                p=(1 + sum(1 for v in null if v >= obs)) / (nperm + 1.0),
                median_nn_km=median([sep_km(pts[i][0], pts[i][1], pts[nn[i][0]][0], pts[nn[i][0]][1])
                                     for i in range(n)]))


def slope(ts, ys):
    n = len(ts)
    if n < 5:
        return dict(n=n)
    mt, my = sum(ts) / n, sum(ys) / n
    sxx = sum((t - mt) ** 2 for t in ts)
    b = sum((t - mt) * (y - my) for t, y in zip(ts, ys)) / sxx
    res = [y - my - b * (t - mt) for t, y in zip(ts, ys)]
    sig = math.sqrt(sum(e * e for e in res) / (n - 2))
    return dict(n=n, b_s_per_yr=b, se=sig / math.sqrt(sxx), sigma=sig)


def part_b(text, nperm=1000):
    ev = []
    for r in csv.DictReader(io.StringIO(text)):
        if str(r.get("ok", "")).strip() != "True":
            continue
        e = dict(id=r["id"], t=fnum(r["t"]), lat=fnum(r["lat"]), lon=fnum(r["lon"]), depth=fnum(r["depth"]),
                 mag=fnum(r["mag"]), res=fnum(r["res"]), res_el=fnum(r.get("res_el")), onset=fnum(r.get("onset_s")),
                 cc=fnum(r.get("cc")), snr=fnum(r.get("snr_df")))
        if all(isf(e[k]) for k in ("t", "lat", "lon", "depth", "res")):
            ev.append(e)
    if len(ev) < 20:
        return dict(n=len(ev), testable=False)
    lat0 = mean([e["lat"] for e in ev])
    lon0 = mean([e["lon"] for e in ev])
    dep0 = mean([e["depth"] for e in ev])
    cl = math.cos(math.radians(lat0))

    def row(e):
        return [1.0, e["lat"] - lat0, (e["lon"] - lon0) * cl, (e["depth"] - dep0) / 100.0]

    fit = [e for e in ev if e["t"] < 2025.0]
    coef = ols([row(e) for e in fit], [e["res"] for e in fit])
    for e in ev:
        pr = sum(c * v for c, v in zip(coef, row(e)))
        e["rc"] = e["res"] - pr
        e["rc_el"] = e["res_el"] - pr if isf(e["res_el"]) else NAN
    mix = em2([e["rc"] for e in ev])
    for e, gi in zip(ev, mix["post"]):
        e["branch"] = 1 if gi > 0.5 else 0
    mid = 0.5 * (mix["mu"][0] + mix["mu"][1])
    out = dict(n=len(ev), testable=True, spatial_coef=coef,
               mixture=dict(mu=mix["mu"], sd=mix["sd"], w=mix["w"], gap_s=mix["gap"],
                            dbic_2_minus_1=mix["bic2"] - mix["bic1"]),
               n_upper=sum(e["branch"] for e in ev))
    pts = [(e["lat"], e["lon"]) for e in ev]
    lab = [e["branch"] for e in ev]
    out["nn1"] = nn_agreement(pts, lab, 1, nperm)
    out["nn5"] = nn_agreement(pts, lab, 5, nperm)
    up = [e for e in ev if e["branch"] == 1]
    lo = [e for e in ev if e["branch"] == 0]
    for k in ("lat", "lon", "depth", "mag", "onset", "cc", "snr", "t"):
        u = [e[k] for e in up if isf(e[k])]
        v = [e[k] for e in lo if isf(e[k])]
        out["median_" + k] = dict(upper=median(u), lower=median(v),
                                  p_upper_lower=mw_less(u, v), p_upper_higher=mw_less(v, u))
    el = [e for e in ev if isf(e["rc_el"])]
    out["element_check"] = dict(n=len(el),
                                frac_abs_diff_gt_0p3=mean([1.0 if abs(e["rc_el"] - e["rc"]) > 0.3 else 0.0 for e in el]),
                                frac_same_branch=mean([1.0 if (e["rc_el"] > mid) == (e["branch"] == 1) else 0.0
                                                       for e in el]))
    periods = (("pre_2008.5", -1e9, 2008.5), ("2008.5-2021.6", 2008.5, 2021.6), ("2021.6-2025", 2021.6, 2025.0),
               ("2025+", 2025.0, 1e9))
    out["frac_upper_by_period"] = {name: dict(n=len([e for e in ev if a <= e["t"] < b]),
                                              frac=mean([float(e["branch"]) for e in ev if a <= e["t"] < b]))
                                   for name, a, b in periods}
    out["slope_lower_pre2025"] = slope([e["t"] for e in lo if e["t"] < 2025], [e["rc"] for e in lo if e["t"] < 2025])
    out["slope_lower_all"] = slope([e["t"] for e in lo], [e["rc"] for e in lo])
    out["slope_upper_all"] = slope([e["t"] for e in up], [e["rc"] for e in up])
    out["slope_all"] = slope([e["t"] for e in ev], [e["rc"] for e in ev])
    out["lower_new_minus_pre2025_median"] = (median([e["rc"] for e in lo if e["t"] >= 2025])
                                             - median([e["rc"] for e in lo if e["t"] < 2025]))
    out["events"] = ev
    return out


# ----------------------------------------------------------------------------- plots
def plots(out_dir, ok, bev):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception:  # noqa: BLE001
        return []
    made = []
    os.makedirs(os.path.join(out_dir, "plots"), exist_ok=True)
    try:
        fig, ax = plt.subplots(1, 2, figsize=(12, 4.5))
        for flag, col, lbl in ((False, "tab:blue", "other"), (True, "tab:red", "2021-22 sequence")):
            rr = [r for r in ok if r["seq21"] == flag]
            ax[0].scatter([r["dt"] for r in rr], [r["D"] for r in rr], s=6, c=col, label=lbl, alpha=0.6)
            ax[1].scatter([r["minsnr"] for r in rr], [r["D"] for r in rr], s=6, c=col, label=lbl, alpha=0.6)
        ax[0].axvline(2.0, color="k", lw=0.6)
        ax[0].axvline(10.0, color="k", lw=0.6)
        ax[0].set_xlabel("|t2 - t1| (yr)")
        ax[1].set_xlabel("min direct SNR")
        ax[1].set_xscale("log")
        for a in ax:
            a.set_yscale("log")
            a.set_ylabel("YKA coda D (1-6 s)")
            a.legend(fontsize=8)
        fig.suptitle("YKA stage B v2 ok pairs (exploratory diagnostics)")
        p = os.path.join(out_dir, "plots", "ykadiag_D.png")
        fig.tight_layout()
        fig.savefig(p, dpi=110)
        plt.close(fig)
        made.append(p)
    except Exception:  # noqa: BLE001
        pass
    try:
        if bev:
            fig, ax = plt.subplots(1, 2, figsize=(12, 4.5))
            for b, col in ((0, "tab:blue"), (1, "tab:orange")):
                ee = [e for e in bev if e["branch"] == b]
                ax[0].scatter([e["lon"] for e in ee], [e["lat"] for e in ee], s=8, c=col, label="branch %d" % b)
                ax[1].scatter([e["t"] for e in ee], [e["rc"] for e in ee], s=8, c=col, label="branch %d" % b)
            ax[0].set_xlabel("lon")
            ax[0].set_ylabel("lat")
            ax[1].set_xlabel("year")
            ax[1].set_ylabel("BC-DF residual minus spatial terms (s)")
            for a in ax:
                a.legend(fontsize=8)
            fig.suptitle("ILAR BC-DF branches (exploratory)")
            p = os.path.join(out_dir, "plots", "bcdf_branches.png")
            fig.tight_layout()
            fig.savefig(p, dpi=110)
            plt.close(fig)
            made.append(p)
    except Exception:  # noqa: BLE001
        pass
    return made


# ----------------------------------------------------------------------------- main
def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-dir", default="out")
    ap.add_argument("--pairs", default="")
    ap.add_argument("--catalog", default="")
    ap.add_argument("--bcdf", default="")
    ap.add_argument("--nboot", type=int, default=500)
    ap.add_argument("--nperm", type=int, default=1000)
    ap.add_argument("--wait", type=int, default=60)
    a = ap.parse_args(argv)
    os.makedirs(a.out_dir, exist_ok=True)
    summary = dict(version=VERSION, exploratory=True)

    rows, ok = [], []
    try:
        rows = read_pairs(get_text(a.pairs, RAW + "results/ykab-v2/yka_pairs.csv", wait=a.wait))
        cat = read_catalog(get_text(a.catalog, USGS, wait=a.wait))
        summary["catalog_events"] = len(cat["items"])
        summary["pairs_located"] = enrich(rows, cat)
        summary["part_y"] = part_y(rows, a.nboot)
        ok = [r for r in rows if r["ok"]]
        cols = ["mult", "ev1", "ev2", "t1", "t2", "ok", "D", "D_c13", "D_c310", "dt", "dB", "sep_km", "ddep_km",
                "minsnr", "dmag", "seq21", "same_gen", "loc"]
        with open(os.path.join(a.out_dir, "ykadiag_pairs.csv"), "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(cols)
            for r in rows:
                w.writerow([r.get(c) for c in cols])
    except Exception as ex:  # noqa: BLE001
        summary["part_y_error"] = repr(ex)

    bev = []
    try:
        pb = part_b(get_text(a.bcdf, RAW + "results/bcdf-v2/bcdf_events.csv", wait=a.wait), a.nperm)
        bev = pb.pop("events", [])
        summary["part_b"] = pb
        cols = ["id", "t", "lat", "lon", "depth", "mag", "res", "rc", "rc_el", "onset", "cc", "snr", "branch"]
        with open(os.path.join(a.out_dir, "bcdf_branches.csv"), "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(cols)
            for e in bev:
                w.writerow([e.get(c) for c in cols])
    except Exception as ex:  # noqa: BLE001
        summary["part_b_error"] = repr(ex)

    summary["plots"] = [os.path.basename(p) for p in plots(a.out_dir, ok, bev)]
    with open(os.path.join(a.out_dir, "icrdiag_summary.json"), "w") as f:
        json.dump(clean(summary), f, indent=1)
    print(json.dumps(clean({k: v for k, v in summary.items() if k != "part_y"}), indent=1)[:4000])
    return summary


if __name__ == "__main__":
    main()

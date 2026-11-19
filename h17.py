#!/usr/bin/env python3
"""h17.py: is the time trend among source-verified YKA pairs a hardware artefact?

Part of inner-core-repeaters by MV2, Munim, Inc. (made with AI agents). Licence: Apache-2.0.
Pre-registered in findings/2026-11-19-prereg17-hardware.md (pre-registration 17), written and
published before this script ran.

Pre-registration 13 (mantlep3.py 0.15.0, job icr-h13-v1) found that among source-verified YKA
pairs (mantle-P source score S <= 1/3) the PKIKP coda difference D grows with time apart
(rho_M 0.158, rho_B 0.127, both family-bootstrap 5th percentiles > 0) while the two clocks could
not be told apart. The YKA records come from two generations (gen "legacy", up to at least
mid-2013, vs "new", from 2015, in the YKA pair table), and pairs far apart in time are more likely to straddle it.
This script repeats the test on pairs recorded by the same hardware generation.

No waveform is fetched. Inputs (public, in the repository):
  results/ykab-v2/yka_pairs.csv at commit f65c966a (the table H13 used; has gen1, gen2)
  results/h13-v1/h13_pairs_scored.csv, mantlep3_pairs.csv, mantlep3_families.csv (on main)
The statistics are copied unchanged from mantlep3.py 0.15.0 / mantlep2.py 0.13.0 and the
script first reproduces H13's published primary numbers.

Writes to --out-dir: h17_summary.json, h17_pairs.csv, h17_bins.csv, plots/h17_overview.png.
"""
import argparse
import csv
import io
import json
import math
import os
import random
import sys
import time
import urllib.request

import numpy as np

VERSION = "0.1.0"
RAW = "https://raw.githubusercontent.com/munimtechnologies/inner-core-repeaters/"
YKA_TABLE = ("f65c966adf2cce3a22ffd9e4629c842959121407", "results/ykab-v2/yka_pairs.csv")
H13_FILES = ("results/h13-v1/h13_pairs_scored.csv", "results/h13-v1/mantlep3_pairs.csv",
             "results/h13-v1/mantlep3_families.csv")
PREREG = "findings/2026-11-19-prereg17-hardware.md"

# H13 primary result as published (h13_summary.json, grade.primary.clock); must be reproduced.
EXPECT = dict(n_pairs=328, n_families=42, n_cross=96, rho_B=0.12651797149019148,
              rho_M=0.1575493725434452, frac_B_gt_M=0.27)

M0324 = frozenset(("usp000bt5x", "usp000fz4d", "usp000hute", "usc000t4wj", "us7000am0e",
                   "us6000h2wz", "us6000ittv", "us6000mk65", "us6000runq"))
TURN, RATIO = 2008.5, 2.5
MIN_STATION_PAIRS = 20
MIN_ST = 2
S_MAX = 1.0 / 3.0
RHO_MIN = 0.10
RHO_NONE_P95 = 0.20
FRAC = 0.90
MIN_V, MIN_FAM, MIN_CROSS = 40, 4, 10
NBOOT = 1000
LOCK_S = 0.3
BINS = [(0.0, 3.0), (3.0, 6.0), (6.0, 10.0), (10.0, 15.0), (15.0, 30.0)]
NAN = float("nan")


# ----------------------------------------------------------------------------- copied helpers
def isf(x):
    try:
        return math.isfinite(float(x))
    except (TypeError, ValueError):
        return False


def fnum(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return NAN


def equiv(t):
    return t if t <= TURN else TURN - (t - TURN) / RATIO


def d_B(t1, t2):
    return abs(equiv(t2) - equiv(t1))


def d_M(t1, t2):
    return abs(t2 - t1)


def is_cross(t1, t2):
    return (t1 < TURN) != (t2 < TURN)


def median(xs):
    xs = sorted(float(x) for x in xs if isf(x))
    n = len(xs)
    if n == 0:
        return NAN
    m = n // 2
    return xs[m] if n % 2 else 0.5 * (xs[m - 1] + xs[m])


def ncc(a, b):
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    if len(a) == 0 or len(a) != len(b):
        return NAN
    a = a - a.mean()
    b = b - b.mean()
    d = math.sqrt(float(np.dot(a, a)) * float(np.dot(b, b)))
    return float(np.dot(a, b) / d) if d > 0 else NAN


def rankdata(v):
    v = np.asarray(v, dtype=float)
    order = np.argsort(v, kind="mergesort")
    r = np.empty(len(v))
    r[order] = np.arange(len(v), dtype=float)
    sv = v[order]
    i = 0
    while i < len(sv):
        j = i
        while j + 1 < len(sv) and sv[j + 1] == sv[i]:
            j += 1
        if j > i:
            r[order[i:j + 1]] = (i + j) / 2.0
        i = j + 1
    return r


def spearman(a, b):
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    m = np.isfinite(a) & np.isfinite(b)
    if m.sum() < 5:
        return NAN
    return ncc(rankdata(a[m]), rankdata(b[m]))


def _pct(xs, q):
    xs = sorted(x for x in xs if isf(x))
    if not xs:
        return NAN
    return xs[int(round(q * (len(xs) - 1)))]


def family_boot(rows, stat, nboot=NBOOT, seed=13):
    by = {}
    for r in rows:
        by.setdefault(r["mult"], []).append(r)
    fams = sorted(by)
    rng = random.Random(seed)
    out = []
    for _ in range(nboot if len(fams) >= 2 else 0):
        sub = [r for f in (rng.choice(fams) for _ in fams) for r in by[f]]
        out.append(stat(sub))
    return out


def decide13(rB, rM, frac_B, frac_M, p95_B, p95_M):
    if not all(isf(x) for x in (rB, rM, frac_B, frac_M, p95_B, p95_M)):
        return "untestable"
    if rB >= RHO_MIN and frac_B >= FRAC:
        return "B"
    if rM >= RHO_MIN and frac_M >= FRAC:
        return "M"
    if max(rB, rM) < RHO_MIN and max(p95_B, p95_M) < RHO_NONE_P95:
        return "none"
    return "inconclusive"


def clock(V, key, nboot=NBOOT, seed=13):
    fams = sorted({r["mult"] for r in V})
    n_cross = sum(1 for r in V if r["cross"])
    out = dict(n_pairs=len(V), n_families=len(fams), n_cross=n_cross,
               testable=bool(len(V) >= MIN_V and len(fams) >= MIN_FAM and n_cross >= MIN_CROSS))
    rB = spearman([r["dB"] for r in V], [r[key] for r in V])
    rM = spearman([r["dM"] for r in V], [r[key] for r in V])

    def st(sub):
        return (spearman([r["dB"] for r in sub], [r[key] for r in sub]),
                spearman([r["dM"] for r in sub], [r[key] for r in sub]))

    boots = [b for b in family_boot(V, st, nboot, seed) if isf(b[0]) and isf(b[1])]
    nb = len(boots)
    frac_B = sum(1 for b, m in boots if b > m) / float(nb) if nb else NAN
    frac_M = sum(1 for b, m in boots if m > b) / float(nb) if nb else NAN
    p95_B = _pct([b for b, _ in boots], 0.95)
    p95_M = _pct([m for _, m in boots], 0.95)
    out.update(rho_B=rB, rho_M=rM, n_boot=nb, frac_B_gt_M=frac_B, frac_M_gt_B=frac_M,
               p05_B=_pct([b for b, _ in boots], 0.05), p95_B=p95_B,
               p05_M=_pct([m for _, m in boots], 0.05), p95_M=p95_M)
    out["decision"] = decide13(rB, rM, frac_B, frac_M, p95_B, p95_M) if out["testable"] else "untestable"
    return out


def source_scores(pair_rows, min_station_pairs=MIN_STATION_PAIRS):
    by, seen = {}, set()
    for r in pair_rows:
        key = tuple(sorted((r["ev1"], r["ev2"])))
        k = (r["station"],) + key
        if k in seen or not isf(r.get("D_late")):
            continue
        seen.add(k)
        by.setdefault(r["station"], []).append((key, float(r["D_late"])))
    pct = {}
    used = []
    for st, lst in sorted(by.items()):
        if len(lst) < min_station_pairs:
            continue
        used.append(st)
        rk = [float(x) + 1.0 for x in rankdata([d for _, d in lst])]
        n = float(len(lst))
        for (key, d), q in zip(lst, rk):
            pct.setdefault(key, []).append(q / n)
    return {k: dict(S=median(v), n_st=len(v)) for k, v in pct.items()}, used


# ----------------------------------------------------------------------------- new in h17
def read_yka(text):
    """ok == True rows of the YKA pair table, both clocks recomputed, with hardware generation."""
    rows = []
    for r in csv.DictReader(io.StringIO(text)):
        if str(r.get("ok", "")).strip() != "True":
            continue
        t1, t2 = fnum(r.get("t1")), fnum(r.get("t2"))
        a, b = (r.get("ev1") or "").strip(), (r.get("ev2") or "").strip()
        if not (isf(t1) and isf(t2) and a and b and a != b):
            continue
        g1, g2 = (r.get("gen1") or "").strip(), (r.get("gen2") or "").strip()
        rows.append(dict(mult=r.get("mult", ""), ev1=a, ev2=b, key=tuple(sorted((a, b))), t1=t1, t2=t2,
                         D=fnum(r.get("D")), D_c13=fnum(r.get("D_c13")), D_c310=fnum(r.get("D_c310")),
                         dB=d_B(t1, t2), dM=d_M(t1, t2), cross=is_cross(t1, t2),
                         m0324=(a in M0324 and b in M0324), gen1=g1, gen2=g2, gclass=gen_class(g1, g2)))
    return rows


def gen_class(g1, g2):
    if g1 not in ("legacy", "new") or g2 not in ("legacy", "new"):
        return "unknown"
    if g1 == g2:
        return "same-" + g1
    return "cross"


def read_scores(text):
    """{sorted (ev1, ev2): dict(S, n_st)} from h13_pairs_scored.csv."""
    out = {}
    for r in csv.DictReader(io.StringIO(text)):
        a, b = (r.get("ev1") or "").strip(), (r.get("ev2") or "").strip()
        S = fnum(r.get("S"))
        if not (a and b and isf(S)):
            continue
        try:
            n_st = int(float(r.get("n_st") or 0))
        except ValueError:
            n_st = 0
        out[tuple(sorted((a, b)))] = dict(S=S, n_st=n_st)
    return out


def verified(rows, scores, key="D", s_max=S_MAX):
    base = [r for r in rows if isf(r.get(key)) and not r["m0324"]]
    scored = [r for r in base if r["key"] in scores and scores[r["key"]]["n_st"] >= MIN_ST]
    V = [r for r in scored if scores[r["key"]]["S"] <= s_max + 1e-12]
    return scored, V


def calibration(scored, scores, key="D", nboot=NBOOT, seed=14):
    rho = spearman([scores[r["key"]]["S"] for r in scored], [r[key] for r in scored])
    cb = family_boot(scored, lambda sub: spearman([scores[r["key"]]["S"] for r in sub],
                                                  [r[key] for r in sub]), nboot, seed)
    return dict(n=len(scored), rho=rho, p05=_pct(cb, 0.05), p95=_pct(cb, 0.95))


def decide17(rho, p05, p95, n, nfam):
    """P17-1: is there a time signal (dM) within same-generation verified pairs?"""
    if n < MIN_V or nfam < MIN_FAM or not all(isf(x) for x in (rho, p05, p95)):
        return "untestable"
    if rho >= RHO_MIN and p05 > 0:
        return "time-signal"
    if rho < RHO_MIN and p95 < RHO_NONE_P95:
        return "no-time-signal"
    return "inconclusive"


def time_signal(V, key="D", nboot=NBOOT, seed=17):
    fams = sorted({r["mult"] for r in V})
    rho = spearman([r["dM"] for r in V], [r[key] for r in V])
    rho_B = spearman([r["dB"] for r in V], [r[key] for r in V])
    bs = [b for b in family_boot(V, lambda sub: spearman([r["dM"] for r in sub], [r[key] for r in sub]),
                                 nboot, seed) if isf(b)]
    p05, p95 = _pct(bs, 0.05), _pct(bs, 0.95)
    return dict(n_pairs=len(V), n_families=len(fams), n_cross_turn=sum(1 for r in V if r["cross"]),
                rho_M=rho, rho_B=rho_B, n_boot=len(bs), p05_M=p05, p95_M=p95,
                decision=decide17(rho, p05, p95, len(V), len(fams)))


def bins_table(V, key="D"):
    out = []
    for lo, hi in BINS:
        sub = [r for r in V if lo <= r["dM"] < hi]
        same = [r[key] for r in sub if r["gclass"].startswith("same-")]
        cr = [r[key] for r in sub if r["gclass"] == "cross"]
        out.append(dict(key=key, dM_lo=lo, dM_hi=hi, n_same=len(same), n_cross_gen=len(cr),
                        median_same=median(same), median_cross_gen=median(cr)))
    return out


def locked_family_stations(fam_rows, lock_s=LOCK_S):
    """(station, mult) whose stack onset sits > lock_s from the reference event's own onset."""
    out = set()
    for r in fam_rows:
        if (r.get("status") or "").strip() != "ok":
            continue
        a, b = fnum(r.get("onset_minus_pred_s")), fnum(r.get("ref_onset_minus_pred_s"))
        if isf(a) and isf(b) and abs(a - b) > lock_s + 1e-9:
            out.add(((r.get("station") or "").strip(), (r.get("mult") or "").strip()))
    return out


def read_csv_text(text):
    return list(csv.DictReader(io.StringIO(text)))


# ----------------------------------------------------------------------------- I/O
def fetch(url, tries=1, wait=60, must=None):
    last = None
    for k in range(tries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "inner-core-repeaters (MV2, Munim, Inc.)"})
            with urllib.request.urlopen(req, timeout=120) as r:
                body = r.read().decode("utf-8", "replace")
            if must is None or must in body:
                return body
            last = "missing marker"
        except Exception as ex:  # noqa: BLE001
            last = repr(ex)[:200]
        if k + 1 < tries:
            time.sleep(wait)
    print("fetch failed", url, last, flush=True)
    return None


def plot(out_dir, V):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    os.makedirs(os.path.join(out_dir, "plots"), exist_ok=True)
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.5))
    col = {"same-legacy": "tab:blue", "same-new": "tab:green", "cross": "tab:red"}
    for g, c in col.items():
        sub = [r for r in V if r["gclass"] == g]
        ax[0].scatter([r["dM"] for r in sub], [r["D"] for r in sub], s=8, c=c, alpha=0.7,
                      label="%s (%d)" % (g, len(sub)))
    ax[0].set_xlabel("time apart dM (yr)")
    ax[0].set_ylabel("YKA core D (source-verified pairs)")
    ax[0].legend(fontsize=8)
    tb = bins_table(V)
    xs = [0.5 * (b["dM_lo"] + b["dM_hi"]) for b in tb]
    ax[1].plot(xs, [b["median_same"] for b in tb], "o-", label="same generation")
    ax[1].plot(xs, [b["median_cross_gen"] for b in tb], "s-", color="tab:red", label="across generations")
    ax[1].set_xlabel("time apart dM (yr, bin centre)")
    ax[1].set_ylabel("median YKA core D")
    ax[1].legend(fontsize=8)
    fig.suptitle("Pre-registration 17: YKA hardware generation (h17 %s)" % VERSION)
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, "plots", "h17_overview.png"), dpi=100)
    plt.close(fig)


def jsonable(o):
    if isinstance(o, dict):
        return {str(k): jsonable(v) for k, v in o.items()}
    if isinstance(o, (list, tuple, set)):
        return [jsonable(v) for v in o]
    if isinstance(o, (np.bool_,)):
        return bool(o)
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating, float)):
        v = float(o)
        return v if math.isfinite(v) else None
    return o


def analyse(yka_rows, scores, mp3_pairs=None, mp3_fams=None, nboot=NBOOT):
    s = {}
    scored, V = verified(yka_rows, scores, "D")
    rep = clock(V, "D", nboot, 13)
    s["reproduction"] = dict(clock=rep, expected=EXPECT, reproduced=bool(
        rep["n_pairs"] == EXPECT["n_pairs"] and rep["n_families"] == EXPECT["n_families"]
        and rep["n_cross"] == EXPECT["n_cross"] and isf(rep["rho_B"]) and isf(rep["rho_M"])
        and abs(rep["rho_B"] - EXPECT["rho_B"]) < 1e-6 and abs(rep["rho_M"] - EXPECT["rho_M"]) < 1e-6))
    counts = {}
    for r in V:
        counts[r["gclass"]] = counts.get(r["gclass"], 0) + 1
    s["verified_by_generation"] = counts
    Vs = [r for r in V if r["gclass"].startswith("same-")]
    s["P17_1_time_signal_same_gen"] = time_signal(Vs, "D", nboot, 17)
    s["P17_2_clock_same_gen"] = clock(Vs, "D", nboot, 17)
    s["headline"] = s["P17_1_time_signal_same_gen"]["decision"]
    sec = {}
    for g in ("same-legacy", "same-new", "cross"):
        sub = [r for r in V if r["gclass"] == g]
        sec["time_signal_" + g] = time_signal(sub, "D", nboot, 17)
    sec["clock_same_legacy"] = clock([r for r in V if r["gclass"] == "same-legacy"], "D", nboot, 17)
    sec["cross_gen_indicator_rho_all_V"] = spearman([1.0 if r["gclass"] == "cross" else 0.0 for r in V],
                                                    [r["D"] for r in V])
    sec["bins_D"] = bins_table(V, "D")
    _, V310 = verified(yka_rows, scores, "D_c310")
    V310s = [r for r in V310 if r["gclass"].startswith("same-")]
    sec["time_signal_same_gen_D_c310"] = time_signal(V310s, "D_c310", nboot, 17)
    sec["clock_same_gen_D_c310"] = clock(V310s, "D_c310", nboot, 17)
    sec["bins_D_c310"] = bins_table(V310, "D_c310")
    s["secondary"] = sec
    if mp3_pairs is not None and mp3_fams is not None:
        locked = locked_family_stations(mp3_fams)
        n_ok = sum(1 for r in mp3_fams if (r.get("status") or "").strip() == "ok")
        base_scores, used0 = source_scores(mp3_pairs)
        kept = [r for r in mp3_pairs
                if ((r.get("station") or "").strip(), (r.get("mult") or "").strip()) not in locked]
        lock_scores, used1 = source_scores(kept)
        out = dict(n_family_stations_ok=n_ok, n_locked=len(locked),
                   locked_by_station={}, stations_used_all=used0, stations_used_unlocked=used1)
        for st, _m in sorted(locked):
            out["locked_by_station"][st] = out["locked_by_station"].get(st, 0) + 1
        for name, sc in (("rescored_all", base_scores), ("rescored_unlocked", lock_scores)):
            sc_scored, sc_V = verified(yka_rows, sc, "D")
            out[name] = dict(calibration=calibration(sc_scored, sc, "D", nboot, 14),
                             clock=clock(sc_V, "D", nboot, 13),
                             time_signal_same_gen=time_signal(
                                 [r for r in sc_V if r["gclass"].startswith("same-")], "D", nboot, 17))
        s["S17_locks"] = out
    return s, V


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-dir", default="/out")
    ap.add_argument("--wait-min", type=float, default=120.0)
    args = ap.parse_args()
    t0 = time.time()
    os.makedirs(args.out_dir, exist_ok=True)
    s = dict(version=VERSION, started=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
             yka_table=YKA_TABLE, h13_files=H13_FILES)
    tries = max(1, int(args.wait_min))
    pre = fetch(RAW + "main/" + PREREG, tries=tries, wait=60, must="Pre-registration 17")
    s["prereg_public_at_start"] = pre is not None
    s["prereg_wait_min"] = round((time.time() - t0) / 60.0, 1)
    if pre is None:   # never analyse before the plan is public
        s["error"] = "pre-registration not public; nothing analysed"
        with open(os.path.join(args.out_dir, "h17_summary.json"), "w") as f:
            json.dump(jsonable(s), f, indent=1)
        return 3
    yka_text = fetch(RAW + YKA_TABLE[0] + "/" + YKA_TABLE[1], tries=10, must="gen1")
    texts = [fetch(RAW + "main/" + p, tries=tries, wait=60, must="ev1" if "families" not in p else "station")
             for p in H13_FILES]
    if yka_text is None or texts[0] is None:
        s["error"] = "inputs missing"
        with open(os.path.join(args.out_dir, "h17_summary.json"), "w") as f:
            json.dump(jsonable(s), f, indent=1)
        return 2
    yka_rows = read_yka(yka_text)
    scores = read_scores(texts[0])
    s["n_yka_ok_rows"], s["n_scores"] = len(yka_rows), len(scores)
    mp3_pairs = read_csv_text(texts[1]) if texts[1] else None
    mp3_fams = read_csv_text(texts[2]) if texts[2] else None
    res, V = analyse(yka_rows, scores, mp3_pairs, mp3_fams)
    s.update(res)
    s["finished"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    s["minutes"] = round((time.time() - t0) / 60.0, 1)
    with open(os.path.join(args.out_dir, "h17_summary.json"), "w") as f:
        json.dump(jsonable(s), f, indent=1)
    with open(os.path.join(args.out_dir, "h17_pairs.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["mult", "ev1", "ev2", "t1", "t2", "dB", "dM", "cross_turn", "gen1", "gen2", "gclass",
                    "S", "D", "D_c310"])
        for r in V:
            w.writerow([r["mult"], r["ev1"], r["ev2"], r["t1"], r["t2"], r["dB"], r["dM"], r["cross"],
                        r["gen1"], r["gen2"], r["gclass"], scores[r["key"]]["S"], r["D"], r["D_c310"]])
    with open(os.path.join(args.out_dir, "h17_bins.csv"), "w", newline="") as f:
        rows = res["secondary"]["bins_D"] + res["secondary"]["bins_D_c310"]
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(jsonable(rows))
    try:
        plot(args.out_dir, V)
    except Exception as ex:  # noqa: BLE001
        print("plot failed", repr(ex)[:200], flush=True)
    print(json.dumps(jsonable(dict(headline=s["headline"], reproduction=s["reproduction"],
                                   by_gen=s["verified_by_generation"],
                                   p17_1=s["P17_1_time_signal_same_gen"])))[:4000], flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())

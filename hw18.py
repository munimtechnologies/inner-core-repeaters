#!/usr/bin/env python3
"""hw18.py: how much does the 2013 YKA hardware change alone change the PKIKP coda?

Part of inner-core-repeaters by MV2, Munim, Inc. (made with AI agents). Licence: Apache-2.0.
Pre-registered in findings/2026-11-19-prereg18-hardware-overlap.md (pre-registration 18),
written and published before this script ran. The script refuses to analyse anything until
that file is public.

Why. Pre-registration 13 (job icr-h13-v1) found a time trend in the PKIKP coda difference D
among source-verified Yellowknife (YKA) repeater pairs. Pairs far apart in time often straddle
the 2013 YKA hardware change, so pre-registration 17 tests pairs within one hardware
generation. This script measures the hardware effect directly. EarthScope station metadata
(checked 2026-11-19) show the legacy elements (CN.YKB*/YKR*, Geotech S-13, 20 sps CNSN
digitiser, 1989-01-25 to 2014-04-02) and the current elements (CN.YKAB*/YKAR*, the same
Geotech S-13 sensors with a new pre-amp and Guralp digitiser, 40 sps, from 2013-06-26) at the
same positions. For about nine months both recorded the same ground motion. For every South
Sandwich Islands earthquake in that overlap, the legacy beam and the current beam of the SAME
event are compared with the YKA pipeline's own alignment, windows and gates (yka.py 0.2.0 at
v0.8.0, onset picked within ak135 +/- 2 s as in the pick amendment). Their D (D_hw) is the
coda change caused by the hardware alone. D_hw is then removed from the cross-generation
pairs of H13's verified set and H13's time trend is tested again.

Inputs (public): yka.py, model.py, phases.py, xcorr.py at commit 7067d4fc; the YKA pair table
results/ykab-v2/yka_pairs.csv at f65c966a; results/h13-v1/h13_pairs_scored.csv on main;
USGS ComCat; waveforms from EarthScope and NRCan FDSN dataselect.

Writes to --out-dir: hw18_summary.json, hw18_events.csv, hw18_elements.csv,
hw18_catalog.csv, plots/hw18_beams.png, plots/hw18_overview.png.
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
import traceback
import urllib.parse
import urllib.request

import numpy as np

VERSION = "0.1.0"
RAW = "https://raw.githubusercontent.com/munimtechnologies/inner-core-repeaters/"
MOD_COMMIT = "7067d4fc2b0a2fdb1378a1c406dd680c8d0b69b1"  # v0.8.0: yka.py 0.2.0 and helpers
MODULES = ("model.py", "phases.py", "xcorr.py", "yka.py")
YKA_TABLE = ("f65c966adf2cce3a22ffd9e4629c842959121407", "results/ykab-v2/yka_pairs.csv")
H13_SCORES = "results/h13-v1/h13_pairs_scored.csv"
PREREG = "findings/2026-11-19-prereg18-hardware-overlap.md"
PREREG_MARK = "Pre-registration 18"

# H13 primary result as published (h13_summary.json, grade.primary.clock); must be reproduced.
EXPECT = dict(n_pairs=328, n_families=42, n_cross=96, rho_B=0.12651797149019148,
              rho_M=0.1575493725434452)

M0324 = frozenset(("usp000bt5x", "usp000fz4d", "usp000hute", "usc000t4wj", "us7000am0e",
                   "us6000h2wz", "us6000ittv", "us6000mk65", "us6000runq"))
TURN, RATIO = 2008.5, 2.5
MIN_ST = 2
S_MAX = 1.0 / 3.0
RHO_MIN = 0.10
RHO_NONE_P95 = 0.20
FRAC = 0.90
MIN_V, MIN_FAM, MIN_CROSS = 40, 4, 10
NBOOT = 1000
BINS = [(0.0, 3.0), (3.0, 6.0), (6.0, 10.0), (10.0, 15.0), (15.0, 30.0)]
NAN = float("nan")

# Pre-registered for H18
USGS = "https://earthquake.usgs.gov/fdsnws/event/1/query"
BOX = dict(minlatitude=-61.5, maxlatitude=-54.0, minlongitude=-31.0, maxlongitude=-22.0)
OVERLAP = ("2013-06-27T00:00:00", "2014-04-02T00:00:00")
MIN_MAG = 4.5
REPICK_S = 2.0
PAIR_THR, MIN_SNR_DIRECT, MIN_SNR_CODA = 0.90, 5.0, 2.0
MIN_COMMON = 5
MIN_EVENTS = 4
SEED_EV = 18
KEYS = ("D", "D_c13", "D_c310", "D_full")


# ----------------------------------------------------------------------------- copied from h17 0.1.0
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


def gen_class(g1, g2):
    if g1 not in ("legacy", "new") or g2 not in ("legacy", "new"):
        return "unknown"
    if g1 == g2:
        return "same-" + g1
    return "cross"


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
                         D_full=fnum(r.get("D_full")),
                         dB=d_B(t1, t2), dM=d_M(t1, t2), cross=is_cross(t1, t2),
                         m0324=(a in M0324 and b in M0324), gen1=g1, gen2=g2, gclass=gen_class(g1, g2),
                         pol_flip=str(r.get("polarity_flipped", "")).strip() == "True"))
    return rows


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


def bins_table(V, key="D"):
    out = []
    for lo, hi in BINS:
        sub = [r for r in V if lo <= r["dM"] < hi]
        same = [r[key] for r in sub if r["gclass"].startswith("same-")]
        cr = [r[key] for r in sub if r["gclass"] == "cross"]
        out.append(dict(key=key, dM_lo=lo, dM_hi=hi, n_same=len(same), n_cross_gen=len(cr),
                        median_same=median(same), median_cross_gen=median(cr),
                        excess=median(cr) - median(same) if same and cr else NAN))
    return out


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


# ----------------------------------------------------------------------------- new in hw18 (pure)
def choose_providers(sets):
    """sets = {provider: {"legacy": set(elements), "new": set(elements)}}.

    For each generation, the provider with the most elements of that generation (ties: the
    alphabetically first provider). None if no provider has any element of it.
    """
    out = {}
    for g in ("legacy", "new"):
        best, nb = None, 0
        for prov in sorted(sets):
            n = len(sets[prov].get(g, ()))
            if n > nb:
                best, nb = prov, n
        out[g] = best
    return out


def correct_D(D, h, mode="mult"):
    """Remove a hardware decorrelation h from a cross-generation D.

    mult (primary): cc_cross = cc_same * cc_hw, so D_corr = 1 - (1 - D) / (1 - h).
    add (sensitivity): D_corr = D - h.
    """
    if not (isf(D) and isf(h)):
        return NAN
    if mode == "add":
        return float(D) - float(h)
    if float(h) >= 1.0:
        return NAN
    return 1.0 - (1.0 - float(D)) / (1.0 - float(h))


def corrected_rows(V, h, key="D", mode="mult"):
    out = []
    for r in V:
        c = dict(r)
        if r["gclass"] == "cross":
            c[key] = correct_D(r[key], h, mode)
        out.append(c)
    return out


def hw_stats(values, nboot=NBOOT, seed=SEED_EV):
    """Median of the per-event hardware D and its event-bootstrap 5th/95th percentiles."""
    xs = [float(x) for x in values if isf(x)]
    out = dict(n=len(xs), median=median(xs), mean=float(np.mean(xs)) if xs else NAN,
               min=min(xs) if xs else NAN, max=max(xs) if xs else NAN, p05=NAN, p95=NAN)
    if len(xs) >= 2:
        rng = random.Random(seed)
        bs = [median([rng.choice(xs) for _ in xs]) for _ in range(nboot)]
        out.update(p05=_pct(bs, 0.05), p95=_pct(bs, 0.95))
    elif len(xs) == 1:
        out.update(p05=xs[0], p95=xs[0])
    return out


def decide18(reproduced, n_events, rho_cons, p05_cons, rho_c, p95_c):
    """H18 primary decision (see the pre-registration)."""
    if not reproduced:
        return "not-reproduced"
    if n_events < MIN_EVENTS or not all(isf(x) for x in (rho_cons, p05_cons, rho_c, p95_c)):
        return "untestable"
    if rho_cons >= RHO_MIN and p05_cons > 0:
        return "trend-survives-hardware"
    if rho_c < RHO_MIN and p95_c < RHO_NONE_P95:
        return "hardware-explains"
    return "inconclusive"


def grade(V, hs, reproduced, key="D", mode="mult", nboot=NBOOT):
    """hs = hw_stats of the per-event hardware D for `key`.

    Conservative correction uses the 95th bootstrap percentile of the median (hs['p95']);
    central correction uses the median.
    """
    n_ev = hs["n"]
    if n_ev < MIN_EVENTS or not isf(hs["median"]) or not isf(hs["p95"]):
        return dict(n_events=n_ev, decision=decide18(reproduced, n_ev, NAN, NAN, NAN, NAN))
    Vc = corrected_rows(V, hs["median"], key, mode)
    Vk = corrected_rows(V, hs["p95"], key, mode)
    cc_ = clock(Vc, key, nboot, 13)
    ck = clock(Vk, key, nboot, 13)
    return dict(n_events=n_ev, h_central=hs["median"], h_conservative=hs["p95"], mode=mode,
                central=cc_, conservative=ck,
                decision=decide18(reproduced, n_ev, ck["rho_M"], ck["p05_M"], cc_["rho_M"], cc_["p95_M"]))


def parse_catalog(text, decimal_year):
    evs = []
    for r in csv.DictReader(io.StringIO(text or "")):
        if r.get("type", "earthquake") != "earthquake":
            continue
        try:
            evs.append(dict(id=r["id"], time=r["time"], lat=float(r["latitude"]), lon=float(r["longitude"]),
                            depth=float(r["depth"]) if r.get("depth") not in (None, "") else 10.0,
                            mag=float(r["mag"]) if r.get("mag") not in (None, "") else NAN,
                            t=decimal_year(r["time"])))
        except (KeyError, ValueError):
            continue
    evs.sort(key=lambda e: e["t"])
    return evs


def pair_ok(cc_ref, D, snr_dir, snr_coda):
    """The YKA pair table's ok gates (yka.py 0.2.0 defaults), applied to a same-event pair."""
    q = list(snr_dir) + list(snr_coda)
    return bool(all(isf(v) for v in q) and isf(D) and isf(cc_ref) and cc_ref >= PAIR_THR
                and min(snr_dir) >= MIN_SNR_DIRECT and min(snr_coda) >= MIN_SNR_CODA)


# ----------------------------------------------------------------------------- waveform part
def load_modules(mod_dir):
    os.makedirs(mod_dir, exist_ok=True)
    for m in MODULES:
        p = os.path.join(mod_dir, m)
        if not os.path.exists(p):
            body = fetch(RAW + MOD_COMMIT + "/" + m, tries=6, wait=20,
                         must='VERSION = "0.2.0"' if m == "yka.py" else None)
            if body is None:
                return None
            with open(p, "w", encoding="utf-8") as f:
                f.write(body)
    if mod_dir not in sys.path:
        sys.path.insert(0, mod_dir)
    import yka  # noqa: E402
    return yka


def fetch_both(yka, session, ev, ph, raw_dir):
    """Both archives for one event, each cached separately. Returns ({prov: Stream}, {prov: status})."""
    from obspy import UTCDateTime, read

    dist = float(yka.gc_deg(ev["lat"], ev["lon"], *yka.YKA_CENTER))
    tdf = ph.t("df", dist, ev["depth"])
    if not isf(tdf):
        return {}, {"all": "no-df-prediction(dist=%.2f,depth=%.1f)" % (dist, ev["depth"])}
    o = UTCDateTime(ev["time"])
    t1, t2 = o + tdf - yka.LONG_PRE, o + tdf + yka.LONG_POST
    streams, status = {}, {}
    for name, url in (("earthscope", yka.EARTHSCOPE), ("nrcan", yka.NRCAN)):
        p = os.path.join(raw_dir, "%s_%s.mseed" % (ev["id"], name))
        if os.path.exists(p):
            try:
                streams[name] = read(p)
                status[name] = "cached"
                continue
            except Exception:  # noqa: BLE001
                os.remove(p)
        if os.path.exists(p + ".nodata"):
            status[name] = "nodata"
            continue
        st, s = yka.dataselect(url, t1, t2, session)
        status[name] = s
        if st is not None:
            streams[name] = st
            try:
                st.write(p, format="MSEED")
            except Exception as ex:  # noqa: BLE001
                print("could not cache", ev["id"], name, repr(ex)[:120], flush=True)
        elif s == "nodata":
            open(p + ".nodata", "w").close()
        time.sleep(0.3)
    return streams, status


def element_sets(yka, streams):
    sets = {}
    for prov, st in streams.items():
        sets[prov] = {"legacy": set(), "new": set()}
        for tr in st:
            c = yka.canonical(tr.stats.station)
            if c is None or tr.stats.channel != "SHZ":
                continue
            sets[prov][yka.generation(tr.stats.station)].add(c)
    return sets


def sub_stream(yka, st, gen, keep):
    from obspy import Stream

    return Stream([tr for tr in st if tr.stats.channel == "SHZ" and yka.generation(tr.stats.station) == gen
                   and yka.canonical(tr.stats.station) in keep])


def process_pair(yka, streams, ev, ph):
    """Legacy and current beams of one event on the same elements. Returns (dL, dN, info)."""
    sets = element_sets(yka, streams)
    prov = choose_providers(sets)
    info = dict(prov_legacy=prov["legacy"], prov_new=prov["new"],
                n_legacy={p: len(s["legacy"]) for p, s in sets.items()},
                n_new={p: len(s["new"]) for p, s in sets.items()})
    if prov["legacy"] is None or prov["new"] is None:
        info["status"] = "no-legacy" if prov["legacy"] is None else "no-new"
        return None, None, info
    common = sets[prov["legacy"]]["legacy"] & sets[prov["new"]]["new"]
    for _ in range(2):
        if len(common) < MIN_COMMON:
            info["status"] = "few-common(%d)" % len(common)
            return None, None, info
        dL, sL = yka.process_event(sub_stream(yka, streams[prov["legacy"]], "legacy", common), ev, ph)
        dN, sN = yka.process_event(sub_stream(yka, streams[prov["new"]], "new", common), ev, ph)
        if dL is None or dN is None:
            info["status"] = "process:" + (sL if dL is None else "ok") + "/" + (sN if dN is None else "ok")
            return None, None, info
        cn = set(str(x) for x in dL["names"]) & set(str(x) for x in dN["names"])
        if cn == set(str(x) for x in dL["names"]) == set(str(x) for x in dN["names"]):
            info.update(status="ok", n_common=len(cn), elements=sorted(cn))
            return dL, dN, info
        common = cn
    info["status"] = "element-sets-unstable"
    return None, None, info


def measure(yka, dL, dN, ev):
    """The pipeline's pair measurement (analyze_multiplet with the two beams as members)."""
    FS, PRE, N = yka.FS, yka.PRE, yka.N
    data = {"legacy": dL, "new": dN}
    root = "legacy" if dL["snr"] >= dN["snr"] else "new"
    other = "new" if root == "legacy" else "legacy"
    a = yka.align_member(data[root]["beam"], data[other]["beam"], FS, PRE, allow_flip=True)
    out = dict(root=root, cc_env=a.get("cc_env"), cc_fine=a.get("cc_fine"), flip=bool(a.get("flip")))
    if not (isf(a["off"]) and isf(a["cc_fine"]) and a["cc_fine"] >= 0.5):
        out["status"] = "alignment-failed"
        return out, []
    off = {root: 0.0, other: float(a["off"])}
    flips = {root: False, other: bool(a["flip"])}
    stack = np.zeros(N)
    for m in ("legacy", "new"):
        stack += yka.fshift(data[m]["beam"], -off[m]) * (-1.0 if flips[m] else 1.0)
    stack /= 2.0
    onset = yka.pick_onset(stack, FS, center=int(PRE * FS), search_s=REPICK_S, sta_s=0.5, lta_s=8.0, thr=6.0)
    if onset is None:
        out["status"] = "no-DF-pick"
        return out, []
    out["onset_minus_pred_s"] = round((onset - PRE * FS) / FS, 3)
    n0, n1 = int((PRE - 35) * FS), int((PRE - 5) * FS)

    def snr(m, w):
        b = data[m]["beam"]
        s = int(round(off[m]))
        lo, hi = onset + int(round(w[0] * FS)) + s, onset + int(round(w[1] * FS)) + s
        if lo < 0 or hi > N:
            return NAN
        nb = np.sqrt(np.mean(b[n0:n1] ** 2))
        return float(np.sqrt(np.mean(b[lo:hi] ** 2)) / nb) if nb > 0 else NAN

    sgn = -1.0 if flips["legacy"] != flips["new"] else 1.0
    xa, xb = dL["beam"], dN["beam"] * sgn
    r = yka.yka_pair(xa, xb, onset, off["legacy"], off["new"], FS)
    if r is None:
        out["status"] = "window-off-trace"
        return out, []
    out.update(cc_ref=r["cc_ref"], cc_df0=r["cc_df0"], cc_df_best=r["cc_df_best"], dt_df_s=r["dt_df_s"],
               lag_new_minus_legacy_s=r["lag_ref_s"], D=1.0 - r["cc_df0"], polarity_flipped=sgn < 0)
    for k, w in yka.CODA_SENS.items():
        rs = yka.yka_pair(xa, xb, onset, off["legacy"], off["new"], FS, window=w)
        out[k] = (1.0 - rs["cc_df0"]) if rs is not None else NAN
    out.update(snr_dir_legacy=snr("legacy", yka.DIRECT), snr_dir_new=snr("new", yka.DIRECT),
               snr_coda_legacy=snr("legacy", yka.CODA), snr_coda_new=snr("new", yka.CODA))
    out["ok"] = pair_ok(out["cc_ref"], out["D"], (out["snr_dir_legacy"], out["snr_dir_new"]),
                        (out["snr_coda_legacy"], out["snr_coda_new"]))
    el_rows = []
    nl, nn = [str(x) for x in dL["names"]], [str(x) for x in dN["names"]]
    XL, XN = np.asarray(dL["X"], dtype=np.float64), np.asarray(dN["X"], dtype=np.float64)
    st_cc = []
    for s in sorted(set(nl) & set(nn)):
        rs = yka.yka_pair(XL[nl.index(s)], XN[nn.index(s)] * sgn, onset, off["legacy"], off["new"], FS)
        if rs is None:
            continue
        el_rows.append(dict(id=ev["id"], element=s, cc_ref=rs["cc_ref"], D=1.0 - rs["cc_df0"],
                            lag_s=rs["lag_ref_s"]))
        if rs["cc_ref"] >= 0.6:
            st_cc.append(rs["cc_df0"])
    out["Ds"] = (1.0 - float(np.median(st_cc))) if len(st_cc) >= 5 else NAN
    out["n_el"] = len(st_cc)
    out["status"] = "ok"
    return out, el_rows, (onset, off, flips)


def write_csv(path, rows, keys):
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(jsonable(r))


def plot_beams(path, yka, shown):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    if not shown:
        return
    FS, PRE = yka.FS, yka.PRE
    fig, axs = plt.subplots(len(shown), 1, figsize=(10, 1.8 * len(shown) + 0.8), squeeze=False)
    for ax, (ev, dL, dN, onset, off, flips, row) in zip(axs[:, 0], shown):
        t = (np.arange(yka.N) - onset) / FS
        for m, d, c in (("legacy", dL, "tab:blue"), ("new", dN, "tab:red")):
            y = yka.fshift(d["beam"], -off[m]) * (-1.0 if flips[m] else 1.0)
            ax.plot(t, y / (np.max(np.abs(y)) or 1.0), color=c, lw=0.8,
                    label="%s (%d el.)" % ("legacy YK*" if m == "legacy" else "current YKA*", len(d["names"])))
        ax.axvspan(yka.CODA[0], yka.CODA[1], color="orange", alpha=0.2)
        ax.set_xlim(-5, 12)
        ax.set_title("%s %s M%.1f   D_hw=%.3f  D_c310=%.3f" % (ev["time"][:10], ev["id"], ev["mag"],
                                                           fnum(row.get("D")), fnum(row.get("D_c310"))),
                     fontsize=8)
        ax.legend(fontsize=6, loc="upper right")
    axs[-1, 0].set_xlabel("s after picked PKIKP onset (1-2 Hz beams; shaded = measured coda 1-6 s)")
    fig.suptitle("Pre-registration 18: same earthquake, legacy vs current YKA hardware (hw18 %s)" % VERSION)
    fig.tight_layout()
    fig.savefig(path, dpi=100)
    plt.close(fig)


def plot_overview(path, ev_rows, V):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    ok = [r for r in ev_rows if r.get("ok")]
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.5))
    bins = np.linspace(0, 1, 41)
    same = [r["D"] for r in V if r["gclass"].startswith("same-") and isf(r["D"])]
    cross = [r["D"] for r in V if r["gclass"] == "cross" and isf(r["D"])]
    ax[0].hist(same, bins=bins, alpha=0.5, density=True, label="H13 verified, same generation (%d)" % len(same))
    ax[0].hist(cross, bins=bins, alpha=0.5, density=True, label="H13 verified, across generations (%d)" % len(cross))
    for r in ok:
        ax[0].axvline(r["D"], color="k", lw=0.8, alpha=0.6)
    ax[0].set_xlabel("D (1-6 s coda); black lines = same event, legacy vs current hardware")
    ax[0].legend(fontsize=7)
    ax[1].scatter([r["snr_coda_legacy"] for r in ok], [r["D"] for r in ok], s=14, label="D (1-6 s)")
    ax[1].scatter([r["snr_coda_legacy"] for r in ok], [r["D_c310"] for r in ok], s=14, marker="s",
                  label="D_c310 (3-10 s)")
    ax[1].set_xscale("log")
    ax[1].set_xlabel("legacy coda SNR")
    ax[1].set_ylabel("hardware-only D")
    ax[1].legend(fontsize=7)
    fig.suptitle("Pre-registration 18: hardware-only coda change vs H13 pair D (hw18 %s)" % VERSION)
    fig.tight_layout()
    fig.savefig(path, dpi=100)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", default="/data")
    ap.add_argument("--out-dir", default="/out")
    ap.add_argument("--wait-min", type=float, default=120.0)
    ap.add_argument("--time-budget-min", type=float, default=150.0)
    args = ap.parse_args()
    t0 = time.time()
    os.makedirs(args.out_dir, exist_ok=True)
    plots = os.path.join(args.out_dir, "plots")
    os.makedirs(plots, exist_ok=True)
    s = dict(version=VERSION, started=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
             rules=dict(overlap=OVERLAP, box=BOX, min_mag=MIN_MAG, repick_s=REPICK_S, pair_thr=PAIR_THR,
                        min_snr_direct=MIN_SNR_DIRECT, min_snr_coda=MIN_SNR_CODA, min_common=MIN_COMMON,
                        min_events=MIN_EVENTS, rho_min=RHO_MIN, rho_none_p95=RHO_NONE_P95, nboot=NBOOT,
                        seed_events=SEED_EV, mod_commit=MOD_COMMIT, yka_table=YKA_TABLE, h13_scores=H13_SCORES))
    summary_path = os.path.join(args.out_dir, "hw18_summary.json")

    def dump():
        with open(summary_path, "w") as f:
            json.dump(jsonable(s), f, indent=1)

    pre = fetch(RAW + "main/" + PREREG, tries=max(1, int(args.wait_min)), wait=60, must=PREREG_MARK)
    s["prereg_public_at_start"] = pre is not None
    s["prereg_wait_min"] = round((time.time() - t0) / 60.0, 1)
    if pre is None:  # never analyse before the plan is public
        s["error"] = "pre-registration not public; nothing analysed"
        dump()
        return 3
    yka_text = fetch(RAW + YKA_TABLE[0] + "/" + YKA_TABLE[1], tries=10, must="gen1")
    sc_text = fetch(RAW + "main/" + H13_SCORES, tries=10, must="ev1")
    if yka_text is None or sc_text is None:
        s["error"] = "inputs missing"
        dump()
        return 2
    yka_rows = read_yka(yka_text)
    scores = read_scores(sc_text)
    _, V = verified(yka_rows, scores, "D")
    rep = clock(V, "D", NBOOT, 13)
    reproduced = bool(rep["n_pairs"] == EXPECT["n_pairs"] and rep["n_families"] == EXPECT["n_families"]
                      and rep["n_cross"] == EXPECT["n_cross"] and isf(rep["rho_B"]) and isf(rep["rho_M"])
                      and abs(rep["rho_B"] - EXPECT["rho_B"]) < 1e-6 and abs(rep["rho_M"] - EXPECT["rho_M"]) < 1e-6)
    s["reproduction"] = dict(clock=rep, expected=EXPECT, reproduced=reproduced)
    counts = {}
    for r in V:
        counts[r["gclass"]] = counts.get(r["gclass"], 0) + 1
    s["verified_by_generation"] = counts
    xg = [r for r in V if r["gclass"] == "cross"]
    s["table_polarity_flipped_frac_cross_gen"] = (sum(1 for r in xg if r["pol_flip"]) / float(len(xg))
                                                  if xg else None)
    s["bins_D"] = bins_table(V, "D")
    s["bins_D_c310"] = bins_table(V, "D_c310")
    dump()

    yka = load_modules("/tmp/hw18_mods")
    if yka is None:
        s["error"] = "pinned modules could not be fetched"
        dump()
        return 2
    s["yka_version"] = getattr(yka, "VERSION", None)
    ph = yka.Phases(yka.load_or_build(os.path.join(args.data_dir, "ak135_yka_grid_v1.npz")))
    q = dict(format="csv", starttime=OVERLAP[0], endtime=OVERLAP[1], minmagnitude=MIN_MAG,
             orderby="time-asc", **BOX)
    cat_text = fetch(USGS + "?" + urllib.parse.urlencode(q), tries=5, wait=30, must="latitude")
    evs = parse_catalog(cat_text, yka.model.decimal_year)
    s["n_catalog"] = len(evs)
    with open(os.path.join(args.out_dir, "hw18_catalog.csv"), "w", newline="") as f:
        f.write(cat_text or "")
    raw_dir = os.path.join(args.data_dir, "hw18_raw")
    os.makedirs(raw_dir, exist_ok=True)
    import requests

    session = requests.Session()
    session.headers["User-Agent"] = "inner-core-repeaters/hw18-%s (MV2, Munim, Inc.; mv2@munimtech.com)" % VERSION
    ev_rows, el_rows, shown = [], [], []
    for k, ev in enumerate(evs):
        row = dict(id=ev["id"], time=ev["time"], mag=ev["mag"], depth=ev["depth"], t=ev["t"])
        if (time.time() - t0) / 60.0 > args.time_budget_min:
            row["status"] = "skipped-time-budget"
            ev_rows.append(row)
            continue
        try:
            streams, fst = fetch_both(yka, session, ev, ph, raw_dir)
            row["fetch"] = ";".join("%s:%s" % kv for kv in sorted(fst.items()))
            dL, dN, info = process_pair(yka, streams, ev, ph) if streams else (None, None, dict(status="no-data"))
            row.update({kk: info.get(kk) for kk in ("prov_legacy", "prov_new", "n_common")})
            row["status"] = info.get("status")
            if dL is not None:
                res = measure(yka, dL, dN, ev)
                m, els = res[0], res[1]
                row.update(m)
                row["snr_legacy"], row["snr_new"] = dL["snr"], dN["snr"]
                el_rows.extend(els)
                if len(res) == 3 and m.get("status") == "ok":
                    onset, off, flips = res[2]
                    shown.append((ev, dL, dN, onset, off, flips, row))
        except Exception as ex:  # noqa: BLE001
            traceback.print_exc()
            row["status"] = "exception:" + repr(ex)[:150]
        row["ok"] = bool(row.get("ok", False))
        ev_rows.append(row)
        print("event %d/%d %s %s ok=%s D=%s" % (k + 1, len(evs), ev["id"], row.get("status"), row["ok"],
                                                row.get("D")), flush=True)
    keys = ["id", "time", "t", "mag", "depth", "fetch", "prov_legacy", "prov_new", "n_common", "status", "ok",
            "root", "cc_env", "cc_fine", "flip", "polarity_flipped", "onset_minus_pred_s", "cc_ref", "cc_df0",
            "cc_df_best", "dt_df_s", "lag_new_minus_legacy_s", "D", "D_c13", "D_c310", "D_full", "Ds", "n_el",
            "snr_legacy", "snr_new", "snr_dir_legacy", "snr_dir_new", "snr_coda_legacy", "snr_coda_new"]
    write_csv(os.path.join(args.out_dir, "hw18_events.csv"), ev_rows, keys)
    write_csv(os.path.join(args.out_dir, "hw18_elements.csv"), el_rows, ["id", "element", "cc_ref", "D", "lag_s"])

    ok = [r for r in ev_rows if r["ok"]]
    s["n_events_measured"] = sum(1 for r in ev_rows if r.get("status") == "ok")
    s["n_events_gated"] = len(ok)
    status_counts = {}
    for r in ev_rows:
        key = str(r.get("status")).split("(")[0]
        status_counts[key] = status_counts.get(key, 0) + 1
    s["status_counts"] = status_counts
    hw = {k: hw_stats([r.get(k) for r in ok]) for k in KEYS}
    s["hardware_D"] = hw
    s["hardware_flip_frac"] = (sum(1 for r in ok if r.get("polarity_flipped")) / float(len(ok))) if ok else None
    s["hardware_lag_s_median"] = median([r.get("lag_new_minus_legacy_s") for r in ok])
    s["hardware_dt_df_s_median"] = median([r.get("dt_df_s") for r in ok])
    s["element_D_median"] = median([r["D"] for r in el_rows if r.get("cc_ref", 0) >= 0.6
                                    and r["id"] in {x["id"] for x in ok}])
    s["P18_1"] = grade(V, hw["D"], reproduced, "D", "mult")
    s["headline"] = s["P18_1"]["decision"]
    s["secondary"] = dict(
        D_additive=grade(V, hw["D"], reproduced, "D", "add"),
        D_c310=grade([r for r in verified(yka_rows, scores, "D_c310")[1]], hw["D_c310"], reproduced, "D_c310", "mult"),
        D_c13=grade([r for r in verified(yka_rows, scores, "D_c13")[1]], hw["D_c13"], reproduced, "D_c13", "mult"),
        excess_vs_hardware=[dict(b, hardware_D_median=hw["D"]["median"]) for b in s["bins_D"]],
    )
    s["finished"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    s["minutes"] = round((time.time() - t0) / 60.0, 1)
    dump()
    try:
        shown.sort(key=lambda x: -min(fnum(x[6].get("snr_legacy")), fnum(x[6].get("snr_new"))))
        plot_beams(os.path.join(plots, "hw18_beams.png"), yka, shown[:6])
    except Exception as ex:  # noqa: BLE001
        print("beam plot failed", repr(ex)[:200], flush=True)
    try:
        plot_overview(os.path.join(plots, "hw18_overview.png"), ev_rows, V)
    except Exception as ex:  # noqa: BLE001
        print("overview plot failed", repr(ex)[:200], flush=True)
    print(json.dumps(jsonable(dict(headline=s["headline"], reproduced=reproduced, n_gated=len(ok),
                                   hardware_D=hw["D"], P18_1=s["P18_1"])))[:4000], flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())

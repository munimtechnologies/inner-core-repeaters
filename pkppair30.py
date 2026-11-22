# pkppair30.py 0.1.0 -- pre-registration 30 of inner-core-repeaters:
# can test-site position alone make the outer-core "change" reported from nuclear-test PKP?
#
# Claim under test: Y. Zhou (2026), JGR Solid Earth 131, e2026JB033798, doi:10.1029/2026JB033798.
# French underground tests (Mururoa/Fangataufa, 1977-1995) recorded at Borovoye (BRVK); 112 pairs
# from 18 tests whose catalogue epicentres agree within 0.05 deg. PKrKP-PP and PKiKP-PP double
# differences are read as 0.1-0.15 s faster in 1982-83 and 1988-90 and 0.15-0.2 s slower in 1995,
# relative to 1977.
#
# PKP and PP leave the source with very different ray parameters (about 2 vs 6 s/deg at 141.6 deg),
# so a shot moved 1 km along the path to BRVK changes PKP-PP by about 35 ms. A 0.05 deg pairing
# tolerance allows 5.5 km. This job predicts, per test, the PKP-PP anomaly that position alone
# makes, under two readings of Zhou's processing, inverts the pair differences for one term per
# year exactly as a pair double-difference inversion does (1977 = 0), and grades the result
# against the reported epoch values with rules fixed in
# findings/2026-11-22-prereg30-pkp-position-plan.md before Table S1 was read by us.
#
#   reading R (raw):        no distance correction; anomaly_i = T(true_i)
#   reading C (catalogue):  pairs corrected with catalogue (Table S1) epicentres in a 1-D model;
#                           anomaly_i = T(true_i) - T(catalogue_i)
#   T = (phase - PP) travel time from TauP iasp91, source depth 0.6 km; "true" = the disclosed
#   shaft position (data/mururoa_shafts_1977_1996.csv).
#
# Inputs fetched at run time from the public repo (main): data/zhou2026_table_s1.csv (Zhou's
# Table S1, transcribed) and data/zhou2026_methods.json (what the paper says about distance
# correction). The job waits for them; without them it ends "no-input".
#
# Part of munimtechnologies/inner-core-repeaters. Licence: Apache-2.0.
# Made with AI agents working in MV2's name (MV2, Munim, Inc.); no human expert has checked it.
import argparse
import calendar
import csv
import hashlib
import io
import json
import math
import os
import re
import sys
import time

import numpy as np

VERSION = "0.1.0"
BRVK = (53.0581, 70.2828)
REF = (-21.835, -138.91)
SRC_DEPTH_KM = 0.6
KM_PER_DEG = 111.195
PAIR_TOL_DEG = 0.05
REF_YEAR = 1977
EPOCHS = {"E1982_83": (1982, 1983), "E1988_90": (1988, 1990), "E1995": (1995, 1995)}
# Zhou's reported epoch values relative to 1977 (s; negative = faster), midpoints of the quoted ranges:
# abstract "0.1-0.15 s faster during 1982-1983 and 1988-1990"; section 3 "0.15 s faster in 1988-1990";
# "0.15-0.2 s slower in 1995".
ZHOU = {"E1982_83": -0.125, "E1988_90": -0.15, "E1995": 0.175}
OBSERVABLES = {"PKiKP_PP": ("PKiKP", "PP"), "PKIKP_PP": ("PKIKP", "PP")}
PRIMARY_OBS = "PKiKP_PP"
SIGN_MIN_S = 0.02          # |pred| below this counts as no sign
EXPLAIN_RMS_S = 0.05       # residual rms for "explains"
PARTLY_FRAC = 0.75         # residual rms below this fraction of rms(Zhou) for "partly"
G0_MIN_MATCHED = 14
REPO_RAW = "https://raw.githubusercontent.com/munimtechnologies/inner-core-repeaters/main/"
S1_PATH = "data/zhou2026_table_s1.csv"
METHODS_PATH = "data/zhou2026_methods.json"
PLAN_PATH = "findings/2026-11-22-prereg30-pkp-position-plan.md"
SHAFTS_PATH = "data/mururoa_shafts_1977_1996.csv"


def log(*a):
    print(time.strftime("%H:%M:%S"), "[pkppair30]", *a, flush=True)


# ---------------------------------------------------------------- pure helpers (unit-tested)
def gc_dist_az(lat1, lon1, lat2, lon2):
    """Spherical great-circle distance (deg) and azimuth at point 1 (deg)."""
    p1, p2 = math.radians(lat1), math.radians(lat2)
    dl = math.radians(lon2 - lon1)
    c = math.sin(p1) * math.sin(p2) + math.cos(p1) * math.cos(p2) * math.cos(dl)
    d = math.degrees(math.acos(max(-1.0, min(1.0, c))))
    az = math.degrees(math.atan2(math.sin(dl) * math.cos(p2),
                                 math.cos(p1) * math.sin(p2) - math.sin(p1) * math.cos(p2) * math.cos(dl)))
    return d, az % 360.0


def along_km(lat1, lon1, lat2, lon2, az):
    """Km from point 1 to point 2 along azimuth az (flat approximation, fine for < 100 km)."""
    dn = (lat2 - lat1) * KM_PER_DEG
    de = ((lon2 - lon1 + 540.0) % 360.0 - 180.0) * KM_PER_DEG * math.cos(math.radians((lat1 + lat2) / 2))
    a = math.radians(az)
    return dn * math.cos(a) + de * math.sin(a)


def parse_coord(s, kind):
    """'21.83S', '-21.83', '21 49.8 S', '138.92W', '-138.92' -> signed decimal degrees.
    French Polynesia is south and west, so an unsigned value without a letter is made negative."""
    s = str(s).strip().replace("°", " ").replace("'", " ").replace("′", " ")
    if not s:
        return None
    hemi = ""
    m = re.search(r"([NSEWnsew])\s*$", s)
    if m:
        hemi = m.group(1).upper()
        s = s[:m.start()].strip()
    parts = [p for p in re.split(r"[\s]+", s) if p]
    try:
        nums = [float(p) for p in parts]
    except ValueError:
        return None
    if not nums:
        return None
    sign = -1.0 if nums[0] < 0 else 1.0
    v = abs(nums[0]) + (nums[1] / 60.0 if len(nums) > 1 else 0.0) + (nums[2] / 3600.0 if len(nums) > 2 else 0.0)
    if hemi in ("S", "W") or sign < 0:
        return -v
    if hemi in ("N", "E"):
        return v
    # unsigned and no letter: this study area is S/W
    if kind == "lat" and 15.0 < v < 30.0:
        return -v
    if kind == "lon" and 130.0 < v < 150.0:
        return -v
    return v


def parse_date(s):
    s = str(s).strip()
    for rx, order in ((r"^(\d{4})[-/.](\d{1,2})[-/.](\d{1,2})", "ymd"),
                      (r"^(\d{1,2})[-/.](\d{1,2})[-/.](\d{4})", "dmy"),
                      (r"^(\d{4})(\d{2})(\d{2})$", "ymd")):
        m = re.match(rx, s)
        if m:
            g = [int(x) for x in m.groups()]
            y, mo, d = (g[0], g[1], g[2]) if order == "ymd" else (g[2], g[1], g[0])
            if 1 <= mo <= 12 and 1 <= d <= 31:
                return "%04d-%02d-%02d" % (y, mo, d)
    return None


def parse_s1(text):
    """Table S1 CSV -> list of dicts (date, time, lat, lon, mag, name). Header names are matched
    loosely (date, time/origin, lat, lon/long, mag/mb, name/event)."""
    rows = []
    rd = csv.DictReader(io.StringIO(text.strip()))
    keys = {k: (k or "").strip().lower() for k in (rd.fieldnames or [])}

    def col(r, *names):
        # exact header first, then header starting with the name, then header containing it
        for test in (lambda low, n: low == n, lambda low, n: low.startswith(n), lambda low, n: n in low):
            for n in names:
                for k, low in keys.items():
                    if test(low, n):
                        v = r.get(k)
                        if v is not None and str(v).strip() != "":
                            return str(v).strip()
        return ""

    for r in rd:
        d = parse_date(col(r, "date"))
        la = parse_coord(col(r, "lat"), "lat")
        lo = parse_coord(col(r, "lon"), "lon")
        if d is None or la is None or lo is None:
            continue
        mag = col(r, "mag", "mb")
        try:
            mag = float(mag)
        except ValueError:
            mag = None
        rows.append(dict(date=d, time=col(r, "time", "origin"), lat=la, lon=lo, mag=mag,
                         name=col(r, "name", "test")))
    return rows


def tod_seconds(t):
    """'HH:MM[:SS.s]' -> seconds of day, or None."""
    m = re.match(r"^\s*(\d{1,2}):(\d{2})(?::(\d{1,2}(?:\.\d*)?))?", str(t or ""))
    if not m:
        return None
    return int(m.group(1)) * 3600 + int(m.group(2)) * 60 + (float(m.group(3)) if m.group(3) else 0.0)


def _day(d):
    y, mo, dd = (int(x) for x in d.split("-"))
    return calendar.timegm((y, mo, dd, 0, 0, 0, 0, 0, 0)) // 86400


def match_events(s1, shafts):
    """Index into shafts for each S1 event (same UT date, else +/-1 day; nearest origin time
    when several), or None."""
    out = []
    sd = [_day(r["date"]) for r in shafts]
    for e in s1:
        de = _day(e["date"])
        best = None
        for off in (0, -1, 1):
            cands = [j for j, x in enumerate(sd) if x == de + off]
            if not cands:
                continue
            if len(cands) == 1:
                best = cands[0]
            else:
                te = tod_seconds(e.get("time"))
                if te is None:
                    best = None
                else:
                    best = min(cands, key=lambda j: abs((tod_seconds(shafts[j].get("time_ut")) or 1e9) - te))
            break
        out.append(best)
    return out


def form_pairs(events, tol_deg=PAIR_TOL_DEG):
    """All (i, j), i < j, whose catalogue epicentres are within tol_deg (Zhou's rule)."""
    pairs = []
    for i in range(len(events)):
        for j in range(i + 1, len(events)):
            d, _ = gc_dist_az(events[i]["lat"], events[i]["lon"], events[j]["lat"], events[j]["lon"])
            if d <= tol_deg + 1e-9:
                pairs.append((i, j))
    return pairs


def _connected(years, pairs, ref):
    adj = {y: set() for y in set(years)}
    for i, j in pairs:
        adj[years[i]].add(years[j])
        adj[years[j]].add(years[i])
    seen, todo = {ref}, [ref]
    while todo:
        y = todo.pop()
        for z in adj.get(y, ()):
            if z not in seen:
                seen.add(z)
                todo.append(z)
    return seen


def year_terms(x, years, pairs, ref=REF_YEAR):
    """Least-squares year terms a_y from pair differences x_j - x_i = a_{y_j} - a_{y_i}, with
    a_ref = 0. Pairs inside one year carry no year information and are skipped. Years not linked to
    ref through pairs are left out. x values that are None drop their pairs."""
    if ref not in years:
        return {}
    conn = _connected(years, [(i, j) for i, j in pairs if x[i] is not None and x[j] is not None], ref)
    free = sorted(y for y in conn if y != ref)
    col = {y: k for k, y in enumerate(free)}
    A, b = [], []
    for i, j in pairs:
        if x[i] is None or x[j] is None or years[i] == years[j]:
            continue
        if years[i] not in conn or years[j] not in conn:
            continue
        row = np.zeros(len(free))
        if years[j] != ref:
            row[col[years[j]]] += 1.0
        if years[i] != ref:
            row[col[years[i]]] -= 1.0
        A.append(row)
        b.append(x[j] - x[i])
    out = {ref: 0.0}
    if free and A:
        sol, *_ = np.linalg.lstsq(np.array(A), np.array(b), rcond=None)
        for y, k in col.items():
            out[y] = float(sol[k])
    return out


def epoch_terms(terms, years):
    """Epoch value = mean of year terms within the epoch, weighted by number of events that year."""
    out = {}
    for k, (a, bb) in EPOCHS.items():
        num = den = 0.0
        for y, v in terms.items():
            if a <= y <= bb:
                w = float(sum(1 for z in years if z == y))
                num += w * v
                den += w
        out[k] = (num / den) if den > 0 else None
    return out


def grade(pred, zhou=ZHOU):
    """Grade predicted epoch terms (s) against Zhou's reported ones (rules of the plan)."""
    if any(pred.get(k) is None for k in zhou):
        return dict(grade="untestable", signs_match=None, rms_res=None, rms_zhou=None)
    sm = sum(1 for k in zhou if abs(pred[k]) >= SIGN_MIN_S and math.copysign(1, pred[k]) == math.copysign(1, zhou[k]))
    res = [zhou[k] - pred[k] for k in zhou]
    rms_res = math.sqrt(sum(r * r for r in res) / len(res))
    rms_z = math.sqrt(sum(zhou[k] ** 2 for k in zhou) / len(zhou))
    if sm == 3 and rms_res <= EXPLAIN_RMS_S:
        g = "explains"
    elif sm >= 2 and rms_res < PARTLY_FRAC * rms_z:
        g = "partly"
    else:
        g = "does-not-explain"
    return dict(grade=g, signs_match=sm, rms_res=rms_res, rms_zhou=rms_z)


def headline(methods, grades):
    """Primary reading from what the paper says: 'none' -> R, 'catalogue-1D' -> C,
    anything else -> both reported, headline names both."""
    dc = (methods or {}).get("distance_correction", "not-stated")
    if dc == "none":
        return dict(reading="R", grade=grades["R"]["grade"], distance_correction=dc)
    if dc == "catalogue-1D":
        return dict(reading="C", grade=grades["C"]["grade"], distance_correction=dc)
    return dict(reading="R+C", grade="R:%s / C:%s" % (grades["R"]["grade"], grades["C"]["grade"]),
                distance_correction=dc)


def load_shafts(path):
    rows = []
    with open(path, encoding="utf-8") as f:
        for r in csv.DictReader(f):
            r["lat"], r["lon"] = float(r["lat"]), float(r["lon"])
            r["approx"] = int(r["approx"] or 0)
            rows.append(r)
    return rows


# ---------------------------------------------------------------- network / TauP parts
def http_text(url, timeout=60):
    import requests
    r = requests.get(url, timeout=timeout, headers={"User-Agent": "inner-core-repeaters/pkppair30 (research)"})
    if r.status_code == 404:
        return None
    r.raise_for_status()
    return r.text


def wait_inputs(wait_min, poll_s=600):
    t0 = time.time()
    while True:
        try:
            s1 = http_text(REPO_RAW + S1_PATH)
            mt = http_text(REPO_RAW + METHODS_PATH)
            if s1 and mt:
                return s1, mt
        except Exception as ex:  # noqa: BLE001
            log("fetch error", repr(ex)[:200])
        if time.time() - t0 > wait_min * 60:
            return None, None
        log("inputs not on main yet; waiting")
        time.sleep(poll_s)


def make_tt():
    from obspy.taup import TauPyModel
    from obspy.geodetics import locations2degrees
    model = TauPyModel("iasp91")
    cache = {}

    def tt(lat, lon):
        d = float(locations2degrees(lat, lon, BRVK[0], BRVK[1]))
        k = round(d, 6)
        if k not in cache:
            first = {}
            for a in model.get_travel_times(source_depth_in_km=SRC_DEPTH_KM, distance_in_degree=d,
                                            phase_list=["PKiKP", "PKIKP", "PP"]):
                first.setdefault(a.name, a.time)
            o = {"dist": d}
            for name, (p1, p2) in OBSERVABLES.items():
                o[name] = (first[p1] - first[p2]) if (p1 in first and p2 in first) else None
            cache[k] = o
        return cache[k]
    return tt


def analyse(s1, shafts, methods, tt):
    out = dict(n_s1=len(s1))
    idx = match_events(s1, shafts)
    out["n_matched"] = sum(1 for j in idx if j is not None)
    years = [int(e["date"][:4]) for e in s1]
    pairs = form_pairs(s1)
    out["n_pairs"] = len(pairs)
    out["pairs_vs_zhou_112"] = len(pairs) - 112
    out["years"] = {str(y): years.count(y) for y in sorted(set(years))}
    _, az = gc_dist_az(REF[0], REF[1], BRVK[0], BRVK[1])
    ev_rows = []
    for e, j, y in zip(s1, idx, years):
        r = dict(date=e["date"], name=e.get("name", ""), year=y, s1_lat=e["lat"], s1_lon=e["lon"],
                 mag=e.get("mag"), shaft_row=None)
        c = tt(e["lat"], e["lon"])
        r["s1_dist"] = c["dist"]
        r["s1_along_km"] = along_km(REF[0], REF[1], e["lat"], e["lon"], az)
        for name in OBSERVABLES:
            r["T_s1_" + name] = c[name]
        if j is not None:
            s = shafts[j]
            t = tt(s["lat"], s["lon"])
            r.update(shaft_row=j, shaft_name=s["name"], shaft=s.get("shaft", ""), atoll=s["atoll"],
                     approx=s["approx"], shaft_lat=s["lat"], shaft_lon=s["lon"],
                     shaft_along_km=along_km(REF[0], REF[1], s["lat"], s["lon"], az),
                     shaft_minus_s1_km=gc_dist_az(e["lat"], e["lon"], s["lat"], s["lon"])[0] * KM_PER_DEG)
            for name in OBSERVABLES:
                r["T_true_" + name] = t[name]
        ev_rows.append(r)
    out["events"] = ev_rows
    gates = dict(G0_matched=out["n_matched"] >= G0_MIN_MATCHED,
                 G1_epochs_present=(REF_YEAR in years) and all(any(a <= y <= b for y in years)
                                                                for a, b in EPOCHS.values()))
    out["gates"] = gates
    res = {}
    for name in OBSERVABLES:
        xr = [(r.get("T_true_" + name)) for r in ev_rows]
        xc = [((r["T_true_" + name] - r["T_s1_" + name]) if (r.get("T_true_" + name) is not None and
                                                           r.get("T_s1_" + name) is not None) else None)
              for r in ev_rows]
        xs = [r.get("T_s1_" + name) for r in ev_rows]
        block = {}
        for reading, x in (("R", xr), ("C", xc)):
            yt = year_terms(x, years, pairs)
            ep = epoch_terms(yt, [y for y, v in zip(years, x) if v is not None])
            g = grade(ep) if all(gates.values()) else dict(grade="untestable")
            block[reading] = dict(year_terms={str(k): v for k, v in sorted(yt.items())}, epochs=ep, **g)
            # report-only: rim shafts only (approx = 0)
            xr_rim = [v if (v is not None and ev_rows[i].get("approx") == 0) else None for i, v in enumerate(x)]
            yt2 = year_terms(xr_rim, years, pairs)
            block[reading]["rim_only_epochs"] = epoch_terms(yt2, [y for y, v in zip(years, xr_rim) if v is not None])
        # report-only: how big a PKP-PP difference the 0.05 deg tolerance itself allows (catalogue positions)
        dd = [abs(xs[j] - xs[i]) for i, j in pairs if xs[i] is not None and xs[j] is not None]
        block["pair_catalogue_dT_abs_s"] = dict(
            n=len(dd), median=(float(np.median(dd)) if dd else None),
            p90=(float(np.percentile(dd, 90)) if dd else None), max=(max(dd) if dd else None))
        grades = {k: block[k] for k in ("R", "C")}
        block["headline"] = headline(methods, grades) if all(gates.values()) else dict(grade="untestable")
        res[name] = block
    out["observables"] = res
    out["primary"] = dict(observable=PRIMARY_OBS, **res[PRIMARY_OBS]["headline"])
    return out


def _clean(o):
    if isinstance(o, float):
        return o if math.isfinite(o) else None
    if isinstance(o, dict):
        return {str(k): _clean(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [_clean(v) for v in o]
    return o


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--shafts", default="mururoa_shafts_1977_1996.csv")
    ap.add_argument("--out-dir", default="out")
    ap.add_argument("--wait-min", type=float, default=0.0)
    a = ap.parse_args(argv)
    os.makedirs(a.out_dir, exist_ok=True)
    t0 = time.time()
    summary = dict(version=VERSION, started=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()), errors=[],
                   zhou_reported=ZHOU, params=dict(pair_tol_deg=PAIR_TOL_DEG, sign_min_s=SIGN_MIN_S,
                                                   explain_rms_s=EXPLAIN_RMS_S, partly_frac=PARTLY_FRAC,
                                                   g0_min_matched=G0_MIN_MATCHED, src_depth_km=SRC_DEPTH_KM,
                                                   ref=REF, brvk=BRVK))
    with open(__file__, "rb") as f:
        summary["code_sha256"] = hashlib.sha256(f.read()).hexdigest()
    try:
        plan = http_text(REPO_RAW + PLAN_PATH)
        summary["plan_public_at_start"] = bool(plan)
        if plan:
            summary["plan_sha256"] = hashlib.sha256(plan.encode("utf-8")).hexdigest()
    except Exception as ex:  # noqa: BLE001
        summary["errors"].append("plan: %r" % ex)
    s1_text, mt_text = wait_inputs(a.wait_min)
    if not s1_text:
        summary["status"] = "no-input"
        summary["primary"] = dict(grade="no-input")
    else:
        summary["s1_sha256"] = hashlib.sha256(s1_text.encode("utf-8")).hexdigest()
        with open(os.path.join(a.out_dir, "zhou2026_table_s1_used.csv"), "w", encoding="utf-8") as f:
            f.write(s1_text)
        try:
            methods = json.loads(mt_text)
        except Exception as ex:  # noqa: BLE001
            methods = {}
            summary["errors"].append("methods json: %r" % ex)
        summary["methods"] = methods
        s1 = parse_s1(s1_text)
        if not os.path.exists(a.shafts):
            txt = http_text(REPO_RAW + SHAFTS_PATH)
            a.shafts = os.path.join(a.out_dir, "mururoa_shafts_1977_1996.csv")
            with open(a.shafts, "w", encoding="utf-8") as f:
                f.write(txt)
        with open(a.shafts, "rb") as f:
            summary["shafts_sha256"] = hashlib.sha256(f.read()).hexdigest()
        shafts = load_shafts(a.shafts)
        try:
            res = analyse(s1, shafts, methods, make_tt())
            events = res.pop("events")
            summary.update(res)
            summary["status"] = "graded"
            keys = []
            for r in events:
                for k in r:
                    if k not in keys:
                        keys.append(k)
            with open(os.path.join(a.out_dir, "pkppair30_events.csv"), "w", newline="", encoding="utf-8") as f:
                w = csv.DictWriter(f, fieldnames=keys)
                w.writeheader()
                for r in events:
                    w.writerow({k: (round(v, 6) if isinstance(v, float) else v) for k, v in r.items()})
        except Exception as ex:  # noqa: BLE001
            summary["status"] = "error"
            summary["errors"].append("analyse: %r" % ex)
    summary["finished"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    summary["runtime_min"] = round((time.time() - t0) / 60.0, 2)
    with open(os.path.join(a.out_dir, "pkppair30_summary.json"), "w") as f:
        json.dump(_clean(summary), f, indent=1)
    log("done", summary.get("status"), summary.get("primary"))
    return 0


if __name__ == "__main__":
    sys.exit(main())

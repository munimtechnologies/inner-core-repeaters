#!/usr/bin/env python3
# Copyright 2026 MV2, Munim, Inc.
# Licensed under the Apache License, Version 2.0 (Apache-2.0). See LICENSE.
# Made with AI agents working in MV2's name.
"""zw35g2.py - gate G2 of pre-registration 35: are the three multiplet-F earthquakes one source?

Pre-registered in findings/2026-11-24-prereg35-g2-source-check.md, written and published before
any mantle-P record of the three multiplet-F events was fetched, and before zw33-v1 (which
measures their core phases at YKA and ILAR) had run.

Pre-registration 35 grades rigid inner-core rotation by "closure": if YKA sees 2005-06-10 and
2022-05-19 as the same after 2012-12-03 looked different, ILAR must also see 2005 = 2022. Its
gate G2 says: "mantle-P source check (rule of pre-registration 11) ... A pair graded S
(different source) makes the grade untestable (source)." Pre-registration 11 ranks group splits
and needs >= 6 events, so multiplet F (3 events) needs a pair-level form of the same rule:

  * Measurement: mantlep2.py v0.13.0 unchanged (pinned commit), on the eight M03-24 events of
    pre-registration 11 plus the three multiplet-F events (11 events) at every IU/II/GT BHZ
    station 30-85 deg from 59.03 S 25.68 W (up to 60, nearest first): 1-2 Hz velocity, aligned
    on the P train, gates SNR_direct >= 5 and SNR_late >= 2, pair D_late = 1 - correlation in
    3-10 s after the stack onset.
  * Reference: at each station, the D_late of every gated M03-24 pair. These eight events are
    relocated repeats (<= 0.5 km apart, Zhang & Wen 2026) and pre-registration 11 graded their
    mantle P "N" (no split), so their D_late is what "same source" looks like at that station.
  * Station score for an F pair (a, c): q = share of reference pairs with D_late < D_late(a, c).
    Counted if both F events are gated and the station has >= 10 gated reference pairs.
    Hit (source differs) if q >= 0.93 (pre-registration 11's hit level p <= 0.07, mirrored).
  * Network decision per pair: pre-registration 11's rule. n counted stations, h hits:
    n < 3 untestable; h >= ceil(n/2) S (different source); h <= max(1, floor(0.2 n)) N (same
    source); otherwise inconclusive.
  * Graded pairs: (2005-06-10, 2012-12-03) and (2005-06-10, 2022-05-19), the two pairs the
    closure grade rests on. (2012-12-03, 2022-05-19) is report-only.
  * Mapping to pre-registration 35 (apply_g2): either graded pair S -> 'untestable (source)';
    both N -> the P35 grade is final and source-checked; otherwise the P35 grade is final but
    reported as 'source check incomplete'.

Top-level imports are standard library only; numpy/obspy/requests are imported inside run().
"""

import argparse
import csv
import json
import math
import os
import sys
import time
import traceback
import urllib.request

VERSION = "0.1.0"
NAN = float("nan")

REPO = "https://raw.githubusercontent.com/munimtechnologies/inner-core-repeaters/"
MANTLEP2_COMMIT = "b4712f0ebbc70b896231620c905ff53a4b8270fe"   # mantlep2.py 0.13.0 (as m34.py)
PREREG_URL = REPO + "main/findings/2026-11-24-prereg35-g2-source-check.md"

OLD8 = ("usp000bt5x", "usp000fz4d", "usp000hute", "usc000t4wj", "us7000am0e", "us6000h2wz",
        "us6000mk65", "us6000runq")
# USGS ComCat (fetched 2026-11-24): (id, origin UTC, lat, lon, depth km, mb)
F_EVENTS = [
    ("usp000dsnt", "2005-06-10T13:08:12.920", -58.706, -25.207, 9.8, 4.7),
    ("usp000jwcf", "2012-12-03T03:11:10.790", -58.864, -25.565, 41.2, 4.9),
    ("us6000hm6g", "2022-05-19T02:24:53.021", -58.796, -25.3186, 35.0, 5.1),
]
F2005, F2012, F2022 = "usp000dsnt", "usp000jwcf", "us6000hm6g"
GRADED_PAIRS = ((F2005, F2012), (F2005, F2022))
REPORT_PAIRS = ((F2012, F2022),)

HIT_Q = 0.93          # station hit if q >= 0.93
MIN_REF_PAIRS = 10    # >= 10 gated M03-24 pairs (>= 5 gated M03-24 events)
MIN_STATIONS = 3


# ----------------------------------------------------------------------------- pure helpers
def isf(x):
    try:
        return math.isfinite(float(x))
    except (TypeError, ValueError):
        return False


def key(a, b):
    return (a, b) if a <= b else (b, a)


def station_pair_score(D, gated, ref_ids, a, c, hit_q=HIT_Q, min_ref=MIN_REF_PAIRS):
    """D: {(e1, e2) sorted: D_late}; gated: {event: bool}. Returns dict with counted, q, hit."""
    out = dict(counted=False, q=None, hit=None, d_pair=None, n_ref=0, ref_median=None, ref_max=None)
    ref = [e for e in ref_ids if gated.get(e)]
    vals = []
    for i in range(len(ref)):
        for j in range(i + 1, len(ref)):
            v = D.get(key(ref[i], ref[j]))
            if isf(v):
                vals.append(float(v))
    out["n_ref"] = len(vals)
    if vals:
        s = sorted(vals)
        n = len(s)
        out["ref_median"] = s[n // 2] if n % 2 else 0.5 * (s[n // 2 - 1] + s[n // 2])
        out["ref_max"] = s[-1]
    if not (gated.get(a) and gated.get(c)):
        out["reason"] = "f_event_not_gated"
        return out
    d = D.get(key(a, c))
    if not isf(d):
        out["reason"] = "no_pair_D"
        return out
    out["d_pair"] = float(d)
    if len(vals) < min_ref:
        out["reason"] = "too_few_reference_pairs"
        return out
    q = sum(1 for v in vals if v < float(d)) / float(len(vals))
    out.update(counted=True, q=q, hit=bool(q >= hit_q - 1e-12))
    return out


def decide(hits, n_min=MIN_STATIONS):
    """Pre-registration 11's network rule (mantlep2.decide), on per-station hit flags."""
    n = len(hits)
    h = sum(1 for x in hits if x)
    if n < n_min:
        return "untestable", h
    if h >= int(math.ceil(n / 2.0)):
        return "S", h
    if h <= max(1, int(math.floor(0.2 * n))):
        return "N", h
    return "inconclusive", h


def apply_g2(p35_grade, pair_decisions):
    """p35_grade: the P35 headline from zw33-v1 (violated / holds / untestable...), or None.
    pair_decisions: {"2005-2012": S|N|inconclusive|untestable, "2005-2022": ...}.
    Returns (final grade, source status)."""
    graded = [pair_decisions.get("2005-2012", "untestable"), pair_decisions.get("2005-2022", "untestable")]
    if any(d == "S" for d in graded):
        return "untestable (source)", "different source"
    status = "source-checked" if all(d == "N" for d in graded) else "source check incomplete"
    return (p35_grade if p35_grade else "pending zw33-v1"), status


def label(pair):
    a, c = pair
    yr = {F2005: "2005", F2012: "2012", F2022: "2022"}
    return "%s-%s" % (yr[a], yr[c])


def grade_pairs(station_scores):
    """station_scores: {station: {pair_label: score dict}} -> {pair_label: decision dict}."""
    out = {}
    for pair in GRADED_PAIRS + REPORT_PAIRS:
        lab = label(pair)
        rows = [(st, sc[lab]) for st, sc in sorted(station_scores.items()) if sc.get(lab, {}).get("counted")]
        dec, h = decide([r["hit"] for _, r in rows])
        qs = sorted(r["q"] for _, r in rows)
        out[lab] = dict(decision=dec, n_stations=len(rows), hits=h, graded=pair in GRADED_PAIRS,
                        stations=[st for st, _ in rows], q=[round(r["q"], 4) for _, r in rows],
                        median_q=(qs[len(qs) // 2] if len(qs) % 2 else 0.5 * (qs[len(qs) // 2 - 1] + qs[len(qs) // 2]))
                        if qs else None)
    return out


# ----------------------------------------------------------------------------- I/O
def fetch(url, tries=4, wait=10):
    for k in range(tries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "inner-core-repeaters (MV2, Munim, Inc.)"})
            with urllib.request.urlopen(req, timeout=60) as r:
                if r.status == 200:
                    return r.read().decode("utf-8", "replace")
        except Exception as ex:  # noqa: BLE001
            print("fetch error", url, repr(ex)[:200], flush=True)
        time.sleep(wait * (k + 1))
    return None


def load_mantlep2(dest):
    os.makedirs(dest, exist_ok=True)
    p = os.path.join(dest, "mantlep2.py")
    if not os.path.exists(p):
        body = fetch(REPO + MANTLEP2_COMMIT + "/mantlep2.py")
        if not body or 'VERSION = "0.13.0"' not in body:
            return None
        compile(body, "mantlep2.py", "exec")
        with open(p, "w") as f:
            f.write(body)
    if dest not in sys.path:
        sys.path.insert(0, dest)
    import mantlep2  # noqa: E402
    return mantlep2


def write_csv(path, rows):
    keys = []
    for r in rows:
        for k in r:
            if k not in keys:
                keys.append(k)
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=keys or ["empty"])
        w.writeheader()
        w.writerows(rows)


def read_p35(path):
    """Optional: P35 headline from zw33-v1's summary (key p35.grade or p35_grade), else None."""
    if not path or not os.path.exists(path):
        return None
    try:
        with open(path, encoding="utf-8") as f:
            s = json.load(f)
    except Exception:  # noqa: BLE001
        return None
    p = s.get("p35")
    if isinstance(p, dict):
        return p.get("grade") or p.get("headline")
    return s.get("p35_grade")


def run(args):
    t0 = time.time()
    os.makedirs(args.out_dir, exist_ok=True)
    s = dict(version=VERSION, started=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
             mantlep2_commit=MANTLEP2_COMMIT, hit_q=HIT_Q, min_ref_pairs=MIN_REF_PAIRS,
             min_stations=MIN_STATIONS, f_events=F_EVENTS)
    s["prereg_public_at_start"] = bool(fetch(PREREG_URL, tries=2, wait=5))
    mp = load_mantlep2("/tmp/zw35g2_pinned")
    if mp is None:
        s["error"] = "mantlep2.py v0.13.0 download failed"
        with open(os.path.join(args.out_dir, "zw35g2_summary.json"), "w") as f:
            json.dump(s, f, indent=1)
        return 2
    if [e[0] for e in mp.EVENTS] != list(OLD8):
        s["error"] = "pinned mantlep2 EVENTS differ from the eight M03-24 events"
        with open(os.path.join(args.out_dir, "zw35g2_summary.json"), "w") as f:
            json.dump(s, f, indent=1)
        return 2
    mp.EVENTS = list(mp.EVENTS) + list(F_EVENTS)
    mp.IDS = [e[0] for e in mp.EVENTS]
    mp.DATES = {e[0]: e[1][:10] for e in mp.EVENTS}

    import requests
    from obspy.taup import TauPyModel
    session = requests.Session()
    session.headers["User-Agent"] = "inner-core-repeaters (MV2, Munim, Inc.)"
    found = None
    for _ in range(3):
        found = mp.discover(session)
        if found:
            break
        time.sleep(30)
    source = "station-service"
    if not found:
        found, source = mp.FALLBACK_STATIONS, "fallback"
    cand = []
    for net, sta, lat, lon in found:
        d = mp.gc_deg(mp.SRC_LAT, mp.SRC_LON, lat, lon)
        if mp.DIST_RANGE[0] <= d <= mp.DIST_RANGE[1]:
            cand.append((d, net, sta, lat, lon))
    cand.sort()
    cand = cand[:mp.MAX_STATIONS]
    s.update(station_source=source, n_candidate_stations=len(cand))

    cache = os.path.join(args.data_dir, "mantlep")
    plots = os.path.join(args.out_dir, "plots")
    os.makedirs(cache, exist_ok=True)
    os.makedirs(plots, exist_ok=True)
    model = TauPyModel("ak135")
    ev_all, pair_all, stations, scores, skipped = [], [], [], {}, []
    for d, net, sta, lat, lon in cand:
        name = "%s.%s" % (net, sta)
        if (time.time() - t0) / 60.0 > args.time_budget_min:
            skipped.append(name)
            continue
        ev_rows, pair_rows = [], []
        try:
            res = mp.run_station(session, model, net, sta, lat, lon, cache, plots, ev_rows, pair_rows, {})
        except Exception as ex:  # noqa: BLE001
            traceback.print_exc()
            res = dict(station=name, status="exception", error=repr(ex)[:300])
        gated = {r["event"]: str(r.get("gated")) == "True" for r in ev_rows}
        D = {key(r["ev1"], r["ev2"]): float(r["D_late"]) for r in pair_rows if isf(r.get("D_late"))}
        sc = {}
        for pair in GRADED_PAIRS + REPORT_PAIRS:
            sc[label(pair)] = station_pair_score(D, gated, OLD8, pair[0], pair[1])
        scores[name] = sc
        stations.append(dict(station=name, dist_deg=round(d, 2), status=res.get("status"),
                             ref=res.get("ref"), n_gated_old=sum(1 for e in OLD8 if gated.get(e)),
                             f_gated={e: bool(gated.get(e)) for e in (F2005, F2012, F2022)},
                             scores=sc))
        ev_all.extend(ev_rows)
        pair_all.extend(pair_rows)
        print(name, json.dumps(sc)[:800], flush=True)

    g = grade_pairs(scores)
    p35 = read_p35(args.p35_summary)
    final, status = apply_g2(p35, {k: v["decision"] for k, v in g.items()})
    s.update(stations=stations, skipped_time_budget=skipped, pairs=g,
             g2=dict(p35_grade_in=p35, final_grade=final, source_status=status),
             finished=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
             minutes=round((time.time() - t0) / 60.0, 1))
    write_csv(os.path.join(args.out_dir, "zw35g2_events.csv"), ev_all)
    write_csv(os.path.join(args.out_dir, "zw35g2_pairs.csv"), pair_all)
    rows = []
    for st, sc in sorted(scores.items()):
        for lab, v in sc.items():
            r = dict(station=st, pair=lab)
            r.update({k: v.get(k) for k in ("counted", "q", "hit", "d_pair", "n_ref", "ref_median", "ref_max", "reason")})
            rows.append(r)
    write_csv(os.path.join(args.out_dir, "zw35g2_station_scores.csv"), rows)
    with open(os.path.join(args.out_dir, "zw35g2_summary.json"), "w") as f:
        json.dump(mp.jsonable(s), f, indent=1)
    print(json.dumps(mp.jsonable(dict(pairs=g, g2=s["g2"])))[:3000], flush=True)
    return 0


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--data-dir", default="/data")
    ap.add_argument("--out-dir", default="/out")
    ap.add_argument("--time-budget-min", type=float, default=100.0)
    ap.add_argument("--p35-summary", default="", help="optional zw33-v1 summary JSON")
    return run(ap.parse_args(argv))


if __name__ == "__main__":
    sys.exit(main())

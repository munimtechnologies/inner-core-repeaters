#!/usr/bin/env python3
# Copyright 2026 MV2, Munim, Inc.
# Licensed under the Apache License, Version 2.0 (Apache-2.0). See LICENSE.
# Made with AI agents working in MV2's name.
"""chilepower.py - power check for a Chile -> Makanchi (MKAR) inner-core clock test.

The Chile pilot (findings/2026-11-20-chile-path-plan.md, graded in
findings/2026-11-20-chile-pilot-results.md) closed the Chile path under its coded rule, partly
because MKAR's "before 2008.5" availability counted 1995-2000 earthquakes that predate the array.
Before any further waveform work, this script asks a catalogue-only question, with its decision
rule fixed in findings/2026-11-20-chile-power-plan.md before it runs:

  Inside MKAR's open-data era (2000-08-28 to 2025-08-15) and distance gate (143-165 deg), are
  there enough candidate repeater pairs that straddle the 2008.5 turn, and do they cover both
  clock-matched ("near") and clock-mismatched ("far") pairings?

Clock: if the inner core retraces its pre-2008.5 motion 2.5 times more slowly (Wang, Vidale et
al. 2024), an event at t2 >= 2008.5 reoccupies the position of t1* = 2008.5 - (t2 - 2008.5) / 2.5.
The clock residual of a pair (t1 < 2008.5 <= t2) is r = t1 - t1*.

Stdlib only. Usage:
  python chilepower.py --catalogue results/chilepilot-v2/chile_catalogue.csv --out-dir out
  python chilepower.py --catalogue-url <raw URL> --plan-url <raw URL> --out-dir /out --wait-min 120
"""

import argparse
import csv
import json
import math
import os
import sys
import time
import urllib.request

VERSION = "0.1.0"
TURN = 2008.5
RATE = 2.5                          # backward motion 2.5x slower (Wang et al. 2024)
MKAR_LAT, MKAR_LON = 46.7698, 82.2983   # mean of KZ MK01-MK09 (fdsnws station, 2026-11-20)
GATE = (143.0, 165.0)               # chilepilot.py distance range
ERA0_ISO = "2000-08-28T00:00:00"    # MKAR SHZ epoch start
ERA1_ISO = "2025-08-15T00:00:00"    # last date a dataselect probe returned MKAR data
EPS = 1e-9

# Candidate pair (pre-registered): both M >= 5.0, epicentres <= 10 km, |dM| <= 0.3, |dz| <= 20 km.
MAIN = dict(name="main", km=10.0, dm=0.3, dz=20.0, min_mag=5.0)
SENS = [dict(name="pilot-proxy", km=5.0, dm=0.2, dz=10.0, min_mag=5.0),
        dict(name="loose-20km", km=20.0, dm=0.3, dz=20.0, min_mag=5.0),
        dict(name="main-m5.5", km=10.0, dm=0.3, dz=20.0, min_mag=5.5)]
NEAR, FAR = 2.0, 4.0                # |r| <= 2 yr near; |r| >= 4 yr far
NEED_PAIRS, NEED_NEAR, NEED_FAR = 60, 15, 15

# chilepilot-v2 close_pairs() on the full catalogue (chile_summary.json), for the reproduction check
PILOT_CLOSE = dict(n_events=2237, pairs=172, straddling=24, straddling_m55=2, events_in_straddling=45)


def log(*a):
    print(time.strftime("%H:%M:%S"), *a, flush=True)


def isf(x):
    return isinstance(x, (int, float)) and not isinstance(x, bool) and math.isfinite(x)


def fnum(x):
    try:
        v = float(x)
    except (TypeError, ValueError):
        return float("nan")
    return v if math.isfinite(v) else float("nan")


def gc_deg(lat1, lon1, lat2, lon2):
    p1, p2 = math.radians(lat1), math.radians(lat2)
    dl = math.radians(lon2 - lon1)
    c = math.sin(p1) * math.sin(p2) + math.cos(p1) * math.cos(p2) * math.cos(dl)
    return math.degrees(math.acos(max(-1.0, min(1.0, c))))


def km_between(lat1, lon1, lat2, lon2):
    return 6371.0 * math.radians(gc_deg(lat1, lon1, lat2, lon2))


def decimal_year(iso):
    from datetime import datetime
    s = str(iso).rstrip("Z")[:19]
    d = datetime.strptime(s, "%Y-%m-%dT%H:%M:%S")
    y0 = datetime(d.year, 1, 1)
    y1 = datetime(d.year + 1, 1, 1)
    return d.year + (d - y0).total_seconds() / (y1 - y0).total_seconds()


ERA0 = decimal_year(ERA0_ISO)
ERA1 = decimal_year(ERA1_ISO)


def mkar_dist(lat, lon):
    return gc_deg(MKAR_LAT, MKAR_LON, lat, lon)


def reoccupied(t2, turn=TURN, rate=RATE):
    """Pre-turn time whose inner-core position an event at t2 (>= turn) reoccupies."""
    return turn - (t2 - turn) / rate


def clock_residual(t1, t2, turn=TURN, rate=RATE):
    return t1 - reoccupied(t2, turn, rate)


def close_pairs(rows, max_km=5.0, max_dm=0.2, max_dz=10.0, min_mag=5.0, turn=TURN):
    """Copied unchanged from chilepilot.py (3df9ac8); used only to reproduce the pilot's count."""
    ev = [r for r in rows if isf(r.get("mag")) and r["mag"] >= min_mag and isf(r.get("lat"))
          and isf(r.get("lon")) and isf(r.get("t"))]
    res = dict(n_events=len(ev), pairs=0, straddling=0, straddling_m55=0, events_in_straddling=0)
    used = set()
    deg = max_km / 111.0
    for i in range(len(ev)):
        a = ev[i]
        for j in range(i + 1, len(ev)):
            b = ev[j]
            if abs(a["lat"] - b["lat"]) > deg or abs(a["mag"] - b["mag"]) > max_dm:
                continue
            za, zb = a.get("depth"), b.get("depth")
            if isf(za) and isf(zb) and abs(za - zb) > max_dz:
                continue
            if km_between(a["lat"], a["lon"], b["lat"], b["lon"]) > max_km:
                continue
            res["pairs"] += 1
            if (a["t"] < turn) != (b["t"] < turn):
                res["straddling"] += 1
                used.add(i)
                used.add(j)
                if min(a["mag"], b["mag"]) >= 5.5:
                    res["straddling_m55"] += 1
    res["events_in_straddling"] = len(used)
    return res


def eligible(rows, era0=ERA0, era1=ERA1, gate=GATE):
    """Events inside MKAR's open-data era and distance gate, with finite t/lat/lon/mag."""
    out = []
    for r in rows:
        if not (isf(r.get("t")) and isf(r.get("lat")) and isf(r.get("lon")) and isf(r.get("mag"))):
            continue
        if r["t"] < era0 or r["t"] > era1:
            continue
        d = mkar_dist(r["lat"], r["lon"])
        if d < gate[0] or d > gate[1]:
            continue
        out.append(dict(r, dist=d))
    return out


def straddling_pairs(ev, km, dm, dz, min_mag, turn=TURN):
    """Every (before, after) pair with both M >= min_mag, epicentres <= km, |dM| <= dm and
    |dz| <= dz (a missing depth passes, as in chilepilot.close_pairs)."""
    sel = [r for r in ev if r["mag"] >= min_mag - EPS]
    before = [r for r in sel if r["t"] < turn]
    after = [r for r in sel if r["t"] >= turn]
    out = []
    for a in before:
        for b in after:
            if abs(a["lat"] - b["lat"]) > km / 111.0 + 0.01:
                continue
            if abs(a["mag"] - b["mag"]) > dm + EPS:
                continue
            za, zb = a.get("depth"), b.get("depth")
            if isf(za) and isf(zb) and abs(za - zb) > dz + EPS:
                continue
            d = km_between(a["lat"], a["lon"], b["lat"], b["lon"])
            if d > km + EPS:
                continue
            out.append(dict(id_before=a["id"], time_before=a["time"], t_before=a["t"], mag_before=a["mag"],
                            depth_before=a.get("depth"), id_after=b["id"], time_after=b["time"],
                            t_after=b["t"], mag_after=b["mag"], depth_after=b.get("depth"),
                            km=d, dm=abs(a["mag"] - b["mag"]),
                            r=clock_residual(a["t"], b["t"], turn)))
    out.sort(key=lambda p: (p["t_after"], p["t_before"]))
    return out


def summarise(pairs):
    near = sum(1 for p in pairs if abs(p["r"]) <= NEAR + EPS)
    far = sum(1 for p in pairs if abs(p["r"]) >= FAR - EPS)
    ids = set()
    for p in pairs:
        ids.add(p["id_before"])
        ids.add(p["id_after"])
    rs = sorted(abs(p["r"]) for p in pairs)
    med = float("nan")
    if rs:
        n = len(rs)
        med = rs[n // 2] if n % 2 else 0.5 * (rs[n // 2 - 1] + rs[n // 2])
    by_year = {}
    for p in pairs:
        y = str(int(math.floor(p["t_after"])))
        by_year[y] = by_year.get(y, 0) + 1
    rbins = {}
    for p in pairs:
        lo = 2 * int(math.floor(p["r"] / 2.0))
        k = "%d..%d" % (lo, lo + 2)
        rbins[k] = rbins.get(k, 0) + 1
    return dict(pairs=len(pairs), near=near, far=far, events=len(ids),
                median_abs_r=med if isf(med) else None,
                after_year_counts=dict(sorted(by_year.items())),
                r_bins=dict(sorted(rbins.items(), key=lambda kv: int(kv[0].split("..")[0]))))


def decide(main):
    ok = main["pairs"] >= NEED_PAIRS and main["near"] >= NEED_NEAR and main["far"] >= NEED_FAR
    reason = "%d candidate pairs (need %d), %d near (need %d), %d far (need %d)" % (
        main["pairs"], NEED_PAIRS, main["near"], NEED_NEAR, main["far"], NEED_FAR)
    return ("mkar-multiplet-pilot" if ok else "close-chile-path"), reason


def analyse(rows):
    res = dict(version=VERSION, turn=TURN, rate=RATE, mkar_centre=[MKAR_LAT, MKAR_LON], gate=list(GATE),
               era=[ERA0_ISO, ERA1_ISO], era_years=[ERA0, ERA1], near_yr=NEAR, far_yr=FAR,
               need=dict(pairs=NEED_PAIRS, near=NEED_NEAR, far=NEED_FAR), n_catalogue=len(rows))
    rep = close_pairs(rows)
    res["reproduce_pilot"] = dict(got=rep, expected=PILOT_CLOSE, ok=(rep == PILOT_CLOSE))
    ev = eligible(rows)
    res["n_eligible"] = len(ev)
    for m in (5.0, 5.5):
        sel = [r for r in ev if r["mag"] >= m - EPS]
        res["eligible_m%.1f" % m] = dict(n=len(sel), before=sum(1 for r in sel if r["t"] < TURN),
                                         after=sum(1 for r in sel if r["t"] >= TURN))
    res["definitions"] = {}
    res["pairs"] = {}
    for d in [MAIN] + SENS:
        p = straddling_pairs(ev, d["km"], d["dm"], d["dz"], d["min_mag"])
        res["definitions"][d["name"]] = d
        res["pairs"][d["name"]] = p
        res[d["name"]] = summarise(p)
    res["decision"], res["decision_reason"] = decide(res["main"])
    m = res["main"]
    res["forecasts"] = dict(G1_pairs_ge_60=m["pairs"] >= NEED_PAIRS,
                            G2_near_ge_15=m["near"] >= NEED_NEAR,
                            G3_decision=res["decision"],
                            G4_median_abs_r_le_4=(m["median_abs_r"] is not None and m["median_abs_r"] <= 4.0),
                            G5_proxy_le_half_of_pilot=res["pilot-proxy"]["pairs"] <= PILOT_CLOSE["straddling"] / 2.0)
    return res


def read_catalogue(path):
    """All rows, as chilepilot.load_catalogue() reads them (no event-type filter)."""
    rows = []
    with open(path, newline="", encoding="utf-8") as f:
        for r in csv.DictReader(f):
            t = r.get("time", "")
            try:
                ty = decimal_year(t)
            except ValueError:
                continue
            rows.append(dict(id=r.get("id", ""), time=t, t=ty, lat=fnum(r.get("latitude")),
                             lon=fnum(r.get("longitude")), depth=fnum(r.get("depth")),
                             mag=fnum(r.get("mag")), magType=r.get("magType", "")))
    return rows


def fetch(url, dest, wait_min, prefix=b"time,"):
    deadline = time.time() + 60.0 * max(0.0, wait_min)
    while True:
        try:
            with urllib.request.urlopen(url, timeout=60) as r:
                data = r.read()
            if data.startswith(prefix):
                with open(dest, "wb") as f:
                    f.write(data)
                return True
            log("not ready yet", len(data), "bytes")
        except Exception as e:  # noqa: BLE001
            log("fetch failed", repr(e)[:200])
        if time.time() > deadline:
            return False
        time.sleep(120)


def write_pairs_csv(res, path):
    cols = ["definition", "id_before", "time_before", "mag_before", "depth_before", "id_after",
            "time_after", "mag_after", "depth_after", "km", "dm", "r"]
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(cols)
        for name, lst in res["pairs"].items():
            for p in lst:
                row = [name]
                for c in cols[1:]:
                    v = p[c]
                    row.append(round(v, 3) if isinstance(v, float) and math.isfinite(v) else v)
                w.writerow(row)


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--catalogue", default="")
    ap.add_argument("--catalogue-url", default="")
    ap.add_argument("--plan-url", default="")
    ap.add_argument("--out-dir", default="out")
    ap.add_argument("--wait-min", type=float, default=0.0)
    a = ap.parse_args(argv)
    os.makedirs(a.out_dir, exist_ok=True)
    plan_public = None
    if a.plan_url:
        plan_public = fetch(a.plan_url, os.path.join(a.out_dir, "plan.md"), a.wait_min, prefix=b"#")
        log("plan public:", plan_public)
    path = a.catalogue
    if not path:
        path = os.path.join(a.out_dir, "chile_catalogue.csv")
        if not fetch(a.catalogue_url, path, a.wait_min):
            log("catalogue unavailable")
            return 2
    rows = read_catalogue(path)
    res = analyse(rows)
    res["catalogue"] = a.catalogue_url or a.catalogue
    res["plan_public_at_start"] = plan_public
    write_pairs_csv(res, os.path.join(a.out_dir, "chilepower_pairs.csv"))
    out = {k: v for k, v in res.items() if k != "pairs"}
    with open(os.path.join(a.out_dir, "chilepower_summary.json"), "w", encoding="utf-8") as f:
        json.dump(out, f, indent=1)
    log("reproduce pilot", json.dumps(res["reproduce_pilot"]))
    log("eligible", res["n_eligible"], json.dumps(res["eligible_m5.0"]), json.dumps(res["eligible_m5.5"]))
    for d in [MAIN] + SENS:
        s = res[d["name"]]
        log(d["name"], "pairs", s["pairs"], "near", s["near"], "far", s["far"], "events", s["events"],
            "median|r|", s["median_abs_r"])
    log("decision", res["decision"], res["decision_reason"])
    log("forecasts", json.dumps(res["forecasts"]))
    return 0


if __name__ == "__main__":
    sys.exit(main())

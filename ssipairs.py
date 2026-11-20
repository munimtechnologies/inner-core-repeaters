#!/usr/bin/env python3
# Copyright 2026 MV2, Munim, Inc.
# Licensed under the Apache License, Version 2.0 (Apache-2.0). See LICENSE.
# Made with AI agents working in MV2's name.
"""ssipairs.py - candidate repeater pairs among large South Sandwich Islands earthquakes.

The SNR pilot (findings/2026-11-19-snr-pilot-plan.md, graded in
findings/2026-11-20-snr-pilot-results.md) returned `testable-from-M6.0`. Its decision rule then
allows a near-antipodal clock test only if the catalogue holds at least 10 candidate pairs with
one event on each side of 2008.5. The plan did not define "candidate pair"; we use the close-pair
proxy fixed in chilepilot.py (commit 3df9ac8, before this pilot's data returned): epicentres
within 5 km, |dM| <= 0.2, |dz| <= 10 km. `close_pairs` below is copied unchanged from
chilepilot.py, except that min_mag is passed in. A looser sensitivity table (10 km, with and
without the magnitude/depth limits) is reported beside it.

Stdlib only. Usage:
  python ssipairs.py --catalogue results/snrpilot-v2/snr_catalogue.csv --out-dir out
  python ssipairs.py --catalogue-url <raw URL> --out-dir /out [--wait-min 60]
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
PAIR_KM, PAIR_DM, PAIR_DZ = 5.0, 0.2, 10.0      # chilepilot.py proxy
NEED = 10                                       # snr-pilot plan: pairs needed at M >= 6.0
MIN_MAGS = (6.0, 5.5)
LOOSE_KM = 10.0


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


def close_pairs(rows, max_km=PAIR_KM, max_dm=PAIR_DM, max_dz=PAIR_DZ, min_mag=5.0, turn=TURN):
    """Copied from chilepilot.py (3df9ac8). Catalogue pairs (both M >= min_mag) within max_km
    epicentral distance, |dM| <= max_dm and |dz| <= max_dz. Returns counts overall and
    straddling the turn, and for M >= 5.5."""
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


def straddling_list(rows, max_km, min_mag, turn=TURN):
    """Every pair (both M >= min_mag, one each side of the turn) within max_km, with no
    magnitude or depth limit; sorted by distance."""
    ev = [r for r in rows if isf(r.get("mag")) and r["mag"] >= min_mag and isf(r.get("lat"))
          and isf(r.get("lon")) and isf(r.get("t"))]
    before = [r for r in ev if r["t"] < turn]
    after = [r for r in ev if r["t"] >= turn]
    out = []
    for a in before:
        for b in after:
            if abs(a["lat"] - b["lat"]) > max_km / 111.0:
                continue
            d = km_between(a["lat"], a["lon"], b["lat"], b["lon"])
            if d > max_km:
                continue
            dz = abs(a["depth"] - b["depth"]) if isf(a.get("depth")) and isf(b.get("depth")) else float("nan")
            out.append(dict(id_before=a["id"], time_before=a["time"], mag_before=a["mag"],
                            depth_before=a.get("depth"), id_after=b["id"], time_after=b["time"],
                            mag_after=b["mag"], depth_after=b.get("depth"), km=round(d, 2),
                            dm=round(abs(a["mag"] - b["mag"]), 2),
                            dz=round(dz, 1) if isf(dz) else "",
                            meets_dm_dz=bool(abs(a["mag"] - b["mag"]) <= PAIR_DM + 1e-9
                                             and (not isf(dz) or dz <= PAIR_DZ + 1e-9))))
    out.sort(key=lambda p: p["km"])
    return out


def read_catalogue(path):
    rows = []
    with open(path, newline="", encoding="utf-8") as f:
        for r in csv.DictReader(f):
            if (r.get("type") or "earthquake") != "earthquake":
                continue
            t = r.get("time", "")
            try:
                ty = decimal_year(t)
            except ValueError:
                continue
            rows.append(dict(id=r.get("id", ""), time=t, t=ty, lat=fnum(r.get("latitude")),
                             lon=fnum(r.get("longitude")), depth=fnum(r.get("depth")),
                             mag=fnum(r.get("mag")), magType=r.get("magType", "")))
    return rows


def analyse(rows):
    res = dict(version=VERSION, turn=TURN, proxy=dict(km=PAIR_KM, dm=PAIR_DM, dz=PAIR_DZ),
               need_pairs_m6=NEED, n_catalogue=len(rows), by_min_mag={})
    for m in MIN_MAGS:
        sel = [r for r in rows if isf(r["mag"]) and r["mag"] >= m]
        lst = straddling_list(rows, LOOSE_KM, m)
        res["by_min_mag"]["%.1f" % m] = dict(
            n_events=len(sel),
            before_turn=sum(1 for r in sel if r["t"] < TURN),
            after_turn=sum(1 for r in sel if r["t"] >= TURN),
            proxy=close_pairs(rows, min_mag=m),
            straddling_within_5km_any=sum(1 for p in lst if p["km"] <= PAIR_KM),
            straddling_within_10km_any=len(lst),
            straddling_within_10km_dm_dz=sum(1 for p in lst if p["meets_dm_dz"]),
            pairs_within_10km=lst)
    s = res["by_min_mag"]["6.0"]["proxy"]["straddling"]
    res["decision"] = ("clock-test" if s >= NEED else "close-near-antipodal-line")
    res["decision_reason"] = "%d straddling M>=6.0 proxy pairs; %d needed" % (s, NEED)
    return res


def fetch(url, dest, wait_min):
    deadline = time.time() + 60.0 * max(0.0, wait_min)
    while True:
        try:
            with urllib.request.urlopen(url, timeout=60) as r:
                data = r.read()
            if data.startswith(b"time,"):
                with open(dest, "wb") as f:
                    f.write(data)
                return True
            log("not a catalogue yet", len(data), "bytes")
        except Exception as e:  # noqa: BLE001
            log("fetch failed", repr(e)[:200])
        if time.time() > deadline:
            return False
        time.sleep(120)


def write_pairs_csv(res, path):
    cols = ["min_mag", "id_before", "time_before", "mag_before", "depth_before", "id_after",
            "time_after", "mag_after", "depth_after", "km", "dm", "dz", "meets_dm_dz"]
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(cols)
        for m, d in res["by_min_mag"].items():
            for p in d["pairs_within_10km"]:
                w.writerow([m] + [p[c] for c in cols[1:]])


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--catalogue", default="")
    ap.add_argument("--catalogue-url", default="")
    ap.add_argument("--out-dir", default="out")
    ap.add_argument("--wait-min", type=float, default=0.0)
    a = ap.parse_args(argv)
    os.makedirs(a.out_dir, exist_ok=True)
    path = a.catalogue
    if not path:
        path = os.path.join(a.out_dir, "snr_catalogue.csv")
        if not fetch(a.catalogue_url, path, a.wait_min):
            log("catalogue unavailable")
            return 2
    rows = read_catalogue(path)
    res = analyse(rows)
    res["catalogue"] = a.catalogue_url or a.catalogue
    with open(os.path.join(a.out_dir, "ssipairs_summary.json"), "w", encoding="utf-8") as f:
        json.dump(res, f, indent=1)
    write_pairs_csv(res, os.path.join(a.out_dir, "ssipairs_pairs.csv"))
    for m, d in res["by_min_mag"].items():
        log("M>=%s" % m, d["n_events"], "events;", "proxy", json.dumps(d["proxy"]),
            "| 5 km any", d["straddling_within_5km_any"], "| 10 km any", d["straddling_within_10km_any"],
            "| 10 km dm/dz", d["straddling_within_10km_dm_dz"])
    log("decision", res["decision"], res["decision_reason"])
    return 0


if __name__ == "__main__":
    sys.exit(main())

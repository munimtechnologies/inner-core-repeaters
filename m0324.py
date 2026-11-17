#!/usr/bin/env python3
"""Patch test: the M03-24 repeating-earthquake patch, 1991-2026, at the ILAR array.

Part of inner-core-repeaters by MV2, Munim, Inc. (made with AI agents). Licence: Apache-2.0.

Zhang & Wen (2026, JGR Solid Earth 131 e2026JB034595; data CC BY 4.0 at
doi:10.5281/zenodo.19625711) mapped inner-core surface change with one South Sandwich
Islands multiplet, "M03-24": six events (2003-03-21, 2008-02-17, 2011-02-20, 2020-07-11,
2022-03-07, 2024-03-21) recorded by ILAR and YKA. They argue the changes need no
inner-core rotation. Wang et al. (2024, Nature 631:340) argue the inner core turned at
2008.5 and is retracing its path 2.5x more slowly, so an event at time t > 2008.5 should
best match a partner from the equivalent date e(t) = 2008.5 - (t - 2008.5)/2.5; for
2025-2026 events that is 2001.5-2001.9.

The same patch (within 35 km of the M03-24 centroid) had USGS events on 2001-01-29 and
2002-06-29 and more than 40 M>=4.3 events in 2025-2026, which no published study has
used. The main ILAR run (run.py) uses M>=5.0 before mid-2024 and M>=4.7 after, so it misses
both 2001-2002 partners, two of the six M03-24 members and 18 of the 46 new events.
This script runs the same objective ILAR pipeline (run.py) on every M>=4.3 event in the
patch and grades the pre-registered tests in findings/2026-11-17-m0324-prereg.md:
  Z0  our pipeline independently re-finds the published M03-24 multiplet
  Z1  new-event best match (model.p3) on the patch: B pick (~2001-2002) vs M pick (latest)
  Z2  out-of-sample rank test (model.grade_set) on the patch's 2025+ pairs
  Z3  calibration (model.p4) on the patch
  Z5  exploratory turn-time fit (model.fit_turn)
"""
import argparse
import csv
import io
import json
import math
import os
import time
import traceback
from types import SimpleNamespace

import numpy as np

import model

VERSION = "0.5.0"
# Zhang & Wen (2026) M03-24 members, with USGS ComCat ids/locations (looked up 2026-11-17).
M0324 = [
    ("usp000bt5x", "2003-03-21T10:47:47.620Z", -59.071, -25.470, 33.0, 4.8),
    ("usp000fz4d", "2008-02-17T19:32:09.100Z", -58.986, -25.697, 66.5, 4.9),
    ("usp000hute", "2011-02-20T17:37:37.320Z", -58.979, -25.857, 10.0, 5.0),
    ("us7000am0e", "2020-07-11T18:58:06.501Z", -59.0899, -25.7164, 57.13, 5.3),
    ("us6000h2wz", "2022-03-07T23:25:40.208Z", -59.0811, -25.7858, 35.0, 5.1),
    ("us6000mk65", "2024-03-21T02:47:21.923Z", -59.0115, -25.6366, 43.89, 5.2),
]
USGS = "https://earthquake.usgs.gov/fdsnws/event/1/query"
PROC_TAG = "proc_m0324_v1"


def centroid(members=M0324):
    """Mean (lat, lon) of the members via unit vectors (safe near the date line)."""
    v = np.zeros(3)
    for _, _, la, lo, _, _ in members:
        a, b = math.radians(la), math.radians(lo)
        v += (math.cos(a) * math.cos(b), math.cos(a) * math.sin(b), math.sin(a))
    v /= np.linalg.norm(v)
    return math.degrees(math.asin(v[2])), math.degrees(math.atan2(v[1], v[0]))


def km_between(la1, lo1, la2, lo2):
    a1, b1, a2, b2 = map(math.radians, (la1, lo1, la2, lo2))
    c = math.sin(a1) * math.sin(a2) + math.cos(a1) * math.cos(a2) * math.cos(b1 - b2)
    return math.degrees(math.acos(max(-1.0, min(1.0, c)))) * 111.195


def member_rows():
    return [dict(id=i, time=t, latitude=str(la), longitude=str(lo), depth=str(d), mag=str(m), magType="mb",
                 type="earthquake") for i, t, la, lo, d, m in M0324]


def merge_catalog(rows, forced):
    """Union of catalogue rows and forced rows by id (catalogue row wins), as event dicts sorted by time."""
    by = {}
    for r in list(forced) + list(rows):
        if r.get("type", "earthquake") != "earthquake":
            continue
        by[r["id"]] = r
    evs = []
    for r in by.values():
        try:
            evs.append(dict(id=r["id"], time=r["time"], lat=float(r["latitude"]), lon=float(r["longitude"]),
                            depth=float(r["depth"]) if r.get("depth") not in (None, "") else 10.0,
                            mag=float(r["mag"]) if r.get("mag") not in (None, "") else float("nan"),
                            magType=r.get("magType", ""), t=model.decimal_year(r["time"])))
        except (KeyError, ValueError, TypeError):
            continue
    evs.sort(key=lambda e: e["t"])
    return evs


def fetch_catalog(data_dir, c_lat, c_lon, radius_km, start, end_iso, min_mag, log=print):
    import requests

    path = os.path.join(data_dir, "catalog_m0324.csv")
    params = dict(format="csv", starttime=start, endtime=end_iso, minmagnitude=min_mag, orderby="time-asc",
                  latitude="%.4f" % c_lat, longitude="%.4f" % c_lon, maxradiuskm=radius_km)
    text = None
    for attempt in range(4):
        try:
            r = requests.get(USGS, params=params, timeout=180)
            if r.status_code == 200:
                text = r.text
                break
            log("USGS status", r.status_code)
        except Exception as ex:
            log("USGS error", repr(ex)[:200])
        time.sleep(10 * (attempt + 1))
    rows = []
    if text is not None:
        rows = list(csv.DictReader(io.StringIO(text)))
        with open(path, "w", newline="") as f:
            f.write(text)
    elif os.path.exists(path):
        with open(path, newline="") as f:
            rows = list(csv.DictReader(f))
    return merge_catalog(rows, member_rows()), text is not None


def z0_membership(comps, member_ids, processed_ids):
    """Z0: does the pipeline re-find the published multiplet? Holds if >= 4 members were
    processed and >= 4 (and >= 2/3 of the processed ones) fall in one component."""
    mem = set(member_ids)
    proc = mem & set(processed_ids)
    best = max((len(mem & set(c)) for c in comps), default=0)
    testable = len(proc) >= 4
    holds = bool(testable and best >= 4 and best >= math.ceil(2 * len(proc) / 3.0))
    return dict(n_members=len(mem), n_processed=len(proc), n_in_largest_shared_component=best,
                testable=testable, holds=holds if testable else None)


def patch_events(evs, c_lat, c_lon, radius_km):
    """Events within radius_km of the centroid, plus all published members whatever their distance."""
    keep = []
    ids = {m[0] for m in M0324}
    for e in evs:
        d = km_between(e["lat"], e["lon"], c_lat, c_lon)
        e["km_from_centroid"] = d
        if d <= radius_km or e["id"] in ids:
            keep.append(e)
    return keep


def grade(rows, comps, processed_ids, nboot=1000):
    ok = [r for r in rows if r.get("ok")]
    outs = [r for r in ok if r["t2"] >= model.NEW_CUT]
    ins = [r for r in ok if r["t2"] < model.NEW_CUT]
    return dict(
        n_rows=len(rows), n_ok=len(ok), n_out_of_sample=len(outs),
        Z0_reproduce_M03_24=z0_membership(comps, [m[0] for m in M0324], processed_ids),
        Z1_new_event_best_match=model.p3(ok),
        Z2_out_of_sample=model.grade_set(outs, min_n=6, min_cross=3, nboot=nboot),
        Z2b_in_sample=model.grade_set(ins, nboot=nboot),
        Z3_calibration=model.p4(ok),
        Z5_fit_all=model.fit_turn(ok) if len(ok) >= 4 else None,
    )


def main():
    import run
    from cluster import components
    from phases import Phases, load_or_build

    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", default="/data")
    ap.add_argument("--out-dir", default="/out")
    ap.add_argument("--end", default=time.strftime("%Y-%m-%dT00:00:00", time.gmtime()))
    ap.add_argument("--start", default="1991-01-01")
    ap.add_argument("--radius-km", type=float, default=35.0)
    ap.add_argument("--min-mag", type=float, default=4.3)
    ap.add_argument("--detect-thr", type=float, default=0.90)
    ap.add_argument("--pair-thr", type=float, default=0.93)
    ap.add_argument("--min-snr-group", type=float, default=3.0)
    ap.add_argument("--min-snr-ref", type=float, default=5.0)
    ap.add_argument("--min-snr-df", type=float, default=3.0)
    ap.add_argument("--max-lag-s", type=float, default=8.0)
    ap.add_argument("--max-station-level", type=int, default=40)
    ap.add_argument("--time-budget-min", type=float, default=150.0)
    args = ap.parse_args()
    os.makedirs(args.data_dir, exist_ok=True)
    os.makedirs(args.out_dir, exist_ok=True)
    plots = os.path.join(args.out_dir, "plots")
    os.makedirs(plots, exist_ok=True)
    proc_dir = os.path.join(args.data_dir, PROC_TAG)
    os.makedirs(proc_dir, exist_ok=True)
    log = run.log
    s = dict(version=VERSION, args=vars(args), started=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()))

    c_lat, c_lon = centroid()
    s["centroid"] = [round(c_lat, 4), round(c_lon, 4)]
    log("M03-24 centroid", s["centroid"])
    ph = Phases(load_or_build(os.path.join(args.data_dir, "ak135_pkp_grid_v1.npz")))
    s["phase_check"] = [dict(dist=d, bc_minus_df=ph.t("bc", d, 40.0) - ph.t("df", d, 40.0),
                             ab_minus_df=ph.t("ab", d, 40.0) - ph.t("df", d, 40.0)) for d in (148.0, 150.0, 152.0)]
    evs_all, cat_ok = fetch_catalog(args.data_dir, c_lat, c_lon, args.radius_km + 10.0, args.start, args.end,
                                    args.min_mag, log=log)
    s["catalog_fetched"] = cat_ok
    evs_list = patch_events(evs_all, c_lat, c_lon, args.radius_km)
    log("patch events:", len(evs_list))
    run.write_csv(os.path.join(args.out_dir, "m0324_catalog.csv"), evs_list)
    evs = {e["id"]: e for e in evs_list}

    # newest and the 1998-2008.5 partners first, so a short budget still covers the test
    def prio(e):
        g = 0 if e["t"] >= 2024.5 else (1 if 1998.0 <= e["t"] < 2008.5 else 2)
        return (g, -(e["mag"] if math.isfinite(e["mag"]) else 0))

    client, coords, data, status = None, None, {}, []
    for e in sorted(evs_list, key=prio):
        pp = os.path.join(proc_dir, e["id"] + ".npz")
        if os.path.exists(pp):
            try:
                data[e["id"]] = run.load_proc(pp)
                status.append(dict(id=e["id"], t=e["t"], mag=e["mag"], status="ok-cached"))
                continue
            except Exception:
                os.remove(pp)
        if run.elapsed_min() > 0.7 * args.time_budget_min:
            status.append(dict(id=e["id"], t=e["t"], mag=e["mag"], status="skipped-time-budget"))
            continue
        try:
            if client is None:
                client = run.get_client()
                coords = run.get_coords(client, args.data_dir)
            st, why = run.fetch_event(client, e, ph, args.data_dir)
            if st is None:
                status.append(dict(id=e["id"], t=e["t"], mag=e["mag"], status=why))
                continue
            d, why = run.process_event(st, e, ph, coords)
            if d is None:
                status.append(dict(id=e["id"], t=e["t"], mag=e["mag"], status=why))
                continue
            run.save_proc(pp, d)
            data[e["id"]] = d
            status.append(dict(id=e["id"], t=e["t"], mag=e["mag"], status="ok", n_sta=len(d["names"]),
                               snr_group=round(d["snr_group"], 2)))
        except Exception as ex:
            status.append(dict(id=e["id"], t=e["t"], mag=e["mag"], status="exception:" + repr(ex)[:150]))
        time.sleep(0.2)
    run.write_csv(os.path.join(args.out_dir, "m0324_event_status.csv"), status,
                  keys=["id", "t", "mag", "status", "n_sta", "snr_group"])
    by_year = {}
    for r in status:
        b = by_year.setdefault(int(r["t"]), dict(tried=0, ok=0))
        b["tried"] += 1
        b["ok"] += int(str(r["status"]).startswith("ok"))
    s["coverage_by_year"] = by_year
    s["n_processed"] = len(data)

    ids = sorted(data, key=lambda i: evs[i]["t"])
    edges, ncomp, nuse = run.detect_pairs(ids, evs, data, 2 * args.radius_km + 20.0, args.detect_thr,
                                          args.max_lag_s, args.min_snr_group)
    s.update(n_compared=ncomp, n_snr_ok=nuse, n_edges=len(edges))
    run.write_csv(os.path.join(args.out_dir, "m0324_edges.csv"),
                  [dict(ev1=a, ev2=b, t1=evs[a]["t"], t2=evs[b]["t"], cc=c, lag_s=l / run.FS)
                   for (a, b), (c, l) in edges.items()])
    comps = sorted([c for c in components(ids, edges) if len(c) >= 2], key=lambda c: -len(c))
    a2 = SimpleNamespace(pair_thr=args.pair_thr, min_snr_ref=args.min_snr_ref, min_snr_df=args.min_snr_df,
                         max_station_level=args.max_station_level)
    rows, infos = [], []
    for k, members in enumerate(comps):
        mid = "Z%03d" % k
        try:
            r, info, fr = run.analyze_multiplet(mid, members, evs, data, edges, a2)
        except Exception as ex:
            traceback.print_exc()
            r, info, fr = [], dict(mult=mid, n=len(members), status="exception:" + repr(ex)[:150]), None
        info["members"] = ";".join("%s(%s)" % (m, evs[m]["time"][:10]) for m in sorted(members, key=lambda m: evs[m]["t"]))
        info["n_published_members"] = len({m[0] for m in M0324} & set(members))
        rows += r
        infos.append(info)
        if fr is not None and k < 12:
            try:
                run.plot_multiplet(os.path.join(plots, "m0324_%s.png" % mid), mid, fr, evs, data, r)
            except Exception:
                traceback.print_exc()
    run.write_csv(os.path.join(args.out_dir, "m0324_multiplets.csv"), infos,
                  keys=["mult", "n", "n_published_members", "status", "root", "pick_minus_pred_s", "df_window_s",
                        "dt_bc", "dt_ab", "members"])
    run.write_csv(os.path.join(args.out_dir, "m0324_pairs.csv"), rows)
    try:
        res = grade(rows, comps, list(data))
    except Exception as ex:
        traceback.print_exc()
        res = dict(error=repr(ex))
    s["tests"] = res
    try:
        run.plot_summary(os.path.join(plots, "m0324_D_vs_time.png"), rows)
    except Exception:
        traceback.print_exc()
    s["finished"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    s["minutes"] = round(run.elapsed_min(), 1)
    with open(os.path.join(args.out_dir, "m0324_summary.json"), "w") as f:
        json.dump(run.jsonable(s), f, indent=1)
    log("done", json.dumps(run.jsonable({k: res.get(k) for k in ("n_ok", "n_out_of_sample")}) if isinstance(res, dict) else {}))


if __name__ == "__main__":
    main()

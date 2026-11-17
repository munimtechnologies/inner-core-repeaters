#!/usr/bin/env python3
"""Forecast-date scan: the backtracking forecast tested on every South Sandwich patch.

Part of inner-core-repeaters by MV2, Munim, Inc. (made with AI agents). Licence: Apache-2.0.

Wang, Vidale et al. (2024, Nature 631:340) argue the inner core turned at T = 2008.5 and is
retracing its path 2.5x more slowly. A repeating earthquake at t > T should then best match
(in its PKIKP waveform at the ILAR array) a partner from the equivalent date
e(t) = T - (t - T)/2.5, which is 2001.1-2001.9 for 2025-2026 events, rather than its most
recent pre-turn partner (what a monotone drift predicts).

m0324.py tests this at one patch (Zhang & Wen's 2026 multiplet M03-24). This script tests it
on every 2025+ South Sandwich Islands event that, in the USGS catalogue, has neighbours
(within r km) both inside the forecast window |t1 - e(t)| <= 1 yr AND later in 2003-2008.5,
so that the backtracking and monotone picks can differ. Waveforms go through the unchanged
ILAR pipeline in run.py. Pre-registered tests: findings/2026-11-17-fdscan-prereg.md.
  F1   best-match vote over all discriminating new events (D), binomial + multiplet bootstrap
  F1x  the same with Dx (D minus reference-window misfit), the robustness check
  F1e  F1 without the M03-24 patch (independent of icr-m0324)
  F3   calibration (model.p4)
  F5   exploratory turn-time fit (model.fit_turn), not a test
"""
import argparse
import csv
import io
import json
import math
import os
import time
import traceback

import numpy as np

import model

VERSION = "0.6.0"
USGS = "https://earthquake.usgs.gov/fdsnws/event/1/query"
BOX = dict(minlatitude=-61.5, maxlatitude=-54.0, minlongitude=-31.0, maxlongitude=-22.0)
PROC_TAG = "proc_fd_v1"
POOL_START = 1999.0
M0324_CENTROID = (-59.036, -25.694)


def km_between(la1, lo1, la2, lo2):
    a1, b1, a2, b2 = map(math.radians, (la1, lo1, la2, lo2))
    c = math.sin(a1) * math.sin(a2) + math.cos(a1) * math.cos(a2) * math.cos(b1 - b2)
    return math.degrees(math.acos(max(-1.0, min(1.0, c)))) * 111.195


def parse_rows(rows):
    """USGS csv rows -> event dicts (earthquakes only, unique ids, sorted by time)."""
    by = {}
    for r in rows:
        if r.get("type", "earthquake") != "earthquake":
            continue
        try:
            by[r["id"]] = dict(id=r["id"], time=r["time"], lat=float(r["latitude"]), lon=float(r["longitude"]),
                               depth=float(r["depth"]) if r.get("depth") not in (None, "") else 10.0,
                               mag=float(r["mag"]) if r.get("mag") not in (None, "") else float("nan"),
                               magType=r.get("magType", ""), t=model.decimal_year(r["time"]))
        except (KeyError, ValueError, TypeError):
            continue
    return sorted(by.values(), key=lambda e: e["t"])


def fetch_catalog(data_dir, end_iso, min_mag_old, min_mag_new, log=print):
    """Pre-turn pool (1999 to mid-2008) and new pool (2025 on), whole SSI box."""
    import requests

    path = os.path.join(data_dir, "catalog_fd.csv")
    rows, complete = [], True
    for s, e, m in [("1999-01-01", "2008-07-02", min_mag_old), ("2025-01-01", end_iso, min_mag_new)]:
        params = dict(format="csv", starttime=s, endtime=e, minmagnitude=m, orderby="time-asc", **BOX)
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
        if text is None:
            complete = False
            continue
        rows += list(csv.DictReader(io.StringIO(text)))
    if not complete and os.path.exists(path):
        with open(path, newline="") as f:
            rows += list(csv.DictReader(f))
    if rows:
        keys = list(rows[0].keys())
        with open(path, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=keys, extrasaction="ignore")
            w.writeheader()
            for r in rows:
                w.writerow(r)
    return parse_rows(rows), complete


def select_cases(evs, r_km=25.0, b_half=1.0, pool_start=POOL_START, T=model.T_TURN, new_cut=model.NEW_CUT):
    """Catalogue-level candidates: a new event (t >= new_cut) with neighbours within r_km both in
    the forecast window |t1 - e(t)| <= b_half and later than it (up to T)."""
    old = [e for e in evs if pool_start <= e["t"] < T]
    cases = []
    for n in evs:
        if n["t"] < new_cut:
            continue
        e2 = model.equiv_date(n["t"])
        nb = [(o, km_between(n["lat"], n["lon"], o["lat"], o["lon"])) for o in old]
        nb = [(o, d) for o, d in nb if d <= r_km]
        b = [o for o, _ in nb if abs(o["t"] - e2) <= b_half]
        m = [o for o, _ in nb if o["t"] > e2 + b_half]
        if b and m:
            cases.append(dict(new=n, e2=e2, b=b, m=m))
    cases.sort(key=lambda c: -(c["new"]["mag"] if math.isfinite(c["new"]["mag"]) else 0.0))
    return cases


def event_order(cases, max_events=0):
    """Events to process, case by case (largest new event first): new, its forecast-window
    neighbours, then its later neighbours, each by magnitude. Duplicates removed."""
    def bym(lst):
        return sorted(lst, key=lambda e: -(e["mag"] if math.isfinite(e["mag"]) else 0.0))

    out, seen = [], set()
    for c in cases:
        for e in [c["new"]] + bym(c["b"]) + bym(c["m"]):
            if e["id"] not in seen:
                seen.add(e["id"])
                out.append(e)
    return out[:max_events] if max_events else out


def binom_sf(k, n):
    """P(X >= k) for X ~ Binomial(n, 1/2)."""
    if n <= 0:
        return 1.0
    return sum(math.comb(n, j) for j in range(max(0, k), n + 1)) / 2.0 ** n


def fd_cases(rows, key="D", new_cut=model.NEW_CUT, T=model.T_TURN):
    """Per new event with >= 2 ok pre-turn partners: best match vs the B pick (partner closest
    to e(t)) and the M pick (latest pre-turn partner). Mirrors model.p3, plus multiplet id,
    chance levels and an SNR control."""
    by = {}
    for r in rows:
        if r.get("ok") and r["t2"] >= new_cut and r["t1"] < T and math.isfinite(r.get(key, float("nan"))):
            by.setdefault(r["ev2"], []).append(r)
    cases = []
    for ev, rs in sorted(by.items()):
        if len(rs) < 2:
            continue
        e2 = model.equiv_date(rs[0]["t2"])
        best = min(rs, key=lambda r: r[key])
        bpk = min(rs, key=lambda r: abs(model.equiv_date(r["t1"]) - e2))
        mpk = max(rs, key=lambda r: r["t1"])
        top = max(rs, key=lambda r: r.get("snr_df1", float("-inf")) if math.isfinite(r.get("snr_df1", float("nan"))) else float("-inf"))
        cases.append(dict(ev2=ev, mult=rs[0].get("mult", ""), t2=rs[0]["t2"], e2=e2, n_partners=len(rs),
                          best_ev1=best["ev1"], best_t1=best["t1"], B_pick=bpk["ev1"], B_pick_t1=bpk["t1"],
                          M_pick=mpk["ev1"], M_pick_t1=mpk["t1"], discriminating=bpk["ev1"] != mpk["ev1"],
                          matches_B=best["ev1"] == bpk["ev1"], matches_M=best["ev1"] == mpk["ev1"],
                          best_is_top_snr=best["ev1"] == top["ev1"],
                          D_best=best[key], D_B=bpk[key], D_M=mpk[key]))
    return cases


def fd_grade(cases, nboot=2000, seed=0, min_cases=6, min_mult=3):
    """Pre-registered decision rule (findings/2026-11-17-fdscan-prereg.md)."""
    dc = [c for c in cases if c["discriminating"]]
    nB = sum(c["matches_B"] for c in dc)
    nM = sum(c["matches_M"] for c in dc)
    n = nB + nM
    mults = sorted({c["mult"] for c in dc})
    chance = sum(1.0 / c["n_partners"] for c in dc)
    pB, pM = binom_sf(nB, n), binom_sf(nM, n)
    frac_B = frac_M = float("nan")
    if len(mults) >= 2:
        groups = {m: [c for c in dc if c["mult"] == m] for m in mults}
        rng = np.random.default_rng(seed)
        d = []
        for _ in range(nboot):
            pick = rng.integers(0, len(mults), len(mults))
            cs = [c for k in pick for c in groups[mults[k]]]
            d.append(sum(c["matches_B"] for c in cs) - sum(c["matches_M"] for c in cs))
        d = np.asarray(d)
        frac_B, frac_M = float((d > 0).mean()), float((d < 0).mean())
    testable = len(dc) >= min_cases and len(mults) >= min_mult
    verdict = "untestable"
    if testable:
        if nB > nM and pB < 0.05 and frac_B >= 0.95:
            verdict = "backtracking-supported"
        elif nM > nB and pM < 0.05 and frac_M >= 0.95:
            verdict = "backtracking-contradicted"
        else:
            verdict = "inconclusive"
    return dict(n_cases=len(cases), n_discriminating=len(dc), n_multiplets=len(mults), n_match_B=nB,
                n_match_M=nM, n_neither=len(dc) - n, chance_expected_B=chance, chance_expected_M=chance,
                p_binom_B=pB, p_binom_M=pM, boot_frac_B_ahead=frac_B, boot_frac_M_ahead=frac_M,
                frac_best_is_top_snr=(sum(c["best_is_top_snr"] for c in dc) / len(dc)) if dc else float("nan"),
                testable=testable, verdict=verdict)


def grade_all(rows, evs, nboot=2000):
    ok = [r for r in rows if r.get("ok")]
    cD = fd_cases(rows, "D")
    cX = fd_cases(rows, "Dx")
    near = {i for i, e in evs.items() if km_between(e["lat"], e["lon"], *M0324_CENTROID) <= 35.0}
    cE = [c for c in cD if c["ev2"] not in near]
    F1, F1x, F1e = fd_grade(cD, nboot), fd_grade(cX, nboot), fd_grade(cE, nboot)
    robust = (F1["verdict"] in ("backtracking-supported", "backtracking-contradicted")
              and F1x["n_discriminating"] > 0
              and np.sign(F1["n_match_B"] - F1["n_match_M"]) == np.sign(F1x["n_match_B"] - F1x["n_match_M"]))
    return dict(n_rows=len(rows), n_ok=len(ok), F1=F1, F1x=F1x, F1e_without_M0324=F1e, F1_robust=bool(robust),
                F3_calibration=model.p4(ok) if ok else None,
                F5_fit_all=model.fit_turn(ok) if len(ok) >= 4 else None,
                p3_crosscheck={k: v for k, v in model.p3(ok).items() if k != "cases"} if ok else None,
                cases=cD)


def main():
    import run
    from cluster import components
    from phases import Phases, load_or_build

    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", default="/data")
    ap.add_argument("--out-dir", default="/out")
    ap.add_argument("--end", default=time.strftime("%Y-%m-%dT00:00:00", time.gmtime()))
    ap.add_argument("--min-mag-old", type=float, default=4.3)
    ap.add_argument("--min-mag-new", type=float, default=4.5)
    ap.add_argument("--radius-km", type=float, default=25.0)
    ap.add_argument("--b-half", type=float, default=1.0)
    ap.add_argument("--detect-thr", type=float, default=0.90)
    ap.add_argument("--pair-thr", type=float, default=0.93)
    ap.add_argument("--min-snr-group", type=float, default=3.0)
    ap.add_argument("--min-snr-ref", type=float, default=5.0)
    ap.add_argument("--min-snr-df", type=float, default=3.0)
    ap.add_argument("--max-lag-s", type=float, default=8.0)
    ap.add_argument("--max-station-level", type=int, default=40)
    ap.add_argument("--max-events", type=int, default=700)
    ap.add_argument("--time-budget-min", type=float, default=200.0)
    ap.add_argument("--fetch-frac", type=float, default=0.7)
    ap.add_argument("--plots", type=int, default=16)
    args = ap.parse_args()
    for d in (args.data_dir, args.out_dir, os.path.join(args.out_dir, "plots"), os.path.join(args.data_dir, PROC_TAG)):
        os.makedirs(d, exist_ok=True)
    plots = os.path.join(args.out_dir, "plots")
    proc_dir = os.path.join(args.data_dir, PROC_TAG)
    log = run.log
    s = dict(version=VERSION, args=vars(args), started=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()))

    ph = Phases(load_or_build(os.path.join(args.data_dir, "ak135_pkp_grid_v1.npz")))
    s["phase_check"] = [dict(dist=d, bc_minus_df=ph.t("bc", d, 30.0) - ph.t("df", d, 30.0),
                             ab_minus_df=ph.t("ab", d, 30.0) - ph.t("df", d, 30.0)) for d in (146.0, 149.0, 152.0)]
    evs_list, cat_ok = fetch_catalog(args.data_dir, args.end, args.min_mag_old, args.min_mag_new, log=log)
    s["catalog_complete"] = cat_ok
    s["n_catalog"] = len(evs_list)
    evs = {e["id"]: e for e in evs_list}
    cases = select_cases(evs_list, args.radius_km, args.b_half)
    s["n_catalog_new"] = sum(1 for e in evs_list if e["t"] >= model.NEW_CUT)
    s["n_catalog_cases"] = len(cases)
    log("catalogue: %d events, %d new, %d candidate cases" % (len(evs_list), s["n_catalog_new"], len(cases)))
    run.write_csv(os.path.join(args.out_dir, "fd_catalog_cases.csv"),
                  [dict(new_id=c["new"]["id"], new_time=c["new"]["time"], new_mag=c["new"]["mag"],
                        e_t=round(c["e2"], 3), n_forecast_window=len(c["b"]), n_later=len(c["m"]),
                        forecast_window_ids=";".join("%s(%s)" % (o["id"], o["time"][:10]) for o in c["b"]),
                        later_ids=";".join("%s(%s)" % (o["id"], o["time"][:10]) for o in c["m"]))
                   for c in cases])
    order = event_order(cases, args.max_events)
    s["n_to_process"] = len(order)

    client, coords, data, status = None, None, {}, []
    stop_at = args.time_budget_min * args.fetch_frac
    for k, e in enumerate(order):
        pp = os.path.join(proc_dir, e["id"] + ".npz")
        if os.path.exists(pp):
            try:
                data[e["id"]] = run.load_proc(pp)
                status.append(dict(id=e["id"], t=e["t"], mag=e["mag"], status="ok-cached"))
                continue
            except Exception:
                os.remove(pp)
        if os.path.exists(os.path.join(args.data_dir, "raw", e["id"] + ".mseed.nodata")):
            status.append(dict(id=e["id"], t=e["t"], mag=e["mag"], status="nodata"))
            continue
        if run.elapsed_min() > stop_at:
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
        if (k + 1) % 50 == 0:
            log("events tried %d/%d, processed %d" % (k + 1, len(order), len(data)))
        time.sleep(0.2)
    run.write_csv(os.path.join(args.out_dir, "fd_event_status.csv"), status,
                  keys=["id", "t", "mag", "status", "n_sta", "snr_group"])
    by_year = {}
    for r in status:
        b = by_year.setdefault(int(r["t"]), dict(tried=0, ok=0))
        b["tried"] += 1
        b["ok"] += int(str(r["status"]).startswith("ok"))
    s["coverage_by_year"] = by_year
    s["n_processed"] = len(data)
    log("processed:", len(data))

    ids = sorted(data, key=lambda i: evs[i]["t"])
    edges, ncomp, nuse = run.detect_pairs(ids, evs, data, 2 * args.radius_km + 10.0, args.detect_thr,
                                          args.max_lag_s, args.min_snr_group)
    s.update(n_compared=ncomp, n_snr_ok=nuse, n_edges=len(edges))
    run.write_csv(os.path.join(args.out_dir, "fd_edges.csv"),
                  [dict(ev1=a, ev2=b, t1=evs[a]["t"], t2=evs[b]["t"], cc=c, lag_s=l / run.FS)
                   for (a, b), (c, l) in edges.items()])
    comps = [c for c in components(ids, edges) if len(c) >= 2
             and any(evs[m]["t"] >= model.NEW_CUT for m in c) and any(evs[m]["t"] < model.T_TURN for m in c)]
    comps.sort(key=lambda c: len(c))  # small first, so a giant chained component cannot eat the budget
    s["n_multiplets_new_and_preturn"] = len(comps)
    a2 = argparse.Namespace(pair_thr=args.pair_thr, min_snr_ref=args.min_snr_ref, min_snr_df=args.min_snr_df,
                            max_station_level=args.max_station_level)
    rows, infos, frames = [], [], []
    for k, members in enumerate(comps):
        mid = "F%03d" % k
        if run.elapsed_min() > args.time_budget_min * 0.95:
            infos.append(dict(mult=mid, n=len(members), status="skipped-time-budget"))
            continue
        try:
            r, info, fr = run.analyze_multiplet(mid, members, evs, data, edges, a2)
        except Exception as ex:
            traceback.print_exc()
            r, info, fr = [], dict(mult=mid, n=len(members), status="exception:" + repr(ex)[:150]), None
        info["members"] = ";".join("%s(%s)" % (m, evs[m]["time"][:10]) for m in sorted(members, key=lambda m: evs[m]["t"]))
        info["centroid"] = "%.3f,%.3f" % (float(np.mean([evs[m]["lat"] for m in members])),
                                          float(np.mean([evs[m]["lon"] for m in members])))
        rows += r
        infos.append(info)
        if fr is not None:
            frames.append((mid, fr, members))
    run.write_csv(os.path.join(args.out_dir, "fd_multiplets.csv"), infos,
                  keys=["mult", "n", "status", "root", "centroid", "pick_minus_pred_s", "df_window_s", "dt_bc",
                        "dt_ab", "members"])
    run.write_csv(os.path.join(args.out_dir, "fd_pairs.csv"), rows)
    try:
        res = grade_all(rows, evs)
    except Exception as ex:
        traceback.print_exc()
        res = dict(error=repr(ex))
    run.write_csv(os.path.join(args.out_dir, "fd_cases.csv"), res.get("cases", []) if isinstance(res, dict) else [])
    s["tests"] = {k: v for k, v in res.items() if k != "cases"} if isinstance(res, dict) else res
    try:
        run.plot_summary(os.path.join(plots, "fd_D_vs_time.png"), rows)
        for mid, fr, members in sorted(frames, key=lambda x: -len(x[2]))[:args.plots]:
            run.plot_multiplet(os.path.join(plots, "fd_%s.png" % mid), mid, fr, evs, data,
                               [r for r in rows if r["mult"] == mid])
    except Exception:
        traceback.print_exc()
    s["finished"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    s["minutes"] = round(run.elapsed_min(), 1)
    with open(os.path.join(args.out_dir, "fd_summary.json"), "w") as f:
        json.dump(run.jsonable(s), f, indent=1)
    t = s["tests"] if isinstance(s["tests"], dict) else {}
    log("done", json.dumps(run.jsonable({k: (t.get(k) or {}).get("verdict") for k in ("F1", "F1x", "F1e_without_M0324")})))


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
# Copyright 2026 MV2, Munim, Inc.
# Licensed under the Apache License, Version 2.0 (Apache-2.0). See LICENSE.
# Made with AI agents working in MV2's name.
"""dprk37probe.py - blind feasibility probe for line X: explosion doublets through the inner core.

Question for a later pre-registration (37): the six underground nuclear tests at Punggye-ri,
North Korea (2006-2017) are near-perfectly repeating sources (same site, explosive source,
relative positions known to about 100 m), and South America and Antarctica lie at PKP
distances (150-180 deg) from them. Do their inner-core waves (PKIKP = DF) change between
tests, relative to their outer-core waves (PKPab = AB), during 2013-2017, the years when the
rival readings of the inner core (one-turn backtracking, near standstill, episodic surface
change) disagree?

This probe is REPORT-ONLY and BLIND. It measures only what a fair pre-registration needs:
which stations recorded which tests, the channel epoch used, and the signal-to-noise ratio of
the DF and AB windows on each single record. It never compares two tests: no cross-
correlation, travel-time difference, amplitude ratio or overlay of two events is computed,
written or plotted. Records are cached under --data-dir for the later test.

0.1.1 (2026-11-25): gc_deg uses the atan2(|a x b|, a.b) form, exact for coincident points
(0.1.0's acos form returned 8.5e-7 deg and failed its own test). Per-record SNRs are written
only as coarse bins (<3, 3-5, 5-10, 10-20, >=20), so the records file cannot be used to form
a per-record DF/AB amplitude ratio before the pre-registration; counts use the exact values.

Top-level imports are stdlib only (the tests run without numpy/obspy).
"""

import argparse
import csv
import json
import math
import os
import sys
import time
import traceback

VERSION = "0.1.1"
NAN = float("nan")
FDSN = "https://service.earthscope.org"
UA = "inner-core-repeaters (MV2, Munim, Inc.)"

SRC_LAT, SRC_LON = 41.30, 129.05
# USGS ComCat (query of 2026-11-25): (id, origin UTC, lat, lon, mb)
EVENTS = [
    ("usp000eurb", "2006-10-09T01:35:28.020", 41.294, 129.094, 4.3),
    ("usp000gxgc", "2009-05-25T00:54:43.120", 41.303, 129.037, 4.7),
    ("usc000f5t0", "2013-02-12T02:57:51.490", 41.299, 129.004, 5.1),
    ("us10004bnm", "2016-01-06T01:30:01.480", 41.2996, 129.0467, 5.1),
    ("us10006n8a", "2016-09-09T00:30:01.440", 41.2869, 129.0783, 5.3),
    ("us2000aert", "2017-09-03T03:30:01.760", 41.3324, 129.0297, 6.3),
]
DIST_RANGE = (150.0, 179.5)
BANDS = ((0.7, 2.0), (1.0, 2.0))
NOISE = (-60.0, -5.0)     # s relative to predicted DF
SIG = (-1.0, 4.0)         # s relative to the predicted phase
SNR_MIN = 5.0
SNR_BIN_EDGES = (3.0, 5.0, 10.0, 20.0)
PRE_S = 90.0              # cached record starts DF - 90 s
POST_S = 120.0            # and ends AB + 120 s (or DF + 180 s without AB)
CHAN_PREF = ("BHZ", "HHZ")
LOC_PREF = ("00", "", "10")
CHUNK = 25


# ----------------------------------------------------------------------------- pure helpers
def isf(x):
    try:
        return math.isfinite(float(x))
    except (TypeError, ValueError):
        return False


def _unit(lat, lon):
    p, l = math.radians(lat), math.radians(lon)
    return (math.cos(p) * math.cos(l), math.cos(p) * math.sin(l), math.sin(p))


def gc_deg(lat1, lon1, lat2, lon2):
    """Great-circle distance in degrees, atan2(|a x b|, a.b): accurate at 0 and 180 deg."""
    a, b = _unit(lat1, lon1), _unit(lat2, lon2)
    cx = a[1] * b[2] - a[2] * b[1]
    cy = a[2] * b[0] - a[0] * b[2]
    cz = a[0] * b[1] - a[1] * b[0]
    dot = a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
    return math.degrees(math.atan2(math.sqrt(cx * cx + cy * cy + cz * cz), dot))


def snr_bin(v):
    """Coarse SNR label for the written records file ('' when not measured)."""
    if not isf(v):
        return ""
    v = float(v)
    labels = ("<3", "3-5", "5-10", "10-20", ">=20")
    for k, edge in enumerate(SNR_BIN_EDGES):
        if v < edge:
            return labels[k]
    return labels[-1]


def blind_rows(rows):
    """Copies of rows with every snr_* value replaced by its coarse bin."""
    out = []
    for r in rows:
        c = dict(r)
        for k in list(c):
            if k.startswith("snr_"):
                c[k] = snr_bin(c[k])
        out.append(c)
    return out


def branches(arrivals, dist):
    """arrivals: [(name, time)] from TauP for PKIKP and PKP. Returns {df, bc, ab} (s or None).
    Two PKP arrivals: the earlier is BC, the later AB. One PKP arrival: AB at >= 150 deg."""
    df = [t for n, t in arrivals if n in ("PKIKP", "PKPdf")]
    pkp = sorted(t for n, t in arrivals if n in ("PKP", "PKPab", "PKPbc"))
    out = dict(df=min(df) if df else None, bc=None, ab=None)
    if len(pkp) >= 2:
        out["bc"], out["ab"] = pkp[0], pkp[-1]
    elif len(pkp) == 1:
        if dist >= 150.0:
            out["ab"] = pkp[0]
        else:
            out["bc"] = pkp[0]
    return out


def windows_clear(b):
    """True when the DF and AB signal windows do not overlap each other (or BC)."""
    if b.get("df") is None or b.get("ab") is None:
        return False
    df1 = b["df"] + SIG[1]
    if b["ab"] + SIG[0] <= df1:
        return False
    if b.get("bc") is not None and b["bc"] + SIG[0] <= df1:
        return False
    return True


def choose_channel(chans):
    """chans: [(loc, cha)] -> the preferred (loc, cha): BHZ before HHZ, loc 00, '', 10, other."""
    def rank(lc):
        loc, cha = lc
        ci = CHAN_PREF.index(cha) if cha in CHAN_PREF else 9
        li = LOC_PREF.index(loc) if loc in LOC_PREF else 5
        return (ci, li, loc)
    c = sorted(set(chans), key=rank)
    if c and c[0][1] in CHAN_PREF:
        return c[0]
    return None


def record_ok(row, band_label, phase):
    v = row.get("snr_%s_%s" % (phase, band_label))
    return isf(v) and float(v) >= SNR_MIN


def pair_counts(rows, band_label):
    """Counts only (no waveform comparison): for each pair of events, how many stations have
    DF SNR >= 5 on both records, also AB SNR >= 5 on both, and the same channel epoch."""
    by_ev = {}
    for r in rows:
        if r.get("status") != "ok":
            continue
        by_ev.setdefault(r["event"], {})[r["station"]] = r
    ids = [e[0] for e in EVENTS]
    out = []
    for i in range(len(ids)):
        for j in range(i + 1, len(ids)):
            a, b = by_ev.get(ids[i], {}), by_ev.get(ids[j], {})
            n_df = n_dfab = n_same = n_same_clear = 0
            for st in sorted(set(a) & set(b)):
                ra, rb = a[st], b[st]
                if not (record_ok(ra, band_label, "df") and record_ok(rb, band_label, "df")):
                    continue
                n_df += 1
                both_ab = record_ok(ra, band_label, "ab") and record_ok(rb, band_label, "ab")
                if both_ab:
                    n_dfab += 1
                same = (ra.get("channel") == rb.get("channel")
                        and ra.get("epoch_start") == rb.get("epoch_start"))
                if same and both_ab:
                    n_same += 1
                    if str(ra.get("windows_clear")) == "True":
                        n_same_clear += 1
            out.append(dict(band=band_label, a=ids[i], b=ids[j],
                            a_date=EVENTS[i][1][:10], b_date=EVENTS[j][1][:10],
                            n_df_both=n_df, n_df_ab_both=n_dfab,
                            n_df_ab_same_epoch=n_same, n_df_ab_same_epoch_clear=n_same_clear))
    return out


def event_summary(rows, band_label):
    out = {}
    for e in EVENTS:
        rr = [r for r in rows if r["event"] == e[0]]
        ok = [r for r in rr if r.get("status") == "ok"]
        df = [r for r in ok if record_ok(r, band_label, "df")]
        dfab = [r for r in df if record_ok(r, band_label, "ab")]
        bins = {}
        for r in df:
            k = "%d-%d" % (int(float(r["dist"]) // 10 * 10), int(float(r["dist"]) // 10 * 10) + 10)
            bins[k] = bins.get(k, 0) + 1
        out[e[1][:10]] = dict(id=e[0], mb=e[4], n_listed=len(rr), n_ok=len(ok), n_df=len(df),
                              n_df_ab=len(dfab), df_by_distance=bins,
                              networks=sorted({r["station"].split(".")[0] for r in df}))
    return out


def write_csv(path, rows):
    keys = []
    for r in rows:
        for k in r:
            if k not in keys:
                keys.append(k)
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=keys or ["empty"])
        w.writeheader()
        for r in rows:
            w.writerow(r)


# ----------------------------------------------------------------------------- numerics
def snr_values(tr, origin, b, band):
    """SNR (rms ratio) of the DF and AB windows against the pre-DF noise, one record only."""
    import numpy as np
    x = tr.copy()
    x.data = x.data.astype("float64")
    x.detrend("demean")
    x.detrend("linear")
    x.taper(0.05)
    x.filter("bandpass", freqmin=band[0], freqmax=band[1], corners=4, zerophase=False)
    fs = float(x.stats.sampling_rate)
    t0 = x.stats.starttime

    def seg(rel0, rel1):
        i0 = int(round(((origin + rel0) - t0) * fs))
        i1 = int(round(((origin + rel1) - t0) * fs))
        if i0 < 0 or i1 > len(x.data) or i1 - i0 < 4:
            return None
        return x.data[i0:i1]

    def rms(v):
        return math.sqrt(float(np.mean(v ** 2))) if v is not None and len(v) else NAN

    nz = rms(seg(b["df"] + NOISE[0], b["df"] + NOISE[1]))
    out = dict(snr_df=NAN, snr_ab=NAN)
    if not (nz > 0):
        return out
    out["snr_df"] = rms(seg(b["df"] + SIG[0], b["df"] + SIG[1])) / nz
    if b.get("ab") is not None:
        out["snr_ab"] = rms(seg(b["ab"] + SIG[0], b["ab"] + SIG[1])) / nz
    return out


# ----------------------------------------------------------------------------- per event
def probe_event(client, model, ev, cache, ttcache, s):
    from obspy import UTCDateTime, read, Stream
    import numpy as np
    eid, t_iso, elat, elon, mb = ev
    origin = UTCDateTime(t_iso)
    rows = []
    try:
        inv = client.get_stations(latitude=SRC_LAT, longitude=SRC_LON, minradius=DIST_RANGE[0],
                                  maxradius=DIST_RANGE[1], channel="BHZ,HHZ",
                                  starttime=origin - 60, endtime=origin + 3600, level="channel")
    except Exception as ex:  # noqa: BLE001
        s.setdefault("station_errors", {})[eid] = repr(ex)[:200]
        return rows
    plan = []
    seen = set()
    for net in inv:
        for sta in net:
            if (net.code, sta.code) in seen:
                continue
            chans = {}
            for ch in sta:
                if ch.code not in CHAN_PREF:
                    continue
                if ch.start_date and ch.start_date > origin:
                    continue
                if ch.end_date and ch.end_date < origin + 3600:
                    continue
                chans[(ch.location_code or "", ch.code)] = ch
            pick = choose_channel(list(chans))
            if pick is None:
                continue
            ch = chans[pick]
            lat = ch.latitude if ch.latitude is not None else sta.latitude
            lon = ch.longitude if ch.longitude is not None else sta.longitude
            dist = gc_deg(elat, elon, lat, lon)
            if not (DIST_RANGE[0] <= dist <= DIST_RANGE[1]):
                continue
            dk = round(dist, 2)
            if dk not in ttcache:
                arr = model.get_travel_times(source_depth_in_km=0.0, distance_in_degree=dk,
                                             phase_list=["PKIKP", "PKP"])
                ttcache[dk] = branches([(a.name, a.time) for a in arr], dk)
            b = ttcache[dk]
            if b["df"] is None:
                continue
            end_rel = (b["ab"] if b["ab"] is not None else b["df"] + 60.0) + POST_S
            seen.add((net.code, sta.code))
            plan.append(dict(net=net.code, sta=sta.code, loc=pick[0], cha=pick[1], dist=dist,
                             b=b, t1=origin + b["df"] - PRE_S, t2=origin + end_rel,
                             epoch=str(ch.start_date)[:19] if ch.start_date else ""))
    edir = os.path.join(cache, eid)
    os.makedirs(edir, exist_ok=True)
    streams = {}
    todo = []
    for p in plan:
        path = os.path.join(edir, "%s.%s.%s.%s.mseed" % (p["net"], p["sta"], p["loc"] or "--",
                                                         p["cha"]))
        p["path"] = path
        if os.path.exists(path):
            try:
                streams[(p["net"], p["sta"])] = read(path)
                continue
            except Exception:  # noqa: BLE001
                pass
        todo.append(p)
    for k in range(0, len(todo), CHUNK):
        chunk = todo[k:k + CHUNK]
        bulk = [(p["net"], p["sta"], p["loc"] or "--", p["cha"], p["t1"], p["t2"]) for p in chunk]
        got = None
        for attempt in range(2):
            try:
                got = client.get_waveforms_bulk(bulk)
                break
            except Exception as ex:  # noqa: BLE001
                msg = repr(ex)[:160]
                if "No data" in msg or "204" in msg:
                    got = Stream()
                    break
                print("bulk error", eid, msg, flush=True)
                time.sleep(20 * (attempt + 1))
        if got is None:
            got = Stream()
            for p in chunk:
                try:
                    got += client.get_waveforms(p["net"], p["sta"], p["loc"] or "--", p["cha"],
                                                p["t1"], p["t2"])
                except Exception:  # noqa: BLE001
                    pass
        for p in chunk:
            sub = got.select(network=p["net"], station=p["sta"], channel=p["cha"])
            sub = Stream([t for t in sub if (t.stats.location or "") == (p["loc"] or "")])
            if len(sub):
                try:
                    sub.write(p["path"], format="MSEED")
                except Exception:  # noqa: BLE001
                    pass
                streams[(p["net"], p["sta"])] = sub
    for p in plan:
        row = dict(event=eid, date=t_iso[:10], station="%s.%s" % (p["net"], p["sta"]),
                   channel="%s.%s" % (p["loc"] or "--", p["cha"]), epoch_start=p["epoch"],
                   dist=round(p["dist"], 3), t_df=round(p["b"]["df"], 2),
                   t_bc=round(p["b"]["bc"], 2) if p["b"]["bc"] is not None else "",
                   t_ab=round(p["b"]["ab"], 2) if p["b"]["ab"] is not None else "",
                   windows_clear=windows_clear(p["b"]))
        st = streams.get((p["net"], p["sta"]))
        if st is None or not len(st):
            row["status"] = "no data"
            rows.append(row)
            continue
        st = st.copy()
        try:
            st.merge(method=0)
        except Exception as ex:  # noqa: BLE001
            row["status"] = "merge error " + repr(ex)[:60]
            rows.append(row)
            continue
        if len(st) != 1:
            row["status"] = "merge left %d traces" % len(st)
            rows.append(row)
            continue
        tr = st[0]
        if np.ma.isMaskedArray(tr.data):
            if np.ma.getmaskarray(tr.data).any():
                row["status"] = "gap"
                rows.append(row)
                continue
            tr.data = np.asarray(tr.data)
        row["fs"] = tr.stats.sampling_rate
        row["covers"] = bool(tr.stats.starttime <= p["t1"] + 5 and tr.stats.endtime >= p["t2"] - 5)
        try:
            for band in BANDS:
                bl = "%.1f-%.1f" % band
                v = snr_values(tr, origin, p["b"], band)
                row["snr_df_" + bl] = round(v["snr_df"], 3) if isf(v["snr_df"]) else ""
                row["snr_ab_" + bl] = round(v["snr_ab"], 3) if isf(v["snr_ab"]) else ""
            row["status"] = "ok"
        except Exception as ex:  # noqa: BLE001
            row["status"] = "processing error " + repr(ex)[:80]
        rows.append(row)
    print("EVENT", t_iso[:10], "listed", len(plan), "ok",
          sum(1 for r in rows if r.get("status") == "ok"), flush=True)
    return rows


def run(args):
    t0 = time.time()
    os.makedirs(args.out_dir, exist_ok=True)
    s = dict(version=VERSION, started=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
             blind="no two events are compared; SNR of single records only, written as bins",
             constants=dict(src=(SRC_LAT, SRC_LON), dist_range=DIST_RANGE, bands=BANDS,
                            noise=NOISE, sig=SIG, snr_min=SNR_MIN, snr_bins=SNR_BIN_EDGES,
                            chan_pref=CHAN_PREF, loc_pref=LOC_PREF),
             events=EVENTS)
    from obspy.clients.fdsn import Client
    from obspy.taup import TauPyModel
    client = Client(base_url=FDSN, _discover_services=False, user_agent=UA, timeout=180)
    model = TauPyModel("ak135")
    cache = os.path.join(args.data_dir, "dprk37")
    os.makedirs(cache, exist_ok=True)
    rows, ttcache = [], {}
    for ev in EVENTS:
        if (time.time() - t0) / 60.0 > args.time_budget_min:
            s.setdefault("skipped_time_budget", []).append(ev[0])
            continue
        try:
            rows.extend(probe_event(client, model, ev, cache, ttcache, s))
        except Exception as ex:  # noqa: BLE001
            traceback.print_exc()
            s.setdefault("event_errors", {})[ev[0]] = repr(ex)[:300]
    write_csv(os.path.join(args.out_dir, "dprk37_records.csv"), blind_rows(rows))
    pairs = []
    for band in BANDS:
        pairs.extend(pair_counts(rows, "%.1f-%.1f" % band))
    write_csv(os.path.join(args.out_dir, "dprk37_pair_counts.csv"), pairs)
    s["events_summary"] = {("%.1f-%.1f" % b): event_summary(rows, "%.1f-%.1f" % b) for b in BANDS}
    s["n_records"] = len(rows)
    s["minutes"] = round((time.time() - t0) / 60.0, 1)
    s["finished"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    with open(os.path.join(args.out_dir, "dprk37_summary.json"), "w", encoding="utf-8") as f:
        json.dump(s, f, indent=1, default=str)
    print(json.dumps(s["events_summary"], default=str)[:3000], flush=True)
    return 0


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--data-dir", default="/data")
    ap.add_argument("--out-dir", default="/out")
    ap.add_argument("--time-budget-min", type=float, default=150.0)
    return run(ap.parse_args(argv))


if __name__ == "__main__":
    sys.exit(main())

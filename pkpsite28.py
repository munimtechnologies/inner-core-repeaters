# pkpsite28.py 0.1.0 -- feasibility scan for re-testing Zhou (2026, JGR Solid Earth,
# doi:10.1029/2026JB033798): "Temporal changes in the Earth's outer core from nuclear-test PKP
# waves" (French Polynesian tests recorded at Borovoye, BRVK, 1977-1995).
#
# This job grades nothing. It measures the numbers a pre-registered test needs:
#   1. geometry and travel-time sensitivity of PKIKP/PKiKP/PKP minus PP (and PPP minus PP) to a
#      shift of the test position along the path to BRVK, and to shot depth (TauP iasp91, ak135);
#   2. the catalogue epicentres of the tests (USGS ComCat and ISC), how many share identical
#      coordinates, how far the two catalogues disagree, and what that disagreement means in
#      double-differential time;
#   3. which tests the public Borovoye archive copy on Zenodo (10.5281/zenodo.16872022) holds,
#      with the recording-system tag, channels, sampling and clipping of each file;
#   4. how many global stations (EarthScope FDSN metadata) had a vertical channel at each test,
#      i.e. whether the tests can be relocated relative to each other with teleseismic P;
#   5. a copy of the LDEO Borovoye archive index pages, for the next run to read.
#
# Part of munimtechnologies/inner-core-repeaters. Licence: Apache-2.0.
# Made with AI agents working in MV2's name (MV2, Munim, Inc.); no human expert has checked it.
import argparse
import csv
import io
import json
import math
import os
import sys
import time

import numpy as np

VERSION = "0.1.0"
BRVK = (53.0581, 70.2828)
MURUROA = (-21.83, -138.90)          # nominal centre of the atoll (lagoon)
FANGATAUFA = (-22.24, -138.74)
BOX = dict(minlatitude=-22.6, maxlatitude=-21.4, minlongitude=-139.4, maxlongitude=-138.3)
T_START, T_END = "1975-01-01", "1997-01-01"
KM_PER_DEG = 111.195
PHASES = ["PKIKP", "PKiKP", "PKP", "PP", "PPP", "PcP", "SKP", "PKKP"]
DIFFS = [("PKIKP", "PP"), ("PKiKP", "PP"), ("PPP", "PP"), ("PKiKP", "PKIKP")]
ZENODO_SUMMARY = ("https://zenodo.org/api/records/16872022/files/"
                  "BRVK_ASCII_file_summary.txt/content")
LDEO_PAGES = [
    "http://www.ldeo.columbia.edu/res/pi/Monitoring/Data/brv_exp_archive.html",
    "http://www.ldeo.columbia.edu/res/pi/Monitoring/Arch/brv_web/brv_archive_web.html",
    "http://www.ldeo.columbia.edu/res/pi/Monitoring/Data/",
]


def log(*a):
    print(time.strftime("%H:%M:%S"), "[pkpsite28]", *a, flush=True)


# ---------------------------------------------------------------- pure helpers (unit-tested)
def gc_dist_az(lat1, lon1, lat2, lon2):
    """Great-circle distance (deg), azimuth at point 1 and back-azimuth at point 2 (deg),
    spherical Earth (geocentric latitude ignored; fine for path geometry at the 0.1-deg level)."""
    p1, p2 = math.radians(lat1), math.radians(lat2)
    dl = math.radians(lon2 - lon1)
    c = math.sin(p1) * math.sin(p2) + math.cos(p1) * math.cos(p2) * math.cos(dl)
    d = math.degrees(math.acos(max(-1.0, min(1.0, c))))
    az = math.degrees(math.atan2(math.sin(dl) * math.cos(p2),
                                 math.cos(p1) * math.sin(p2) - math.sin(p1) * math.cos(p2) * math.cos(dl)))
    baz = math.degrees(math.atan2(-math.sin(dl) * math.cos(p1),
                                  math.cos(p2) * math.sin(p1) - math.sin(p2) * math.cos(p1) * math.cos(dl)))
    return d, az % 360.0, baz % 360.0


def offset_point(lat, lon, az, dist_km):
    """Point reached from (lat, lon) after dist_km along azimuth az (spherical)."""
    d = dist_km / 6371.0
    p1, l1, a = math.radians(lat), math.radians(lon), math.radians(az)
    p2 = math.asin(math.sin(p1) * math.cos(d) + math.cos(p1) * math.sin(d) * math.cos(a))
    l2 = l1 + math.atan2(math.sin(a) * math.sin(d) * math.cos(p1),
                         math.cos(d) - math.sin(p1) * math.sin(p2))
    return math.degrees(p2), (math.degrees(l2) + 540.0) % 360.0 - 180.0


def along_across_km(lat1, lon1, lat2, lon2, az):
    """Components (km) of the move from point 1 to point 2 along azimuth az and across it
    (positive to the right of az). Local flat approximation, good to <1% for < 50 km."""
    dn = (lat2 - lat1) * KM_PER_DEG
    de = ((lon2 - lon1 + 540.0) % 360.0 - 180.0) * KM_PER_DEG * math.cos(math.radians((lat1 + lat2) / 2))
    a = math.radians(az)
    along = dn * math.cos(a) + de * math.sin(a)
    across = -dn * math.sin(a) + de * math.cos(a)
    return along, across


def slope(xs, ys):
    """Least-squares slope of y on x, ignoring non-finite pairs; nan if < 2 points."""
    pts = [(x, y) for x, y in zip(xs, ys) if x is not None and y is not None
           and math.isfinite(x) and math.isfinite(y)]
    if len(pts) < 2:
        return float("nan")
    x = np.array([p[0] for p in pts], float)
    y = np.array([p[1] for p in pts], float)
    xm = x - x.mean()
    den = float(np.dot(xm, xm))
    return float(np.dot(xm, y - y.mean()) / den) if den > 0 else float("nan")


def parse_iso(s):
    """'1977-11-24T17:00:00.5Z' or with space -> POSIX seconds (UTC)."""
    s = s.strip().replace("Z", "").replace(" ", "T")
    if "." in s:
        main, frac = s.split(".", 1)
        frac = float("0." + "".join(ch for ch in frac if ch.isdigit()) or "0")
    else:
        main, frac = s, 0.0
    tm = time.strptime(main, "%Y-%m-%dT%H:%M:%S")
    return float(_timegm(tm)) + frac


def _timegm(tm):
    import calendar
    return calendar.timegm(tm)


def parse_zenodo_line(line):
    """One line of BRVK_ASCII_file_summary.txt -> dict, or None.
    Fields: T0 DT NSAMP EPOCH NEG% OK% POS% CHANNEL FILENAME; tag = token after '.BRV.' or
    '.brvk.' when it is alphabetic (SKM, TSG, SKD, KODB, KODM, ...), else ''."""
    parts = line.split()
    if len(parts) < 9 or "T" not in parts[0]:
        return None
    try:
        rec = dict(t0=parts[0], t0_s=parse_iso(parts[0]), dt=float(parts[1]), nsamp=int(parts[2]),
                   neg=float(parts[4]), ok=float(parts[5]), pos=float(parts[6]),
                   channel=parts[7], file=parts[8])
    except ValueError:
        return None
    toks = rec["file"].split(".")
    tag = ""
    for i, tok in enumerate(toks):
        if tok.lower() in ("brv", "brvk") and i + 1 < len(toks):
            nxt = toks[i + 1]
            if nxt.isalpha() and nxt.upper() == nxt:
                tag = nxt
            break
    rec["tag"] = tag
    return rec


def match_by_time(t_list, t, tol):
    """Indices i with |t_list[i] - t| <= tol."""
    return [i for i, x in enumerate(t_list) if abs(x - t) <= tol]


def close_pairs(events, max_deg):
    """Pairs (i, j), i < j, whose epicentres differ by <= max_deg in both latitude and
    longitude-on-the-sphere (great-circle separation <= max_deg)."""
    out = []
    for i in range(len(events)):
        for j in range(i + 1, len(events)):
            d, _, _ = gc_dist_az(events[i]["lat"], events[i]["lon"], events[j]["lat"], events[j]["lon"])
            if d <= max_deg + 1e-12:
                out.append((i, j, d))
    return out


# ---------------------------------------------------------------- network parts
def http_get(url, timeout=60, verify=True):
    import requests
    r = requests.get(url, timeout=timeout, verify=verify,
                     headers={"User-Agent": "inner-core-repeaters/pkpsite28 (research)"})
    return r


def part_geometry(out, summary):
    from obspy.taup import TauPyModel
    d0, az, baz = gc_dist_az(MURUROA[0], MURUROA[1], BRVK[0], BRVK[1])
    dF, azF, _ = gc_dist_az(FANGATAUFA[0], FANGATAUFA[1], BRVK[0], BRVK[1])
    geo = dict(dist_deg=d0, az_to_brvk=az, baz_at_brvk=baz, fangataufa_dist_deg=dF,
               fangataufa_az=azF)
    rows, sens = [], {}
    for mname in ("iasp91", "ak135"):
        model = TauPyModel(mname)
        # arrivals at the nominal distance
        arrs = model.get_travel_times(source_depth_in_km=0.6, distance_in_degree=d0, phase_list=PHASES)
        geo["arrivals_" + mname] = [dict(phase=a.name, time=a.time, p_s_per_deg=a.ray_param_sec_degree,
                                         incident=a.incident_angle) for a in arrs]
        # along-path shifts (positive = toward BRVK = smaller distance)
        shifts = list(np.arange(-15.0, 15.01, 1.0))
        tab = {ph: [] for ph in ("PKIKP", "PKiKP", "PP", "PPP")}
        for s in shifts:
            dd = d0 - s / KM_PER_DEG
            aa = model.get_travel_times(source_depth_in_km=0.6, distance_in_degree=dd,
                                        phase_list=list(tab))
            first = {}
            for a in aa:
                first.setdefault(a.name, a.time)
            for ph in tab:
                tab[ph].append(first.get(ph, float("nan")))
            rows.append(dict(model=mname, kind="along", shift_km=s, dist_deg=dd,
                             **{ph: first.get(ph, float("nan")) for ph in tab}))
        for a_, b_ in DIFFS:
            y = [x - z for x, z in zip(tab[a_], tab[b_])]
            sens["%s_%s-%s_ms_per_km_along" % (mname, a_, b_)] = 1000.0 * slope(shifts, y)
        # depth (0.0 to 2.0 km)
        depths = [0.0, 0.5, 1.0, 1.5, 2.0]
        dtab = {ph: [] for ph in ("PKIKP", "PKiKP", "PP", "PPP")}
        for z in depths:
            aa = model.get_travel_times(source_depth_in_km=z, distance_in_degree=d0, phase_list=list(dtab))
            first = {}
            for a in aa:
                first.setdefault(a.name, a.time)
            for ph in dtab:
                dtab[ph].append(first.get(ph, float("nan")))
            rows.append(dict(model=mname, kind="depth", shift_km=z, dist_deg=d0,
                             **{ph: first.get(ph, float("nan")) for ph in dtab}))
        for a_, b_ in DIFFS:
            y = [x - z for x, z in zip(dtab[a_], dtab[b_])]
            sens["%s_%s-%s_ms_per_km_depth" % (mname, a_, b_)] = 1000.0 * slope(depths, y)
        # where PKP (ab/bc) starts near the B caustic
        pk = []
        for dd in np.arange(140.0, 147.01, 0.25):
            aa = model.get_travel_times(source_depth_in_km=0.6, distance_in_degree=float(dd),
                                        phase_list=["PKP", "PKIKP", "PKiKP"])
            pk.append(dict(dist=float(dd), phases=sorted({a.name + ("@%.2f" % a.time) for a in aa})))
        geo["pkp_near_caustic_" + mname] = pk
    with open(os.path.join(out, "pkp28_sensitivity.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    geo["sensitivity"] = sens
    summary["geometry"] = geo
    return geo


def fetch_usgs():
    q = ("https://earthquake.usgs.gov/fdsnws/event/1/query?format=csv&starttime=%s&endtime=%s"
         "&minlatitude=%s&maxlatitude=%s&minlongitude=%s&maxlongitude=%s&minmagnitude=3.0"
         "&orderby=time-asc&limit=2000") % (T_START, T_END, BOX["minlatitude"], BOX["maxlatitude"],
                                          BOX["minlongitude"], BOX["maxlongitude"])
    r = http_get(q)
    r.raise_for_status()
    ev = []
    for row in csv.DictReader(io.StringIO(r.text)):
        try:
            ev.append(dict(src="usgs", id=row.get("id", ""), time=row["time"], t=parse_iso(row["time"]),
                           lat=float(row["latitude"]), lon=float(row["longitude"]),
                           depth=float(row["depth"]) if row.get("depth") else float("nan"),
                           mag=float(row["mag"]) if row.get("mag") else float("nan"),
                           magtype=row.get("magType", ""), type=row.get("type", ""),
                           net=row.get("net", ""), place=row.get("place", "")))
        except (KeyError, ValueError):
            continue
    return ev, q


def fetch_isc():
    from obspy.clients.fdsn import Client
    from obspy import UTCDateTime
    c = Client("ISC", timeout=180)
    cat = c.get_events(starttime=UTCDateTime(T_START), endtime=UTCDateTime(T_END),
                       minmagnitude=3.0, **BOX)
    ev = []
    for e in cat:
        o = e.preferred_origin() or (e.origins[0] if e.origins else None)
        if o is None:
            continue
        m = e.preferred_magnitude() or (e.magnitudes[0] if e.magnitudes else None)
        ev.append(dict(src="isc", id=str(e.resource_id).split("=")[-1], time=str(o.time),
                       t=float(o.time.timestamp), lat=float(o.latitude), lon=float(o.longitude),
                       depth=(float(o.depth) / 1000.0) if o.depth is not None else float("nan"),
                       mag=float(m.mag) if m is not None and m.mag is not None else float("nan"),
                       magtype=(m.magnitude_type or "") if m is not None else "",
                       type=str(e.event_type or ""), net="", place="",
                       lat_err_km=_unc(o, "latitude"), lon_err_km=_unc(o, "longitude")))
    return sorted(ev, key=lambda x: x["t"])


def _unc(o, which):
    try:
        u = getattr(o, which + "_errors").uncertainty
        if u is None:
            return float("nan")
        return float(u) * KM_PER_DEG if which == "latitude" else float(u) * KM_PER_DEG * math.cos(math.radians(o.latitude))
    except Exception:  # noqa: BLE001
        return float("nan")


def part_catalogues(out, summary, geo):
    az = geo["az_to_brvk"]
    sens = geo["sensitivity"]
    k_ikp = sens.get("iasp91_PKiKP-PP_ms_per_km_along", float("nan"))
    res = dict(errors=[])
    usgs, isc = [], []
    try:
        usgs, q = fetch_usgs()
        res["usgs_query"] = q
    except Exception as ex:  # noqa: BLE001
        res["errors"].append("usgs: " + repr(ex)[:200])
    try:
        isc = fetch_isc()
    except Exception as ex:  # noqa: BLE001
        res["errors"].append("isc: " + repr(ex)[:200])
    res["n_usgs"], res["n_isc"] = len(usgs), len(isc)
    for name, ev in (("usgs", usgs), ("isc", isc)):
        coords = {}
        for e in ev:
            coords.setdefault((round(e["lat"], 3), round(e["lon"], 3)), []).append(e["time"][:10])
        res[name + "_distinct_coords"] = len(coords)
        res[name + "_top_shared_coords"] = sorted(([list(k), len(v)] for k, v in coords.items()),
                                                  key=lambda x: -x[1])[:10]
        pairs = close_pairs(ev, 0.05)
        res[name + "_pairs_within_0p05deg"] = len(pairs)
        al = [abs(along_across_km(ev[i]["lat"], ev[i]["lon"], ev[j]["lat"], ev[j]["lon"], az)[0])
              for i, j, _ in pairs]
        if al:
            res[name + "_pair_along_km_median_max"] = [float(np.median(al)), float(np.max(al))]
            res[name + "_pair_implied_PKiKP-PP_ms_median_max"] = [float(np.median(al)) * abs(k_ikp),
                                                                  float(np.max(al)) * abs(k_ikp)]
    # cross-match the two catalogues
    rows = []
    ti = [e["t"] for e in isc]
    for e in usgs:
        m = match_by_time(ti, e["t"], 30.0)
        row = dict(time=e["time"], usgs_id=e["id"], usgs_lat=e["lat"], usgs_lon=e["lon"],
                   usgs_mag=e["mag"], usgs_magtype=e["magtype"], usgs_type=e["type"],
                   site=("fangataufa" if gc_dist_az(e["lat"], e["lon"], *FANGATAUFA)[0]
                         < gc_dist_az(e["lat"], e["lon"], *MURUROA)[0] else "mururoa"))
        if m:
            j = min(m, key=lambda k: abs(ti[k] - e["t"]))
            o = isc[j]
            al, ac = along_across_km(e["lat"], e["lon"], o["lat"], o["lon"], az)
            row.update(isc_id=o["id"], isc_time=o["time"], isc_lat=o["lat"], isc_lon=o["lon"],
                       isc_mag=o["mag"], isc_magtype=o["magtype"], isc_type=o["type"],
                       isc_lat_err_km=o.get("lat_err_km"), isc_lon_err_km=o.get("lon_err_km"),
                       dt_s=o["t"] - e["t"], along_km=al, across_km=ac,
                       implied_PKiKP_PP_ms=al * k_ikp)
        rows.append(row)
    matched = [r for r in rows if "along_km" in r]
    res["n_matched"] = len(matched)
    if matched:
        a = np.array([abs(r["along_km"]) for r in matched])
        res["usgs_vs_isc_along_km_median_p90_max"] = [float(np.median(a)), float(np.percentile(a, 90)),
                                                      float(a.max())]
    keys = []
    for r in rows:
        for k in r:
            if k not in keys:
                keys.append(k)
    with open(os.path.join(out, "pkp28_events.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys or ["time"])
        w.writeheader()
        for r in rows:
            w.writerow(r)
    with open(os.path.join(out, "pkp28_isc_events.csv"), "w", newline="") as f:
        if isc:
            w = csv.DictWriter(f, fieldnames=list(isc[0].keys()))
            w.writeheader()
            w.writerows(isc)
    summary["catalogues"] = res
    return usgs, isc


def part_brv_zenodo(out, summary, usgs, isc):
    res = dict(errors=[])
    try:
        r = http_get(ZENODO_SUMMARY, timeout=120)
        r.raise_for_status()
        recs = [x for x in (parse_zenodo_line(l) for l in r.text.splitlines()) if x]
    except Exception as ex:  # noqa: BLE001
        res["errors"].append("zenodo: " + repr(ex)[:200])
        recs = []
    res["n_files"] = len(recs)
    if recs:
        res["first_t0"], res["last_t0"] = recs[0]["t0"], recs[-1]["t0"]
    events = usgs if usgs else isc
    t0s = [x["t0_s"] for x in recs]
    rows = []
    for e in events:
        for i in match_by_time(t0s, e["t"], 600.0):
            x = recs[i]
            rows.append(dict(event_time=e["time"], event_mag=e["mag"], file_t0=x["t0"],
                             t0_minus_origin_s=x["t0_s"] - e["t"], dt=x["dt"], nsamp=x["nsamp"],
                             len_s=x["dt"] * x["nsamp"], tag=x["tag"], channel=x["channel"],
                             clip_neg=x["neg"], ok=x["ok"], clip_pos=x["pos"], file=x["file"]))
    res["n_files_matched"] = len(rows)
    res["n_events_with_files"] = len({r["event_time"] for r in rows})
    by_year = {}
    for r in rows:
        y = r["event_time"][:4]
        by_year.setdefault(y, set()).add(r["event_time"])
    res["events_with_files_by_year"] = {k: len(v) for k, v in sorted(by_year.items())}
    tags = {}
    for r in rows:
        tags.setdefault(r["tag"] or "untagged", set()).add(r["event_time"][:4])
    res["tags_years"] = {k: sorted(v) for k, v in tags.items()}
    with open(os.path.join(out, "pkp28_brv_zenodo.csv"), "w", newline="") as f:
        fn = ["event_time", "event_mag", "file_t0", "t0_minus_origin_s", "dt", "nsamp", "len_s", "tag",
              "channel", "clip_neg", "ok", "clip_pos", "file"]
        w = csv.DictWriter(f, fieldnames=fn)
        w.writeheader()
        w.writerows(rows)
    summary["brv_zenodo"] = res


def part_ldeo(out, summary):
    res = []
    for u in LDEO_PAGES:
        rec = dict(url=u)
        for verify in (True, False):
            try:
                r = http_get(u, timeout=60, verify=verify)
                rec.update(status=r.status_code, final_url=r.url, verify=verify, bytes=len(r.content))
                fn = "ldeo_" + str(len(res)) + ".html"
                with open(os.path.join(out, fn), "wb") as f:
                    f.write(r.content[:1500000])
                rec["saved"] = fn
                break
            except Exception as ex:  # noqa: BLE001
                rec["error_verify_%s" % verify] = repr(ex)[:200]
        res.append(rec)
    summary["ldeo"] = res


def part_stations(out, summary, usgs, isc, budget_s, t_start):
    events = usgs if usgs else isc
    rows, res = [], dict(errors=[])
    base = "https://service.earthscope.org/fdsnws/station/1/query"
    for e in events:
        if time.time() - t_start > budget_s:
            res["errors"].append("station scan stopped by time budget")
            break
        from obspy import UTCDateTime
        t = UTCDateTime(e["t"])
        q = ("%s?level=channel&format=text&channel=?HZ,?LZ&starttime=%s&endtime=%s&latitude=%.3f"
             "&longitude=%.3f&minradius=25&maxradius=100") % (
            base, (t).strftime("%Y-%m-%dT%H:%M:%S"), (t + 3600).strftime("%Y-%m-%dT%H:%M:%S"),
            e["lat"], e["lon"])
        row = dict(event_time=e["time"], mag=e["mag"])
        try:
            r = http_get(q, timeout=120)
            if r.status_code == 204:
                stas, chans = set(), set()
            else:
                r.raise_for_status()
                stas, chans, nets = set(), set(), set()
                for line in r.text.splitlines():
                    if line.startswith("#") or "|" not in line:
                        continue
                    p = line.split("|")
                    stas.add(p[0] + "." + p[1])
                    chans.add(p[3])
                    nets.add(p[0])
                row["networks"] = " ".join(sorted(nets))
            row["n_stations_25_100deg"] = len(stas)
            row["channels"] = " ".join(sorted(chans))
        except Exception as ex:  # noqa: BLE001
            row["error"] = repr(ex)[:160]
        try:
            q2 = ("%s?level=channel&format=text&station=BRVK&starttime=%s&endtime=%s") % (
                base, t.strftime("%Y-%m-%dT%H:%M:%S"), (t + 3600).strftime("%Y-%m-%dT%H:%M:%S"))
            r2 = http_get(q2, timeout=60)
            if r2.status_code == 200:
                row["brvk_fdsn_channels"] = " ".join(sorted({l.split("|")[0] + "." + l.split("|")[3]
                                                             for l in r2.text.splitlines()
                                                             if "|" in l and not l.startswith("#")}))
            else:
                row["brvk_fdsn_channels"] = ""
        except Exception as ex:  # noqa: BLE001
            row["brvk_error"] = repr(ex)[:160]
        rows.append(row)
        time.sleep(0.5)
    keys = []
    for r in rows:
        for k in r:
            if k not in keys:
                keys.append(k)
    with open(os.path.join(out, "pkp28_stations.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys or ["event_time"])
        w.writeheader()
        w.writerows(rows)
    by_year = {}
    for r in rows:
        by_year.setdefault(r["event_time"][:4], []).append(r.get("n_stations_25_100deg", 0) or 0)
    res["median_stations_by_year"] = {k: float(np.median(v)) for k, v in sorted(by_year.items())}
    res["n_events_scanned"] = len(rows)
    summary["stations"] = res


def _clean(o):
    if isinstance(o, float):
        return o if math.isfinite(o) else None
    if isinstance(o, dict):
        return {k: _clean(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [_clean(v) for v in o]
    return o


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-dir", default="out")
    ap.add_argument("--time-budget-min", type=float, default=90.0)
    a = ap.parse_args(argv)
    os.makedirs(a.out_dir, exist_ok=True)
    t_start = time.time()
    summary = dict(version=VERSION, started=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
                   errors=[], note="feasibility scan; grades nothing")

    def save():
        with open(os.path.join(a.out_dir, "pkp28_summary.json"), "w") as f:
            json.dump(_clean(summary), f, indent=1)

    geo, usgs, isc = None, [], []
    try:
        geo = part_geometry(a.out_dir, summary)
        log("geometry", round(geo["dist_deg"], 3), round(geo["az_to_brvk"], 2), geo["sensitivity"])
    except Exception as ex:  # noqa: BLE001
        summary["errors"].append("geometry: " + repr(ex)[:300])
    save()
    if geo is None:
        geo = dict(az_to_brvk=gc_dist_az(MURUROA[0], MURUROA[1], *BRVK)[1], sensitivity={})
    try:
        usgs, isc = part_catalogues(a.out_dir, summary, geo)
        log("catalogues", len(usgs), len(isc))
    except Exception as ex:  # noqa: BLE001
        summary["errors"].append("catalogues: " + repr(ex)[:300])
    save()
    try:
        part_brv_zenodo(a.out_dir, summary, usgs, isc)
        log("zenodo", summary["brv_zenodo"].get("n_events_with_files"))
    except Exception as ex:  # noqa: BLE001
        summary["errors"].append("zenodo: " + repr(ex)[:300])
    save()
    try:
        part_ldeo(a.out_dir, summary)
        log("ldeo", [x.get("status") for x in summary["ldeo"]])
    except Exception as ex:  # noqa: BLE001
        summary["errors"].append("ldeo: " + repr(ex)[:300])
    save()
    try:
        part_stations(a.out_dir, summary, usgs, isc, a.time_budget_min * 60.0, t_start)
        log("stations", summary["stations"].get("median_stations_by_year"))
    except Exception as ex:  # noqa: BLE001
        summary["errors"].append("stations: " + repr(ex)[:300])
    summary["finished"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    summary["runtime_min"] = round((time.time() - t_start) / 60.0, 2)
    save()
    log("done", summary["runtime_min"], "min; errors:", summary["errors"])
    return 0


if __name__ == "__main__":
    sys.exit(main())

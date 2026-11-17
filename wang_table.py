#!/usr/bin/env python3
"""Wang et al. (2024) repeater list -> USGS ComCat, plus new candidates at the same places.

Part of inner-core-repeaters by MV2, Munim, Inc. (made with AI agents). Licence: Apache-2.0.

Wang, Vidale et al. 2024 (Nature 631:340, doi:10.1038/s41586-024-07536-4) list their 121
South Sandwich Islands (SSI) events "in 42 locations" in Supplementary Table 1 (an xlsx).
This script

  1. downloads that workbook and dumps every sheet to CSV, cell for cell, so the table can
     be read and cited without Excel;
  2. finds each event in USGS ComCat (by origin time, then location), so the pipeline in
     run.py can be checked for recall (which of Wang's repeaters it finds) and so their
     ComCat magnitudes show what a magnitude cut loses;
  3. groups the matched events into locations (the table's own label column if it has one,
     else single-linkage on ComCat epicentres) and lists every ComCat event since
     2023-02-01 near each location: candidates for new, out-of-sample members of the
     published multiplets.

Nothing here measures waveforms; it only prepares the event lists.
"""
import argparse
import csv
import datetime as dt
import io
import json
import math
import os
import re
import sys
import time

USGS = "https://earthquake.usgs.gov/fdsnws/event/1/query"
XLSX_URLS = [
    "https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-024-07536-4/MediaObjects/"
    "41586_2024_7536_MOESM1_ESM.xlsx",
    "https://pmc.ncbi.nlm.nih.gov/articles/instance/11236701/bin/41586_2024_7536_MOESM1_ESM.xlsx",
]
BOX = dict(minlatitude=-62.5, maxlatitude=-53.0, minlongitude=-33.0, maxlongitude=-19.0)
LAT_RANGE = (-62.5, -53.0)
LON_RANGE = (-33.0, -19.0)
UA = "inner-core-repeaters/wang_table (MV2, Munim, Inc.; mv2@munimtech.com)"
# run.py v0.1.0 catalogue cuts (to flag what the first ILAR run could not see)
V1_NEW_FROM = 2024.5
V1_MIN_OLD, V1_MIN_NEW = 5.0, 4.7


def log(*a):
    print(*a, flush=True)


# ----------------------------------------------------------------------------- parsing helpers
_ISO = re.compile(r"((?:19|20)\d\d)[-/.](\d{1,2})[-/.](\d{1,2})(?:[ T_,]+(\d{1,2}):(\d{1,2})(?::(\d{1,2}(?:\.\d*)?))?)?")
_JUL = re.compile(r"\b((?:19|20)\d\d)[-/.:](\d{3})(?:[ T_,/:-]+(\d{1,2}):?(\d{2}):?(\d{2}(?:\.\d*)?)?)?\b")
_COMPACT = re.compile(r"(?<!\d)((?:19|20)\d\d)(\d\d)(\d\d)(?:[ T_.-]?(\d\d)(\d\d)(\d\d(?:\.\d+)?)?)?(?!\d)")
_HMS = re.compile(r"\b(\d{1,2}):(\d{2})(?::(\d{2}(?:\.\d*)?))?\b")


def _mk(y, mo, d, h=0, mi=0, s=0.0):
    try:
        y, mo, d, h, mi = int(y), int(mo), int(d), int(h or 0), int(mi or 0)
        s = float(s or 0.0)
        if not (1985 <= y <= 2027 and 1 <= mo <= 12 and 1 <= d <= 31 and 0 <= h < 24 and 0 <= mi < 60
                and 0 <= s < 61):
            return None
        return dt.datetime(y, mo, d, h, mi, tzinfo=dt.timezone.utc) + dt.timedelta(seconds=s)
    except (ValueError, TypeError):
        return None


def parse_time_text(s):
    """Return (datetime, has_time) from free text, or (None, False)."""
    s = str(s or "").strip()
    if not s:
        return None, False
    m = _ISO.search(s)
    if m:
        t = _mk(m.group(1), m.group(2), m.group(3), m.group(4), m.group(5), m.group(6))
        if t is not None:
            return t, m.group(4) is not None
    m = _COMPACT.search(s)
    if m:
        t = _mk(m.group(1), m.group(2), m.group(3), m.group(4), m.group(5), m.group(6))
        if t is not None:
            return t, m.group(4) is not None
    m = _JUL.search(s)
    if m:
        try:
            y, j = int(m.group(1)), int(m.group(2))
            if 1 <= j <= 366:
                base = dt.datetime(y, 1, 1, tzinfo=dt.timezone.utc) + dt.timedelta(days=j - 1)
                t = _mk(base.year, base.month, base.day, m.group(3), m.group(4), m.group(5))
                if t is not None:
                    return t, m.group(3) is not None
        except ValueError:
            pass
    return None, False


def _num(v):
    if isinstance(v, bool):
        return None
    if isinstance(v, (int, float)):
        return float(v) if math.isfinite(float(v)) else None
    try:
        return float(str(v).strip())
    except (ValueError, TypeError):
        return None


def row_event(cells, header=None):
    """Best-effort (time, has_time, lat, lon, mag, label) from one spreadsheet row.

    cells: list of cell values (datetime, date, time, float, int, str or None).
    header: list of lower-case header strings (same length) or None.
    """
    header = header or [""] * len(cells)
    t, has_time = None, False
    tcell = None
    # 1. real datetime/date cells
    for v in cells:
        if isinstance(v, dt.datetime) and 1985 <= v.year <= 2027:
            t = v.replace(tzinfo=dt.timezone.utc) if v.tzinfo is None else v.astimezone(dt.timezone.utc)
            has_time = (v.hour, v.minute, v.second, v.microsecond) != (0, 0, 0, 0)
            break
        if isinstance(v, dt.date) and not isinstance(v, dt.datetime) and 1985 <= v.year <= 2027:
            t = dt.datetime(v.year, v.month, v.day, tzinfo=dt.timezone.utc)
            break
    for v in cells:
        if isinstance(v, dt.time):
            tcell = v
            break
    # 2. split Y / M / D / h / m / s columns
    if t is None:
        idx = {}
        for k, h in enumerate(header):
            hh = h.strip().lower()
            for key, names in (("y", ("year", "yr", "yyyy")), ("mo", ("month", "mon", "mm")),
                               ("d", ("day", "dd")), ("h", ("hour", "hr", "hh")),
                               ("mi", ("minute", "min")), ("s", ("second", "sec", "ss"))):
                if hh in names and key not in idx:
                    idx[key] = k
        if all(k in idx for k in ("y", "mo", "d")):
            g = {k: _num(cells[i]) for k, i in idx.items()}
            if all(g.get(k) is not None for k in ("y", "mo", "d")):
                t = _mk(g["y"], g["mo"], g["d"], g.get("h") or 0, g.get("mi") or 0, g.get("s") or 0.0)
                has_time = t is not None and "h" in idx
    # 3. free text, cell by cell, then the joined row
    if t is None:
        for v in cells:
            if isinstance(v, str):
                t, has_time = parse_time_text(v)
                if t is not None:
                    break
    if t is None:
        t, has_time = parse_time_text(" ".join(str(v) for v in cells if v is not None and not isinstance(v, float)))
    if t is not None and not has_time:
        tt = None
        if tcell is not None:
            tt = (tcell.hour, tcell.minute, tcell.second + tcell.microsecond * 1e-6)
        else:
            for v in cells:
                if isinstance(v, str):
                    m = _HMS.search(v)
                    if m and not _ISO.search(v):
                        tt = (int(m.group(1)), int(m.group(2)), float(m.group(3) or 0))
                        break
        if tt and 0 <= tt[0] < 24 and 0 <= tt[1] < 60:
            t = t.replace(hour=0, minute=0, second=0, microsecond=0) + dt.timedelta(
                hours=tt[0], minutes=tt[1], seconds=tt[2])
            has_time = True
    # lat / lon / mag by value range and header
    lat = lon = mag = None
    label = None
    for k, v in enumerate(cells):
        h = header[k].strip().lower() if k < len(header) else ""
        x = _num(v)
        if (any(w in h for w in ("multiplet", "location", "group", "cluster", "family", "doublet"))
                or h in ("loc", "loc.", "site")) and v not in (None, "") and label is None:
            label = str(v).strip()
        if x is None:
            continue
        if lat is None and LAT_RANGE[0] <= x <= LAT_RANGE[1] and "lon" not in h:
            lat = x
        elif lon is None and (LON_RANGE[0] <= x <= LON_RANGE[1] or 360 + LON_RANGE[0] <= x <= 360 + LON_RANGE[1]):
            lon = x - 360 if x > 180 else x
        elif mag is None and ("mag" in h or h in ("mb", "ms", "mw", "m")) and 3.0 <= x <= 8.5:
            mag = x
    return t, has_time, lat, lon, mag, label


def gc_km(lat1, lon1, lat2, lon2):
    p1, p2 = math.radians(lat1), math.radians(lat2)
    dp, dl = p2 - p1, math.radians(lon2 - lon1)
    h = math.sin(dp / 2) ** 2 + math.cos(p1) * math.cos(p2) * math.sin(dl / 2) ** 2
    return 6371.0 * 2 * math.asin(min(1.0, math.sqrt(h)))


def decimal_year(t):
    a = dt.datetime(t.year, 1, 1, tzinfo=dt.timezone.utc)
    b = dt.datetime(t.year + 1, 1, 1, tzinfo=dt.timezone.utc)
    return t.year + (t - a).total_seconds() / (b - a).total_seconds()


def match_event(ev, catalog, max_dt_s=90.0):
    """ev: dict(time, has_time, lat, lon). catalog: list of dict(time, lat, lon, mag, id).

    With an origin time: the nearest ComCat event within max_dt_s (ties broken by distance).
    With a date only: the events that day; the nearest one if lat/lon known, the largest if not.
    Returns (catalog row or None, how)."""
    t = ev["time"]
    if t is None:
        return None, "no-time"
    if ev["has_time"]:
        c = [(abs((r["time"] - t).total_seconds()), r) for r in catalog
             if abs((r["time"] - t).total_seconds()) <= max_dt_s]
        if not c:
            return None, "no-comcat-within-%ds" % int(max_dt_s)
        if ev.get("lat") is not None and ev.get("lon") is not None:
            c.sort(key=lambda x: (x[0] > 20.0, gc_km(ev["lat"], ev["lon"], x[1]["lat"], x[1]["lon"]), x[0]))
        else:
            c.sort(key=lambda x: x[0])
        return c[0][1], "time(dt=%.1fs)" % c[0][0]
    day = [r for r in catalog if r["time"].date() == t.date()]
    if not day:
        return None, "no-comcat-that-day"
    if ev.get("lat") is not None and ev.get("lon") is not None:
        day.sort(key=lambda r: gc_km(ev["lat"], ev["lon"], r["lat"], r["lon"]))
        return day[0], "date+location(n=%d)" % len(day)
    day.sort(key=lambda r: -(r["mag"] if r["mag"] is not None else 0))
    return day[0], "date-only-largest(n=%d)" % len(day)


def single_link(points, km):
    """points: list of (lat, lon). Returns a cluster index per point."""
    n = len(points)
    parent = list(range(n))

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    for i in range(n):
        for j in range(i + 1, n):
            if gc_km(points[i][0], points[i][1], points[j][0], points[j][1]) <= km:
                ri, rj = find(i), find(j)
                if ri != rj:
                    parent[rj] = ri
    roots = {}
    return [roots.setdefault(find(i), len(roots)) for i in range(n)]


def in_v1_catalog(t_dec, mag):
    if mag is None or not math.isfinite(mag):
        return False
    return mag >= (V1_MIN_NEW if t_dec >= V1_NEW_FROM else V1_MIN_OLD)


# ----------------------------------------------------------------------------- I/O
def http_get(url, params=None, timeout=180, tries=4):
    import requests

    last = None
    for k in range(tries):
        try:
            r = requests.get(url, params=params, timeout=timeout, headers={"User-Agent": UA})
            if r.status_code == 200 and r.content:
                return r
            last = "status %d" % r.status_code
        except Exception as ex:  # network hiccup
            last = repr(ex)[:200]
        time.sleep(5 * (k + 1))
    raise RuntimeError("GET failed %s: %s" % (url, last))


def comcat(start, end, minmag):
    r = http_get(USGS, dict(format="csv", starttime=start, endtime=end, minmagnitude=minmag,
                            orderby="time-asc", **BOX))
    out = []
    for row in csv.DictReader(io.StringIO(r.text)):
        if row.get("type", "earthquake") != "earthquake":
            continue
        try:
            t = dt.datetime.fromisoformat(row["time"].replace("Z", "+00:00"))
            out.append(dict(id=row["id"], time=t, lat=float(row["latitude"]), lon=float(row["longitude"]),
                            depth=_num(row.get("depth")), mag=_num(row.get("mag")),
                            magType=row.get("magType", "")))
        except (KeyError, ValueError):
            continue
    return out


def write_csv(path, rows, keys=None):
    with open(path, "w", newline="", encoding="utf-8") as f:
        if not rows:
            return
        keys = keys or list(rows[0].keys())
        w = csv.DictWriter(f, fieldnames=keys, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)


def cell_str(v):
    if v is None:
        return ""
    if isinstance(v, (dt.datetime, dt.date, dt.time)):
        return v.isoformat()
    return str(v)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", default="/data")
    ap.add_argument("--out-dir", default="/out")
    ap.add_argument("--xlsx", default="", help="local copy of the workbook (else downloaded)")
    ap.add_argument("--cand-from", default="2023-02-01")
    ap.add_argument("--cand-min-mag", type=float, default=3.5)
    ap.add_argument("--end", default=time.strftime("%Y-%m-%d", time.gmtime()))
    ap.add_argument("--link-km", type=float, default=20.0)
    args = ap.parse_args()
    os.makedirs(args.data_dir, exist_ok=True)
    os.makedirs(args.out_dir, exist_ok=True)
    summary = dict(started=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()))

    # 1. workbook
    path = args.xlsx or os.path.join(args.data_dir, "wang2024_supp.xlsx")
    if not os.path.exists(path):
        for u in XLSX_URLS:
            try:
                r = http_get(u)
                if r.content[:2] == b"PK":
                    with open(path, "wb") as f:
                        f.write(r.content)
                    summary["xlsx_url"] = u
                    break
                log("not an xlsx from", u, r.content[:80])
            except Exception as ex:
                log("download failed", u, repr(ex)[:200])
    if not os.path.exists(path):
        log("could not get the workbook")
        summary["error"] = "workbook unavailable"
        json.dump(summary, open(os.path.join(args.out_dir, "wang_summary.json"), "w"), indent=1)
        sys.exit(1)
    import openpyxl

    wb = openpyxl.load_workbook(path, data_only=True, read_only=True)
    sheets = {}
    for k, ws in enumerate(wb.worksheets):
        rows = [list(r) for r in ws.iter_rows(values_only=True)]
        sheets[ws.title] = rows
        safe = re.sub(r"[^A-Za-z0-9]+", "_", ws.title).strip("_") or "sheet"
        with open(os.path.join(args.out_dir, "wang_sheet%d_%s.csv" % (k + 1, safe)), "w", newline="",
                  encoding="utf-8") as f:
            w = csv.writer(f)
            for r in rows:
                w.writerow([cell_str(v) for v in r])
        log("sheet %d %r: %d rows x %d cols" % (k + 1, ws.title, len(rows), max((len(r) for r in rows), default=0)))
    summary["sheets"] = {name: len(rows) for name, rows in sheets.items()}

    # 2. events: the sheet that yields the most rows with a time and an SSI lat/lon (or time only)
    best = None
    for name, rows in sheets.items():
        header = None
        evs = []
        for r in rows:
            strs = [v for v in r if isinstance(v, str) and v.strip()]
            t, has_time, lat, lon, mag, label = row_event(r, header)
            if t is None and len(strs) >= 3:
                header = [cell_str(v).lower() for v in r]
                continue
            if t is not None:
                evs.append(dict(sheet=name, raw=" | ".join(cell_str(v) for v in r if v is not None)[:300],
                                time=t, has_time=has_time, lat=lat, lon=lon, mag_table=mag, label=label))
        # de-duplicate (an event can appear in several rows if the sheet lists pairs)
        seen, uniq = set(), []
        for e in evs:
            key = (e["time"].strftime("%Y%m%d%H%M%S") if e["has_time"]
                   else e["time"].strftime("%Y%m%d") + str(e["lat"]))
            if key not in seen:
                seen.add(key)
                uniq.append(e)
        evs = uniq
        score = sum(1 for e in evs if e["has_time"]) * 2 + len(evs)
        log("sheet %r: %d distinct events (%d with origin time)" % (name, len(evs),
                                                                   sum(e["has_time"] for e in evs)))
        if evs and (best is None or score > best[0]):
            best = (score, name, evs)
    if best is None:
        log("no event rows recognised; see the sheet CSVs")
        summary["error"] = "no event rows recognised"
        json.dump(summary, open(os.path.join(args.out_dir, "wang_summary.json"), "w"), indent=1)
        return
    _, ev_sheet, wevs = best
    summary["event_sheet"] = ev_sheet
    summary["n_table_events"] = len(wevs)
    log("event sheet %r: %d distinct events" % (ev_sheet, len(wevs)))

    cat_old = comcat("1989-01-01", "2023-03-01", 4.0)  # stays well under ComCat's 20,000-row limit
    log("ComCat 1989-2023 M>=4.0 in box:", len(cat_old))
    for e in wevs:
        r, how = match_event(e, cat_old)
        e["match"] = how
        if r is not None:
            e.update(comcat_id=r["id"], comcat_time=r["time"].isoformat(), comcat_lat=r["lat"],
                     comcat_lon=r["lon"], comcat_depth=r["depth"], comcat_mag=r["mag"],
                     comcat_magType=r["magType"],
                     sep_km=round(gc_km(e["lat"], e["lon"], r["lat"], r["lon"]), 1)
                     if e["lat"] is not None and e["lon"] is not None else None,
                     t_dec=round(decimal_year(r["time"]), 4),
                     in_v1_catalog=in_v1_catalog(decimal_year(r["time"]), r["mag"]))
    matched = [e for e in wevs if e.get("comcat_id")]
    mags = sorted(e["comcat_mag"] for e in matched if e.get("comcat_mag") is not None)
    summary.update(n_matched=len(matched),
                   n_in_v1_catalog=sum(1 for e in matched if e.get("in_v1_catalog")),
                   comcat_mag_quartiles=[mags[int(q * (len(mags) - 1))] for q in (0, 0.25, 0.5, 0.75, 1)]
                   if mags else None,
                   years=sorted({e["time"].year for e in wevs}))
    log("matched %d/%d; inside the v1 catalogue cuts: %d" % (len(matched), len(wevs), summary["n_in_v1_catalog"]))

    # 3. locations
    labels = [e.get("label") for e in matched]
    use_labels = (sum(1 for x in labels if x) >= 0.8 * max(1, len(matched)) and len(set(labels)) >= 5
                  and len(set(labels)) <= 0.6 * len(matched))  # a per-event ID column is not a location
    if use_labels:
        for e in matched:
            e["loc"] = "T:" + str(e["label"])
    else:
        cl = single_link([(e["comcat_lat"], e["comcat_lon"]) for e in matched], args.link_km)
        for e, c in zip(matched, cl):
            e["loc"] = "C%02d" % c
    summary["location_source"] = "table label column" if use_labels else "single-link %.0f km on ComCat" % args.link_km
    locs = {}
    for e in matched:
        locs.setdefault(e["loc"], []).append(e)
    summary["n_locations"] = len(locs)
    log("locations:", len(locs), "(%s)" % summary["location_source"])

    keys = ["sheet", "time", "has_time", "lat", "lon", "mag_table", "label", "loc", "match", "comcat_id",
            "comcat_time", "comcat_lat", "comcat_lon", "comcat_depth", "comcat_mag", "comcat_magType", "sep_km",
            "t_dec", "in_v1_catalog", "raw"]
    write_csv(os.path.join(args.out_dir, "wang_events_comcat.csv"),
              [dict(e, time=e["time"].isoformat()) for e in wevs], keys=keys)
    with open(os.path.join(args.out_dir, "wang_comcat_ids.txt"), "w") as f:
        for e in matched:
            f.write(e["comcat_id"] + "\n")

    # 4. candidates since cand-from at each location
    cat_new = comcat(args.cand_from, args.end, args.cand_min_mag)
    log("ComCat %s..%s M>=%.1f in box: %d" % (args.cand_from, args.end, args.cand_min_mag, len(cat_new)))
    cand, loc_rows = [], []
    for name, mem in sorted(locs.items()):
        la = sum(e["comcat_lat"] for e in mem) / len(mem)
        lo = sum(e["comcat_lon"] for e in mem) / len(mem)
        spread = max(gc_km(la, lo, e["comcat_lat"], e["comcat_lon"]) for e in mem)
        rad = min(50.0, max(25.0, 1.5 * spread))
        last = max(e["comcat_time"] for e in mem)
        near = []
        for r in cat_new:
            d = gc_km(la, lo, r["lat"], r["lon"])
            if d <= rad and r["time"].isoformat() > last:
                near.append((d, r))
        loc_rows.append(dict(loc=name, n=len(mem), lat=round(la, 3), lon=round(lo, 3), spread_km=round(spread, 1),
                             radius_km=rad, first=min(e["comcat_time"] for e in mem)[:10], last=last[:10],
                             members=";".join(e["comcat_id"] for e in sorted(mem, key=lambda x: x["comcat_time"])),
                             n_candidates=len(near)))
        for d, r in sorted(near, key=lambda x: x[1]["time"]):
            td = decimal_year(r["time"])
            cand.append(dict(loc=name, loc_n=len(mem), id=r["id"], time=r["time"].isoformat(), lat=r["lat"],
                             lon=r["lon"], depth=r["depth"], mag=r["mag"], magType=r["magType"],
                             dist_km=round(d, 1), t_dec=round(td, 4), in_v1_catalog=in_v1_catalog(td, r["mag"])))
    write_csv(os.path.join(args.out_dir, "wang_locations.csv"), loc_rows)
    write_csv(os.path.join(args.out_dir, "candidates_since_2023.csv"), cand)
    summary.update(n_candidates=len(cand), n_candidates_2025on=sum(1 for c in cand if c["t_dec"] >= 2025.0),
                   n_candidates_in_v1=sum(1 for c in cand if c["in_v1_catalog"]),
                   n_locations_with_2025on=len({c["loc"] for c in cand if c["t_dec"] >= 2025.0}))
    summary["finished"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    with open(os.path.join(args.out_dir, "wang_summary.json"), "w") as f:
        json.dump(summary, f, indent=1, default=str)
    log(json.dumps(summary, indent=1, default=str))


if __name__ == "__main__":
    main()

# pkploc29.py 0.1.0 -- where were the bombs? Position-only predictions of PKP-PP at Borovoye
# for the French underground tests behind Zhou (2026, JGR Solid Earth, doi:10.1029/2026JB033798,
# "Temporal changes in the Earth's outer core from nuclear-test PKP waves").
#
# Zhou pairs tests whose catalogue epicentres (Yang et al. 2002 / CMR database) lie within
# 0.05 deg and reads 0.1-0.2 s changes of PKrKP-PP and PKiKP-PP as outer-core change. PKP and
# PP leave the source with very different ray parameters, so any real difference in shot
# position along the path to BRVK changes PKP-PP by tens of ms per km.
#
# This job grades nothing. It builds the inputs of pre-registration 28:
#   1. the shaft-site table (mururoa_shafts_1977_1996.csv, transcribed from the English
#      Wikipedia series pages, which compile French Ministry of Defence / IAEA site maps)
#      checked coordinate by coordinate against the live wikitext;
#   2. the Yang, North & Romney "Worldwide nuclear explosions" appendix (LDEO copy) rows for
#      French Polynesia, as text, parsed where possible;
#   3. ISC and USGS catalogue epicentres for the same tests (matched by origin time);
#   4. the official French list (CIVEN) as text, for zone cross-checks;
#   5. an OCR text of Blacknest AWE Report O2/95 (UK array epicentres of French tests);
#   6. per test and per location set: distance to BRVK and exact TauP (iasp91) PKIKP-PP,
#      PKiKP-PP and PPP-PP relative to a fixed reference point, and per-epoch means for
#      Zhou's epochs (1977, 1982-83, 1988-90, 1995).
#
# Part of munimtechnologies/inner-core-repeaters. Licence: Apache-2.0 (code).
# mururoa_shafts_1977_1996.csv is derived from Wikipedia (CC BY-SA 4.0).
# Made with AI agents working in MV2's name (MV2, Munim, Inc.); no human expert has checked it.
import argparse
import csv
import json
import math
import os
import re
import subprocess
import sys
import time

import numpy as np

VERSION = "0.1.0"
BRVK = (53.0581, 70.2828)
REF = (-21.835, -138.91)              # Blacknest AWE O2/95 reference point: centre of Mururoa
KM_PER_DEG = 111.195
SRC_DEPTH_KM = 0.6
PHASES = ["PKIKP", "PKiKP", "PP", "PPP"]
EPOCHS = {"E1977": (1977, 1977), "E1982_83": (1982, 1983), "E1988_90": (1988, 1990),
          "E1995": (1995, 1995)}
WIKI_PAGES = ["1975–78_French_nuclear_tests", "1979–80_French_nuclear_tests",
              "1981–82_French_nuclear_tests", "1983–85_French_nuclear_tests",
              "1986–88_French_nuclear_tests", "1989–91_French_nuclear_tests",
              "1995–96_French_nuclear_tests"]
YANG_PDF = "https://www.ldeo.columbia.edu/~richards/my_papers/WW_nuclear_tests_IASPEI_HB.pdf"
CIVEN_PDF = "https://www.civen.fr/sites/civen/files/documents/liste-des-essais-nucleaires-francais_0.pdf"
AVEN_PDF = "https://www.aven.org/wp-content/uploads/2025/05/liste-complete-des-essais-francais.pdf"
BLACKNEST_PDF = "https://www.blacknest.gov.uk/publications/oreports/AWE_Report_O2_95.pdf"
MAX_TXT = 1800000


def log(*a):
    print(time.strftime("%H:%M:%S"), "[pkploc29]", *a, flush=True)


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


def along_across_km(lat1, lon1, lat2, lon2, az):
    """Move from point 1 to point 2 split into km along azimuth az and across it (right +)."""
    dn = (lat2 - lat1) * KM_PER_DEG
    de = ((lon2 - lon1 + 540.0) % 360.0 - 180.0) * KM_PER_DEG * math.cos(math.radians((lat1 + lat2) / 2))
    a = math.radians(az)
    return dn * math.cos(a) + de * math.sin(a), -dn * math.sin(a) + de * math.cos(a)


def dms_to_dec(d, m, s, hemi):
    v = float(d) + float(m or 0) / 60.0 + float(s or 0) / 3600.0
    return -v if hemi.upper() in ("S", "W") else v


def parse_coord_templates(wikitext):
    """All {{coord|...}} templates in wikitext -> list of (lat, lon) decimal degrees.
    Handles decimal (|21.87|S|138.86|W), DMS (|21|51|10|S|138|53|39|W), DM and signed decimal."""
    out = []
    for m in re.finditer(r"\{\{\s*[Cc]oord\s*\|([^{}]*)\}\}", wikitext):
        parts = [p.strip() for p in m.group(1).split("|")]
        vals = [p for p in parts if "=" not in p]
        hemis = [i for i, p in enumerate(vals) if p.upper() in ("N", "S", "E", "W")]
        try:
            if len(hemis) >= 2:
                i1, i2 = hemis[0], hemis[1]
                lat = dms_to_dec(*(vals[:i1] + ["0", "0"])[:3], vals[i1])
                lon_parts = vals[i1 + 1:i2]
                lon = dms_to_dec(*(lon_parts + ["0", "0"])[:3], vals[i2])
            else:
                lat, lon = float(vals[0]), float(vals[1])
        except (ValueError, IndexError, TypeError):
            continue
        out.append((lat, lon))
    return out


def verify_rows(rows, coords, tol=0.0006):
    """For each row (dict with lat, lon) say whether some wikitext coordinate is within tol."""
    res = []
    for r in rows:
        ok = any(abs(r["lat"] - a) <= tol and abs(r["lon"] - b) <= tol for a, b in coords)
        res.append(ok)
    return res


_DATE_RES = [
    re.compile(r"\b(19[789]\d)[/\-\.](\d{1,2})[/\-\.](\d{1,2})\b"),
    re.compile(r"\b(\d{1,2})[/\-\.](\d{1,2})[/\-\.](19[789]\d)\b"),
    re.compile(r"\b(19[789]\d)(\d\d)(\d\d)\b"),
    re.compile(r"\b(19[789]\d)\s+(\d{1,2})\s+(\d{1,2})\b"),
]


def find_date(line):
    """First plausible date (1970-1999) in a text line -> 'YYYY-MM-DD' or None."""
    for k, rx in enumerate(_DATE_RES):
        m = rx.search(line)
        if not m:
            continue
        g = m.groups()
        if k == 1:
            d, mo, y = int(g[0]), int(g[1]), int(g[2])
        else:
            y, mo, d = int(g[0]), int(g[1]), int(g[2])
        if 1 <= mo <= 12 and 1 <= d <= 31:
            return "%04d-%02d-%02d" % (y, mo, d)
    return None


def find_latlon(line):
    """(lat, lon) for French Polynesia in a text line, or None. Accepts signed decimals
    (-21.83 -138.92) or unsigned with S/W letters (21.83S 138.92W)."""
    lat = lon = None
    for m in re.finditer(r"(-?\d{2,3}\.\d+)\s*([NSEW])?", line):
        v, h = float(m.group(1)), (m.group(2) or "")
        if h == "S" or (not h and -23.5 < v < -21.0):
            if 21.0 < abs(v) < 23.5 and lat is None:
                lat = -abs(v)
                continue
        if h == "W" or (not h and -140.5 < v < -138.0):
            if 138.0 < abs(v) < 140.5 and lon is None:
                lon = -abs(v)
    if lat is None or lon is None:
        return None
    return lat, lon


def epoch_of(year):
    for k, (a, b) in EPOCHS.items():
        if a <= year <= b:
            return k
    return None


def epoch_means(rows, key, min_yield=None, atoll="Mururoa"):
    """Mean of rows[key] per Zhou epoch (finite values only), and each epoch minus E1977."""
    groups = {k: [] for k in EPOCHS}
    for r in rows:
        if atoll and r.get("atoll") != atoll:
            continue
        if min_yield is not None and not (r.get("yield_kt") is not None and r["yield_kt"] >= min_yield):
            continue
        v = r.get(key)
        e = epoch_of(int(r["date"][:4]))
        if e and v is not None and math.isfinite(v):
            groups[e].append(v)
    out = {}
    for k, v in groups.items():
        out[k] = dict(n=len(v), mean=(float(np.mean(v)) if v else None),
                      sd=(float(np.std(v, ddof=1)) if len(v) > 1 else None))
    base = out["E1977"]["mean"]
    for k in out:
        m = out[k]["mean"]
        out[k]["minus_E1977"] = (m - base) if (m is not None and base is not None) else None
    return out


def load_shafts(path):
    rows = []
    with open(path, encoding="utf-8") as f:
        for r in csv.DictReader(f):
            r["lat"], r["lon"] = float(r["lat"]), float(r["lon"])
            r["approx"] = int(r["approx"] or 0)
            r["yield_kt"] = float(r["yield_kt"]) if r["yield_kt"] else None
            rows.append(r)
    return rows


def origin_seconds(date, time_ut):
    """POSIX seconds of date + 'HH:MM[:SS.s]' (UTC); None if time is empty."""
    import calendar
    if not time_ut:
        return None
    p = time_ut.split(":")
    h, mi = int(p[0]), int(p[1])
    s = float(p[2]) if len(p) > 2 and p[2] else 0.0
    y, mo, d = (int(x) for x in date.split("-"))
    return calendar.timegm((y, mo, d, h, mi, 0, 0, 0, 0)) + s


# ---------------------------------------------------------------- network / tool parts
def http_get(url, timeout=120, verify=True):
    import requests
    return requests.get(url, timeout=timeout, verify=verify,
                        headers={"User-Agent": "inner-core-repeaters/pkploc29 (research)"})


def fetch_bytes(url, timeout=180):
    import urllib3
    urllib3.disable_warnings()
    last = None
    for verify in (True, False):
        try:
            r = http_get(url, timeout=timeout, verify=verify)
            r.raise_for_status()
            return r.content, verify
        except Exception as ex:  # noqa: BLE001
            last = ex
    raise last


def pdf_text(pdf_path):
    r = subprocess.run(["pdftotext", "-layout", pdf_path, "-"], capture_output=True, timeout=600)
    return r.stdout.decode("utf-8", "replace")


def part_wiki(out, summary, shafts):
    import urllib.parse
    res = dict(pages={}, errors=[])
    coords = []
    for p in WIKI_PAGES:
        url = "https://en.wikipedia.org/w/index.php?title=%s&action=raw" % urllib.parse.quote(p)
        try:
            r = http_get(url, timeout=60)
            r.raise_for_status()
            c = parse_coord_templates(r.text)
            coords.extend(c)
            res["pages"][p] = dict(bytes=len(r.text), n_coords=len(c))
            with open(os.path.join(out, "wiki_%s.txt" % p.split("_")[0].replace("–", "-")), "w",
                      encoding="utf-8") as f:
                f.write(r.text[:MAX_TXT])
        except Exception as ex:  # noqa: BLE001
            res["errors"].append("%s: %r" % (p, ex))
        time.sleep(1.0)
    ok = verify_rows(shafts, coords)
    for r, v in zip(shafts, ok):
        r["wiki_verified"] = int(v)
    res["n_rows"], res["n_verified"] = len(shafts), int(sum(ok))
    res["unverified"] = [r["name"] for r, v in zip(shafts, ok) if not v]
    summary["wiki"] = res


def part_yang(out, summary, shafts):
    res = dict(errors=[])
    try:
        data, verified = fetch_bytes(YANG_PDF)
        p = os.path.join(out, "_yang.pdf")
        with open(p, "wb") as f:
            f.write(data)
        txt = pdf_text(p)
        os.remove(p)
        res.update(bytes=len(data), tls_verified=verified, text_chars=len(txt))
        lines = [l for l in txt.splitlines() if find_latlon(l)]
        with open(os.path.join(out, "pkp29_yang_fp_lines.txt"), "w", encoding="utf-8") as f:
            f.write("\n".join(lines))
        res["n_fp_lines"] = len(lines)
        by_date = {}
        for l in lines:
            d = find_date(l)
            ll = find_latlon(l)
            if d and ll:
                by_date.setdefault(d, []).append(ll)
        res["n_fp_dated"] = len(by_date)
        n = 0
        for r in shafts:
            v = by_date.get(r["date"])
            if v:
                r["yang_lat"], r["yang_lon"] = v[0]
                n += 1
        res["n_matched_to_shafts"] = n
        if not lines:
            with open(os.path.join(out, "pkp29_yang_text_head.txt"), "w", encoding="utf-8") as f:
                f.write(txt[:200000])
    except Exception as ex:  # noqa: BLE001
        res["errors"].append(repr(ex)[:300])
    summary["yang"] = res


def part_catalogues(out, summary, shafts):
    res = dict(errors=[])
    ts = [origin_seconds(r["date"], r["time_ut"]) for r in shafts]
    # USGS
    try:
        q = ("https://earthquake.usgs.gov/fdsnws/event/1/query?format=csv&starttime=1977-01-01"
             "&endtime=1996-02-01&minlatitude=-22.6&maxlatitude=-21.4&minlongitude=-139.4"
             "&maxlongitude=-138.3&minmagnitude=3&orderby=time-asc&limit=2000")
        r = http_get(q)
        r.raise_for_status()
        import io
        ev = list(csv.DictReader(io.StringIO(r.text)))
        res["n_usgs"] = len(ev)
        for e in ev:
            t = _iso_s(e["time"])
            j = _nearest(ts, t, 90.0, shafts, e["time"][:10])
            if j is not None:
                shafts[j].update(usgs_lat=float(e["latitude"]), usgs_lon=float(e["longitude"]),
                                 usgs_mag=e.get("mag"), usgs_dt_s=t - (ts[j] or t))
    except Exception as ex:  # noqa: BLE001
        res["errors"].append("usgs: " + repr(ex)[:200])
    # ISC
    try:
        from obspy.clients.fdsn import Client
        from obspy import UTCDateTime
        c = Client("ISC", timeout=240)
        cat = c.get_events(starttime=UTCDateTime("1977-01-01"), endtime=UTCDateTime("1996-02-01"),
                           minlatitude=-22.6, maxlatitude=-21.4, minlongitude=-139.4,
                           maxlongitude=-138.3, minmagnitude=3.0)
        res["n_isc"] = len(cat)
        for e in cat:
            o = e.preferred_origin() or (e.origins[0] if e.origins else None)
            if o is None:
                continue
            t = float(o.time.timestamp)
            j = _nearest(ts, t, 90.0, shafts, str(o.time)[:10])
            if j is not None:
                shafts[j].update(isc_lat=float(o.latitude), isc_lon=float(o.longitude),
                                 isc_dt_s=t - (ts[j] or t))
    except Exception as ex:  # noqa: BLE001
        res["errors"].append("isc: " + repr(ex)[:200])
    for k in ("usgs", "isc"):
        res["n_matched_" + k] = sum(1 for r in shafts if (k + "_lat") in r)
    summary["catalogues"] = res


def _iso_s(s):
    import calendar
    s = s.replace("Z", "")
    main, _, frac = s.partition(".")
    tm = time.strptime(main, "%Y-%m-%dT%H:%M:%S")
    return calendar.timegm(tm) + (float("0." + frac) if frac else 0.0)


def _nearest(ts, t, tol, shafts, date):
    best, bd = None, tol
    for j, x in enumerate(ts):
        if x is None:
            if shafts[j]["date"] == date and best is None:
                best = j
            continue
        if abs(x - t) <= bd:
            best, bd = j, abs(x - t)
    return best


def part_pdf_text(out, summary, key, url, fname, ocr=False, budget_s=3600):
    res = dict(url=url, errors=[])
    try:
        data, verified = fetch_bytes(url)
        p = os.path.join(out, "_%s.pdf" % key)
        with open(p, "wb") as f:
            f.write(data)
        res.update(bytes=len(data), tls_verified=verified)
        txt = pdf_text(p)
        res["text_chars"] = len(txt.strip())
        if ocr and len(txt.strip()) < 2000:
            txt = ocr_pdf(p, os.path.join(out, "_ocr_" + key), budget_s, res)
        os.remove(p)
        with open(os.path.join(out, fname), "w", encoding="utf-8") as f:
            f.write(txt[:MAX_TXT])
        rows = [l for l in txt.splitlines() if find_date(l)]
        res["n_dated_lines"] = len(rows)
        res["n_fp_latlon_lines"] = sum(1 for l in txt.splitlines() if find_latlon(l))
    except Exception as ex:  # noqa: BLE001
        res["errors"].append(repr(ex)[:300])
    summary[key] = res


def ocr_pdf(pdf, tmp, budget_s, res):
    os.makedirs(tmp, exist_ok=True)
    t0 = time.time()
    info = subprocess.run(["pdfinfo", pdf], capture_output=True).stdout.decode("utf-8", "replace")
    m = re.search(r"Pages:\s+(\d+)", info)
    npages = int(m.group(1)) if m else 0
    res["pages"] = npages
    texts = []
    for i in range(1, npages + 1):
        if time.time() - t0 > budget_s:
            res["errors"].append("ocr stopped at page %d by time budget" % i)
            break
        base = os.path.join(tmp, "p")
        subprocess.run(["pdftoppm", "-r", "200", "-gray", "-f", str(i), "-l", str(i), "-png", pdf, base],
                       capture_output=True, timeout=300)
        pngs = sorted(x for x in os.listdir(tmp) if x.endswith(".png"))
        for png in pngs:
            pp = os.path.join(tmp, png)
            r = subprocess.run(["tesseract", pp, "-", "--psm", "6"], capture_output=True, timeout=300)
            texts.append("=== page %d ===\n%s" % (i, r.stdout.decode("utf-8", "replace")))
            os.remove(pp)
    res["ocr_pages_done"] = len(texts)
    return "\n".join(texts)


def part_predict(out, summary, shafts):
    from obspy.taup import TauPyModel
    from obspy.geodetics import locations2degrees
    model = TauPyModel("iasp91")
    cache = {}

    def tt(lat, lon):
        d = float(locations2degrees(lat, lon, BRVK[0], BRVK[1]))
        k = round(d, 5)
        if k not in cache:
            first = {}
            for a in model.get_travel_times(source_depth_in_km=SRC_DEPTH_KM, distance_in_degree=d,
                                            phase_list=PHASES):
                first.setdefault(a.name, (a.time, a.ray_param_sec_degree))
            cache[k] = (d, first)
        return cache[k]

    d0, f0 = tt(*REF)
    _, az = gc_dist_az(REF[0], REF[1], BRVK[0], BRVK[1])
    summary["reference"] = dict(lat=REF[0], lon=REF[1], dist_deg=d0, az_to_brvk=az,
                                ray_params_s_per_deg={k: v[1] for k, v in f0.items()},
                                times_s={k: v[0] for k, v in f0.items()})

    def diffs(lat, lon):
        d, f = tt(lat, lon)
        o = dict(dist_deg=d)
        for a, b in (("PKIKP", "PP"), ("PKiKP", "PP"), ("PPP", "PP")):
            if a in f and b in f and a in f0 and b in f0:
                o["%s_%s_ms" % (a, b)] = 1000.0 * ((f[a][0] - f[b][0]) - (f0[a][0] - f0[b][0]))
        return o

    sets = ["shaft", "yang", "usgs", "isc"]
    for r in shafts:
        for s in sets:
            la = r["lat"] if s == "shaft" else r.get(s + "_lat")
            lo = r["lon"] if s == "shaft" else r.get(s + "_lon")
            if la is None or lo is None:
                continue
            al, ac = along_across_km(REF[0], REF[1], la, lo, az)
            r[s + "_along_km"], r[s + "_across_km"] = al, ac
            for k, v in diffs(la, lo).items():
                r["%s_%s" % (s, k)] = v
        for s in sets[1:]:
            if (s + "_PKiKP_PP_ms") in r:
                r[s + "_minus_shaft_PKiKP_PP_ms"] = r[s + "_PKiKP_PP_ms"] - r["shaft_PKiKP_PP_ms"]
                r[s + "_minus_shaft_along_km"] = r[s + "_along_km"] - r["shaft_along_km"]
    # sensitivity from the shaft set
    xs = [r["shaft_along_km"] for r in shafts if r["atoll"] == "Mururoa"]
    ep = {}
    for key in ("PKIKP_PP_ms", "PKiKP_PP_ms", "PPP_PP_ms"):
        ys = [r.get("shaft_" + key) for r in shafts if r["atoll"] == "Mururoa"]
        pts = [(x, y) for x, y in zip(xs, ys) if y is not None]
        if len(pts) > 2:
            a = np.polyfit([p[0] for p in pts], [p[1] for p in pts], 1)
            ep["slope_" + key + "_per_km_along"] = float(a[0])
    summary["sensitivity"] = ep
    # epoch means
    em = {}
    for key in ["shaft_along_km", "shaft_PKiKP_PP_ms", "shaft_PKIKP_PP_ms", "shaft_PPP_PP_ms",
                "yang_minus_shaft_PKiKP_PP_ms", "isc_minus_shaft_PKiKP_PP_ms",
                "usgs_minus_shaft_PKiKP_PP_ms"]:
        em[key] = dict(all=epoch_means(shafts, key), ge20kt=epoch_means(shafts, key, min_yield=20.0))
    summary["epochs"] = em
    keys = []
    for r in shafts:
        for k in r:
            if k not in keys:
                keys.append(k)
    with open(os.path.join(out, "pkp29_tests.csv"), "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        for r in shafts:
            w.writerow({k: (round(v, 5) if isinstance(v, float) else v) for k, v in r.items()})


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
    ap.add_argument("--shafts", default="mururoa_shafts_1977_1996.csv")
    ap.add_argument("--out-dir", default="out")
    ap.add_argument("--ocr-budget-min", type=float, default=60.0)
    a = ap.parse_args(argv)
    os.makedirs(a.out_dir, exist_ok=True)
    t0 = time.time()
    summary = dict(version=VERSION, started=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
                   note="scan for pre-registration 28; grades nothing", errors=[])
    shafts = load_shafts(a.shafts)
    summary["n_tests"] = len(shafts)

    def save():
        with open(os.path.join(a.out_dir, "pkp29_summary.json"), "w") as f:
            json.dump(_clean(summary), f, indent=1)

    steps = [
        ("wiki", lambda: part_wiki(a.out_dir, summary, shafts)),
        ("yang", lambda: part_yang(a.out_dir, summary, shafts)),
        ("catalogues", lambda: part_catalogues(a.out_dir, summary, shafts)),
        ("predict", lambda: part_predict(a.out_dir, summary, shafts)),
        ("civen", lambda: part_pdf_text(a.out_dir, summary, "civen", CIVEN_PDF, "pkp29_civen.txt")),
        ("aven", lambda: part_pdf_text(a.out_dir, summary, "aven", AVEN_PDF, "pkp29_aven.txt")),
        ("blacknest", lambda: part_pdf_text(a.out_dir, summary, "blacknest", BLACKNEST_PDF,
                                            "pkp29_blacknest_ocr.txt", ocr=True,
                                            budget_s=a.ocr_budget_min * 60.0)),
    ]
    for name, fn in steps:
        try:
            fn()
            log(name, "ok")
        except Exception as ex:  # noqa: BLE001
            summary["errors"].append("%s: %r" % (name, ex))
            log(name, "FAILED", repr(ex)[:200])
        save()
    summary["finished"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    summary["runtime_min"] = round((time.time() - t0) / 60.0, 2)
    save()
    log("done", summary["runtime_min"], "min; errors:", summary["errors"])
    return 0


if __name__ == "__main__":
    sys.exit(main())

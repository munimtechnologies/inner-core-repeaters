#!/usr/bin/env python3
"""geom.py: where do this project's PKIKP paths cross Earth's inner-core boundary (ICB)?

Part of inner-core-repeaters by MV2, Munim, Inc. (made with AI agents). Licence: Apache-2.0.
Companion to pre-registration 16 (near-antipodal clock test, nea.py). The reading rules that use
this geometry are fixed in findings/2026-11-19-prereg16-addendum-geometry.md, published together
with this file and before any geometry or any P16 result exists.

Why it matters. A rigid turn-back of the whole inner core (Wang, Vidale et al. 2024) predicts the
same "equivalent date" clock on every path. A local change of the inner-core surface (Wen 2006;
Zhang & Wen 2026) predicts the clock only where a path crosses the changing patch. So a
backtracking clock on the near-antipodal paths speaks against a single local patch only if those
paths cross the ICB away from where the ILAR and YKA paths cross it. This script measures that.

For the South Sandwich Islands (SSI) multiplet centroid and every candidate station of nea.py
(BHZ, networks IU II IC G GE GT CU, 155-179.5 deg) plus the ILAR and YKA arrays, it computes with
ObsPy TauP (ak135): the source-side ICB entry point, the receiver-side exit point, the bottoming
point, the inner-core path length, the share of it within 300 km of the ICB, and the distance of
each path's crossing points from the nearest ILAR/YKA crossing point (km, measured on the ICB).

Modes:
  python geom.py --out-dir /out                      geometry (needs network + obspy)
  python geom.py --grade --geom-csv G --nea-dir N --out-dir O   pre-registered reading of P16
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
NAN = float("nan")

REPO = "https://raw.githubusercontent.com/munimtechnologies/inner-core-repeaters/"
INPUT_COMMIT = "d0186eece9f9bae982a55a3302a0773cf8627d98"   # v0.14.6, same inputs as nea.py
CENTROID = (-59.0365, -25.6939)                            # results/m0324-v1 centroid
ANTIPODE = (-CENTROID[0], CENTROID[1] + 180.0 if CENTROID[1] < 0 else CENTROID[1] - 180.0)
MIN_DIST, MAX_DIST = 155.0, 179.5
NETWORKS = "IU,II,IC,G,GE,GT,CU"
MAX_STATIONS = 10
MIN_COVER = 20
LOC_PREF = ("00", "", "10", "20", "60")
STATION_BASES = ("https://service.earthscope.org", "https://service.iris.edu")
# Reference arrays (array centres; a few km of error here is far below the thresholds at the ICB).
REFERENCES = {"IM.ILAR": (64.7714, -146.8866), "CN.YKA": (62.4931, -114.6053)}

R_EARTH = 6371.0
DEFAULT_DEPTH_KM = 20.0
# Pre-registered constants (addendum to pre-registration 16)
THR_KM = 250.0            # primary: crossing points >= 250 km apart on the ICB count as separate
THR_ALT = (150.0, 400.0)  # sensitivity, reported only
F_SEP_MIN = 0.75          # weighted share of "separate" stations needed for reading R1 / R3-separate
UPPER_IC_KM = 300.0       # "upper inner core" = within 300 km below the ICB


def log(*a):
    print(time.strftime("[%H:%M:%S]"), *a, flush=True)


# ----------------------------------------------------------------------------- pure helpers
def isf(x):
    try:
        return math.isfinite(float(x))
    except (TypeError, ValueError):
        return False


def gc_deg(lat1, lon1, lat2, lon2):
    """Great-circle distance in degrees on a sphere (haversine)."""
    p1, p2 = math.radians(lat1), math.radians(lat2)
    dp, dl = p2 - p1, math.radians(lon2 - lon1)
    a = math.sin(dp / 2) ** 2 + math.cos(p1) * math.cos(p2) * math.sin(dl / 2) ** 2
    return math.degrees(2.0 * math.asin(min(1.0, math.sqrt(a))))


def azimuth_deg(lat1, lon1, lat2, lon2):
    """Initial azimuth (deg clockwise from north) of the great circle from point 1 to point 2."""
    p1, p2 = math.radians(lat1), math.radians(lat2)
    dl = math.radians(lon2 - lon1)
    y = math.sin(dl) * math.cos(p2)
    x = math.cos(p1) * math.sin(p2) - math.sin(p1) * math.cos(p2) * math.cos(dl)
    return math.degrees(math.atan2(y, x)) % 360.0


def destination(lat, lon, az, ang):
    """Point reached from (lat, lon) after `ang` degrees of arc along initial azimuth `az`."""
    p1, l1 = math.radians(lat), math.radians(lon)
    t, d = math.radians(az), math.radians(ang)
    sp2 = math.sin(p1) * math.cos(d) + math.cos(p1) * math.sin(d) * math.cos(t)
    p2 = math.asin(max(-1.0, min(1.0, sp2)))
    l2 = l1 + math.atan2(math.sin(t) * math.sin(d) * math.cos(p1), math.cos(d) - math.sin(p1) * sp2)
    lon2 = (math.degrees(l2) + 180.0) % 360.0 - 180.0
    return math.degrees(p2), lon2


def ic_metrics(dists_rad, depths_km, icb_depth, src, az, r_earth=R_EARTH, upper_km=UPPER_IC_KM, tol=0.5):
    """Inner-core part of one ray path given as plain lists (epicentral angle in radians, depth in
    km). Returns None if the path never reaches the ICB."""
    idx = [i for i, z in enumerate(depths_km) if z >= icb_depth - tol]
    if not idx:
        return None
    i0, i1 = idx[0], idx[-1]
    ib = max(range(len(depths_km)), key=lambda i: depths_km[i])
    r_icb = r_earth - icb_depth
    e_ang, x_ang, b_ang = (math.degrees(dists_rad[i]) for i in (i0, i1, ib))
    e_ll, x_ll, b_ll = (destination(src[0], src[1], az, a) for a in (e_ang, x_ang, b_ang))
    length = upper = lever_w = 0.0
    for i in range(i0, i1):
        ra, rb = r_earth - depths_km[i], r_earth - depths_km[i + 1]
        da, db = dists_rad[i], dists_rad[i + 1]
        seg = math.hypot(rb * math.cos(db) - ra * math.cos(da), rb * math.sin(db) - ra * math.sin(da))
        zm, dm = 0.5 * (depths_km[i] + depths_km[i + 1]), 0.5 * (da + db)
        lat_m, _ = destination(src[0], src[1], az, math.degrees(dm))
        length += seg
        if zm - icb_depth <= upper_km:
            upper += seg
        lever_w += seg * (r_earth - zm) * math.cos(math.radians(lat_m))
    return dict(entry_ang=e_ang, exit_ang=x_ang, bottom_ang=b_ang,
                entry_lat=e_ll[0], entry_lon=e_ll[1], exit_lat=x_ll[0], exit_lon=x_ll[1],
                bottom_lat=b_ll[0], bottom_lon=b_ll[1], bottom_below_icb_km=max(depths_km) - icb_depth,
                ic_path_km=length, upper_ic_frac=(upper / length) if length > 0 else NAN,
                mean_axis_dist_km=(lever_w / length) if length > 0 else NAN,
                entry_km_per_deg_rot=r_icb * math.cos(math.radians(e_ll[0])) * math.pi / 180.0,
                exit_km_per_deg_rot=r_icb * math.cos(math.radians(x_ll[0])) * math.pi / 180.0)


def icb_km(p, q, r_icb):
    """Distance on the ICB (km) between two crossing points given as (lat, lon)."""
    return math.radians(gc_deg(p[0], p[1], q[0], q[1])) * r_icb


def nearest_ref(pt, refs, r_icb):
    """refs: [(label, (lat, lon))]. Returns (km, label) of the closest reference crossing point."""
    best = (NAN, "")
    for lab, q in refs:
        d = icb_km(pt, q, r_icb)
        if not isf(best[0]) or d < best[0]:
            best = (d, lab)
    return best


def classify(entry_sep_km, exit_sep_km, thr=THR_KM):
    """Path class from the distances of its two ICB crossings to the nearest reference crossing."""
    if not (isf(entry_sep_km) and isf(exit_sep_km)):
        return "unknown"
    e, x = entry_sep_km < thr, exit_sep_km < thr
    if not e and not x:
        return "separate"
    if e and x:
        return "shared-both"
    return "shared-entry" if e else "shared-exit"


def f_separate(classes, weights):
    """Weighted share of 'separate' paths. classes, weights: {station: value}. NaN if no weight."""
    tot = sum(float(w) for s, w in weights.items() if s in classes and float(w) > 0)
    if tot <= 0:
        return NAN
    sep = sum(float(w) for s, w in weights.items() if classes.get(s) == "separate" and float(w) > 0)
    return sep / tot


def reading(headline, fsep, f_min=F_SEP_MIN):
    """Pre-registered reading of the P16 headline given the weighted share of separate paths."""
    far = isf(fsep) and fsep >= f_min
    if headline == "B-supported":
        return "R1" if far else "R2"
    if headline == "M-supported":
        return "R3-separate" if far else "R3-shared"
    if headline == "no-time-signal":
        return "R4"
    return "none"


READINGS = {
    "R1": "Backtracking clock on ICB crossings separate from the ILAR/YKA crossings: favours a "
          "motion of the whole inner core over one local surface patch.",
    "R2": "Backtracking clock, but the paths share ICB crossings with ILAR/YKA: compatible with "
          "both a whole-inner-core motion and one local surface patch.",
    "R3-separate": "Time signal on separate crossings follows plain time order, not the "
                   "backtracking clock: the 2024 clock does not describe these parts of the inner core.",
    "R3-shared": "Plain time order on paths that share crossings with ILAR/YKA: conflicts with the "
                 "backtracking clock seen on those same crossings.",
    "R4": "No time signal on near-antipodal paths: read together with their upper-inner-core share "
          "(upper_ic_frac) against ILAR's.",
    "none": "No reading beyond the P16 headline.",
}


def parse_members(s):
    out = []
    for part in str(s).split(";"):
        part = part.strip()
        if part:
            out.append(part.split("(")[0])
    return out


def choose_loc(chans, iso, pref=LOC_PREF):
    """Same rule as nea.py: preferred location code whose BHZ epoch covers `iso`, or None."""
    s = str(iso)[:19]
    cover = {loc for loc, a, b in chans if str(a)[:19] <= s < (str(b)[:19] if b else "9999")}
    for loc in pref:
        if loc in cover:
            return loc
    return sorted(cover)[0] if cover else None


def parse_station_text(body):
    """FDSN station text (level=channel) -> {(net, sta): dict(lat, lon, chans=[(loc, start, end)])}."""
    out = {}
    for line in str(body).splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        f = line.split("|")
        if len(f) < 17:
            continue
        net, sta, loc = f[0].strip(), f[1].strip(), f[2].strip()
        try:
            lat, lon = float(f[4]), float(f[5])
        except ValueError:
            continue
        start = f[15].strip()[:19]
        end = f[16].strip()[:19] or "9999"
        d = out.setdefault((net, sta), dict(lat=lat, lon=lon, chans=[]))
        d["lat"], d["lon"] = lat, lon
        d["chans"].append((loc, start, end))
    return out


def select_candidates(stations, event_times, centroid=CENTROID):
    """nea.py's station rule: distance window, cover = events inside a BHZ epoch, rank by cover,
    keep at most MAX_STATIONS with >= MIN_COVER (one per station code)."""
    found = []
    for (net, sta), d in stations.items():
        dist = gc_deg(centroid[0], centroid[1], d["lat"], d["lon"])
        if not (MIN_DIST <= dist <= MAX_DIST):
            continue
        cover = sum(1 for t in event_times if choose_loc(d["chans"], str(t).rstrip("Z")) is not None)
        found.append(dict(net=net, sta=sta, lat=d["lat"], lon=d["lon"], dist=dist, cover=cover))
    found.sort(key=lambda s: (-s["cover"], s["net"], s["sta"]))
    seen, chosen = set(), []
    for s in found:
        if s["cover"] < MIN_COVER or s["sta"] in seen:
            continue
        seen.add(s["sta"])
        chosen.append("%s.%s" % (s["net"], s["sta"]))
        if len(chosen) >= MAX_STATIONS:
            break
    return found, chosen


def median(xs):
    v = sorted(float(x) for x in xs if isf(x))
    n = len(v)
    if not n:
        return NAN
    return v[n // 2] if n % 2 else 0.5 * (v[n // 2 - 1] + v[n // 2])


# ----------------------------------------------------------------------------- IO
def fetch_text(url, tries=4, wait=10, timeout=90):
    for k in range(tries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "inner-core-repeaters (MV2, Munim, Inc.)"})
            with urllib.request.urlopen(req, timeout=timeout) as r:
                if r.status == 200:
                    return r.read().decode("utf-8")
        except Exception as ex:  # noqa: BLE001
            log("fetch error", url, repr(ex)[:160])
        time.sleep(wait * (k + 1))
    return None


def read_csv_text(body):
    return list(csv.DictReader(body.splitlines()))


def load_events():
    """Multiplet event times and depths from results/m0324-v1 at the pinned commit."""
    mb = fetch_text(REPO + INPUT_COMMIT + "/results/m0324-v1/m0324_multiplets.csv")
    cb = fetch_text(REPO + INPUT_COMMIT + "/results/m0324-v1/m0324_catalog.csv")
    if not mb or not cb:
        return [], []
    cat = {r["id"]: r for r in read_csv_text(cb) if r.get("id")}
    ids = set()
    for r in read_csv_text(mb):
        if r.get("status") == "ok":
            m = [i for i in parse_members(r.get("members", "")) if i in cat]
            if len(m) >= 2:
                ids.update(m)
    times = [cat[i]["time"] for i in sorted(ids)]
    depths = [cat[i].get("depth") for i in sorted(ids)]
    return times, depths


def write_csv(path, rows):
    keys = []
    for r in rows:
        for k in r:
            if k not in keys:
                keys.append(k)
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        for r in rows:
            w.writerow({k: ("" if (isinstance(v, float) and not math.isfinite(v)) else v) for k, v in r.items()})


def jsonable(o):
    if isinstance(o, dict):
        return {str(k): jsonable(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [jsonable(v) for v in o]
    if isinstance(o, float):
        return o if math.isfinite(o) else None
    return o


# ----------------------------------------------------------------------------- geometry run
def run_geometry(out_dir):
    from obspy.taup import TauPyModel

    os.makedirs(os.path.join(out_dir, "plots"), exist_ok=True)
    s = dict(version=VERSION, started=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()), model="ak135",
             input_commit=INPUT_COMMIT, centroid=CENTROID,
             constants=dict(thr_km=THR_KM, thr_alt=THR_ALT, f_sep_min=F_SEP_MIN, upper_ic_km=UPPER_IC_KM))
    times, depths = load_events()
    s["n_events"] = len(times)
    depth = median(depths)
    s["source_depth_km"] = depth if isf(depth) else DEFAULT_DEPTH_KM
    s["source_depth_from"] = "median catalogue depth" if isf(depth) else "default"
    model = TauPyModel("ak135")
    tm = model.model
    icb_depth = float(getattr(tm, "iocb_depth", 5153.5))
    r_icb = R_EARTH - icb_depth
    s["icb_depth_km"], s["icb_radius_km"] = icb_depth, r_icb

    body, used = None, ""
    q = ("/fdsnws/station/1/query?network=%s&channel=BHZ&level=channel&format=text"
         "&latitude=%.4f&longitude=%.4f&maxradius=%.1f" % (NETWORKS, ANTIPODE[0], ANTIPODE[1], 180.0 - MIN_DIST))
    for base in STATION_BASES:
        body = fetch_text(base + q)
        if body:
            used = base
            break
    s["station_service"] = used
    stations = parse_station_text(body or "")
    found, chosen = select_candidates(stations, times)
    s["n_candidates"], s["predicted_chosen"] = len(found), chosen

    def geom(lat, lon):
        dist = gc_deg(CENTROID[0], CENTROID[1], lat, lon)
        az = azimuth_deg(CENTROID[0], CENTROID[1], lat, lon)
        arr = model.get_ray_paths(source_depth_in_km=s["source_depth_km"], distance_in_degree=dist,
                                  phase_list=["PKIKP"])
        if not arr:
            return dict(dist=dist, az=az, status="no-PKIKP")
        a = arr[0]
        m = ic_metrics([float(v) for v in a.path["dist"]], [float(v) for v in a.path["depth"]], icb_depth,
                       CENTROID, az)
        if m is None:
            return dict(dist=dist, az=az, status="no-ICB")
        return dict(dist=dist, az=az, status="ok", p_s_per_deg=float(a.ray_param_sec_degree), **m)

    rows, refs = [], []
    for name, (lat, lon) in REFERENCES.items():
        g = geom(lat, lon)
        rows.append(dict(station=name, role="reference", lat=lat, lon=lon, cover="", **g))
        if g.get("status") == "ok":
            refs += [(name + ":entry", (g["entry_lat"], g["entry_lon"])), (name + ":exit", (g["exit_lat"], g["exit_lon"]))]
    for f in found:
        name = "%s.%s" % (f["net"], f["sta"])
        g = geom(f["lat"], f["lon"])
        row = dict(station=name, role="predicted-chosen" if name in chosen else "candidate", lat=f["lat"],
                   lon=f["lon"], cover=f["cover"], **g)
        if g.get("status") == "ok" and refs:
            ek, el = nearest_ref((g["entry_lat"], g["entry_lon"]), refs, r_icb)
            xk, xl = nearest_ref((g["exit_lat"], g["exit_lon"]), refs, r_icb)
            row.update(entry_sep_km=ek, entry_nearest=el, exit_sep_km=xk, exit_nearest=xl,
                       path_class=classify(ek, xk))
            for t in THR_ALT:
                row["path_class_%d" % int(t)] = classify(ek, xk, t)
        rows.append(row)
    write_csv(os.path.join(out_dir, "geom_paths.csv"), rows)

    pc = [r for r in rows if r["role"] == "predicted-chosen"]
    s["predicted_chosen_classes"] = {r["station"]: r.get("path_class", "unknown") for r in pc}
    s["f_sep_predicted_unweighted"] = f_separate(s["predicted_chosen_classes"], {r["station"]: 1.0 for r in pc})
    s["references"] = {r["station"]: {k: r.get(k) for k in ("dist", "az", "entry_lat", "entry_lon", "exit_lat",
                                                              "exit_lon", "bottom_below_icb_km", "ic_path_km",
                                                              "upper_ic_frac", "mean_axis_dist_km")}
                       for r in rows if r["role"] == "reference"}
    for key in ("upper_ic_frac", "bottom_below_icb_km", "mean_axis_dist_km"):
        s["predicted_chosen_median_" + key] = median([r.get(key) for r in pc])
    try:
        plot(os.path.join(out_dir, "plots", "geom_icb_crossings.png"), rows)
    except Exception as ex:  # noqa: BLE001
        log("plot failed", repr(ex)[:160])
    s["finished"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    with open(os.path.join(out_dir, "geom_summary.json"), "w") as f:
        json.dump(jsonable(s), f, indent=1)
    log(json.dumps(jsonable({k: s[k] for k in ("n_candidates", "predicted_chosen_classes",
                                                "f_sep_predicted_unweighted")}))[:2000])
    return 0


def plot(path, rows):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(10, 5.5))
    for r in rows:
        if r.get("status") != "ok":
            continue
        ref = r["role"] == "reference"
        col = "C3" if ref else ("C0" if r["role"] == "predicted-chosen" else "0.6")
        ax.plot([r["entry_lon"]], [r["entry_lat"]], "o", color=col, ms=6 if ref else 4)
        ax.plot([r["exit_lon"]], [r["exit_lat"]], "s", color=col, ms=6 if ref else 4)
        ax.plot([r["lon"]], [r["lat"]], "^", color=col, ms=6 if ref else 3, alpha=0.6)
        if ref or r["role"] == "predicted-chosen":
            ax.text(r["exit_lon"] + 1, r["exit_lat"] + 1, r["station"].split(".")[-1], fontsize=6, color=col)
    ax.plot([CENTROID[1]], [CENTROID[0]], "*", color="k", ms=10)
    ax.set_xlim(-180, 180)
    ax.set_ylim(-90, 90)
    ax.set_xlabel("longitude")
    ax.set_ylabel("latitude")
    ax.grid(alpha=0.3)
    ax.set_title("ak135 PKIKP inner-core-boundary crossings from the SSI centroid (star): circle = entry, "
                 "square = exit, triangle = station\nred = ILAR/YKA, blue = predicted P16 stations, grey = "
                 "other candidates", fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=110)
    plt.close(fig)


# ----------------------------------------------------------------------------- grading mode
def run_grade(geom_csv, nea_dir, out_dir):
    """Reading of the P16 headline with the geometry, weights = ok pairs per station in nea_pairs.csv."""
    with open(geom_csv, newline="") as f:
        g = list(csv.DictReader(f))
    with open(os.path.join(nea_dir, "nea_summary.json")) as f:
        nea = json.load(f)
    weights = {}
    with open(os.path.join(nea_dir, "nea_pairs.csv"), newline="") as f:
        for r in csv.DictReader(f):
            if str(r.get("ok")).strip() == "True":
                weights[r["station"]] = weights.get(r["station"], 0.0) + 1.0
    classes = {r["station"]: r.get("path_class") or "unknown" for r in g if r.get("role") != "reference"}
    head = (nea.get("grade") or {}).get("headline", "untestable")
    fs = f_separate(classes, weights)
    rd = reading(head, fs)
    res = dict(version=VERSION, headline=head, weights=weights, classes={s: classes.get(s, "unknown") for s in weights},
               f_sep_weighted=fs, reading=rd, reading_text=READINGS[rd],
               stations_chosen=nea.get("stations_chosen"))
    os.makedirs(out_dir, exist_ok=True)
    with open(os.path.join(out_dir, "geom_reading.json"), "w") as f:
        json.dump(jsonable(res), f, indent=1)
    log(json.dumps(jsonable(res))[:2000])
    return 0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-dir", default="/out")
    ap.add_argument("--grade", action="store_true")
    ap.add_argument("--geom-csv", default="")
    ap.add_argument("--nea-dir", default="")
    a = ap.parse_args()
    if a.grade:
        return run_grade(a.geom_csv, a.nea_dir, a.out_dir)
    return run_geometry(a.out_dir)


if __name__ == "__main__":
    sys.exit(main())

#!/usr/bin/env python3
"""m27.py: blind YKA test of the 2026 South Sandwich Islands events in ILAR multiplets Z006,
Z007 and Z009 (pre-registration 15).

Part of inner-core-repeaters by MV2, Munim, Inc. (made with AI agents). Licence: Apache-2.0.
Pre-registered in findings/2026-11-19-prereg15-z006-2026-event.md, written and published
before any YKA or mantle-P waveform of these families was fetched.

The ILAR pipeline (icr-m0324-v1, results/m0324-v1, run 2026-11-17) put the 2026-11-13 event
us7000ty3q in multiplet Z006 with events from 2000-10-08, 2007-11-21, 2012-08-23, 2019-05-01,
2021-09-04 and 2023-10-15, not in M03-24 (Z004). On the backtracking clock of Wang, Vidale
et al. (2024; turn 2008.5, back-motion 2.5x slower) 2026.87 maps to an effective 2001.16, so
the best-matching partner should be the 2000-10-08 event (effective distance 0.39 yr). If the
inner core has simply kept changing, or stayed quiet since 2023, the 2023-10-15 event should
match best. This compares the two at YKA, an array whose path differs from ILAR's except near
the source and in the inner core, with a mantle-P source control at the stations of
pre-registration 13.

Parts:
  Y  m0324_check.run_yka (v0.8.0, pinned) with yka_repick's +/-2 s picker, exactly as in
     icr-m26-v1, on each family (Z006 primary; Z007 and Z009 reported only).
  M  mantlep2.py v0.13.0 (pinned commit) with its event list replaced by the Z006 events, at
     the eight stations of pre-registration 13.
"""
import argparse
import calendar
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

# USGS ComCat (id, origin UTC, lat, lon, depth km, mag) as in results/m0324-v1/m0324_catalog.csv.
Z006 = [
    ("usp000a1fh", "2000-10-08T14:14:45.800Z", -58.902, -25.607, 33.0, 4.3),
    ("usp000ft0f", "2007-11-21T03:24:34.320Z", -58.983, -25.663, 62.0, 4.5),
    ("usp000jqn5", "2012-08-23T15:54:50.880Z", -58.895, -25.452, 35.0, 5.2),
    ("us70003jaz", "2019-05-01T21:08:56.344Z", -59.0194, -25.6465, 65.71, 4.7),
    ("us7000f7ug", "2021-09-04T07:44:29.044Z", -58.9817, -25.5482, 44.8, 4.9),
    ("us6000lfpq", "2023-10-15T08:58:34.423Z", -58.9574, -25.4713, 35.0, 4.8),
    ("us7000ty3q", "2026-11-13T19:04:07.703Z", -58.962, -25.7253, 63.81, 5.1),
]
Z007 = [
    ("usp000dzku", "2005-09-17T11:56:22.430Z", -59.152, -25.828, 39.7, 4.7),
    ("usp000hsps", "2011-01-07T09:45:17.020Z", -59.184, -25.948, 30.4, 4.7),
    ("us100044ec", "2015-12-06T08:53:28.450Z", -59.1301, -25.7922, 50.04, 4.8),
    ("us6000ddkg", "2021-02-01T05:19:59.928Z", -59.1499, -25.9477, 71.8, 5.1),
    ("us7000fpnn", "2021-10-29T06:38:20.317Z", -59.1418, -25.7541, 48.06, 5.0),
    ("us7000iss3", "2022-11-27T18:07:47.127Z", -59.143, -25.9197, 68.248, 4.9),
    ("us6000ty0j", "2026-09-27T21:13:35.077Z", -59.1168, -25.9907, 117.745, 5.0),
]
Z009 = [
    ("usp000h6wk", "2010-02-03T02:18:27.950Z", -58.941, -25.738, 49.8, 5.2),
    ("us7000gf0k", "2022-01-17T06:42:13.520Z", -59.1926, -25.5922, 35.0, 4.4),
    ("us7000txqr", "2026-11-12T03:45:56.644Z", -59.0326, -26.0288, 121.557, 5.4),
]
FAMILIES = {
    "Z006": dict(events=Z006, new="us7000ty3q", primary=True, b_partner="usp000a1fh", m_partner="us6000lfpq"),
    "Z007": dict(events=Z007, new="us6000ty0j", primary=False),
    "Z009": dict(events=Z009, new="us7000txqr", primary=False),
}
H13_STATIONS = ("IU.QSPA", "GT.VNDA", "IU.LCO", "GT.DBIC", "GT.PLCA", "GT.LPAZ", "GT.LBTB", "IU.LSZ")

# Pre-registered constants (pre-registration 15)
KEY = "D"                 # YKA PKIKP coda 1-6 s (yka.py's measured window); D_c310, D_c13 secondary
RATIO = 1.5               # key-contrast ratio needed for a B or M label
MIN_PARTNERS_RHO = 4      # ok partners needed for the clock correlations
MIN_SRC_STATIONS = 2      # mantle-P stations needed for the source control
MIN_OLD_GATED = 3         # for the family-membership score q (secondary)
TURN, BACK_RATIO = 2008.5, 2.5

REPO = "https://raw.githubusercontent.com/munimtechnologies/inner-core-repeaters/"
MANTLEP2_COMMIT = "b4712f0ebbc70b896231620c905ff53a4b8270fe"   # v0.13.0
PREREG_URL = REPO + "main/findings/2026-11-19-prereg15-z006-2026-event.md"


# ----------------------------------------------------------------------------- pure helpers
def isf(x):
    try:
        return math.isfinite(float(x))
    except (TypeError, ValueError):
        return False


def fnum(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return NAN


def median(xs):
    v = sorted(float(x) for x in xs if isf(x))
    n = len(v)
    if not n:
        return NAN
    return v[n // 2] if n % 2 else 0.5 * (v[n // 2 - 1] + v[n // 2])


def truthy(x):
    return x is True or str(x).strip() == "True"


def decimal_year(iso):
    """Same convention as z004_events.csv: year + (seconds since 1 Jan) / (days * 86400)."""
    s = iso.strip().rstrip("Z")
    date, hms = s.split("T")
    y, mo, d = (int(v) for v in date.split("-"))
    h, mi, sec = hms.split(":")
    t0 = calendar.timegm((y, 1, 1, 0, 0, 0, 0, 0, 0))
    t = calendar.timegm((y, mo, d, int(h), int(mi), 0, 0, 0, 0)) + float(sec)
    days = 366 if calendar.isleap(y) else 365
    return y + (t - t0) / 86400.0 / days


def event_row(e):
    eid, t, la, lo, dep, mag = e
    return dict(id=eid, time=t, lat=la, lon=lo, depth=dep, mag=mag, t=decimal_year(t))


def equiv(t):
    return t if t <= TURN else TURN - (t - TURN) / BACK_RATIO


def ranks(v):
    idx = sorted(range(len(v)), key=lambda i: v[i])
    r = [0.0] * len(v)
    i = 0
    while i < len(idx):
        j = i
        while j + 1 < len(idx) and v[idx[j + 1]] == v[idx[i]]:
            j += 1
        for k in range(i, j + 1):
            r[idx[k]] = 0.5 * (i + j)
        i = j + 1
    return r


def spearman(a, b):
    pr = [(x, y) for x, y in zip(a, b) if isf(x) and isf(y)]
    if len(pr) < 3:
        return NAN
    ra, rb = ranks([p[0] for p in pr]), ranks([p[1] for p in pr])
    ma, mb = sum(ra) / len(ra), sum(rb) / len(rb)
    sab = sum((x - ma) * (y - mb) for x, y in zip(ra, rb))
    saa = sum((x - ma) ** 2 for x in ra)
    sbb = sum((y - mb) ** 2 for y in rb)
    return sab / math.sqrt(saa * sbb) if saa > 0 and sbb > 0 else NAN


def partner_values(rows, new_id, key=KEY, ok_only=True):
    """{partner id: value} for rows pairing new_id with any other event."""
    out = {}
    for r in rows:
        a, b = r.get("ev1"), r.get("ev2")
        if new_id not in (a, b):
            continue
        if ok_only and not truthy(r.get("ok")):
            continue
        v = fnum(r.get(key))
        if isf(v):
            out[b if a == new_id else a] = v
    return out


def key_contrast(pv, b_partner, m_partner, ratio=RATIO):
    """P15-1: K = D(new, M partner) / D(new, B partner). B if K >= ratio, M if K <= 1/ratio."""
    out = dict(D_B=pv.get(b_partner), D_M=pv.get(m_partner), K=None, label="untestable")
    db, dm = fnum(pv.get(b_partner)), fnum(pv.get(m_partner))
    if not (isf(db) and isf(dm)) or db < 0 or dm < 0 or (db == 0 and dm == 0):
        return out
    k = float("inf") if db == 0 else dm / db
    out["K"] = k if math.isfinite(k) else None
    if k >= ratio:
        out["label"] = "B"
    elif k <= 1.0 / ratio:
        out["label"] = "M"
    else:
        out["label"] = "inconclusive"
    return out


def clock_test(pv, t_new, times, min_n=MIN_PARTNERS_RHO):
    """S15-a: Spearman of D(new, x) against the backtracking and the plain time distance."""
    xs = sorted(e for e in pv if e in times)
    d = [pv[e] for e in xs]
    gb = [abs(equiv(t_new) - equiv(times[e])) for e in xs]
    gm = [abs(t_new - times[e]) for e in xs]
    best = min(xs, key=lambda e: pv[e]) if xs else None
    b_pick = min(xs, key=lambda e: abs(equiv(t_new) - equiv(times[e]))) if xs else None
    m_pick = min(xs, key=lambda e: abs(t_new - times[e])) if xs else None
    ok = len(xs) >= min_n
    return dict(n=len(xs), partners=xs, D=d, dB=gb, dM=gm, rho_B=spearman(gb, d) if ok else NAN,
                rho_M=spearman(gm, d) if ok else NAN, best_match=best, b_pick=b_pick, m_pick=m_pick)


def source_ratios(D, gated, new_id, b_partner, m_partner):
    """Mantle-P ratio Km = D_late(new, M partner) / D_late(new, B partner) at one station, or NaN.
    D: {(a, b) sorted tuple: D_late}; counted only if all three events are gated."""
    if not (gated.get(new_id) and gated.get(b_partner) and gated.get(m_partner)):
        return NAN
    db = D.get(tuple(sorted((new_id, b_partner))), NAN)
    dm = D.get(tuple(sorted((new_id, m_partner))), NAN)
    if not (isf(db) and isf(dm)) or db <= 0:
        return NAN
    return dm / db


def source_score(D, gated, old_ids, new_id, min_old=MIN_OLD_GATED):
    """Secondary (as m26): q = share of gated old events whose median D_late to the other gated
    events is below the new event's."""
    old = [e for e in old_ids if gated.get(e)]
    out = dict(n_old_gated=len(old), new_gated=bool(gated.get(new_id)), counted=False, q=None)
    if not out["new_gated"] or len(old) < min_old:
        return out
    members = old + [new_id]
    m = {}
    for e in members:
        vals = [D.get(tuple(sorted((e, f))), NAN) for f in members if f != e]
        vals = [v for v in vals if isf(v)]
        m[e] = median(vals) if len(vals) >= 2 else NAN
    olds = [m[e] for e in old if isf(m[e])]
    if not isf(m[new_id]) or len(olds) < min_old:
        return out
    q = sum(1 for v in olds if v < m[new_id]) / float(len(olds))
    out.update(counted=True, q=q, m_new=m[new_id], m_old_median=median(olds))
    return out


def headline(core_label, kms, ratio=RATIO, n_min=MIN_SRC_STATIONS):
    """Pre-registered headline of pre-registration 15."""
    if core_label == "untestable":
        return "untestable"
    if core_label not in ("B", "M"):
        return core_label
    kms = [k for k in kms if isf(k)]
    if len(kms) < n_min:
        return "source-unverified:" + core_label
    km = median(kms)
    if core_label == "B":
        return "source-confounded:B" if km >= ratio else "B-supported"
    return "source-confounded:M" if km <= 1.0 / ratio else "M-supported"


def jsonable(o):
    if isinstance(o, dict):
        return {str(k): jsonable(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [jsonable(v) for v in o]
    if isinstance(o, float):
        return o if math.isfinite(o) else None
    if hasattr(o, "item") and not isinstance(o, (str, bytes)):
        try:
            return jsonable(o.item())
        except Exception:  # noqa: BLE001
            return str(o)
    return o


def write_csv(path, rows):
    keys = []
    for r in rows:
        for k in r:
            if k not in keys:
                keys.append(k)
    with open(path, "w", newline="") as f:
        if keys:
            w = csv.DictWriter(f, fieldnames=keys)
            w.writeheader()
            w.writerows(rows)


def fetch(url, tries=6, wait=20, timeout=90):
    for k in range(tries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "inner-core-repeaters (MV2, Munim, Inc.)"})
            with urllib.request.urlopen(req, timeout=timeout) as r:
                if r.status == 200:
                    return r.read().decode("utf-8")
        except Exception as ex:  # noqa: BLE001
            print("fetch error", url, repr(ex)[:160], flush=True)
        time.sleep(wait * (k + 1))
    return None


# ----------------------------------------------------------------------------- part Y (YKA)
def run_yka_part(mods, data_dir, out_dir, s):
    if mods not in sys.path:
        sys.path.insert(0, mods)
    import m0324_check as mc  # noqa: E402  (pinned v0.8.0)
    import xcorr  # noqa: E402
    import yka_repick  # noqa: E402  (v0.9.2, as published)

    yka = mc.load_yka_module("/tmp/icr_mods")
    if yka is None:
        s["yka_error"] = "yka.py download failed"
        return {}
    out = {}
    for fam, spec in FAMILIES.items():
        evs = {e[0]: event_row(e) for e in spec["events"]}
        pn = os.path.join(data_dir, "yka_raw", spec["new"] + ".mseed.nodata")
        if os.path.exists(pn):                 # a stale marker must not block a retry
            os.remove(pn)
            s.setdefault("removed_nodata_markers", []).append(spec["new"])
        sub = os.path.join(out_dir, "yka_" + fam)
        plots = os.path.join(sub, "plots")
        os.makedirs(plots, exist_ok=True)
        rec = []
        yka.pick_onset = yka_repick.constrained_picker(xcorr.pick_onset, yka_repick.REPICK_S, record=rec)
        sr = {}
        try:
            mc.run_yka(evs, list(evs), data_dir, sub, plots, sr)   # files are named yka_z004_* by run_yka
        except Exception as ex:  # noqa: BLE001
            traceback.print_exc()
            sr["yka"] = dict(status="exception", error=repr(ex)[:300])
        rows, st = [], []
        p = os.path.join(sub, "yka_z004_pairs.csv")
        if os.path.exists(p) and os.path.getsize(p) > 0:
            with open(p, newline="") as f:
                rows = list(csv.DictReader(f))
        p2 = os.path.join(sub, "yka_z004_status.csv")
        if os.path.exists(p2) and os.path.getsize(p2) > 0:
            with open(p2, newline="") as f:
                st = list(csv.DictReader(f))
        y = sr.get("yka", {})
        info = y.get("multiplet", {}) if isinstance(y, dict) else {}
        out[fam] = dict(rows=rows, status=st, evs=evs)
        s["yka_" + fam] = dict(n_events=len(evs), n_pairs=len(rows), n_ok=sum(1 for r in rows if truthy(r.get("ok"))),
                               multiplet=info, pick_valid=yka_repick.pick_valid(info.get("pick_minus_pred_s")),
                               status=st or None, pick_diagnostics=rec[:3])
    return out


# ----------------------------------------------------------------------------- part M (mantle P)
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


def run_mantle_part(data_dir, out_dir, s, deadline):
    spec = FAMILIES["Z006"]
    mp = load_mantlep2(os.path.join(out_dir, "_pinned"))
    if mp is None:
        s["mantle_error"] = "mantlep2.py v0.13.0 download failed"
        return [], [], {}
    mp.EVENTS = [(e[0], e[1].rstrip("Z"), e[2], e[3], e[4], e[5]) for e in spec["events"]]
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
    coords = {"%s.%s" % (n, st): (la, lo) for n, st, la, lo in (found or [])}
    cache = os.path.join(data_dir, "mantlep")
    plots = os.path.join(out_dir, "mantle", "plots")
    os.makedirs(cache, exist_ok=True)
    os.makedirs(plots, exist_ok=True)
    model = TauPyModel("ak135")
    new, bp, mpart = spec["new"], spec["b_partner"], spec["m_partner"]
    olds = [e[0] for e in spec["events"] if e[0] != new]
    ev_all, pair_all, per = [], [], {}
    for name in H13_STATIONS:
        if name not in coords:
            per[name] = dict(status="not_in_station_service")
            continue
        if time.time() > deadline:
            per[name] = dict(status="skipped_time_budget")
            continue
        net, sta = name.split(".")
        ev_rows, pair_rows = [], []
        try:
            res = mp.run_station(session, model, net, sta, coords[name][0], coords[name][1], cache, plots,
                                 ev_rows, pair_rows, {})
        except Exception as ex:  # noqa: BLE001
            traceback.print_exc()
            per[name] = dict(status="exception", error=repr(ex)[:200])
            continue
        gated = {r["event"]: truthy(r.get("gated")) for r in ev_rows}
        D = {tuple(sorted((r["ev1"], r["ev2"]))): fnum(r["D_late"]) for r in pair_rows}
        km = source_ratios(D, gated, new, bp, mpart)
        sc = source_score(D, gated, olds, new)
        sc.update(status=res.get("status"), ref=res.get("ref"), onset_minus_pred_s=res.get("onset_minus_pred_s"),
                  Km=km, gated=gated,
                  D_late_new_B=D.get(tuple(sorted((new, bp)))), D_late_new_M=D.get(tuple(sorted((new, mpart)))))
        per[name] = sc
        ev_all += ev_rows
        pair_all += pair_rows
        print(name, json.dumps(jsonable(sc))[:600], flush=True)
    return ev_all, pair_all, per


# ----------------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mods", required=True, help="folder with the pinned v0.8.0 modules (fetch_pinned.py)")
    ap.add_argument("--data-dir", default="/data")
    ap.add_argument("--out-dir", default="/out")
    ap.add_argument("--time-budget-min", type=float, default=120.0)
    args = ap.parse_args()
    t0 = time.time()
    deadline = t0 + 60.0 * args.time_budget_min
    os.makedirs(args.out_dir, exist_ok=True)
    mods = os.path.abspath(args.mods)
    s = dict(version=VERSION, started=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
             families={k: dict(new=v["new"], events=[e[0] for e in v["events"]], primary=v["primary"])
                       for k, v in FAMILIES.items()},
             constants=dict(key=KEY, ratio=RATIO, min_partners_rho=MIN_PARTNERS_RHO,
                            min_src_stations=MIN_SRC_STATIONS, min_old_gated=MIN_OLD_GATED))
    s["prereg_public_at_start"] = bool(fetch(PREREG_URL, tries=2, wait=5))

    yk = {}
    try:
        yk = run_yka_part(mods, args.data_dir, args.out_dir, s)
    except Exception as ex:  # noqa: BLE001
        traceback.print_exc()
        s["yka_error"] = repr(ex)[:300]
    grade = {}
    for fam, spec in FAMILIES.items():
        if fam not in yk:
            continue
        rows = yk[fam]["rows"]
        times = {e: v["t"] for e, v in yk[fam]["evs"].items()}
        t_new = times[spec["new"]]
        g = dict(t_new=t_new, equiv_new=equiv(t_new))
        for key in ("D", "D_c310", "D_c13"):
            pv = partner_values(rows, spec["new"], key, True)
            g["clock_" + key] = clock_test(pv, t_new, times)
            g["ungated_clock_" + key] = clock_test(partner_values(rows, spec["new"], key, False), t_new, times)
            if spec["primary"]:
                g["contrast_" + key] = key_contrast(pv, spec["b_partner"], spec["m_partner"])
                g["ungated_contrast_" + key] = key_contrast(partner_values(rows, spec["new"], key, False),
                                                            spec["b_partner"], spec["m_partner"])
        grade[fam] = g

    ev_all, pair_all, per = [], [], {}
    try:
        ev_all, pair_all, per = run_mantle_part(args.data_dir, args.out_dir, s, deadline)
    except Exception as ex:  # noqa: BLE001
        traceback.print_exc()
        s["mantle_error"] = repr(ex)[:300]
    s["mantle_stations"] = per
    kms = [v.get("Km") for v in per.values() if isf(v.get("Km"))]
    qs = [v["q"] for v in per.values() if v.get("counted")]
    core = grade.get("Z006", {}).get("contrast_" + KEY, {}).get("label", "untestable")
    grade["P15_0_source"] = dict(Km=kms, median_Km=median(kms), n_stations=len(kms),
                                 membership_q=qs, membership_median_q=median(qs))
    grade["P15_1_core"] = core
    grade["headline"] = headline(core, kms)
    s["grade"] = grade
    s["finished"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    s["minutes"] = round((time.time() - t0) / 60.0, 1)
    write_csv(os.path.join(args.out_dir, "m27_mantle_events.csv"), ev_all)
    write_csv(os.path.join(args.out_dir, "m27_mantle_pairs.csv"), pair_all)
    for fam in yk:
        write_csv(os.path.join(args.out_dir, "m27_yka_pairs_%s.csv" % fam), yk[fam]["rows"])
    with open(os.path.join(args.out_dir, "m27_summary.json"), "w") as f:
        json.dump(jsonable(s), f, indent=1)
    print(json.dumps(jsonable(grade))[:4000], flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())

#!/usr/bin/env python3
"""m26.py: prospective test of a new 2026 South Sandwich Islands earthquake at multiplet M03-24.

Part of inner-core-repeaters by MV2, Munim, Inc. (made with AI agents). Licence: Apache-2.0.
Pre-registered in findings/2026-11-19-prereg14-m0324-2026-event.md (pre-registration 14),
written and published before any waveform of the new event was fetched or looked at.

USGS ComCat lists us7000ty3q (2026-11-13 19:04:07.7 UTC, -58.962, -25.7253, 63.8 km, mb 5.1)
inside the catalogue footprint of multiplet M03-24 (Z004). At YKA, M03-24's PKIKP late coda
(3-10 s, "D_c310") splits into X = {2003, 2020, 2022-03, 2024, 2025-12} and Y = {2008, 2011,
2014} (found exploratorily in icr-m0324chk-v2); mantlep2-v1 showed the split is not in the
mantle P (decision N), i.e. it is on the core path. No published study uses a 2026 event.

Two parts, both with code already published in this repo:
  M  mantle-P source check: mantlep2.py v0.13.0 (pinned commit) with the new event added, at
     the eight stations of pre-registration 13. Score q per station = share of the old gated
     events whose median mantle-P D_late to the others is below the new event's.
  Y  core-coda state at YKA: m0324_check.run_yka (v0.8.0, pinned) with yka_repick's +/-2 s
     picker, exactly as icr-m0324chk-v2, on the nine Z004 events plus the new event.
     Group of the new event from median D_c310 to the X and Y members (ok pairs only).
Three YKA runs: "old9" (reproduction of the published values), "primary" (+ us7000ty3q) and
"extra" (+ us7000ty3q + us6000thwz, 2026-08-04 mb 4.5; exploratory).
"""
import argparse
import calendar
import csv
import io
import json
import math
import os
import sys
import time
import traceback
import urllib.request

VERSION = "0.1.0"
NAN = float("nan")

NEW_ID = "us7000ty3q"
NEW_EVENT = dict(id=NEW_ID, time="2026-11-13T19:04:07.703Z", lat=-58.962, lon=-25.7253, depth=63.81, mag=5.1)
EXTRA_EVENTS = [dict(id="us6000thwz", time="2026-08-04T02:28:42.201Z", lat=-58.9198, lon=-25.596,
                     depth=35.0, mag=4.5)]

X_GROUP = ("usp000bt5x", "us7000am0e", "us6000h2wz", "us6000mk65", "us6000runq")
Y_GROUP = ("usp000fz4d", "usp000hute", "usc000t4wj")
EXCLUDED = ("us6000ittv",)        # 2022-10-15, SNR 6 at YKA, left out of the split as before
OLD8 = ("usp000bt5x", "usp000fz4d", "usp000hute", "usc000t4wj", "us7000am0e", "us6000h2wz",
        "us6000mk65", "us6000runq")
H13_STATIONS = ("IU.QSPA", "GT.VNDA", "IU.LCO", "GT.DBIC", "GT.PLCA", "GT.LPAZ", "GT.LBTB", "IU.LSZ")

# Pre-registered thresholds (pre-registration 14)
NEAR = 0.05          # D_c310 at most this: "alike" (published within-group max 0.033)
FAR = 0.07           # D_c310 at least this: "different" (published between-group min 0.098)
RATIO = 3.0
MIN_X, MIN_Y = 3, 2  # ok partners needed in each group
SRC_SAME_Q = 0.75
SRC_DIFF_Q = 0.999
MIN_OLD_GATED = 5
MIN_SRC_STATIONS = 3
TURN, BACK_RATIO = 2008.5, 2.5   # Wang et al. 2024 backtracking clock

REPO = "https://raw.githubusercontent.com/munimtechnologies/inner-core-repeaters/"
MANTLEP2_COMMIT = "b4712f0ebbc70b896231620c905ff53a4b8270fe"   # v0.13.0, mantlep2.py 0.13.0
MANTLEP2_PAIRS = REPO + "e4c618a566be2324fa96dfa94f9c326539c90846/results/mantlep2-v1/mantlep2_pairs.csv"
YKA_V2_PAIRS = REPO + "main/results/m0324chk-v2/yka_z004_pairs.csv"
PREREG_URL = REPO + "main/findings/2026-11-19-prereg14-m0324-2026-event.md"


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
    """Same convention as z004_events.csv: year + (day_of_year - 1 + fraction of day) / days."""
    s = iso.strip().rstrip("Z")
    date, hms = s.split("T")
    y, mo, d = (int(v) for v in date.split("-"))
    h, mi, sec = hms.split(":")
    t0 = calendar.timegm((y, 1, 1, 0, 0, 0, 0, 0, 0))
    t = calendar.timegm((y, mo, d, int(h), int(mi), 0, 0, 0, 0)) + float(sec)
    days = 366 if calendar.isleap(y) else 365
    return y + (t - t0) / 86400.0 / days


def event_row(ev):
    r = dict(ev)
    r["t"] = decimal_year(ev["time"])
    return r


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


def partner_values(rows, new_id, key="D_c310", ok_only=True):
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


def classify_state(dX, dY, nX, nY, near=NEAR, far=FAR, ratio=RATIO):
    """Pre-registered P14-1 label of the new event's late-coda state."""
    if nX < MIN_X or nY < MIN_Y or not (isf(dX) and isf(dY)):
        return "untestable"
    if dX <= near and dY >= far and dY >= ratio * dX:
        return "X"
    if dY <= near and dX >= far and dX >= ratio * dY:
        return "Y"
    if min(dX, dY) > near:
        return "new-state"
    return "ambiguous"


def state_test(rows, new_id, key="D_c310", ok_only=True):
    pv = partner_values(rows, new_id, key, ok_only)
    xs = {e: v for e, v in pv.items() if e in X_GROUP}
    ys = {e: v for e, v in pv.items() if e in Y_GROUP}
    dX, dY = median(xs.values()), median(ys.values())
    return dict(key=key, ok_only=ok_only, nX=len(xs), nY=len(ys), dX=dX, dY=dY,
                label=classify_state(dX, dY, len(xs), len(ys)), partners=pv)


def backtracking_order(rows, new_id, t_new, times, key="D_c310"):
    """S14-a: Spearman of D(new, x) against |equiv(t_new) - equiv(t_x)| over ok X partners."""
    pv = partner_values(rows, new_id, key, True)
    xs = [e for e in X_GROUP if e in pv and e in times]
    d = [pv[e] for e in xs]
    g = [abs(equiv(t_new) - equiv(times[e])) for e in xs]
    return dict(n=len(xs), rho=spearman(g, d), partners=xs, dB=g, D=d)


def source_score(D, gated, old_ids, new_id, min_old=MIN_OLD_GATED):
    """q = share of gated old events whose median mantle-P D_late to the other gated events
    (old + new) is below the new event's. D: {(a, b) sorted tuple: value}."""
    old = [e for e in old_ids if gated.get(e)]
    out = dict(n_old_gated=len(old), new_gated=bool(gated.get(new_id)), counted=False, q=None)
    if not out["new_gated"] or len(old) < min_old:
        return out
    members = old + [new_id]

    def d(a, b):
        return D.get(tuple(sorted((a, b))), NAN)

    m = {}
    for e in members:
        vals = [d(e, f) for f in members if f != e]
        vals = [v for v in vals if isf(v)]
        m[e] = median(vals) if len(vals) >= 3 else NAN
    olds = [m[e] for e in old if isf(m[e])]
    if not isf(m[new_id]) or len(olds) < min_old:
        return out
    q = sum(1 for v in olds if v < m[new_id]) / float(len(olds))
    out.update(counted=True, q=q, m_new=m[new_id], m_old_median=median(olds), m_old_max=max(olds),
               n_old_scored=len(olds))
    return out


def decide_source(qs, n_min=MIN_SRC_STATIONS):
    qs = [q for q in qs if isf(q)]
    if len(qs) < n_min:
        return "untestable"
    m = median(qs)
    if m <= SRC_SAME_Q:
        return "same"
    if m >= SRC_DIFF_Q:
        return "different"
    return "inconclusive"


def headline(source, state):
    if state == "untestable":
        return "untestable"
    if source == "different":
        return "not-a-repeat"
    if source != "same":
        return "source-unverified:" + state
    return {"X": "X-confirmed", "Y": "change-detected", "new-state": "change-detected"}.get(state, "ambiguous")


def max_abs_diff(rows_a, rows_b, key, stations=False):
    """Largest |difference| of `key` over pairs present in both row lists."""
    def k(r):
        p = tuple(sorted((r["ev1"], r["ev2"])))
        return (r.get("station"),) + p if stations else p
    b = {k(r): fnum(r.get(key)) for r in rows_b}
    diffs = [abs(fnum(r.get(key)) - b[k(r)]) for r in rows_a if k(r) in b and isf(r.get(key)) and isf(b[k(r)])]
    return dict(n=len(diffs), max_abs_diff=max(diffs) if diffs else None)


def read_csv_text(text):
    return list(csv.DictReader(io.StringIO(text))) if text else []


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


# ----------------------------------------------------------------------------- I/O
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
    import yka_repick  # noqa: E402  (this folder, v0.9.2 as published)

    yka = mc.load_yka_module("/tmp/icr_mods")
    if yka is None:
        s["yka_error"] = "yka.py download failed"
        return {}
    for ev in [NEW_EVENT] + EXTRA_EVENTS:      # a stale "nodata" marker must not block a retry
        pn = os.path.join(data_dir, "yka_raw", ev["id"] + ".mseed.nodata")
        if os.path.exists(pn):
            os.remove(pn)
            s.setdefault("removed_nodata_markers", []).append(ev["id"])
    old = mc.read_events(os.path.join(mods, "z004_events.csv"))
    runs = {"old9": [], "primary": [NEW_EVENT], "extra": [NEW_EVENT] + EXTRA_EVENTS}
    out = {}
    for name, extra in runs.items():
        evs = dict(old)
        for ev in extra:
            evs[ev["id"]] = event_row(ev)
        sub = os.path.join(out_dir, "yka_" + name)
        plots = os.path.join(sub, "plots")
        os.makedirs(plots, exist_ok=True)
        rec = []
        yka.pick_onset = yka_repick.constrained_picker(xcorr.pick_onset, yka_repick.REPICK_S, record=rec)
        sr = {}
        try:
            mc.run_yka(evs, list(evs), data_dir, sub, plots, sr)
        except Exception as ex:  # noqa: BLE001
            traceback.print_exc()
            sr["yka"] = dict(status="exception", error=repr(ex)[:300])
        rows = []
        p = os.path.join(sub, "yka_z004_pairs.csv")
        if os.path.exists(p) and os.path.getsize(p) > 0:
            with open(p, newline="") as f:
                rows = list(csv.DictReader(f))
        st = []
        p2 = os.path.join(sub, "yka_z004_status.csv")
        if os.path.exists(p2) and os.path.getsize(p2) > 0:
            with open(p2, newline="") as f:
                st = list(csv.DictReader(f))
        y = sr.get("yka", {})
        info = y.get("multiplet", {}) if isinstance(y, dict) else {}
        out[name] = dict(rows=rows, status=st, evs=evs)
        s["yka_" + name] = dict(n_events=len(evs), n_pairs=len(rows), n_ok=sum(1 for r in rows if truthy(r["ok"])),
                                multiplet=info, pick_valid=yka_repick.pick_valid(info.get("pick_minus_pred_s")),
                                status=[r for r in st if r.get("id") in [e["id"] for e in extra]] or None,
                                pick_diagnostics=rec[:3])
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
    mp = load_mantlep2(os.path.join(out_dir, "_pinned"))
    if mp is None:
        s["mantle_error"] = "mantlep2.py v0.13.0 download failed"
        return [], [], {}
    t = NEW_EVENT["time"].rstrip("Z")
    mp.EVENTS = list(mp.EVENTS) + [(NEW_ID, t, NEW_EVENT["lat"], NEW_EVENT["lon"], NEW_EVENT["depth"],
                                    NEW_EVENT["mag"])]
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
        sc = source_score(D, gated, list(OLD8), NEW_ID)
        new_row = [r for r in ev_rows if r["event"] == NEW_ID]
        sc.update(status=res.get("status"), ref=res.get("ref"), onset_minus_pred_s=res.get("onset_minus_pred_s"),
                  new_event=new_row[0] if new_row else None,
                  cc_direct_new_median=median([r["cc_direct"] for r in pair_rows if NEW_ID in (r["ev1"], r["ev2"])
                                               and gated.get(r["ev1"]) and gated.get(r["ev2"])]))
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
             new_event=NEW_EVENT, extra_events=EXTRA_EVENTS, X=X_GROUP, Y=Y_GROUP, excluded=EXCLUDED,
             thresholds=dict(near=NEAR, far=FAR, ratio=RATIO, min_x=MIN_X, min_y=MIN_Y, src_same_q=SRC_SAME_Q,
                             src_diff_q=SRC_DIFF_Q, min_old_gated=MIN_OLD_GATED, min_src_stations=MIN_SRC_STATIONS))
    s["prereg_public_at_start"] = bool(fetch(PREREG_URL, tries=2, wait=5))

    # ---- part Y
    yk = {}
    try:
        yk = run_yka_part(mods, args.data_dir, args.out_dir, s)
    except Exception as ex:  # noqa: BLE001
        traceback.print_exc()
        s["yka_error"] = repr(ex)[:300]
    pub = read_csv_text(fetch(YKA_V2_PAIRS, tries=3))
    for name in ("old9", "primary", "extra"):
        if name in yk:
            s["R14_yka_" + name] = {k: max_abs_diff(yk[name]["rows"], pub, k) for k in ("D", "D_c310")}
    grade = {}
    if "primary" in yk:
        rows = yk["primary"]["rows"]
        times = {e: v["t"] for e, v in yk["primary"]["evs"].items()}
        t_new = event_row(NEW_EVENT)["t"]
        grade["P14_1_state"] = state_test(rows, NEW_ID, "D_c310", True)
        grade["S14_ungated_state"] = state_test(rows, NEW_ID, "D_c310", False)
        grade["S14_D_1to6_values"] = partner_values(rows, NEW_ID, "D", True)
        grade["S14_D_c13_values"] = partner_values(rows, NEW_ID, "D_c13", True)
        grade["S14a_backtracking_order"] = backtracking_order(rows, NEW_ID, t_new, times)
        grade["t_new"] = t_new
        grade["equiv_new"] = equiv(t_new)
    if "extra" in yk:
        for ev in EXTRA_EVENTS:
            grade["S14b_extra_" + ev["id"]] = state_test(yk["extra"]["rows"], ev["id"], "D_c310", True)
        grade["S14b_primary_event_in_extra_run"] = state_test(yk["extra"]["rows"], NEW_ID, "D_c310", True)

    # ---- part M
    ev_all, pair_all, per = [], [], {}
    try:
        ev_all, pair_all, per = run_mantle_part(args.data_dir, args.out_dir, s, deadline)
    except Exception as ex:  # noqa: BLE001
        traceback.print_exc()
        s["mantle_error"] = repr(ex)[:300]
    s["mantle_stations"] = per
    qs = [v["q"] for v in per.values() if v.get("counted")]
    grade["P14_0_source"] = dict(decision=decide_source(qs), n_stations=len(qs), q=qs, median_q=median(qs))
    pubm = read_csv_text(fetch(MANTLEP2_PAIRS, tries=3))
    old_pairs = [r for r in pair_all if NEW_ID not in (r["ev1"], r["ev2"])]
    s["R14_mantle_old_pairs"] = max_abs_diff(old_pairs, pubm, "D_late", stations=True)

    st = grade.get("P14_1_state", {}).get("label", "untestable")
    grade["headline"] = headline(grade["P14_0_source"]["decision"], st)
    s["grade"] = grade
    s["finished"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    s["minutes"] = round((time.time() - t0) / 60.0, 1)
    write_csv(os.path.join(args.out_dir, "m26_mantle_events.csv"), ev_all)
    write_csv(os.path.join(args.out_dir, "m26_mantle_pairs.csv"), pair_all)
    for name in yk:
        write_csv(os.path.join(args.out_dir, "m26_yka_pairs_%s.csv" % name), yk[name]["rows"])
    with open(os.path.join(args.out_dir, "m26_summary.json"), "w") as f:
        json.dump(jsonable(s), f, indent=1)
    print(json.dumps(jsonable(grade))[:4000], flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())

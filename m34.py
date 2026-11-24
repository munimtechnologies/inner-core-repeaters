#!/usr/bin/env python3
"""m34.py: mirror test of inner-core backtracking with unused members of multiplet M03-24.

Part of inner-core-repeaters by MV2, Munim, Inc. (made with AI agents). Licence: Apache-2.0.
Pre-registered in findings/2026-11-24-prereg34-m0324-mirror-plan.md (pre-registration 34),
written and published before any waveform of the 31 candidates below was fetched or looked at.

Background. At YKA the late PKIKP coda (3-10 s, D_c310) of the South Sandwich multiplet M03-24
splits into X = {2003-03, 2020-07, 2022-03, 2024-03, 2025-12} and Y = {2008-02, 2011-02, 2014-12};
mantle P shows no split (pre-registration 11, grade N), so the split is on the core path.
A single-turn rotation clock (inner core turns back once, Wang et al. 2024) puts every member on
one effective-time axis, so the states must form contiguous blocks there: X, then (optionally)
intermediate, then Y. Local, episodic change of the inner-core surface (Zhang & Wen 2026) does not
require this. Members between 2003 and 2008 and between 2014-12 and 2020-07 are what can break it,
and no published study of this multiplet uses any.

For each of 31 USGS ComCat candidates (within 25 km of 59.03 S 25.68 W, mb >= 4.6, 2003-03-22..
2008-02-16 and 2014-12-10..2020-07-10):
  M  membership: mantlep2.py v0.13.0 (pinned) at the eight stations of pre-registration 13 with
     the old eight members plus this one candidate; per-station q as in pre-registration 14.
  Y  state at YKA: m0324_check.run_yka (v0.8.0, pinned) with yka_repick's +/-2 s picker on the
     nine Z004 events plus this one candidate; label from median D_c310 to X and to Y (ok pairs).
Then the mirror tests R1 (Wang clock) and R2 (any single-turn clock on a grid) on the confirmed
members, and the recovery bracket.
"""
import argparse
import calendar
import csv
import io
import json
import math
import os
import shutil
import sys
import time
import traceback
import urllib.request

VERSION = "0.1.0"
NAN = float("nan")

X_GROUP = ("usp000bt5x", "us7000am0e", "us6000h2wz", "us6000mk65", "us6000runq")
Y_GROUP = ("usp000fz4d", "usp000hute", "usc000t4wj")
OLD8 = ("usp000bt5x", "usp000fz4d", "usp000hute", "usc000t4wj", "us7000am0e", "us6000h2wz",
        "us6000mk65", "us6000runq")
OLD_TIMES = {  # decimal years, z004_events.csv convention
    "usp000bt5x": "2003-03-21T10:47:47.620Z", "usp000fz4d": "2008-02-17T19:32:09.100Z",
    "usp000hute": "2011-02-20T17:37:37.320Z", "usc000t4wj": "2014-12-09T19:23:36.400Z",
    "us7000am0e": "2020-07-11T18:58:06.501Z", "us6000h2wz": "2022-03-07T23:25:40.208Z",
    "us6000mk65": "2024-03-21T02:47:21.923Z", "us6000runq": "2025-12-13T14:04:33.132Z"}
H13_STATIONS = ("IU.QSPA", "GT.VNDA", "IU.LCO", "GT.DBIC", "GT.PLCA", "GT.LPAZ", "GT.LBTB", "IU.LSZ")

# USGS ComCat query (2026-11-24): latitude=-59.03, longitude=-25.68, maxradiuskm=25,
# minmagnitude=4.6; every event in the two windows. (id, time, lat, lon, depth km, mag)
CANDIDATES = [
    ("usp000bxap", "2003-05-16T05:55:03.270Z", -59.156, -25.831, 33.0, 4.8),
    ("usp000by6r", "2003-05-26T11:47:54.500Z", -59.016, -25.372, 33.0, 4.8),
    ("usp000c3yx", "2003-08-03T05:38:57.550Z", -58.915, -25.561, 33.0, 5.0),
    ("usp000cf4p", "2003-12-13T00:27:06.550Z", -59.02, -25.452, 33.0, 4.6),
    ("usp000dxs7", "2005-08-21T07:12:18.520Z", -59.243, -25.675, 58.9, 5.1),
    ("usp000dzku", "2005-09-17T11:56:22.430Z", -59.152, -25.828, 39.7, 4.7),
    ("usp000e7e5", "2006-01-02T14:45:26.860Z", -59.221, -25.871, 10.0, 4.6),
    ("usp000eu1v", "2006-09-29T05:50:49.400Z", -59.001, -25.757, 70.4, 4.8),
    ("usp000f0ed", "2006-12-17T19:05:42.570Z", -59.022, -25.661, 37.8, 4.6),
    ("usp000fq8f", "2007-10-10T10:08:43.750Z", -58.944, -25.561, 46.4, 4.7),
    ("usp000fwhh", "2008-01-12T05:53:54.990Z", -58.999, -25.681, 48.0, 4.6),
    ("us1000328u", "2015-08-14T09:37:03.230Z", -58.9719, -25.6942, 10.0, 4.9),
    ("us10003y36", "2015-11-13T14:04:05.970Z", -59.1712, -25.489, 35.0, 4.9),
    ("us1000419q", "2015-11-22T07:45:34.960Z", -58.9358, -25.5681, 35.0, 4.8),
    ("us100044ec", "2015-12-06T08:53:28.450Z", -59.1301, -25.7922, 50.04, 4.8),
    ("us10004bgp", "2016-01-05T10:33:40.340Z", -59.0472, -25.5758, 53.43, 4.8),
    ("us10004fdy", "2016-01-12T18:02:13.900Z", -59.2136, -25.8201, 62.88, 5.0),
    ("us20004vf6", "2016-01-29T04:30:46.290Z", -58.9775, -25.8249, 97.67, 5.1),
    ("us10005gaq", "2016-04-30T12:50:43.560Z", -58.9991, -25.9011, 97.73, 4.7),
    ("us10006d1x", "2016-08-11T18:06:06.410Z", -58.9588, -25.3369, 43.87, 5.0),
    ("us2000asmb", "2017-09-15T08:34:36.770Z", -59.006, -25.3238, 10.0, 4.8),
    ("us2000au57", "2017-09-25T06:18:03.290Z", -59.0431, -25.6394, 10.0, 5.0),
    ("us1000dnn1", "2018-04-14T18:02:50.880Z", -59.2251, -25.5066, 35.0, 4.6),
    ("us1000dsql", "2018-04-22T22:30:10.200Z", -58.9399, -25.3884, 71.89, 4.6),
    ("us1000i46d", "2018-12-02T16:35:28.770Z", -58.9263, -25.6575, 35.0, 4.7),
    ("us1000is3e", "2019-01-13T12:08:00.670Z", -59.1232, -25.6036, 63.71, 4.6),
    ("us1000jgkm", "2019-03-07T21:29:48.910Z", -59.1053, -26.0556, 35.0, 4.6),
    ("us70003jaz", "2019-05-01T21:08:56.344Z", -59.0194, -25.6465, 65.71, 4.7),
    ("us70006w5y", "2020-01-08T02:27:39.193Z", -59.025, -25.6708, 38.72, 4.8),
    ("us70007pxg", "2020-02-14T09:55:28.590Z", -58.905, -25.5422, 35.0, 5.1),
    ("us7000aj2c", "2020-07-06T20:53:15.378Z", -58.823, -25.6731, 60.12, 5.2),
]

# ---- pre-registered thresholds (pre-registration 34)
X_NEAR = 0.07        # a group is "near" if its median D_c310 is at most this
FAR = 0.07           # the other group must be at least this far ...
RATIO = 3.0          # ... and at least RATIO times the near value
BETWEEN_MIN = 0.098  # smallest published X-Y D_c310: below it for both groups -> intermediate
MIN_X, MIN_Y = 3, 2
SRC_SAME_Q = 0.75
SRC_DIFF_Q = 0.999
MIN_OLD_GATED = 5
MIN_SRC_STATIONS = 3
WANG_TURN, WANG_RATIO = 2008.5, 2.5
GRID_TURNS = [round(2007.0 + 0.1 * k, 1) for k in range(51)]     # 2007.0 .. 2012.0
GRID_RATIOS = [round(0.5 + 0.1 * k, 1) for k in range(46)]       # 0.5 .. 5.0
ORDER = {"X": 0, "I": 1, "Y": 2}
STATE_CODE = {"X": "X", "Y": "Y", "intermediate": "I"}

REPO = "https://raw.githubusercontent.com/munimtechnologies/inner-core-repeaters/"
MANTLEP2_COMMIT = "b4712f0ebbc70b896231620c905ff53a4b8270fe"   # v0.13.0, mantlep2.py 0.13.0
YKA_V2_PAIRS = REPO + "main/results/m0324chk-v2/yka_z004_pairs.csv"
PREREG_URL = REPO + "main/findings/2026-11-24-prereg34-m0324-mirror-plan.md"


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


def cand_dict(c):
    cid, t, lat, lon, dep, mag = c
    return dict(id=cid, time=t, lat=lat, lon=lon, depth=dep, mag=mag, t=decimal_year(t))


def equiv(t, turn=WANG_TURN, ratio=WANG_RATIO):
    """Effective date on a single-turn clock: forward until `turn`, then back `ratio` times slower."""
    return t if t <= turn else turn - (t - turn) / ratio


def old_times():
    return {e: decimal_year(v) for e, v in OLD_TIMES.items()}


def wang_boundary():
    """(largest effective date of an old X member, smallest of an old Y member) on the Wang clock."""
    ot = old_times()
    return max(equiv(ot[e]) for e in X_GROUP), min(equiv(ot[e]) for e in Y_GROUP)


def wang_prediction(t):
    """What the Wang clock and the old eight members require of a member at time t."""
    lo, hi = wang_boundary()
    e = equiv(t)
    if e < lo:
        return "X"
    if e > hi:
        return "Y"
    return "free"


def classify34(dX, dY, nX, nY):
    """Pre-registered P34 state label from median D_c310 to the X and Y members."""
    if nX < MIN_X or nY < MIN_Y or not (isf(dX) and isf(dY)):
        return "untestable"
    if dX <= X_NEAR and dY >= FAR and dY >= RATIO * dX:
        return "X"
    if dY <= X_NEAR and dX >= FAR and dX >= RATIO * dY:
        return "Y"
    if max(dX, dY) < BETWEEN_MIN:
        return "intermediate"
    if min(dX, dY) >= BETWEEN_MIN:
        return "new-state"
    return "ambiguous"


def partner_values(rows, new_id, key="D_c310", ok_only=True):
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


def state_test(rows, new_id, key="D_c310", ok_only=True):
    pv = partner_values(rows, new_id, key, ok_only)
    xs = {e: v for e, v in pv.items() if e in X_GROUP}
    ys = {e: v for e, v in pv.items() if e in Y_GROUP}
    dX, dY = median(xs.values()), median(ys.values())
    s = dX / (dX + dY) if isf(dX) and isf(dY) and (dX + dY) > 0 else NAN
    return dict(key=key, ok_only=ok_only, nX=len(xs), nY=len(ys), dX=dX, dY=dY, s=s,
                label=classify34(dX, dY, len(xs), len(ys)), partners=pv)


def source_score(D, gated, old_ids, new_id, min_old=MIN_OLD_GATED):
    """Pre-registration 14's q: share of gated old events whose median mantle-P D_late to the
    other gated events (old + new) is below the new event's."""
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


def violations(members, turn, ratio):
    """members: {id: (t, code)} with code in X/I/Y. A violation is a pair whose state order
    (X < I < Y) disagrees with its effective-date order on the clock (turn, ratio)."""
    items = sorted(members.items())
    out = []
    for i in range(len(items)):
        for j in range(i + 1, len(items)):
            (a, (ta, ca)), (b, (tb, cb)) = items[i], items[j]
            if ORDER[ca] == ORDER[cb]:
                continue
            lo, hi = (a, b) if ORDER[ca] < ORDER[cb] else (b, a)
            tl, th = members[lo][0], members[hi][0]
            if equiv(tl, turn, ratio) >= equiv(th, turn, ratio):
                out.append((lo, hi))
    return out


def old_members():
    ot = old_times()
    m = {e: (ot[e], "X") for e in X_GROUP}
    m.update({e: (ot[e], "Y") for e in Y_GROUP})
    return m


def compatible_clocks(members, turns=GRID_TURNS, ratios=GRID_RATIOS):
    return [(T, r) for T in turns for r in ratios if not violations(members, T, r)]


def mirror_tests(new):
    """new: {id: (t, code)} confirmed new members with code X/I/Y. Returns R1 and R2."""
    old = old_members()
    if not new:
        return dict(R1=dict(label="untestable", n_new=0), R2=dict(label="untestable", n_new=0),
                    baseline_compatible=len(compatible_clocks(old)))
    allm = dict(old)
    allm.update(new)
    v = violations(allm, WANG_TURN, WANG_RATIO)
    if not v:
        r1 = "consistent"
    else:
        robust = True
        for k in new:
            sub = {e: x for e, x in allm.items() if e != k}
            if not violations(sub, WANG_TURN, WANG_RATIO):
                robust = False
                break
        r1 = "violated" if robust else "violated-fragile"
    comp = compatible_clocks(allm)
    if comp:
        r2 = "compatible"
    else:
        r2 = "contradicted"
        for k in new:
            sub = {e: x for e, x in allm.items() if e != k}
            if compatible_clocks(sub):
                r2 = "contradicted-fragile"
                break
    turns = sorted(set(T for T, _ in comp))
    rats = sorted(set(r for _, r in comp))
    return dict(R1=dict(label=r1, n_new=len(new), violations=v),
                R2=dict(label=r2, n_new=len(new), n_compatible=len(comp),
                        turn_range=[turns[0], turns[-1]] if turns else None,
                        ratio_range=[rats[0], rats[-1]] if rats else None,
                        wang_in_set=(WANG_TURN, WANG_RATIO) in comp),
                baseline_compatible=len(compatible_clocks(old)))


def brackets(new):
    """Outbound bracket (last X before first Y, t < 2011) and return bracket (last Y before first
    X, t > 2011) over old and confirmed new members with code X or Y."""
    allm = old_members()
    allm.update(new)
    out = {}
    for leg, (a, b), sel in (("outbound", ("X", "Y"), lambda t: t < 2011.0),
                             ("return", ("Y", "X"), lambda t: t > 2011.0)):
        ta = [t for t, c in allm.values() if sel(t) and c == a]
        tb = [t for t, c in allm.values() if sel(t) and c == b]
        if not ta or not tb:
            out[leg] = dict(label="untestable")
            continue
        last_a, first_b = max(ta), min(tb)
        flips = sum(1 for t in ta if t > first_b) + sum(1 for t in tb if t < last_a)
        out[leg] = dict(last_before=round(last_a, 4), first_after=round(first_b, 4),
                        width_yr=round(first_b - last_a, 3), label="clean" if last_a < first_b else "overlap",
                        n_out_of_order=flips,
                        intermediates=sorted(round(t, 4) for t, c in allm.values() if sel(t) and c == "I"))
    return out


def max_abs_diff(rows_a, rows_b, key):
    def k(r):
        return tuple(sorted((r["ev1"], r["ev2"])))
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
def read_rows(path):
    if os.path.exists(path) and os.path.getsize(path) > 0:
        with open(path, newline="") as f:
            return list(csv.DictReader(f))
    return []


def run_yka_part(mods, data_dir, out_dir, s, deadline):
    if mods not in sys.path:
        sys.path.insert(0, mods)
    import m0324_check as mc  # noqa: E402  (pinned v0.8.0)
    import xcorr  # noqa: E402
    import yka_repick  # noqa: E402  (this folder, v0.9.2 as published)

    yka = mc.load_yka_module("/tmp/icr_mods")
    if yka is None:
        s["yka_error"] = "yka.py download failed"
        return {}, []
    old = mc.read_events(os.path.join(mods, "z004_events.csv"))
    runs = [("old9", None)] + [(c[0], cand_dict(c)) for c in CANDIDATES]
    results, rows_all = {}, []
    for name, cand in runs:
        if time.time() > deadline:
            results[name] = dict(status="skipped_time_budget")
            continue
        if cand is not None:
            pn = os.path.join(data_dir, "yka_raw", cand["id"] + ".mseed.nodata")
            if os.path.exists(pn):
                os.remove(pn)
        evs = dict(old)
        if cand is not None:
            evs[cand["id"]] = cand
        sub = os.path.join(out_dir, "yka_runs", name)
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
        rows = read_rows(os.path.join(sub, "yka_z004_pairs.csv"))
        st = read_rows(os.path.join(sub, "yka_z004_status.csv"))
        y = sr.get("yka", {})
        info = y.get("multiplet", {}) if isinstance(y, dict) else {}
        res = dict(status="ran" if rows else "no_pairs", n_pairs=len(rows),
                   n_ok=sum(1 for r in rows if truthy(r.get("ok"))),
                   pick_minus_pred_s=info.get("pick_minus_pred_s"),
                   pick_valid=yka_repick.pick_valid(info.get("pick_minus_pred_s")),
                   cand_status=[r for r in st if cand is not None and r.get("id") == cand["id"]] or None)
        if cand is not None:
            res["state"] = state_test(rows, cand["id"], "D_c310", True)
            res["state_ungated"] = state_test(rows, cand["id"], "D_c310", False)
            res["D_c13"] = partner_values(rows, cand["id"], "D_c13", True)
            for r in rows:
                if cand["id"] in (r.get("ev1"), r.get("ev2")):
                    rr = dict(r)
                    rr["candidate"] = cand["id"]
                    rows_all.append(rr)
        else:
            res["rows"] = rows
        results[name] = res
        print("YKA", name, json.dumps(jsonable({k: v for k, v in res.items() if k != "rows"}))[:500], flush=True)
    return results, rows_all


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


def run_mantle_part(data_dir, out_dir, s, deadline, order):
    mp = load_mantlep2("/tmp/m34_pinned")
    if mp is None:
        s["mantle_error"] = "mantlep2.py v0.13.0 download failed"
        return {}, [], []
    base = list(mp.EVENTS)
    if [e[0] for e in base] != list(OLD8):
        s["mantle_error"] = "pinned mantlep2 EVENTS differ from OLD8"
        return {}, [], []
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
    s["mantle_stations_found"] = sorted(n for n in H13_STATIONS if n in coords)
    cache = os.path.join(data_dir, "mantlep")
    os.makedirs(cache, exist_ok=True)
    model = TauPyModel("ak135")
    out, ev_all, pair_all = {}, [], []
    for cid in order:
        c = [x for x in CANDIDATES if x[0] == cid][0]
        if time.time() > deadline:
            out[cid] = dict(decision="untestable", reason="skipped_time_budget")
            continue
        mp.EVENTS = base + [(cid, c[1].rstrip("Z"), c[2], c[3], c[4], c[5])]
        mp.IDS = [e[0] for e in mp.EVENTS]
        mp.DATES = {e[0]: e[1][:10] for e in mp.EVENTS}
        tmp_plots = os.path.join("/tmp/m34_plots", cid)
        os.makedirs(tmp_plots, exist_ok=True)
        per = {}
        for name in H13_STATIONS:
            if name not in coords:
                per[name] = dict(status="not_in_station_service")
                continue
            net, sta = name.split(".")
            ev_rows, pair_rows = [], []
            try:
                res = mp.run_station(session, model, net, sta, coords[name][0], coords[name][1], cache,
                                     tmp_plots, ev_rows, pair_rows, {})
            except Exception as ex:  # noqa: BLE001
                traceback.print_exc()
                per[name] = dict(status="exception", error=repr(ex)[:200])
                continue
            gated = {r["event"]: truthy(r.get("gated")) for r in ev_rows}
            D = {tuple(sorted((r["ev1"], r["ev2"]))): fnum(r["D_late"]) for r in pair_rows}
            sc = source_score(D, gated, list(OLD8), cid)
            sc.update(status=res.get("status"), ref=res.get("ref"))
            per[name] = sc
            for r in ev_rows:
                if r["event"] == cid:
                    r = dict(r)
                    r["candidate"] = cid
                    ev_all.append(r)
            for r in pair_rows:
                if cid in (r["ev1"], r["ev2"]):
                    r = dict(r)
                    r["candidate"] = cid
                    pair_all.append(r)
        qs = [v["q"] for v in per.values() if v.get("counted")]
        dec = decide_source(qs)
        out[cid] = dict(decision=dec, n_stations=len(qs), q=qs, median_q=median(qs), stations=per)
        if dec == "same":
            dst = os.path.join(out_dir, "mantle_plots", cid)
            try:
                shutil.copytree(tmp_plots, dst, dirs_exist_ok=True)
            except Exception:  # noqa: BLE001
                pass
        print("MANTLE", cid, dec, qs, flush=True)
    return out, ev_all, pair_all


# ----------------------------------------------------------------------------- grading
def grade(yres, mres):
    table, new = [], {}
    for c in CANDIDATES:
        cid = c[0]
        cd = cand_dict(c)
        y = yres.get(cid, {})
        st = y.get("state", {}) if isinstance(y, dict) else {}
        m = mres.get(cid, {})
        label = st.get("label", "untestable")
        dec = m.get("decision", "untestable")
        member = dec == "same"
        row = dict(id=cid, time=c[1], t=round(cd["t"], 5), equiv_wang=round(equiv(cd["t"]), 4),
                   leg="outbound" if cd["t"] < 2011.0 else "return", wang_prediction=wang_prediction(cd["t"]),
                   mag=c[5], depth=c[4], source=dec, n_src_stations=m.get("n_stations", 0),
                   median_q=m.get("median_q"), yka_status=y.get("status", "not_run"), nX=st.get("nX"),
                   nY=st.get("nY"), dX=st.get("dX"), dY=st.get("dY"), s=st.get("s"), state=label,
                   confirmed_member=member)
        row["used_in_mirror"] = bool(member and label in STATE_CODE)
        if row["used_in_mirror"]:
            new[cid] = (cd["t"], STATE_CODE[label])
            p = row["wang_prediction"]
            row["wang_check"] = "n/a" if p == "free" else ("agrees" if STATE_CODE[label] == p else "disagrees")
        table.append(row)
    g = mirror_tests(new)
    g["brackets"] = brackets(new)
    g["n_candidates"] = len(CANDIDATES)
    g["n_confirmed"] = sum(1 for r in table if r["confirmed_member"])
    g["n_used"] = len(new)
    g["n_different"] = sum(1 for r in table if r["source"] == "different")
    g["state_counts_confirmed"] = {}
    for r in table:
        if r["confirmed_member"]:
            g["state_counts_confirmed"][r["state"]] = g["state_counts_confirmed"].get(r["state"], 0) + 1
    return table, g


# ----------------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mods", required=True, help="folder with the pinned v0.8.0 modules (fetch_pinned.py)")
    ap.add_argument("--data-dir", default="/data")
    ap.add_argument("--out-dir", default="/out")
    ap.add_argument("--time-budget-min", type=float, default=200.0)
    args = ap.parse_args()
    t0 = time.time()
    total = 60.0 * args.time_budget_min
    os.makedirs(args.out_dir, exist_ok=True)
    mods = os.path.abspath(args.mods)
    s = dict(version=VERSION, started=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
             n_candidates=len(CANDIDATES), wang_boundary=wang_boundary(),
             thresholds=dict(x_near=X_NEAR, far=FAR, ratio=RATIO, between_min=BETWEEN_MIN, min_x=MIN_X,
                             min_y=MIN_Y, src_same_q=SRC_SAME_Q, src_diff_q=SRC_DIFF_Q,
                             min_old_gated=MIN_OLD_GATED, min_src_stations=MIN_SRC_STATIONS,
                             grid_turns=[GRID_TURNS[0], GRID_TURNS[-1]], grid_ratios=[GRID_RATIOS[0], GRID_RATIOS[-1]]))
    s["prereg_public_at_start"] = bool(fetch(PREREG_URL, tries=2, wait=5))

    yres, yrows = {}, []
    try:
        yres, yrows = run_yka_part(mods, args.data_dir, args.out_dir, s, t0 + 0.45 * total)
    except Exception as ex:  # noqa: BLE001
        traceback.print_exc()
        s["yka_error"] = repr(ex)[:300]
    pub = read_csv_text(fetch(YKA_V2_PAIRS, tries=3))
    if "old9" in yres and "rows" in yres["old9"]:
        s["R34_yka_old9"] = {k: max_abs_diff(yres["old9"]["rows"], pub, k) for k in ("D", "D_c310")}
        write_csv(os.path.join(args.out_dir, "m34_yka_old9_pairs.csv"), yres["old9"].pop("rows"))

    # pre-registered order: candidates with a testable YKA state first, then the rest
    def testable(cid):
        return yres.get(cid, {}).get("state", {}).get("label", "untestable") != "untestable"
    order = [c[0] for c in CANDIDATES if testable(c[0])] + [c[0] for c in CANDIDATES if not testable(c[0])]
    mres, mev, mpairs = {}, [], []
    try:
        mres, mev, mpairs = run_mantle_part(args.data_dir, args.out_dir, s, t0 + total, order)
    except Exception as ex:  # noqa: BLE001
        traceback.print_exc()
        s["mantle_error"] = repr(ex)[:300]

    table, g = grade(yres, mres)
    s["grade"] = g
    s["yka"] = yres
    s["mantle"] = mres
    s["finished"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    s["minutes"] = round((time.time() - t0) / 60.0, 1)
    write_csv(os.path.join(args.out_dir, "m34_candidates.csv"), table)
    write_csv(os.path.join(args.out_dir, "m34_yka_pairs.csv"), yrows)
    write_csv(os.path.join(args.out_dir, "m34_mantle_events.csv"), mev)
    write_csv(os.path.join(args.out_dir, "m34_mantle_pairs.csv"), mpairs)
    with open(os.path.join(args.out_dir, "m34_summary.json"), "w") as f:
        json.dump(jsonable(s), f, indent=1)
    print(json.dumps(jsonable(g))[:4000], flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())

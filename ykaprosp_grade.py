#!/usr/bin/env python3
"""ykaprosp_grade.py v0.1.0: the grader for pre-registration 12 (prospective YKA test).

Part of inner-core-repeaters by MV2, Munim, Inc. (made with AI agents). Licence: Apache-2.0.

Pre-registration: findings/2026-11-18-yka-prospective-prereg.md. It promised that this grader is
published before the interim look (on or after 2027-11-19) and before any YKA waveform of an
event from 2026-11-19T00:00Z onward is measured. It was written on 2026-11-19 (UTC), before any
such measurement. The choices the pre-registration left open are fixed here and listed in
findings/2026-11-19-prereg12-grader.md.

Input: the pairs table (yka_pairs.csv) that yka_b.py 0.10.1 --repick writes when run on the
extended South Sandwich Islands catalogue (families merged by ILAR correlation, unchanged
run.py threshold), plus a USGS ComCat CSV for event locations (fetched if not given).

What it grades, in order (all thresholds come from the pre-registration):
  selection  ok pair (direct cc >= 0.90, direct SNR >= 5, coda SNR >= 2, set by yka_b), ComCat
             epicentral separation <= 14.7 km, min(YKA direct SNR of the two events) >= 11.7.
  P4-prosp   selected pairs with a new event. near: partner on/after 2024-11-19. far: partner
             before 2016-11-19 and > 5 years from the new event on the backtracking clock.
             Testable if >= 5 near and >= 5 far; holds if median D_c13(near) < median D_c13(far).
             If untestable or failed, nothing downstream is interpreted.
  F1-prosp   a case is a new event with >= 2 selected partners before 2008.5. B pick: partner
             closest to e(t) = 2008.5 - (t - 2008.5)/2.5; M pick: latest pre-turn partner; best
             match: lowest D_c13. Testable with >= 6 discriminating cases in >= 3 families.
             supported: n_B > n_M, one-sided binomial p < 0.025, n_B > n_M in >= 95% of 2,000
             whole-family bootstraps. contradicted: the same with B and M swapped. Else
             inconclusive. At the interim look an untestable vote is not graded.
  extras     (ungraded) the same with D (+1 to +6 s), and with 2025-2026 partners removed.

Stdlib only.
"""
import argparse
import bisect
import csv
import io
import json
import math
import os
import random
import sys
import time
from datetime import datetime, timezone

VERSION = "0.1.0"
PREREG = "findings/2026-11-18-yka-prospective-prereg.md"
NAN = float("nan")
_EPOCH0 = datetime(1970, 1, 1, tzinfo=timezone.utc)

T_TURN, RATIO = 2008.5, 2.5
MAX_SEP_KM = 14.7
MIN_SNR = 11.7
FAR_CLOCK_YEARS = 5.0
MIN_NEAR = MIN_FAR = 5
MIN_CASES, MIN_FAMILIES = 6, 3
ALPHA = 0.025
BOOT_FRAC = 0.95
NBOOT = 2000
PRIMARY_KEY = "D_c13"

NEW_START_ISO = "2026-11-19T00:00:00Z"
NEAR_START_ISO = "2024-11-19T00:00:00Z"
FAR_BEFORE_ISO = "2016-11-19T00:00:00Z"
POOL_DROP_ISO = ("2025-01-01T00:00:00Z", NEW_START_ISO)
LOOKS = {"interim": "2027-11-19T00:00:00Z", "final": "2028-11-19T00:00:00Z"}

USGS = ("https://earthquake.usgs.gov/fdsnws/event/1/query?format=csv&starttime=1990-01-01"
        "&endtime=%s&minlatitude=-61.5&maxlatitude=-54&minlongitude=-31&maxlongitude=-22"
        "&minmagnitude=4.3&orderby=time-asc&limit=20000")


# ----------------------------------------------------------------------------- time
def epoch(s):
    """Seconds since 1970 for an ISO-8601 UTC time such as 2001-05-29T06:22:55.100Z."""
    s = str(s).strip().replace(" ", "T")
    if s.endswith("Z"):
        s = s[:-1]
    elif len(s) > 19 and s[-6] in "+-" and s[-3] == ":":
        if s[-6:] not in ("+00:00", "-00:00"):
            raise ValueError("not UTC: " + s)
        s = s[:-6]
    date, _, tm = s.partition("T")
    y, mo, d = (int(x) for x in date.split("-"))
    hh = mm = 0
    ss = 0.0
    if tm:
        p = tm.split(":")
        hh = int(p[0])
        mm = int(p[1]) if len(p) > 1 else 0
        ss = float(p[2]) if len(p) > 2 else 0.0
    return (datetime(y, mo, d, hh, mm, tzinfo=timezone.utc) - _EPOCH0).total_seconds() + ss


def iso(ep):
    return datetime.fromtimestamp(ep, timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def decyear(ep):
    """Decimal year (same convention as yka_b.py: fraction of the calendar year elapsed)."""
    y = datetime.fromtimestamp(ep, timezone.utc).year
    a = (datetime(y, 1, 1, tzinfo=timezone.utc) - _EPOCH0).total_seconds()
    b = (datetime(y + 1, 1, 1, tzinfo=timezone.utc) - _EPOCH0).total_seconds()
    return y + (ep - a) / (b - a)


def equiv(t):
    """Backtracking clock: the pre-turn date a post-turn date maps to (Wang et al. 2024)."""
    t = float(t)
    return t if t <= T_TURN else T_TURN - (t - T_TURN) / RATIO


NEW_START = epoch(NEW_START_ISO)
NEAR_START = epoch(NEAR_START_ISO)
FAR_BEFORE = epoch(FAR_BEFORE_ISO)
POOL_DROP = (epoch(POOL_DROP_ISO[0]), epoch(POOL_DROP_ISO[1]))


def look_open(look, now_ep):
    return now_ep >= epoch(LOOKS[look])


# ----------------------------------------------------------------------------- helpers
def fnum(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return NAN


def finite(x):
    return isinstance(x, (int, float)) and math.isfinite(x)


def median(xs):
    xs = sorted(x for x in xs if finite(x))
    n = len(xs)
    if n == 0:
        return NAN
    return xs[n // 2] if n % 2 else 0.5 * (xs[n // 2 - 1] + xs[n // 2])


def sep_km(lat1, lon1, lat2, lon2):
    p1, p2 = math.radians(lat1), math.radians(lat2)
    dp, dl = p2 - p1, math.radians(lon2 - lon1)
    h = math.sin(dp / 2) ** 2 + math.cos(p1) * math.cos(p2) * math.sin(dl / 2) ** 2
    return 2 * 6371.0 * math.asin(min(1.0, math.sqrt(h)))


def binom_tail(k, n):
    """One-sided P(X >= k), X ~ Binomial(n, 1/2)."""
    if n <= 0:
        return 1.0
    return sum(math.comb(n, j) for j in range(k, n + 1)) / 2.0 ** n


# ----------------------------------------------------------------------------- tables
def read_catalog(text):
    by_id, items = {}, []
    for r in csv.DictReader(io.StringIO(text)):
        try:
            e = dict(lat=float(r["latitude"]), lon=float(r["longitude"]), epoch=epoch(r["time"]))
        except (KeyError, ValueError, TypeError):
            continue
        by_id[r.get("id", "")] = e
        items.append(e)
    items.sort(key=lambda e: e["epoch"])
    return dict(by_id=by_id, items=items, keys=[e["epoch"] for e in items])


def locate(ev_id, when, cat, tol=5.0):
    """ComCat entry by event id, else by origin time within tol seconds."""
    if ev_id in cat["by_id"]:
        return cat["by_id"][ev_id]
    try:
        t = epoch(when)
    except (ValueError, TypeError):
        return None
    i = bisect.bisect_left(cat["keys"], t)
    best = None
    for j in (i - 1, i):
        if 0 <= j < len(cat["items"]):
            d = abs(cat["keys"][j] - t)
            if d <= tol and (best is None or d < best[0]):
                best = (d, cat["items"][j])
    return best[1] if best else None


def read_pairs(path):
    rows = []
    with open(path, newline="") as f:
        for r in csv.DictReader(f):
            rr = dict(r)
            for k in ("t1", "t2", "D", "D_c13", "D_c310", "snr_dir1", "snr_dir2"):
                rr[k] = fnum(r.get(k))
            rr["ok"] = str(r.get("ok")).strip() == "True"
            rows.append(rr)
    return rows


def orient(r):
    """A pair that includes an event from 2026-11-19 on, as (new event, partner); else None.
    If both events are new, the later one is the new event and the earlier one the partner."""
    try:
        e1, e2 = epoch(r["time1"]), epoch(r["time2"])
    except (KeyError, ValueError, TypeError):
        return None
    n1, n2 = e1 >= NEW_START, e2 >= NEW_START
    if not (n1 or n2):
        return None
    if n1 and n2:
        new, part = (2, 1) if e2 >= e1 else (1, 2)
    else:
        new, part = (2, 1) if n2 else (1, 2)
    ep = {1: e1, 2: e2}
    t = {}
    for i in (1, 2):
        v = r.get("t%d" % i)
        t[i] = v if finite(v) else decyear(ep[i])
    s1, s2 = r.get("snr_dir1"), r.get("snr_dir2")
    minsnr = min(s1, s2) if finite(s1) and finite(s2) else NAN
    return dict(mult=r.get("mult", ""), new_ev=r.get("ev%d" % new, ""), partner=r.get("ev%d" % part, ""),
                new_time=r["time%d" % new], partner_time=r["time%d" % part], new_ep=ep[new],
                partner_ep=ep[part], t_new=t[new], t_partner=t[part], both_new=n1 and n2, ok=bool(r.get("ok")),
                minsnr=minsnr, D=r.get("D", NAN), D_c13=r.get("D_c13", NAN), D_c310=r.get("D_c310", NAN))


def prospective_pairs(rows, cat):
    out = []
    for r in rows:
        p = orient(r)
        if p is None:
            continue
        a = locate(p["new_ev"], p["new_time"], cat)
        b = locate(p["partner"], p["partner_time"], cat)
        p["sep_km"] = sep_km(a["lat"], a["lon"], b["lat"], b["lon"]) if a and b else NAN
        out.append(p)
    return out


def exclusion(p, key):
    """Why a prospective pair is not selected ('' if it is)."""
    if not p["ok"]:
        return "not_ok"
    if not finite(p.get(key)):
        return "no_" + key
    if not finite(p["sep_km"]):
        return "no_location"
    if p["sep_km"] > MAX_SEP_KM:
        return "separation"
    if not (finite(p["minsnr"]) and p["minsnr"] >= MIN_SNR):
        return "snr"
    return ""


# ----------------------------------------------------------------------------- tests
def is_near(p):
    return p["partner_ep"] >= NEAR_START


def is_far(p):
    return p["partner_ep"] < FAR_BEFORE and abs(equiv(p["t_new"]) - equiv(p["t_partner"])) > FAR_CLOCK_YEARS


def p4(sel, key=PRIMARY_KEY):
    near = [p[key] for p in sel if is_near(p)]
    far = [p[key] for p in sel if is_far(p)]
    mn, mf = median(near), median(far)
    testable = len(near) >= MIN_NEAR and len(far) >= MIN_FAR
    holds = bool(testable and mn < mf)
    return dict(n_near=len(near), n_far=len(far), median_near=mn, median_far=mf, testable=testable,
                holds=holds, verdict=("holds" if holds else "fails") if testable else "untestable")


def cases(sel, key=PRIMARY_KEY):
    by = {}
    for p in sel:
        if p["t_partner"] < T_TURN:
            by.setdefault((p["mult"], p["new_ev"]), []).append(p)
    out = []
    for (mult, ev), ps in sorted(by.items()):
        if len(ps) < 2:
            continue
        e = equiv(ps[0]["t_new"])
        best = min(ps, key=lambda p: (p[key], p["t_partner"]))
        b = min(ps, key=lambda p: (abs(p["t_partner"] - e), p["t_partner"]))
        m = max(ps, key=lambda p: p["t_partner"])
        out.append(dict(
            mult=mult, new_ev=ev, t_new=round(ps[0]["t_new"], 3), e_t_new=round(e, 3), n_partners=len(ps),
            partners=";".join("%s@%.2f:%s=%.3f" % (p["partner"], p["t_partner"], key, p[key])
                              for p in sorted(ps, key=lambda p: p["t_partner"])),
            best=best["partner"], best_t=round(best["t_partner"], 3), best_value=round(best[key], 4),
            B_pick=b["partner"], B_t=round(b["t_partner"], 3), M_pick=m["partner"], M_t=round(m["t_partner"], 3),
            discriminating=b["partner"] != m["partner"], matches_B=best["partner"] == b["partner"],
            matches_M=best["partner"] == m["partner"], chance=1.0 / len(ps)))
    return out


def grade_vote(cs, nboot=NBOOT, seed=0, alpha=ALPHA, boot_frac=BOOT_FRAC):
    dc = [c for c in cs if c["discriminating"]]
    nB = sum(1 for c in dc if c["matches_B"])
    nM = sum(1 for c in dc if c["matches_M"])
    fams = sorted({c["mult"] for c in dc})
    testable = len(dc) >= MIN_CASES and len(fams) >= MIN_FAMILIES
    fb = fm = NAN
    if len(fams) >= 2:
        rng = random.Random(seed)
        groups = {f: [c for c in dc if c["mult"] == f] for f in fams}
        b_ahead = m_ahead = 0
        for _ in range(nboot):
            pick = [groups[fams[rng.randrange(len(fams))]] for _ in fams]
            xb = sum(1 for g in pick for c in g if c["matches_B"])
            xm = sum(1 for g in pick for c in g if c["matches_M"])
            b_ahead += xb > xm
            m_ahead += xm > xb
        fb, fm = b_ahead / nboot, m_ahead / nboot
    pB, pM = binom_tail(nB, nB + nM), binom_tail(nM, nB + nM)
    verdict = "untestable"
    if testable:
        if nB > nM and pB < alpha and fb >= boot_frac:
            verdict = "supported"
        elif nM > nB and pM < alpha and fm >= boot_frac:
            verdict = "contradicted"
        else:
            verdict = "inconclusive"
    return dict(n_cases=len(cs), n_discriminating=len(dc), n_families=len(fams), n_match_B=nB, n_match_M=nM,
                n_neither=len(dc) - nB - nM, chance_expected_B=sum(c["chance"] for c in dc), p_binom_B=pB,
                p_binom_M=pM, boot_frac_B_ahead=fb, boot_frac_M_ahead=fm, nboot=nboot, alpha=alpha,
                testable=testable, verdict=verdict)


def decide(p4r, vote, look):
    if not p4r["testable"]:
        return "not interpreted (P4-prosp untestable)"
    if not p4r["holds"]:
        return "not interpreted (P4-prosp fails)"
    if not vote["testable"]:
        return "not graded at interim (vote untestable)" if look == "interim" else "untestable"
    return vote["verdict"]


def grade_key(pp, key, look, nboot=NBOOT, drop_pool=False):
    pool = [p for p in pp if not (drop_pool and POOL_DROP[0] <= p["partner_ep"] < POOL_DROP[1])]
    reasons, sel = {}, []
    for p in pool:
        why = exclusion(p, key)
        if why:
            reasons[why] = reasons.get(why, 0) + 1
        else:
            sel.append(p)
    p4r = p4(sel, key)
    cs = cases(sel, key)
    vote = grade_vote(cs, nboot=nboot)
    res = dict(key=key, drop_2025_2026_partners=drop_pool, n_prospective=len(pool), n_selected=len(sel),
               n_new_events=len({p["new_ev"] for p in sel}), excluded=reasons, P4=p4r, F1=vote,
               decision=decide(p4r, vote, look))
    return res, sel, cs


def grade_all(rows, cat, look, nboot=NBOOT):
    pp = prospective_pairs(rows, cat)
    prim, sel, cs = grade_key(pp, PRIMARY_KEY, look, nboot)
    prim["graded"] = True
    ex_d, _, _ = grade_key(pp, "D", look, nboot)
    ex_pool, _, _ = grade_key(pp, PRIMARY_KEY, look, nboot, drop_pool=True)
    ex_d["graded"] = ex_pool["graded"] = False
    return dict(primary=prim, extras=dict(D=ex_d, pool_drop=ex_pool), selected=sel, cases=cs)


# ----------------------------------------------------------------------------- io
def clean(o):
    if isinstance(o, dict):
        return {str(k): clean(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [clean(v) for v in o]
    if isinstance(o, float) and not math.isfinite(o):
        return None
    return o


def fetch_catalog(now_ep, tries=8, wait_s=60):
    import urllib.request

    url = USGS % datetime.fromtimestamp(now_ep + 86400, timezone.utc).strftime("%Y-%m-%d")
    for _ in range(tries):
        try:
            with urllib.request.urlopen(url, timeout=120) as r:
                text = r.read().decode("utf-8")
            if "latitude" in text:
                return text
        except Exception:
            pass
        time.sleep(wait_s)
    raise RuntimeError("could not fetch the ComCat catalogue")


def write_csv(path, rows, cols):
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(cols)
        for r in rows:
            w.writerow([r.get(c) for c in cols])


def main(argv=None):
    ap = argparse.ArgumentParser(description="Grade pre-registration 12 (prospective YKA test).")
    ap.add_argument("--pairs", required=True, help="yka_pairs.csv from yka_b.py 0.10.1 --repick")
    ap.add_argument("--catalog", default="", help="ComCat CSV (fetched from USGS if empty)")
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--look", choices=sorted(LOOKS), required=True)
    ap.add_argument("--now", default="", help="ISO time to use as 'now' (default: the clock)")
    ap.add_argument("--nboot", type=int, default=NBOOT)
    a = ap.parse_args(argv)
    os.makedirs(a.out_dir, exist_ok=True)
    now = epoch(a.now) if a.now else time.time()
    head = dict(version=VERSION, prereg=PREREG, look=a.look, look_opens=LOOKS[a.look], now=iso(now))
    if not look_open(a.look, now):
        head["status"] = "too-early: the %s look opens %s; nothing graded" % (a.look, LOOKS[a.look])
        with open(os.path.join(a.out_dir, "ykaprosp_grade.json"), "w") as f:
            json.dump(head, f, indent=1)
        print(head["status"])
        return 2
    if a.catalog:
        with open(a.catalog, encoding="utf-8") as f:
            text = f.read()
    else:
        text = fetch_catalog(now)
        with open(os.path.join(a.out_dir, "comcat_ssi.csv"), "w", encoding="utf-8", newline="") as f:
            f.write(text)
    res = grade_all(read_pairs(a.pairs), read_catalog(text), a.look, nboot=a.nboot)
    for p in res["selected"]:
        p["near"], p["far"] = is_near(p), is_far(p)
    write_csv(os.path.join(a.out_dir, "ykaprosp_pairs.csv"), res["selected"],
              ["mult", "new_ev", "partner", "new_time", "partner_time", "t_new", "t_partner", "both_new", "sep_km",
               "minsnr", "D", "D_c13", "D_c310", "near", "far"])
    cols = ["mult", "new_ev", "t_new", "e_t_new", "n_partners", "partners", "best", "best_t", "best_value",
            "B_pick", "B_t", "M_pick", "M_t", "discriminating", "matches_B", "matches_M", "chance"]
    write_csv(os.path.join(a.out_dir, "ykaprosp_cases.csv"), res["cases"], cols)
    head["status"] = "graded"
    head.update(primary=res["primary"], extras=res["extras"])
    with open(os.path.join(a.out_dir, "ykaprosp_grade.json"), "w") as f:
        json.dump(clean(head), f, indent=1)
    pr = res["primary"]
    print("prereg 12 %s look: P4-prosp %s (near %d, far %d); vote %d discriminating cases in %d families; "
          "decision: %s" % (a.look, pr["P4"]["verdict"], pr["P4"]["n_near"], pr["P4"]["n_far"],
                            pr["F1"]["n_discriminating"], pr["F1"]["n_families"], pr["decision"]))
    return 0


if __name__ == "__main__":
    sys.exit(main())

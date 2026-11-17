#!/usr/bin/env python3
"""yka_b.py: YKA stage B on the merged ILAR families (pre-registration 8).

ILAR's outer-core phases (PKPbc+PKPab) define the repeating-earthquake families; the
Yellowknife array (YKA), where PKIKP is strong, measures the inner core. This script

  1. fetches yka.py v0.2.0 and its helpers (model, phases, xcorr) and the published ILAR
     family tables from the public repo at a pinned commit;
  2. merges the ilar-v1, m0324-v1 and fdscan-v1 families (connected components of the
     union of their repeater graphs);
  3. runs yka.py stage B unchanged on the merged families (YKA catalogue down to M4.3,
     only family members fetched);
  4. grades the new-event vote F1-Y / F1x-Y of findings/2026-11-17-yka-stageb-prereg.md.

By MV2, Munim, Inc. (made with AI agents). Licence: Apache-2.0.
"""
import argparse
import csv
import json
import math
import os
import sys
import time

VERSION = "0.9.0"
COMMIT = "7067d4fc2b0a2fdb1378a1c406dd680c8d0b69b1"  # v0.8.0 of munimtechnologies/inner-core-repeaters
RAW = "https://raw.githubusercontent.com/munimtechnologies/inner-core-repeaters/%s/" % COMMIT
MODULES = ("model.py", "phases.py", "xcorr.py", "yka.py")
REMOTE_TABLES = {
    "ilar_v1_multiplets.csv": "results/ilar-v1/multiplets.csv",
    "ilar_v1_pairs.csv": "results/ilar-v1/pairs.csv",
    "m0324_multiplets.csv": "results/m0324-v1/m0324_multiplets.csv",
    "m0324_pairs.csv": "results/m0324-v1/m0324_pairs.csv",
}
# (prefix, multiplets file, pairs file). fd_multiplets.csv is icr-fdscan-v1/out/fd_multiplets.csv
# unchanged; fd_pairs_ok.csv holds the 16 ok rows of its fd_pairs.csv (mult, ev1, ev2, D, ok).
SOURCES = (("I", "ilar_v1_multiplets.csv", "ilar_v1_pairs.csv"),
           ("Z", "m0324_multiplets.csv", "m0324_pairs.csv"),
           ("F", "fd_multiplets.csv", "fd_pairs_ok.csv"))

T_TURN, RATIO, NEW_CUT = 2008.5, 2.5, 2025.0
ANOMALY_EVENT = "us6000runq"  # 2025-12-13, pre-registration 7
T0 = time.time()


def log(*a):
    print("[%7.1f min]" % ((time.time() - T0) / 60.0), *a, flush=True)


def equiv(t):
    t = float(t)
    return t if t <= T_TURN else T_TURN - (t - T_TURN) / RATIO


def parse_members(s):
    out = []
    for part in str(s or "").split(";"):
        part = part.strip()
        if part:
            out.append(part.split("(")[0].strip())
    return [p for p in out if p]


def member_dates(s):
    """{'id': '(YYYY-MM-DD)'} so merged rows keep the readable dates."""
    d = {}
    for part in str(s or "").split(";"):
        part = part.strip()
        if part:
            i = part.split("(")[0].strip()
            d[i] = part[len(i):]
    return d


# ----------------------------------------------------------------------------- downloads
def fetch(url, path, must_contain=None, tries=6):
    import requests

    for attempt in range(tries):
        try:
            r = requests.get(url, timeout=90)
            if r.status_code == 200 and (must_contain is None or must_contain in r.text):
                with open(path, "w", encoding="utf-8", newline="") as f:
                    f.write(r.text)
                return True
            log("fetch", url, "status", r.status_code)
        except Exception as ex:
            log("fetch error", url, repr(ex)[:200])
        time.sleep(10 * (attempt + 1))
    return False


def get_inputs(mod_dir, tab_dir):
    os.makedirs(mod_dir, exist_ok=True)
    for m in MODULES:
        p = os.path.join(mod_dir, m)
        if not os.path.exists(p):
            need = 'VERSION = "0.2.0"' if m == "yka.py" else None
            if not fetch(RAW + m, p, must_contain=need):
                raise SystemExit("could not fetch %s" % m)
    for local, remote in REMOTE_TABLES.items():
        p = os.path.join(tab_dir, local)
        if not os.path.exists(p):
            if not fetch(RAW + remote, p, must_contain="mult,"):
                raise SystemExit("could not fetch %s" % remote)


# ----------------------------------------------------------------------------- merge
def read_families(path, prefix):
    fams = []
    with open(path, newline="") as f:
        for r in csv.DictReader(f):
            ids = parse_members(r.get("members"))
            if len(ids) >= 2:
                fams.append(dict(name=prefix + r["mult"], ids=ids, dates=member_dates(r.get("members"))))
    return fams


def merge_families(fams):
    """Connected components of the union of families (families sharing an event merge)."""
    parent = {}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[max(ra, rb)] = min(ra, rb)

    dates = {}
    for fam in fams:
        for i in fam["ids"]:
            parent.setdefault(i, i)
            dates.setdefault(i, fam["dates"].get(i, ""))
        for i in fam["ids"][1:]:
            union(fam["ids"][0], i)
    comps = {}
    for i in parent:
        comps.setdefault(find(i), set()).add(i)
    srcs = {}
    for fam in fams:
        srcs.setdefault(find(fam["ids"][0]), []).append(fam["name"])
    out = []
    for root, ids in comps.items():
        ids = sorted(ids, key=lambda i: (dates.get(i, ""), i))
        out.append(dict(ids=ids, sources=sorted(srcs.get(root, [])), dates=dates))
    out.sort(key=lambda c: (-len(c["ids"]), c["ids"][0]))
    for k, c in enumerate(out):
        c["mult"] = "Y%03d" % k
    return out


def write_merged(path, merged):
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["mult", "n", "sources", "members"])
        for c in merged:
            w.writerow([c["mult"], len(c["ids"]), "+".join(c["sources"]),
                        ";".join(i + c["dates"].get(i, "") for i in c["ids"])])


def merge_ok_pairs(paths, out_path):
    """ILAR ok pairs (ev1 earlier) for the cross-array comparison; first source wins."""
    seen, rows = set(), []
    for p in paths:
        if not os.path.exists(p):
            continue
        with open(p, newline="") as f:
            for r in csv.DictReader(f):
                if str(r.get("ok")) != "True":
                    continue
                try:
                    d = float(r["D"])
                except (TypeError, ValueError):
                    continue
                k = (r["ev1"], r["ev2"])
                if k in seen or not math.isfinite(d):
                    continue
                seen.add(k)
                rows.append(dict(ev1=r["ev1"], ev2=r["ev2"], D=d, ok="True"))
    with open(out_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["ev1", "ev2", "D", "ok"])
        w.writeheader()
        for r in rows:
            w.writerow(r)
    return len(rows)


# ----------------------------------------------------------------------------- the vote
def fnum(x):
    try:
        v = float(x)
    except (TypeError, ValueError):
        return float("nan")
    return v


def read_pairs(path):
    rows = []
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        return rows
    with open(path, newline="") as f:
        for r in csv.DictReader(f):
            rr = dict(r)
            for k in ("t1", "t2", "D", "Dx", "D_c13", "D_c310", "D_full", "cc_ref", "snr_dir1", "snr_dir2",
                      "snr_coda1", "snr_coda2"):
                rr[k] = fnum(r.get(k))
            rr["ok"] = str(r.get("ok")) == "True"
            rows.append(rr)
    return rows


def cases_from_pairs(rows, key="D"):
    by = {}
    for r in rows:
        if not r["ok"] or not math.isfinite(r.get(key, float("nan"))):
            continue
        if r["t2"] >= NEW_CUT and r["t1"] < T_TURN:
            by.setdefault((r["mult"], r["ev2"]), []).append(r)
    cases = []
    for (mult, ev), rs in sorted(by.items()):
        if len(rs) < 2:
            continue
        t2 = rs[0]["t2"]
        e2 = equiv(t2)
        best = min(rs, key=lambda r: r[key])
        b = min(rs, key=lambda r: (abs(r["t1"] - e2), r["t1"]))
        m = max(rs, key=lambda r: r["t1"])
        top = max(rs, key=lambda r: r["snr_dir1"] if math.isfinite(r["snr_dir1"]) else -1.0)
        cases.append(dict(
            mult=mult, ev2=ev, t2=round(t2, 3), e_t2=round(e2, 3), n_partners=len(rs),
            partners=";".join("%s@%.2f:D=%.3f" % (r["ev1"], r["t1"], r[key]) for r in sorted(rs, key=lambda r: r["t1"])),
            best_ev1=best["ev1"], best_t1=round(best["t1"], 3), best_D=round(best[key], 4),
            B_pick=b["ev1"], B_t1=round(b["t1"], 3), M_pick=m["ev1"], M_t1=round(m["t1"], 3),
            discriminating=b["ev1"] != m["ev1"], matches_B=best["ev1"] == b["ev1"],
            matches_M=best["ev1"] == m["ev1"], best_is_top_snr=best["ev1"] == top["ev1"],
            chance=1.0 / len(rs)))
    return cases


def binom_tail(k, n):
    """One-sided P(X >= k), X ~ Binomial(n, 1/2)."""
    if n <= 0:
        return 1.0
    return sum(math.comb(n, j) for j in range(k, n + 1)) / 2.0 ** n


def grade_vote(cases, nboot=2000, seed=0, min_cases=6, min_mults=3):
    dc = [c for c in cases if c["discriminating"]]
    nB = sum(1 for c in dc if c["matches_B"])
    nM = sum(1 for c in dc if c["matches_M"])
    mults = sorted({c["mult"] for c in dc})
    testable = len(dc) >= min_cases and len(mults) >= min_mults
    fb = fm = float("nan")
    if len(mults) >= 2:
        import random

        rng = random.Random(seed)
        groups = {m: [c for c in dc if c["mult"] == m] for m in mults}
        b_ahead = m_ahead = 0
        for _ in range(nboot):
            pick = [groups[mults[rng.randrange(len(mults))]] for _ in mults]
            xb = sum(1 for g in pick for c in g if c["matches_B"])
            xm = sum(1 for g in pick for c in g if c["matches_M"])
            b_ahead += xb > xm
            m_ahead += xm > xb
        fb, fm = b_ahead / nboot, m_ahead / nboot
    pB, pM = binom_tail(nB, nB + nM), binom_tail(nM, nB + nM)
    verdict = "untestable"
    if testable:
        if nB > nM and pB < 0.05 and fb >= 0.95:
            verdict = "supported"
        elif nM > nB and pM < 0.05 and fm >= 0.95:
            verdict = "contradicted"
        else:
            verdict = "inconclusive"
    return dict(n_cases=len(cases), n_discriminating=len(dc), n_multiplets=len(mults), n_match_B=nB,
                n_match_M=nM, n_neither=len(dc) - nB - nM,
                chance_expected_B=sum(c["chance"] for c in dc), p_binom_B=pB, p_binom_M=pM,
                boot_frac_B_ahead=fb, boot_frac_M_ahead=fm,
                frac_best_is_top_snr=(sum(c["best_is_top_snr"] for c in cases) / len(cases)) if cases else None,
                testable=testable, verdict=verdict)


def clean(o):
    if isinstance(o, dict):
        return {str(k): clean(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [clean(v) for v in o]
    if isinstance(o, float) and not math.isfinite(o):
        return None
    return o


def grade_all(pairs_path, nboot=2000):
    rows = read_pairs(pairs_path)
    res = dict(n_rows=len(rows), n_ok=sum(r["ok"] for r in rows),
               n_ok_new=sum(1 for r in rows if r["ok"] and r["t2"] >= NEW_CUT),
               n_ok_new_preturn=sum(1 for r in rows if r["ok"] and r["t2"] >= NEW_CUT and r["t1"] < T_TURN))
    cases = cases_from_pairs(rows, "D")
    res["F1_Y"] = grade_vote(cases, nboot=nboot)
    casesx = cases_from_pairs(rows, "Dx")
    res["F1x_Y"] = grade_vote(casesx, nboot=nboot)
    f1, f1x = res["F1_Y"]["verdict"], res["F1x_Y"]["verdict"]
    opposite = {("supported", "contradicted"), ("contradicted", "supported")}
    res["F1_Y_robust"] = f1 in ("supported", "contradicted") and (f1, f1x) not in opposite
    res["sensitivity"] = {k: grade_vote(cases_from_pairs(rows, k), nboot=300) for k in ("D_c13", "D_c310", "D_full")}
    res["anomaly_event_pairs"] = [
        dict(ev1=r["ev1"], ev2=r["ev2"], t1=round(r["t1"], 3), t2=round(r["t2"], 3), D=r["D"], Dx=r["Dx"],
             cc_direct=r["cc_ref"], snr_dir1=r["snr_dir1"], snr_dir2=r["snr_dir2"], ok=r["ok"])
        for r in rows if ANOMALY_EVENT in (r["ev1"], r["ev2"])]
    return res, cases


def write_cases(path, cases):
    keys = ["mult", "ev2", "t2", "e_t2", "n_partners", "partners", "best_ev1", "best_t1", "best_D", "B_pick", "B_t1",
            "M_pick", "M_t1", "discriminating", "matches_B", "matches_M", "best_is_top_snr", "chance"]
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        for c in cases:
            w.writerow(c)


# ----------------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", default="/data")
    ap.add_argument("--out-dir", default="/out")
    ap.add_argument("--tables-dir", default=".")
    ap.add_argument("--nboot", type=int, default=2000)
    ap.add_argument("--time-budget-min", type=float, default=150.0)
    ap.add_argument("--grade-only", default="", help="grade an existing yka_pairs.csv and stop")
    args = ap.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)
    if args.grade_only:
        res, cases = grade_all(args.grade_only, nboot=args.nboot)
        print(json.dumps(clean(res), indent=1))
        return
    mod_dir = os.path.join(args.data_dir, "mods_" + COMMIT[:7])
    fetched = os.path.join(args.out_dir, "inputs")  # kept with the outputs for reproducibility
    os.makedirs(fetched, exist_ok=True)
    get_inputs(mod_dir, fetched)

    def locate(name):
        p = os.path.join(args.tables_dir, name)
        return p if os.path.exists(p) else os.path.join(fetched, name)

    fams = []
    for prefix, mp, _ in SOURCES:
        fams += read_families(locate(mp), prefix)
    merged = merge_families(fams)
    mpath = os.path.join(args.out_dir, "merged_multiplets.csv")
    write_merged(mpath, merged)
    ppath = os.path.join(args.out_dir, "merged_ilar_ok_pairs.csv")
    n_il = merge_ok_pairs([locate(p) for _, _, p in SOURCES], ppath)
    members = {i for c in merged for i in c["ids"]}
    with_new = [c for c in merged if any(c["dates"].get(i, "")[1:5] in ("2025", "2026") for i in c["ids"])]
    log("families read", len(fams), "merged", len(merged), "members", len(members), "with 2025+ member",
        len(with_new), "ILAR ok pairs", n_il)

    sys.path.insert(0, mod_dir)
    import yka  # noqa: E402  (v0.2.0, pinned commit)

    sys.argv = ["yka.py", "--data-dir", args.data_dir, "--out-dir", args.out_dir,
                "--min-mag-old", "4.3", "--min-mag-new", "4.3", "--new-from", "2024-07-01",
                "--ilar-multiplets", mpath, "--ilar-pairs", ppath, "--max-events", str(len(members)),
                "--time-budget-min", str(args.time_budget_min), "--fetch-frac", "0.8", "--plots", "30",
                "--quicklooks", "6"]
    yka.main()

    res, cases = grade_all(os.path.join(args.out_dir, "yka_pairs.csv"), nboot=args.nboot)
    res.update(version=VERSION, commit=COMMIT, n_families_read=len(fams), n_merged=len(merged),
               n_members=len(members), n_merged_with_new=len(with_new), n_ilar_ok_pairs=n_il)
    try:
        with open(os.path.join(args.out_dir, "yka_predictions.json")) as f:
            res["yka_predictions"] = json.load(f)
    except Exception as ex:
        res["yka_predictions"] = dict(error=repr(ex))
    write_cases(os.path.join(args.out_dir, "yka_b_cases.csv"), cases)
    with open(os.path.join(args.out_dir, "yka_b_summary.json"), "w") as f:
        json.dump(clean(res), f, indent=1)
    log("F1-Y:", json.dumps(clean(res["F1_Y"])))
    log("F1x-Y:", json.dumps(clean(res["F1x_Y"])))
    log("done")


if __name__ == "__main__":
    main()

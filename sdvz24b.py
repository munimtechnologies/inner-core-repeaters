#!/usr/bin/env python3
# Copyright 2026 MV2, Munim, Inc.
# Licensed under the Apache License, Version 2.0 (Apache-2.0). See LICENSE.
# Made with AI agents working in MV2's name.
"""sdvz24b.py - pre-registration 24, amendment 1: run sdvz24.py 0.1.0 (pinned commit) with
the processing sets grouped as the plan says ("same file name apart from the component
code") and paired like with like.

Why: in job sdv-zframe24-v1, sdvz24.set_key 0.1.0 fell back to folder|rate|start for the
authors' file names (2009.R.sac, 2009.R.sac.bp ...). For 2009 both files are 100 Hz, so the
plain and the band-passed ('.bp') files fell into one set and the '.bp' traces overwrote the
plain ones; for 2000 the '.bp' files are 1000 Hz and stayed separate. The graded primary pair
was therefore plain 2000 against band-passed 2009.

What this does (plan: findings/2026-11-22-prereg24-amendment-1.md):
 1. PRIMARY: sdvz24.run with set_key_v2 (component token replaced by '?'), so the primary pair
    is the authors' plain traces 2000.?.sac / 2009.?.sac. Grades exactly as in the plan.
 2. SECONDARY: the same with the '.bp' sets sorted first, i.e. the authors' band-passed pair.
 3. Combined grade with pre-registration 23's long-period delta (+0.55 deg, its gate G5 passed),
    as the plan prescribes (0.1.0 only used the Rayleigh-Z delta, which was unavailable).
 4. REPORTED: how the authors' '.bp' traces relate to their plain traces in each event (lag
    and signed correlation after the common 0.4-2 Hz zero-phase band-pass), per window.

Top-level imports are stdlib + numpy; the pinned sdv23.py / sdvz24.py are fetched at run time.
"""

import argparse
import hashlib
import importlib
import json
import math
import os
import sys
import tempfile
import time
import types
import urllib.request

import numpy as np

VERSION = "0.1.0"
REPO_RAW = "https://raw.githubusercontent.com/munimtechnologies/inner-core-repeaters/"
CODE_COMMIT = "dbbd91d3f26f9dc2165cbf59a5c05f7d82d8521d"   # release v0.21.0
MUST = {"sdv23.py": b'VERSION = "0.1.1"', "sdvz24.py": b'VERSION = "0.1.0"'}
AMEND_PATH = "findings/2026-11-22-prereg24-amendment-1.md"
P23_DELTA = 0.5500000000004341   # results/sdv23-v2/sdv23_summary.json, top-level "delta"
P23_G5 = True                    # findings/2026-11-21-prereg23-results.md, gate G5 passed
COMP_TOKENS = ("R", "T", "Z", "N", "E", "BHR", "BHT", "BHZ", "BHN", "BHE", "HHR", "HHT")
EXPECTED_SUFFIX = {"primary": ".sac", "secondary": ".sac.bp"}
DIAG = dict(sr=200.0, maxlag_s=2.0, sks_win=(-5.0, 10.0), s_win=(-5.0, 15.0))


def log(*a):
    print(time.strftime("%H:%M:%S"), "[sdvz24b]", *a, flush=True)


# ----------------------------------------------------------------- pure helpers

def set_key_v2(rel_path, cha, sr, start_iso):
    """Same file name apart from the component code. A channel code inside the name is
    replaced as in 0.1.0; otherwise the last dot-separated token that is a component code
    ('R', 'T', 'Z', 'BHR', ...) becomes '?'. Falls back to folder|rate|start only when the
    name has no component token."""
    base_name = rel_path.replace("\\", "/")
    d, f = os.path.split(base_name)
    if cha and cha in f:
        return d + "/" + f.replace(cha, cha[:-1] + "?")
    parts = f.split(".")
    for i in range(len(parts) - 1, -1, -1):
        if parts[i].upper() in COMP_TOKENS:
            parts[i] = "?"
            return d + "/" + ".".join(parts)
    return "%s|%.3f|%s" % (d, float(sr), start_iso[:21])


def set_suffix(key):
    """'data_IU.SDV/2000.?.sac.bp' -> '.sac.bp'; '' if the key has no '?'."""
    f = os.path.basename((key or "").replace("!", ""))
    return f.split("?", 1)[1] if "?" in f else ""


def bp_first(key):
    """Key variant that makes '.bp' sets sort before every other set of the folder."""
    if set_suffix(key).endswith(".bp"):
        d, f = os.path.split(key)
        return d + "/!" + f
    return key


def wrap180(x):
    return (float(x) + 180.0) % 360.0 - 180.0


def grade_combined(rho, delta, consistent_tol=1.0, misrot_min=2.0):
    """As sdvz24.grade_combined: eps = rho - delta (wrapped)."""
    if rho is None or delta is None or not np.isfinite(rho) or not np.isfinite(delta):
        return "unresolved", float("nan")
    eps = wrap180(rho - delta)
    if abs(eps) <= consistent_tol:
        return "consistent", eps
    if abs(eps) >= misrot_min:
        return "misrotated", eps
    return "unclear", eps


def xcorr_signed(a, b, maxlag):
    """Normalised cross-correlation of b against a for lags -maxlag..maxlag (samples;
    > 0 = b later). The central part a[m:n-m] is compared with b shifted by k.
    Returns dict(lag_abs, cc_abs) at the largest |cc| (cc keeps its sign) and
    dict(lag_pos, cc_pos) at the largest positive cc; lags parabolically refined."""
    a = np.nan_to_num(np.asarray(a, float))
    b = np.nan_to_num(np.asarray(b, float))
    n, m = len(a), int(maxlag)
    if len(b) != n or n <= 2 * m + 2:
        raise ValueError("xcorr_signed: bad lengths")
    seg = a[m:n - m]
    na = math.sqrt(float(np.dot(seg, seg)))
    cc = np.full(2 * m + 1, np.nan)
    for j, k in enumerate(range(-m, m + 1)):
        s = b[m + k:n - m + k]
        nb = math.sqrt(float(np.dot(s, s)))
        if na > 0 and nb > 0:
            cc[j] = float(np.dot(seg, s)) / (na * nb)
    if not np.any(np.isfinite(cc)):
        return dict(lag_abs=float("nan"), cc_abs=float("nan"), lag_pos=float("nan"),
                    cc_pos=float("nan"))

    def refine(j, y):
        if 0 < j < len(y) - 1 and np.all(np.isfinite(y[j - 1:j + 2])):
            d = y[j - 1] - 2 * y[j] + y[j + 1]
            if d != 0:
                return j + 0.5 * (y[j - 1] - y[j + 1]) / d
        return float(j)
    c0 = np.nan_to_num(cc, nan=0.0)
    ja = int(np.argmax(np.abs(c0)))
    jp = int(np.argmax(c0))
    sgn = 1.0 if c0[ja] >= 0 else -1.0
    return dict(lag_abs=refine(ja, sgn * c0) - m, cc_abs=float(c0[ja]),
                lag_pos=refine(jp, c0) - m, cc_pos=float(c0[jp]))


def sha256(b):
    return hashlib.sha256(b).hexdigest()


def fetch(url, tries=3, timeout=120):
    err = ""
    for k in range(tries):
        try:
            with urllib.request.urlopen(url, timeout=timeout) as r:
                return r.read()
        except Exception as ex:  # noqa: BLE001
            err = repr(ex)[:160]
            time.sleep(10 * (k + 1))
    raise RuntimeError("fetch failed %s: %s" % (url, err))


def jsonable(o):
    if isinstance(o, np.bool_):
        return bool(o)
    if isinstance(o, np.integer):
        return int(o)
    if isinstance(o, np.floating):
        return float(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    return str(o)


# ----------------------------------------------------------------- run

def load_pinned(dest):
    """Fetch sdv23.py and sdvz24.py at CODE_COMMIT into dest and import sdvz24."""
    shas = {}
    for name, must in MUST.items():
        body = fetch(REPO_RAW + CODE_COMMIT + "/" + name)
        if must not in body:
            raise RuntimeError("unexpected %s at %s" % (name, CODE_COMMIT))
        with open(os.path.join(dest, name), "wb") as f:
            f.write(body)
        shas[name] = sha256(body)
    sys.path.insert(0, dest)
    for mod in ("sdvz24", "sdv23"):
        sys.modules.pop(mod, None)
    return importlib.import_module("sdvz24"), shas


def run_variant(m, keyfun, data_dir, out_dir):
    m.set_key = keyfun
    os.makedirs(out_dir, exist_ok=True)
    return m.run(types.SimpleNamespace(data_dir=data_dir, out_dir=out_dir, wait_min=0.0))


def grades_of(S, expected_suffix):
    prim = S.get("primary_set", {})
    like = all(set_suffix(v) == expected_suffix for v in prim.values()) and len(prim) == 2
    fr = S.get("frame", {})
    rho = fr.get("rho", float("nan"))
    rho = float(rho) if rho is not None else float("nan")
    if fr.get("grade") not in (None, "untestable") and P23_G5:
        comb, eps = grade_combined(rho, P23_DELTA)
    else:
        comb, eps = "unresolved", float("nan")
    me = S.get("measurement", {})
    return dict(primary_sets=prim, like_with_like=bool(like), frame=fr.get("grade"),
                rho=rho, frame_gates=fr.get("gates"), combined_p23=comb, eps_p23=eps,
                delta_p23=P23_DELTA, measurement=me.get("grade"), x0=me.get("x0"),
                measurement_gates=me.get("gates"), lags0=me.get("lags0"),
                slope_s_per_deg=me.get("slope_s_per_deg"), x_at=me.get("x_at"),
                ratios=me.get("ratios"))


def filter_diag(m, data_dir, geometry):
    """Authors' '.bp' traces against their plain traces, per event, component, window."""
    from obspy import UTCDateTime
    P = m.PARAMS
    root = os.path.join(data_dir, "sdvz24", "x")
    sr = DIAG["sr"]
    g0, g1 = geometry["tSKS"] - 40.0, geometry["tS"] + 60.0
    wins = {"all": (g0 + 5.0, g1 - 5.0),
            "SKS": (geometry["tSKS"] + DIAG["sks_win"][0], geometry["tSKS"] + DIAG["sks_win"][1]),
            "S": (geometry["tS"] + DIAG["s_win"][0], geometry["tS"] + DIAG["s_win"][1])}
    rows = []
    for eid, (tiso, _, _, _) in m.EVENTS.items():
        yr = tiso[:4]
        org = UTCDateTime(tiso)
        for comp in ("R", "T", "Z"):
            p1 = "data_IU.SDV/%s.%s.sac" % (yr, comp)
            p2 = p1 + ".bp"
            if not (os.path.exists(os.path.join(root, p1)) and os.path.exists(os.path.join(root, p2))):
                rows.append(dict(event=eid, comp=comp, status="missing"))
                continue
            a = m.their_trace(root, dict(path=p1, index=0), org, P)
            b = m.their_trace(root, dict(path=p2, index=0), org, P)
            for wn, (t1, t2) in wins.items():
                ml = DIAG["maxlag_s"]
                ga = m.on_grid(a[0], a[1], a[2], t1 - ml, t2 + ml, sr)
                gb = m.on_grid(b[0], b[1], b[2], t1 - ml, t2 + ml, sr)
                r = xcorr_signed(ga, gb, int(round(ml * sr)))
                rows.append(dict(event=eid, comp=comp, window=wn, status="ok",
                                 lag_abs_s=r["lag_abs"] / sr, cc_abs=r["cc_abs"],
                                 lag_pos_s=r["lag_pos"] / sr, cc_pos=r["cc_pos"],
                                 sr_plain=a[1], sr_bp=b[1]))
    diff = {}
    for comp in ("R", "T", "Z"):
        for wn in ("all", "SKS", "S"):
            v = {r["event"]: r for r in rows if r.get("status") == "ok" and r["comp"] == comp
                 and r["window"] == wn}
            if len(v) == 2:
                diff["%s_%s" % (comp, wn)] = v["E2009a"]["lag_abs_s"] - v["E2000a"]["lag_abs_s"]
    return rows, diff


def run(args):
    t0 = time.time()
    os.makedirs(args.out_dir, exist_ok=True)
    S = dict(version=VERSION, code_commit=CODE_COMMIT, errors=[],
             started=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()))
    with open(os.path.abspath(__file__), "rb") as f:
        S["code_sha256"] = sha256(f.read())
    tmp = tempfile.mkdtemp(prefix="z24b_")
    m, shas = load_pinned(tmp)
    S["pinned_sha256"] = shas
    ok, psha = m.base.check_plan(args.wait_min, AMEND_PATH)
    S["amendment_public_at_start"], S["amendment_sha256"] = ok, psha
    log("amendment public:", ok, psha[:16])

    variants = {}
    for name, keyfun, sub in (("primary", set_key_v2, "plain"),
                              ("secondary", lambda *a: bp_first(set_key_v2(*a)), "bp")):
        try:
            Sv = run_variant(m, keyfun, args.data_dir, os.path.join(args.out_dir, sub))
            variants[name] = grades_of(Sv, EXPECTED_SUFFIX[name])
            if name == "primary":
                S["geometry"] = Sv.get("geometry")
            log(name, variants[name]["primary_sets"], variants[name]["frame"],
                variants[name]["rho"], variants[name]["combined_p23"],
                variants[name]["measurement"], variants[name]["x0"])
        except Exception as ex:  # noqa: BLE001
            S["errors"].append("%s: %s" % (name, repr(ex)[:200]))
            variants[name] = dict(frame="untestable", measurement="untestable",
                                  combined_p23="unresolved")
    S["variants"] = variants
    pr = variants.get("primary", {})
    S["graded"] = dict(
        gate_A1_like_with_like=bool(pr.get("like_with_like")),
        frame=pr.get("frame") if pr.get("like_with_like") else "untestable",
        combined=pr.get("combined_p23") if pr.get("like_with_like") else "unresolved",
        measurement=pr.get("measurement") if pr.get("like_with_like") else "untestable")
    log("GRADED", S["graded"])

    try:
        if S.get("geometry"):
            rows, diff = filter_diag(m, args.data_dir, S["geometry"])
            S["bp_vs_plain"] = rows
            S["bp_delay_2009_minus_2000_s"] = diff
            log("bp vs plain differential delays", diff)
    except Exception as ex:  # noqa: BLE001
        S["errors"].append("filter_diag: " + repr(ex)[:200])

    S["finished"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    S["runtime_min"] = round((time.time() - t0) / 60.0, 1)
    with open(os.path.join(args.out_dir, "sdvz24b_summary.json"), "w") as f:
        json.dump(S, f, indent=1, default=jsonable)
    return S


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", default="/data")
    ap.add_argument("--out-dir", default="/out")
    ap.add_argument("--wait-min", type=float, default=120.0)
    run(ap.parse_args(argv))
    return 0


if __name__ == "__main__":
    sys.exit(main())

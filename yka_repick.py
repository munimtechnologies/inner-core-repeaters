#!/usr/bin/env python3
"""yka_repick.py: the M03-24 (Z004) check at YKA with the PKIKP onset picked near ak135.

Part of inner-core-repeaters by MV2, Munim, Inc. (made with AI agents). Licence: Apache-2.0.
Pre-registered in findings/2026-11-18-yka-pick-amendment.md (written before this ran).

Why: icr-m0324chk-v1 ran yka.py v0.2.0 on the nine Z004 events at YKA. Its STA/LTA picker
takes the FIRST trigger within ak135 +/- 6 s. On the very quiet 9-event stack it triggered
4.8 s before ak135 on weak coherent energy, while the strong, impulsive PKIKP at YKA (about
139 deg from the South Sandwich Islands) arrives within about 1 s of ak135 (stage A quicklooks
2011-02-20 and 2025-12-15). So the "direct pulse" window held pre-PKIKP energy, the direct
SNR gate failed for every pair, and the "coda" window held only the first ~1 s of PKIKP.

What this does, and nothing else:
  * fetches nothing itself; it imports the v0.8.0 modules (m0324_check, model, phases, xcorr)
    from --mods (see fetch_pinned.py) and yka.py v0.2.0 the way m0324_check does;
  * replaces yka's onset picker by the SAME picker restricted to ak135 +/- REPICK_S (2.0 s);
  * runs m0324_check.run_yka unchanged (same alignment, windows, gates, grade_y rules);
  * records, as diagnostics only, the original +/-6 s pick and an envelope onset.
"""
import argparse
import json
import math
import os
import sys
import time
import traceback

import numpy as np

VERSION = "0.9.2"
REPICK_S = 2.0       # pre-registered: onset searched within ak135 +/- 2 s
VALID_PICK_S = 2.0   # pre-registered: a multiplet pick is valid if |pick - ak135| <= 2 s
ENV_FRAC = 0.25      # diagnostic envelope onset: first sample >= 25% of the max envelope
ORIG_SEARCH_S = 6.0  # what yka.py v0.2.0 uses


def constrained_picker(base_pick, search_s=REPICK_S, record=None, fs_default=40.0):
    """Wrap pick_onset so any caller's search_s is replaced by search_s.

    If record is a list, each call appends a dict of diagnostics (all in seconds relative to
    the predicted sample `center`): the constrained pick, the original +/-6 s pick and the
    envelope onset.
    """
    def pick(x, fs, center=None, **kw):
        kw.pop("search_s", None)
        out = base_pick(x, fs, center, search_s=search_s, **kw)
        if record is not None:
            try:
                orig = base_pick(x, fs, center, search_s=ORIG_SEARCH_S, **kw)
                env = envelope_onset(x, fs, center)
                rel = lambda i: None if i is None else round((i - center) / float(fs), 3)
                record.append(dict(pick_constrained_s=rel(out), pick_original_s=rel(orig),
                                   envelope_onset_s=rel(env[0]), envelope_max_s=rel(env[1]),
                                   search_s=search_s))
            except Exception as ex:  # diagnostics must never change the measurement
                record.append(dict(error=repr(ex)[:200]))
        return out
    return pick


def envelope_onset(x, fs, center, half_s=ORIG_SEARCH_S, frac=ENV_FRAC, smooth_s=0.1):
    """(onset, argmax) sample indices of a smoothed |x| envelope within center +/- half_s.

    onset = first sample in the range whose envelope >= frac * the range maximum.
    """
    x = np.asarray(x, dtype=float)
    w = max(1, int(round(smooth_s * fs)))
    e = np.sqrt(np.convolve(x ** 2, np.ones(w) / w, mode="same"))
    a = max(0, int(center - half_s * fs))
    b = min(len(x), int(center + half_s * fs))
    if b <= a:
        return None, None
    seg = e[a:b]
    k = int(np.argmax(seg))
    if not seg[k] > 0:
        return None, None
    hit = np.nonzero(seg >= frac * seg[k])[0]
    return int(a + hit[0]), int(a + k)


def pick_valid(pick_minus_pred_s, lim=VALID_PICK_S):
    """Pre-registered validity of a multiplet-stack pick relative to ak135."""
    try:
        v = float(pick_minus_pred_s)
    except (TypeError, ValueError):
        return False
    return math.isfinite(v) and abs(v) <= lim


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mods", required=True, help="folder with the pinned v0.8.0 modules")
    ap.add_argument("--data-dir", default="/data")
    ap.add_argument("--out-dir", default="/out/repick")
    args = ap.parse_args()
    t0 = time.time()
    os.makedirs(args.out_dir, exist_ok=True)
    plots = os.path.join(args.out_dir, "plots")
    os.makedirs(plots, exist_ok=True)
    mods = os.path.abspath(args.mods)
    if mods not in sys.path:
        sys.path.insert(0, mods)
    import m0324_check as mc  # noqa: E402
    import xcorr  # noqa: E402

    s = dict(version=VERSION, rule=dict(repick_search_s=REPICK_S, valid_pick_s=VALID_PICK_S,
                                        original_search_s=ORIG_SEARCH_S),
             started=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()))
    yka = mc.load_yka_module("/tmp/icr_mods")
    if yka is None:
        s["error"] = "yka.py download failed"
    else:
        rec = []
        yka.pick_onset = constrained_picker(xcorr.pick_onset, REPICK_S, record=rec)
        evs = mc.read_events(os.path.join(mods, "z004_events.csv"))
        try:
            mc.run_yka(evs, list(evs), args.data_dir, args.out_dir, plots, s)
        except Exception as ex:
            traceback.print_exc()
            s["yka"] = dict(status="exception", error=repr(ex)[:300])
        s["pick_diagnostics"] = rec
        y = s.get("yka", {})
        info = y.get("multiplet", {}) if isinstance(y, dict) else {}
        s["pick_valid"] = pick_valid(info.get("pick_minus_pred_s"))
        tests = y.get("tests") if isinstance(y, dict) else None
        if tests:
            s["Y1R"] = tests.get("Y1")
            s["Y0R"] = tests.get("Y0")
            if not s["pick_valid"]:
                s["Y1R_note"] = "pick invalid: Y1R/Y0R reported but not interpreted"
    s["finished"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    s["minutes"] = round((time.time() - t0) / 60.0, 1)
    with open(os.path.join(args.out_dir, "repick_summary.json"), "w") as f:
        json.dump(mc.jsonable(s), f, indent=1)
    print(json.dumps({k: s.get(k) for k in ("pick_valid", "Y1R", "Y0R", "pick_diagnostics")}, default=str)[:3000],
          flush=True)


if __name__ == "__main__":
    main()

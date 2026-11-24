#!/usr/bin/env python3
# Copyright 2026 MV2, Munim, Inc.
# Licensed under the Apache License, Version 2.0 (Apache-2.0). See LICENSE.
# Made with AI agents working in MV2's name.
"""closure35.py - pre-registration 35: the rotation-closure test for repeating earthquakes.

Idea. If the only thing that changes between two repeats of an earthquake is the inner
core's orientation (rigid differential rotation; source and mantle fixed), every
array's waveform is a function of ONE angle phi(t). So when one array A sees a full
"recovery" (event c looks like an earlier event a, after an intermediate event b looked
different), the simplest rotation reading is phi(c) = phi(a), and then EVERY array must
see a and c as the same. An array B that sees a and c as different breaks this closure:
rotation can only survive it through a coincidence (A's response returning to the same
waveform at a different angle). Local surface change has no such constraint.

Input: a CSV with columns channel,event_i,event_j and either label (same / changed / nd)
or D plus D_floor (dissimilarity and its noise floor; labels then come from the
pre-registered ratio rule in label_from_ratio). Events are ISO dates (YYYY-MM-DD).

Output (--out-dir): closure35_triples.csv (one row per recovery triple x other channel)
and closure35_summary.json (counts and the headline grade).

Stdlib only.
"""

import argparse
import csv
import json
import os
import sys

VERSION = "0.1.0"
SAME_MAX_RATIO = 2.0      # D / D_floor <= 2.0  -> same
CHANGED_MIN_RATIO = 3.0   # D / D_floor >= 3.0  -> changed ; in between -> nd
FLOOR_MIN = 0.005         # D_floor is never taken below this


def label_from_ratio(d, d_floor):
    """Pre-registered labelling rule. Returns 'same', 'changed' or 'nd'.
    NaN, None or negative inputs give 'nd'."""
    try:
        d = float(d)
        f = float(d_floor)
    except (TypeError, ValueError):
        return "nd"
    if d != d or f != f or d < 0 or f < 0:
        return "nd"
    r = d / max(f, FLOOR_MIN)
    if r <= SAME_MAX_RATIO:
        return "same"
    if r >= CHANGED_MIN_RATIO:
        return "changed"
    return "nd"


def _key(a, b):
    return (a, b) if a <= b else (b, a)


def build_labels(rows):
    """rows: dicts with channel, event_i, event_j and label or (D, D_floor).
    Returns {channel: {(early, late): label}}. A pair given twice with different labels
    becomes 'nd'."""
    out = {}
    for r in rows:
        ch = (r.get("channel") or "").strip()
        a = (r.get("event_i") or "").strip()
        b = (r.get("event_j") or "").strip()
        if not ch or not a or not b or a == b:
            continue
        lab = (r.get("label") or "").strip().lower()
        if lab not in ("same", "changed", "nd"):
            lab = label_from_ratio(r.get("D"), r.get("D_floor"))
        k = _key(a, b)
        chl = out.setdefault(ch, {})
        if k in chl and chl[k] != lab:
            chl[k] = "nd"
        else:
            chl[k] = lab
    return out


def get(labels, ch, a, b):
    if a == b:
        return "same"
    return labels.get(ch, {}).get(_key(a, b), "missing")


def implied(labels, ch, a, b, c):
    """Label of (a, c) at channel ch implied by transitivity through b, or None.
    a~b and b!~c -> changed; a!~b and b~c -> changed; a~b and b~c -> same."""
    ab = get(labels, ch, a, b)
    bc = get(labels, ch, b, c)
    if ab == "same" and bc == "changed":
        return "changed"
    if ab == "changed" and bc == "same":
        return "changed"
    if ab == "same" and bc == "same":
        return "same"
    return None


def recovery_triples(labels, ch):
    """All (a, b, c) with a < b < c (dates) such that at channel ch: changed(a, b) and
    same(a, c). (changed(b, c) is then implied and is not required.)"""
    evs = sorted({e for k in labels.get(ch, {}) for e in k})
    out = []
    for i, a in enumerate(evs):
        for j in range(i + 1, len(evs)):
            b = evs[j]
            if get(labels, ch, a, b) != "changed":
                continue
            for c in evs[j + 1:]:
                if get(labels, ch, a, c) == "same":
                    out.append((a, b, c))
    return out


def grade(labels):
    """For each recovery triple at a channel A and each other channel B, the closure
    grade of B:
      holds              B measured same(a, c)
      violated           B measured changed(a, c)
      violated-implied   B has no usable (a, c) label but transitivity through b gives changed
      holds-implied      transitivity through b gives same
      untestable         otherwise
    Headline over measured grades only: 'violated' if any triple has a measured violation
    and no channel measured 'holds' for that same triple; 'holds' if some triple has
    measured holds and none is violated; 'mixed' if both; 'untestable' if no measured
    grade."""
    rows = []
    chans = sorted(labels)
    for A in chans:
        for (a, b, c) in recovery_triples(labels, A):
            for B in chans:
                if B == A:
                    continue
                m = get(labels, B, a, c)
                if m == "same":
                    g = "holds"
                elif m == "changed":
                    g = "violated"
                else:
                    imp = implied(labels, B, a, b, c)
                    g = {"changed": "violated-implied", "same": "holds-implied"}.get(imp,
                                                                                  "untestable")
                rows.append(dict(recovery_channel=A, a=a, b=b, c=c, other_channel=B,
                                 other_ac=m, other_ab=get(labels, B, a, b),
                                 other_bc=get(labels, B, b, c), grade=g))
    by_triple = {}
    for r in rows:
        by_triple.setdefault((r["recovery_channel"], r["a"], r["b"], r["c"]), []).append(
            r["grade"])
    n_viol = sum(1 for g in by_triple.values() if "violated" in g and "holds" not in g)
    n_hold = sum(1 for g in by_triple.values() if "holds" in g and "violated" not in g)
    n_both = sum(1 for g in by_triple.values() if "holds" in g and "violated" in g)
    if n_viol and not n_hold and not n_both:
        head = "violated"
    elif n_hold and not n_viol and not n_both:
        head = "holds"
    elif n_viol or n_hold or n_both:
        head = "mixed"
    else:
        head = "untestable"
    n_imp = sum(1 for g in by_triple.values() if "violated-implied" in g)
    summary = dict(version=VERSION, n_channels=len(chans), n_triples=len(by_triple),
                   triples_violated=n_viol, triples_holds=n_hold, triples_mixed=n_both,
                   triples_violated_implied=n_imp, headline=head,
                   rule=dict(same_max_ratio=SAME_MAX_RATIO,
                             changed_min_ratio=CHANGED_MIN_RATIO, floor_min=FLOOR_MIN))
    return rows, summary


def read_csv(path):
    with open(path, newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def run(args):
    os.makedirs(args.out_dir, exist_ok=True)
    labels = build_labels(read_csv(args.pairs))
    rows, summary = grade(labels)
    summary["input"] = os.path.basename(args.pairs)
    keys = ["recovery_channel", "a", "b", "c", "other_channel", "other_ac", "other_ab",
            "other_bc", "grade"]
    with open(os.path.join(args.out_dir, "closure35_triples.csv"), "w", newline="",
              encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        for r in rows:
            w.writerow(r)
    with open(os.path.join(args.out_dir, "closure35_summary.json"), "w",
              encoding="utf-8") as f:
        json.dump(summary, f, indent=1)
    print(json.dumps(summary))
    return 0


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--pairs", required=True)
    ap.add_argument("--out-dir", default="/out")
    return run(ap.parse_args(argv))


if __name__ == "__main__":
    sys.exit(main())

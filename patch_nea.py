#!/usr/bin/env python3
"""Addendum T to pre-registration 16: fix the TauP phase name in nea.py 0.1.0.

nea.py 0.1.0 (tag v0.14.7) asked ObsPy's TauP for the phase "PKPab". ObsPy does not accept branch
suffixes, so every travel-time call raised ValueError and compete job icr-nea-v2 analysed no
waveform at all (1,339 of 1,340 event-station rows: taup-error; the other row had no epoch).

This script changes exactly two expressions and nothing else:
  phase_list=["PKIKP", "PKPab"]                         -> phase_list=["PKIKP", "PKP"]
  next((a.time for a in arr if a.name == "PKPab"), NAN) -> max((a.time for a in arr if a.name == "PKP"), default=NAN)
"PKP" returns every PKP ray that turns in the outer core (the AB and, where it exists, BC branch).
At 155-179.5 deg the AB ray is the latest-arriving PKP, so the AB time is the maximum PKP time.
No window, gate, threshold, station rule or grading rule changes.

It refuses to patch (exit 3) unless the addendum is public on the repo's main branch, and it
writes nea_patch.json (sha256 before/after) to --out-dir.

Part of inner-core-repeaters by MV2, Munim, Inc. (made with AI agents). Licence: Apache-2.0.
"""
import argparse
import hashlib
import json
import os
import re
import sys
import time
import urllib.request

ADDENDUM_URL = ("https://raw.githubusercontent.com/munimtechnologies/inner-core-repeaters/main/"
                "findings/2026-11-19-prereg16-addendum-taup-fix.md")
MARK = "ADDENDUM-T-TAUP-FIX"
PATCH_ID = "prereg16-addendum-T"

NEW_PHASES = 'phase_list=["PKIKP", "PKP"]'
NEW_AB = 'max((a.time for a in arr if a.name == "PKP"), default=NAN)'
RULES = [
    (re.compile(r'phase_list=\[\s*([\'"])PKIKP\1\s*,\s*([\'"])PKPab\2\s*\]'), NEW_PHASES),
    (re.compile(r'next\(\(\s*a\.time\s+for\s+a\s+in\s+arr\s+if\s+a\.name\s*==\s*([\'"])PKPab\1\s*\)\s*,\s*NAN\s*\)'),
     NEW_AB),
]


def ab_time(arrivals, nan=float("nan")):
    """Same rule as NEW_AB: the latest-arriving ray named "PKP" (the AB branch at 155-179.5 deg)."""
    return max((a.time for a in arrivals if a.name == "PKP"), default=nan)


def apply_patch(src):
    """Return the patched source. Raises ValueError unless each rule matches exactly once."""
    out = src
    for rx, new in RULES:
        out, n = rx.subn(lambda _m, _new=new: _new, out)
        if n != 1:
            raise ValueError("pattern %r matched %d times (need exactly 1)" % (rx.pattern[:40], n))
    if re.search(r'name\s*==\s*[\'"]PKPab[\'"]|[\'"]PKPab[\'"]\s*\]', out):
        raise ValueError("a PKPab phase request is still present after patching")
    return out


def addendum_public(wait_min):
    deadline = time.time() + 60.0 * wait_min
    while True:
        try:
            req = urllib.request.Request(ADDENDUM_URL + "?t=%d" % int(time.time()),
                                         headers={"User-Agent": "inner-core-repeaters (MV2, Munim, Inc.)"})
            with urllib.request.urlopen(req, timeout=60) as r:
                if r.status == 200 and MARK in r.read().decode("utf-8", "replace"):
                    return True
        except Exception as ex:  # noqa: BLE001
            print("addendum fetch failed", repr(ex)[:160], flush=True)
        if time.time() > deadline:
            return False
        time.sleep(300)


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--file", default="nea.py")
    ap.add_argument("--out-dir", default="/out")
    ap.add_argument("--wait-min", type=float, default=120.0)
    a = ap.parse_args(argv)
    if not addendum_public(a.wait_min):
        print("addendum T is not public on main; refusing to patch", flush=True)
        return 3
    with open(a.file, "r", encoding="utf-8") as f:
        src = f.read()
    if 'VERSION = "0.1.0"' not in src:
        print("nea.py is not version 0.1.0; refusing to patch", flush=True)
        return 2
    new = apply_patch(src)
    with open(a.file, "w", encoding="utf-8", newline="") as f:
        f.write(new)
    rec = dict(patch=PATCH_ID, addendum=ADDENDUM_URL,
               sha256_before=hashlib.sha256(src.encode("utf-8")).hexdigest(),
               sha256_after=hashlib.sha256(new.encode("utf-8")).hexdigest(),
               replaced=[NEW_PHASES, NEW_AB])
    os.makedirs(a.out_dir, exist_ok=True)
    with open(os.path.join(a.out_dir, "nea_patch.json"), "w") as f:
        json.dump(rec, f, indent=1)
    print(json.dumps(rec), flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())

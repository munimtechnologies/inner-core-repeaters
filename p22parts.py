#!/usr/bin/env python3
# Copyright 2026 MV2, Munim, Inc.
# Licensed under the Apache License, Version 2.0 (Apache-2.0). See LICENSE.
# Made with AI agents working in MV2's name.
"""p22parts.py - run pre-registration 22 (dtdiff22.py, unchanged, pinned commit 2b162f4) on
pre-registration 21's pair table when that table is published in parts.

Why: results/mkar21-v1/mkar20_pairs.csv (12,569 lines) is over the 2 MB single-file limit of
our publishing route, so it is committed as results/mkar21-v1/mkar20_pairs.part1.csv ...
plus a manifest results/mkar21-v1/mkar20_pairs.split.json (sha256 of the whole table and of
each part). This wrapper replaces only dtdiff22.wait_for for the pair-table path: it polls
main for the manifest and the parts, checks every sha256, joins the parts (part 1 whole, later
parts without their header line) and hands dtdiff22 the exact original bytes. Nothing in the
measurement or the decision rules changes. Standard library only at top level.

Usage (after fetch22.py has placed dtdiff22.py next to this file):
    python p22parts.py --data-dir /data --out-dir /out --wait-min 240 --time-budget-min 240
"""

import argparse
import hashlib
import json
import os
import sys
import time

VERSION = "0.1.0"
MANIFEST = "results/mkar21-v1/mkar20_pairs.split.json"
POLL_S = 120


class PartsError(ValueError):
    pass


def sha(b):
    return hashlib.sha256(b).hexdigest()


def join_parts(manifest, bodies):
    """manifest: dict from mkar20_pairs.split.json; bodies: {part name: bytes}.
    Returns the joined original bytes, or raises PartsError if any check fails."""
    parts = manifest.get("parts") or []
    if not parts:
        raise PartsError("manifest lists no parts")
    header = None
    chunks = []
    for i, p in enumerate(parts):
        name = p.get("name")
        if name not in bodies:
            raise PartsError("missing part %s" % name)
        b = bodies[name]
        if sha(b) != p.get("sha256"):
            raise PartsError("sha256 mismatch for %s" % name)
        nl = b.find(b"\n")
        if nl < 0:
            raise PartsError("part %s has no header line" % name)
        h = b[:nl + 1]
        if header is None:
            header = h
            chunks.append(b)
        else:
            if h != header:
                raise PartsError("header of %s differs from part 1" % name)
            chunks.append(b[len(h):])
    full = b"".join(chunks)
    if sha(full) != manifest.get("sha256"):
        raise PartsError("sha256 of the joined table does not match the manifest")
    if "bytes" in manifest and len(full) != int(manifest["bytes"]):
        raise PartsError("joined size differs from the manifest")
    return full


def part_path(name):
    return os.path.dirname(MANIFEST) + "/" + name


def make_wait_for(orig, fetch, pairs_path, record, sleep=time.sleep, clock=time.time):
    """Return a wait_for(path, must, wait_min) that serves pairs_path from the verified parts and
    passes every other path to orig. fetch(repo_path) -> bytes (raises if absent).
    record (dict) receives the manifest sha256, part names and the last error."""

    def wait_for(path, must, wait_min):
        if path != pairs_path:
            return orig(path, must, wait_min)
        t0 = clock()
        while True:
            try:
                manifest = json.loads(fetch(MANIFEST).decode("utf-8"))
                bodies = {p["name"]: fetch(part_path(p["name"])) for p in manifest.get("parts", [])}
                full = join_parts(manifest, bodies)
                if must in full.lower():
                    record.update(manifest_sha256=manifest.get("sha256"),
                                  parts=[p["name"] for p in manifest["parts"]], error="")
                    return full, round((clock() - t0) / 60.0, 1)
                record["error"] = "joined table lacks the expected header"
            except Exception as ex:  # noqa: BLE001
                record["error"] = repr(ex)[:300]
            if clock() - t0 > wait_min * 60:
                return None, round((clock() - t0) / 60.0, 1)
            sleep(POLL_S)

    return wait_for


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", default="/data")
    ap.add_argument("--out-dir", default="/out")
    ap.add_argument("--wait-min", type=float, default=240.0)
    ap.add_argument("--time-budget-min", type=float, default=240.0)
    a = ap.parse_args(argv)
    here = os.path.dirname(os.path.abspath(__file__))
    if here not in sys.path:
        sys.path.insert(0, here)
    os.makedirs(a.out_dir, exist_ok=True)
    import dtdiff22  # placed here by fetch22.py at the pinned commit

    record = dict(version=VERSION, manifest=MANIFEST, dtdiff22_version=dtdiff22.VERSION)

    def fetch(path):
        return dtdiff22.fetch_bytes(dtdiff22.RAW + "main/" + path, tries=1, wait=0, timeout=60)

    dtdiff22.wait_for = make_wait_for(dtdiff22.wait_for, fetch, dtdiff22.PAIRS, record)
    rc = dtdiff22.main(["--data-dir", a.data_dir, "--out-dir", a.out_dir, "--wait-min", str(a.wait_min),
                        "--time-budget-min", str(a.time_budget_min)])
    with open(os.path.join(a.out_dir, "p22parts.json"), "w", encoding="utf-8") as f:
        json.dump(record, f, indent=1)
    return rc


if __name__ == "__main__":
    sys.exit(main())

#!/usr/bin/env python3
"""fetch_pinned.py: download the published inner-core-repeaters files at a pinned commit.

Part of inner-core-repeaters by MV2, Munim, Inc. (made with AI agents). Licence: Apache-2.0.

Compete jobs run on whichever host is free, and dataset volumes are per host. This puts the
exact v0.8.0 code (commit 7067d4f) next to a job without copying files by hand, so a job can
rebuild its own cache with the published run.py / m0324.py before measuring.
"""
import argparse
import os
import sys
import time
import urllib.request

COMMIT = "7067d4fc2b0a2fdb1378a1c406dd680c8d0b69b1"  # v0.8.0
RAW = "https://raw.githubusercontent.com/munimtechnologies/inner-core-repeaters/%s/" % COMMIT
FILES = ("model.py", "phases.py", "xcorr.py", "cluster.py", "run.py", "m0324.py", "m0324_check.py",
         "bcdf.py", "test_bcdf.py", "test_m0324_check.py",
         "z004_events.csv", "z004_edges.csv", "z004_published_D.csv")


def get(url, tries=6):
    for k in range(tries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "inner-core-repeaters (MV2, Munim, Inc.)"})
            with urllib.request.urlopen(req, timeout=60) as r:
                if r.status == 200:
                    return r.read()
        except Exception as ex:
            print("fetch error", url, repr(ex)[:200], flush=True)
        time.sleep(5 * (k + 1))
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dest", required=True)
    ap.add_argument("--files", nargs="*", default=list(FILES))
    args = ap.parse_args()
    os.makedirs(args.dest, exist_ok=True)
    bad = []
    for name in args.files:
        body = get(RAW + name)
        if not body:
            bad.append(name)
            continue
        if name.endswith(".py"):
            try:
                compile(body.decode("utf-8"), name, "exec")
            except SyntaxError as ex:
                bad.append("%s (%s)" % (name, ex))
                continue
        with open(os.path.join(args.dest, name), "wb") as f:
            f.write(body)
        print("ok", name, len(body), "bytes", flush=True)
    if bad:
        print("FAILED:", bad, flush=True)
        sys.exit(1)


if __name__ == "__main__":
    main()

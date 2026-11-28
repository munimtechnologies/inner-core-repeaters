#!/usr/bin/env python3
# Copyright 2026 MV2, Munim, Inc.
# Licensed under the Apache License, Version 2.0 (Apache-2.0). See LICENSE.
# Made with AI agents working in MV2's name; no human expert has checked it.
"""get_frozen.py: put the frozen arm M front end (m2measure 0.1.0) next to a job, byte for byte.

Part of inner-core-repeaters by MV2, Munim, Inc. (mv2@munimtech.com). Licence: Apache-2.0.
Downloads m2measure/m2measure.py and m2measure/test_m2measure.py at the commit that froze them
(a7b58bd, 2026-11-27, findings/2026-11-27-prereg38-armM-m2measure-frozen.md) and refuses to continue
unless each file's git blob SHA-1 equals the one GitHub lists for that commit. Nothing is edited.
Standard library only.
"""
import argparse
import hashlib
import os
import sys
import time
import urllib.request

COMMIT = "a7b58bdb1d27cfe0fce09fd54ced519d557e5346"
RAW = "https://raw.githubusercontent.com/munimtechnologies/inner-core-repeaters/%s/" % COMMIT
FILES = {
    "m2measure/m2measure.py": "eb92f4cfdfd2f218b0aed2bffb6cc9e2e5cfd87c",
    "m2measure/test_m2measure.py": "61e7144c2979107669beb807c52033e79264876d",
}


def blob_sha(body):
    """git hash-object: SHA-1 of b'blob <len>\\0' + body."""
    h = hashlib.sha1()
    h.update(b"blob %d\x00" % len(body))
    h.update(body)
    return h.hexdigest()


def get(url, tries=6):
    for k in range(tries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "inner-core-repeaters (MV2, Munim, Inc.)"})
            with urllib.request.urlopen(req, timeout=60) as r:
                if r.status == 200:
                    return r.read()
        except Exception as ex:  # noqa: BLE001
            print("fetch error", url, repr(ex)[:200], flush=True)
        time.sleep(5 * (k + 1))
    return None


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--dest", default=".")
    args = ap.parse_args(argv)
    bad = []
    for path, want in FILES.items():
        body = get(RAW + path)
        if body is None:
            bad.append("%s (download failed)" % path)
            continue
        got = blob_sha(body)
        if got != want:
            bad.append("%s (blob %s, expected %s)" % (path, got, want))
            continue
        out = os.path.join(args.dest, os.path.basename(path))
        with open(out, "wb") as f:
            f.write(body)
        print("ok", path, len(body), "bytes, blob", got, flush=True)
    if bad:
        print("FAILED:", bad, flush=True)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())

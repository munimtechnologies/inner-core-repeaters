# star30.py 0.1.0 -- report-only companion to pre-registration 30 of inner-core-repeaters
# (amendment 1, findings/2026-11-22-prereg30-amendment-1.md).
#
# Zhou (2026) Table S1 lists 23 French tests; stars mark the 18 that formed her 112 pairs. The
# pre-registered job pkppair30.py 0.1.0 (release v0.29.0) pairs ALL rows of the table, as the plan
# says. This wrapper runs exactly the same code on the 18 starred rows only, so the pair set is
# Zhou's own. Its grade is report-only; the headline stays the one of pkppair30.py.
#
# It does not change pkppair30.py: it downloads the file pinned at tag v0.29.0, replaces only the
# table parser with one that first drops unstarred rows, and calls the original main().
#
# Part of munimtechnologies/inner-core-repeaters. Licence: Apache-2.0.
# Made with AI agents working in MV2's name (MV2, Munim, Inc.); no human expert has checked it.
import argparse
import hashlib
import json
import os
import sys

VERSION = "0.1.0"
PINNED = "https://raw.githubusercontent.com/munimtechnologies/inner-core-repeaters/v0.29.0/pkppair30.py"


def starred_only(text):
    """Keep the header line and the rows whose first field (the date) ends with '*'."""
    lines = [ln for ln in str(text).strip().splitlines() if ln.strip()]
    if not lines:
        return ""
    head, body = lines[0], lines[1:]
    keep = [ln for ln in body if ln.split(",")[0].strip().endswith("*")]
    return "\n".join([head] + keep) + "\n"


def fetch_pinned(path="pkppair30.py"):
    if os.path.exists(path):
        return path
    import requests
    r = requests.get(PINNED, timeout=60, headers={"User-Agent": "inner-core-repeaters/star30 (research)"})
    r.raise_for_status()
    with open(path, "w", encoding="utf-8") as f:
        f.write(r.text)
    return path


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-dir", default="out")
    a, _ = ap.parse_known_args(argv)
    fetch_pinned()
    sys.path.insert(0, os.getcwd())
    import pkppair30 as P
    orig = P.parse_s1
    P.parse_s1 = lambda text: orig(starred_only(text))
    rc = P.main(argv)
    sp = os.path.join(a.out_dir, "pkppair30_summary.json")
    if os.path.exists(sp):
        with open(sp) as f:
            s = json.load(f)
        s["variant"] = "starred-18 (report-only, amendment 1)"
        s["star30_version"] = VERSION
        with open(__file__, "rb") as f:
            s["star30_sha256"] = hashlib.sha256(f.read()).hexdigest()
        with open(P.__file__, "rb") as f:
            s["pkppair30_file_sha256"] = hashlib.sha256(f.read()).hexdigest()
        with open(os.path.join(a.out_dir, "star30_summary.json"), "w") as f:
            json.dump(s, f, indent=1)
    return rc


if __name__ == "__main__":
    sys.exit(main())

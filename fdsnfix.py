#!/usr/bin/env python3
# Copyright 2026 MV2, Munim, Inc.
# Licensed under the Apache License, Version 2.0 (Apache-2.0). See LICENSE.
# Made with AI agents working in MV2's name.
"""fdsnfix.py - run an already-published pilot unchanged, with ObsPy's FDSN service
discovery switched off.

Why (findings/2026-11-20-fdsn-discovery-addendum.md): by 2026-11-20 EarthScope's
dataselect application.wadl lists the query resource with path "" instead of "query".
ObsPy's Client then drops "dataselect" from its discovered services, and every
get_waveforms / get_waveforms_bulk call raises FDSNNoServiceException ("The current
client does not have a dataselect service"), although
https://service.earthscope.org/fdsnws/dataselect/1/query still serves data. icr-snrpilot-v1
lost its four hours to this.

What this does, and nothing more:
1. downloads the listed files from the public repo at one pinned commit, writes them to the
   current directory and records their sha256 in <out>/fdsnfix.json;
2. makes obspy.clients.fdsn.Client default to _discover_services=False (ObsPy then assumes
   the standard dataselect / station / event paths);
3. checks with one 10-second request that waveforms really arrive (fails fast if not);
4. runs the pinned script as __main__ with the given arguments.
The analysis code, constants and decision rules of the pinned script are not touched.

Usage:
  python fdsnfix.py --sha <commit> --files a.py test_a.py --run a.py --out-dir /out -- <args>
"""

import argparse
import hashlib
import json
import os
import runpy
import sys
import time
import urllib.request

VERSION = "0.1.0"
RAW = "https://raw.githubusercontent.com/munimtechnologies/inner-core-repeaters/%s/%s"
PROBE = ("IU", "MAJO", "00", "BHZ", "2024-05-15T20:30:00", 10.0)
BASE = "https://service.earthscope.org"


def log(*a):
    print(time.strftime("%H:%M:%S"), "[fdsnfix]", *a, flush=True)


def split_args(argv):
    """Split argv at the first '--': (own args, args passed to the pinned script)."""
    argv = list(argv)
    if "--" in argv:
        i = argv.index("--")
        return argv[:i], argv[i + 1:]
    return argv, []


def raw_url(sha, path):
    if not sha or any(c not in "0123456789abcdef" for c in sha.lower()) or len(sha) < 7:
        raise ValueError("commit must be a hex sha: %r" % sha)
    if path.startswith("/") or ".." in path.split("/"):
        raise ValueError("bad repo path: %r" % path)
    return RAW % (sha, path)


def fetch(url, tries=5, wait=10, timeout=90):
    err = ""
    for k in range(tries):
        try:
            with urllib.request.urlopen(url, timeout=timeout) as r:
                return r.read()
        except Exception as ex:  # noqa: BLE001
            err = repr(ex)[:160]
            time.sleep(wait * (k + 1))
    raise RuntimeError("fetch failed %s: %s" % (url, err))


def sha256(data):
    return hashlib.sha256(data).hexdigest()


def default_no_discovery(cls):
    """Make cls(...) default to _discover_services=False. Idempotent; returns cls."""
    if getattr(cls, "_fdsnfix_patched", False):
        return cls
    orig = cls.__init__

    def init(self, *a, **kw):
        kw.setdefault("_discover_services", False)
        orig(self, *a, **kw)

    init.__doc__ = orig.__doc__
    cls.__init__ = init
    cls._fdsnfix_patched = True
    return cls


def patch_obspy():
    from obspy.clients.fdsn import client as fdsn_client
    import obspy.clients.fdsn as fdsn
    default_no_discovery(fdsn_client.Client)
    assert fdsn.Client is fdsn_client.Client
    return fdsn_client.Client


def probe(Client, tries=3):
    """One short waveform request through the patched Client. (ok, message)."""
    from obspy import UTCDateTime
    net, sta, loc, cha, t0, dur = PROBE
    msg = ""
    for k in range(tries):
        try:
            c = Client(BASE, timeout=120)
            if "dataselect" not in getattr(c, "services", {}):
                return False, "patched client still has no dataselect service"
            st = c.get_waveforms(net, sta, loc, cha, UTCDateTime(t0), UTCDateTime(t0) + dur)
            if len(st) > 0 and sum(tr.stats.npts for tr in st) > 0:
                return True, "%d trace(s), %d samples" % (len(st), sum(tr.stats.npts for tr in st))
            msg = "empty stream"
        except Exception as ex:  # noqa: BLE001
            msg = repr(ex)[:200]
        time.sleep(15 * (k + 1))
    return False, msg


def main(argv=None):
    own, rest = split_args(sys.argv[1:] if argv is None else argv)
    ap = argparse.ArgumentParser()
    ap.add_argument("--sha", required=True)
    ap.add_argument("--files", nargs="+", required=True)
    ap.add_argument("--run", required=True)
    ap.add_argument("--out-dir", default="/out")
    ap.add_argument("--no-probe", action="store_true")
    a = ap.parse_args(own)
    if a.run not in a.files:
        raise SystemExit("--run must be one of --files")
    os.makedirs(a.out_dir, exist_ok=True)
    rec = dict(version=VERSION, sha=a.sha, files={}, args=rest,
               started=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()))
    for p in a.files:
        data = fetch(raw_url(a.sha, p))
        with open(os.path.basename(p), "wb") as f:
            f.write(data)
        rec["files"][p] = dict(sha256=sha256(data), bytes=len(data))
        log("fetched", p, rec["files"][p]["sha256"][:16], len(data), "bytes")
    Client = patch_obspy()
    if not a.no_probe:
        ok, msg = probe(Client)
        rec["probe"] = dict(ok=ok, message=msg)
        log("probe", ok, msg)
    with open(os.path.join(a.out_dir, "fdsnfix.json"), "w") as f:
        json.dump(rec, f, indent=1)
    if not a.no_probe and not rec["probe"]["ok"]:
        raise SystemExit("probe failed: dataselect still unreachable; not running " + a.run)
    sys.argv = [os.path.basename(a.run)] + rest
    sys.path.insert(0, os.getcwd())
    log("running", a.run, " ".join(rest))
    runpy.run_path(os.path.basename(a.run), run_name="__main__")
    return 0


if __name__ == "__main__":
    sys.exit(main())

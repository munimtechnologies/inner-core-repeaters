#!/usr/bin/env python3
"""Fetch mantlep2.py (pinned commit) next to mantlep3.py inside the job container.

Part of inner-core-repeaters by MV2, Munim, Inc. (made with AI agents). Licence: Apache-2.0.
In the repository itself mantlep2.py already sits next to mantlep3.py and this step is a no-op.
"""
import hashlib
import os
import sys
import time
import urllib.request

COMMIT = "b4712f0ebbc70b896231620c905ff53a4b8270fe"   # v0.13.0, mantlep2.py as run by icr-mantlep2-v1
URL = "https://raw.githubusercontent.com/munimtechnologies/inner-core-repeaters/%s/mantlep2.py" % COMMIT


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    dest = os.path.join(here, "mantlep2.py")
    if os.path.exists(dest):
        print("mantlep2.py already present")
        return 0
    last = None
    for k in range(10):
        try:
            req = urllib.request.Request(URL, headers={"User-Agent": "inner-core-repeaters (MV2, Munim, Inc.)"})
            with urllib.request.urlopen(req, timeout=60) as r:
                body = r.read()
            if b'VERSION = "0.13.0"' not in body:
                raise RuntimeError("unexpected mantlep2.py content")
            with open(dest, "wb") as f:
                f.write(body)
            print("fetched mantlep2.py sha256", hashlib.sha256(body).hexdigest())
            return 0
        except Exception as ex:  # noqa: BLE001
            last = ex
            time.sleep(30)
    print("could not fetch mantlep2.py:", repr(last))
    return 1


if __name__ == "__main__":
    sys.exit(main())

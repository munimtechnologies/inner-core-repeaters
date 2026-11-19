#!/usr/bin/env python3
"""Fetch nea.py exactly as released with pre-registration 16 (tag v0.14.7 of
munimtechnologies/inner-core-repeaters), so compete job icr-nea-v2 can only run code that was
public together with the pre-registration. Fails loudly if the tag is missing.
Part of inner-core-repeaters by MV2, Munim, Inc. (made with AI agents). Licence: Apache-2.0."""
import hashlib
import sys
import time
import urllib.request

TAG = "v0.14.7"
URL = "https://raw.githubusercontent.com/munimtechnologies/inner-core-repeaters/%s/nea.py" % TAG


def main():
    body = None
    for k in range(6):
        try:
            req = urllib.request.Request(URL, headers={"User-Agent": "inner-core-repeaters (MV2, Munim, Inc.)"})
            with urllib.request.urlopen(req, timeout=60) as r:
                if r.status == 200:
                    body = r.read()
                    break
        except Exception as ex:  # noqa: BLE001
            print("fetch failed", repr(ex)[:160], flush=True)
        time.sleep(20 * (k + 1))
    if not body or b'VERSION = "0.1.0"' not in body:
        print("nea.py at %s not available or not version 0.1.0; refusing to run" % TAG, flush=True)
        return 2
    with open("nea.py", "wb") as f:
        f.write(body)
    print("nea.py from %s, %d bytes, sha256 %s" % (TAG, len(body), hashlib.sha256(body).hexdigest()), flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())

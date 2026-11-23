# brvpkp31.py 0.1.1 -- pre-registration 31 of inner-core-repeaters (Borovoye waveforms, PKP-PP and PPP-PP).
#
# All of the measurement code is brvpkp31_core.py, the 0.1.0 code written before any waveform was looked at,
# committed unchanged. This file changes exactly one function, despike(), because the repository's own unit
# test (test_brvpkp31.test_despike) failed on 0.1.0 when it was first run (2026-11-22 23:45 UTC), before any
# waveform was downloaded:
#   1. Scale collapse. On a smooth or integer-quantised trace a 7-sample running median equals the centre
#      sample almost everywhere, so the residual is mostly exact zeros, its MAD is 0, and 0.1.0 silently
#      removed nothing. 0.1.1 uses max(1.4826 MAD of the residual, 1.4826 MAD of first differences / sqrt 2):
#      a glitch must exceed 20 robust sample-to-sample noise steps.
#   2. Replacement. Replacing a flagged run by the running median is biased on steep slopes when the run is
#      two or three samples long (the median is pulled toward the glitch). 0.1.1 replaces flagged samples by
#      linear interpolation between the nearest unflagged samples.
# Everything else (windows, gates, grades, pairing, inputs) is unchanged. Plan, with this change recorded:
# findings/2026-11-23-prereg31-brv-waveform-plan.md.
#
# Part of munimtechnologies/inner-core-repeaters. Licence: Apache-2.0.
# Made with AI agents working in MV2's name (MV2, Munim, Inc.); no human expert has checked it.
import hashlib
import json
import os
import sys

import numpy as np

import brvpkp31_core as _core
from brvpkp31_core import *  # noqa: F401,F403  (public API, unchanged)
from brvpkp31_core import _clean, _sd, _wf_types  # noqa: F401  (used by tests)

VERSION = "0.1.1"
_core.VERSION = VERSION
PARAMS = _core.PARAMS


def despike(x, k=PARAMS["spike_k"], width=PARAMS["spike_width"], max_run=PARAMS["spike_max_run"]):
    """Replace short glitches by linear interpolation. A sample is flagged when its residual from a
    width-sample running median exceeds k robust sigmas, sigma = max(1.4826 MAD(residual),
    1.4826 MAD(first differences)/sqrt(2)); only runs of <= max_run flagged samples are replaced.
    -> (cleaned copy, number of samples replaced)."""
    x = np.asarray(x, float).copy()
    n = len(x)
    if n < max(width, 3):
        return x, 0
    h = width // 2
    pad = np.pad(x, h, mode="edge")
    med = np.median(np.lib.stride_tricks.sliding_window_view(pad, width), axis=1)
    r = x - med
    s_res = 1.4826 * float(np.median(np.abs(r - np.median(r))))
    d = np.diff(x)
    s_dif = 1.4826 * float(np.median(np.abs(d - np.median(d)))) / np.sqrt(2.0)
    scale = max(s_res, s_dif)
    if not scale > 0:
        return x, 0
    idx = np.flatnonzero(np.abs(r) > k * scale)
    bad = np.zeros(n, bool)
    if idx.size:
        for run in np.split(idx, np.flatnonzero(np.diff(idx) > 1) + 1):
            if len(run) <= max_run:
                bad[run] = True
    count = int(bad.sum())
    if count:
        good = np.flatnonzero(~bad)
        x[bad] = np.interp(np.flatnonzero(bad), good, x[good])
    return x, count


# The pipeline in brvpkp31_core calls despike() through its own module globals.
_core.despike = despike


def main(argv=None):
    rc = _core.main(argv)
    # Record this wrapper's own hash next to the core's (the core records only its own file).
    try:
        out_dir = "out"
        args = list(sys.argv[1:] if argv is None else argv)
        if "--out-dir" in args:
            out_dir = args[args.index("--out-dir") + 1]
        p = os.path.join(out_dir, "brv31_summary.json")
        if os.path.isfile(p):
            with open(p) as f:
                s = json.load(f)
            with open(__file__, "rb") as f:
                s["wrapper_sha256"] = hashlib.sha256(f.read()).hexdigest()
            s["wrapper_version"] = VERSION
            with open(p, "w") as f:
                json.dump(s, f, indent=1)
    except Exception as ex:  # noqa: BLE001
        print("[brvpkp31] could not add wrapper hash: %r" % ex, flush=True)
    return rc


if __name__ == "__main__":
    sys.exit(main())

#!/usr/bin/env python3
# Copyright 2026 MV2, Munim, Inc.
# Licensed under the Apache License, Version 2.0 (Apache-2.0). See LICENSE.
# Made with AI agents working in MV2's name.
"""zw33.py - pre-registration 33 (E1) and the pre-registration 35 measurement.

Pre-registered in findings/2026-11-24-prereg33-zw26-e1-closure-plan.md, written and
published before any waveform of Zhang & Wen's Zenodo bundles was read.

Data: Zhang, X. & Wen, L. (2026), JGR Solid Earth, doi 10.1029/2026JB034595; their raw
array records on Zenodo, doi 10.5281/zenodo.19625711 (CC BY 4.0): Data_YKAR.zip and
Data_ILAR.zip (multiplet M03-24, six events each) and Data_multiF.zip (multiplet F:
2005-06-10, 2012-12-03, 2022-05-19, both arrays).

E1   With one fixed dissimilarity D and floor D_floor, do the authors' own records
     reproduce their M03-24 change / quiescence / recovery labels at each array?
P35  The decisive rotation-closure pair: is ILAR 2005-06-10 vs 2022-05-19 'same' or
     'changed', given that YKA recovers (gate G1)? Graded with closure35.py 0.1.0.
P35b Secondary: the same closure rule on the M03-24 labels measured here.

Top-level imports are stdlib only (tests run without numpy/obspy); the numerical parts
import numpy, scipy and obspy inside functions.
"""

import argparse
import calendar
import csv
import hashlib
import io
import json
import math
import os
import re
import sys
import time
import traceback
import urllib.request
import zipfile

VERSION = "0.1.0"
# closure35 0.1.0 labelling rule, repeated here so the pure helpers need no import;
# main() imports closure35 and refuses to run if its constants differ.
SAME_MAX_RATIO = 2.0
CHANGED_MIN_RATIO = 3.0
FLOOR_MIN = 0.005
NAN = float("nan")

ZENODO_API = "https://zenodo.org/api/records/19625711"
ZIPS = ("Data_YKAR.zip", "Data_ILAR.zip", "Data_multiF.zip")
COMCAT = ("https://earthquake.usgs.gov/fdsnws/event/1/query?format=csv&starttime=%s&endtime=%s"
          "&minlatitude=-64&maxlatitude=-53&minlongitude=-36&maxlongitude=-14&minmagnitude=4")
PREREG_URL = ("https://raw.githubusercontent.com/munimtechnologies/inner-core-repeaters/main/"
              "findings/2026-11-24-prereg33-zw26-e1-closure-plan.md")
UA = {"User-Agent": "inner-core-repeaters (MV2, Munim, Inc.)"}

# ---- pre-registered measurement constants (pre-registration 33)
FS = 40.0                      # resampling rate (Hz)
PRE = 40.0                     # s before predicted PKIKP (ak135) in each element cut
POST = 40.0                    # s after
N = int(round((PRE + POST) * FS))
BAND = (1.0, 2.0)              # primary band (Hz); 4-pole Butterworth, causal
BANDS_REPORT = ((0.5, 2.0),)   # report-only band
DIRECT = (-0.3, 0.9)           # direct PKIKP window around the onset: alignment only
WINDOWS = (("D_c310", (3.0, 10.0)), ("D_c16", (1.0, 6.0)), ("D_c13", (1.0, 3.0)))
PRIMARY = "D_c310"             # the graded window
NOISE = (-35.0, -5.0)          # noise window, s relative to predicted PKIKP
PAIR_CC_MIN = 0.90             # direct-window correlation needed for a usable pair
SNR_DIRECT_MIN = 5.0
SNR_CODA_MIN = 2.0
MIN_COMMON = 6                 # elements present in both events after rejection/statics
STATIC_MAX_S = 0.5
STATIC_CC_MIN = 0.5
COARSE_LAG_S = 5.0
FINE_LAG_S = 0.4
REF_LAG_S = 0.3
RECOVERY_RATIO_MAX = 0.67      # recovery: D(a,c) <= 0.67 D(a,b)
RECOVERY_RATIO_FAIL = 1.0      # no recovery: D(a,c) >= D(a,b)
ARRAY_CENTER = {"YKA": (62.4931, -114.6053), "ILAR": (64.7714, -146.8866)}
DEFAULT_EVENT = (-59.03, -25.68, 35.0)

# ---- paper labels used by E1 (Zhang & Wen 2026, Results and Fig. 6, as copied 2026-11-24)
E1_LABELS = (
    ("YKA", "2003-03-21", "2008-02-17", "changed"),
    ("YKA", "2008-02-17", "2011-02-20", "changed"),
    ("YKA", "2011-02-20", "2022-03-07", "changed"),
    ("ILAR", "2003-03-21", "2008-02-17", "changed"),
    ("ILAR", "2008-02-17", "2011-02-20", "same"),
    ("ILAR", "2011-02-20", "2022-03-07", "changed"),
)
E1_RECOVERY = (
    ("YKA", "2003-03-21", "2011-02-20", "2022-03-07"),
    ("ILAR", "2003-03-21", "2011-02-20", "2022-03-07"),
)
E1_REPORT_ONLY = (
    ("YKA", "2022-03-07", "2024-03-21", "same"),
    ("ILAR", "2022-03-07", "2024-03-21", "same"),
)
F_EVENTS = ("2005-06-10", "2012-12-03", "2022-05-19")


# ----------------------------------------------------------------------------- pure helpers
def isf(x):
    try:
        return math.isfinite(float(x))
    except (TypeError, ValueError):
        return False


def parse_event_folder(name):
    m = re.match(r"^Event_(\d{4})_(\d{2})_(\d{2})_(\d{2})_(\d{2})_(\d{2})$", name)
    if not m:
        return None
    y, mo, d, h, mi, s = m.groups()
    return dict(iso="%s-%s-%sT%s:%s:%s" % (y, mo, d, h, mi, s), date="%s-%s-%s" % (y, mo, d))


def canonical_element(sta):
    """Short-period element name, unified across the 2013 YKA renaming (YKAB0 -> YKB0,
    YKAR1 -> YKR1). ILAR IL01..IL19 unchanged. Broadband or other stations -> None."""
    s = (sta or "").strip().upper()
    m = re.match(r"^YKA?([BR])(\d)$", s)
    if m:
        return "YK%s%s" % m.groups()
    m = re.match(r"^IL(\d{2})$", s)
    if m and 1 <= int(m.group(1)) <= 19:
        return s
    return None


def element_array(el):
    if not el:
        return None
    if el.startswith("YK"):
        return "YKA"
    if el.startswith("IL"):
        return "ILAR"
    return None


def parse_filename(fn):
    parts = fn.split(".")
    if len(parts) != 5 or parts[4].lower() != "sac":
        return None
    return dict(net=parts[0], sta=parts[1], loc=parts[2], cha=parts[3])


def iso_to_epoch(s):
    s = s.strip().rstrip("Z")
    frac = 0.0
    if "." in s:
        s, f = s.split(".", 1)
        frac = float("0." + f)
    return calendar.timegm(time.strptime(s, "%Y-%m-%dT%H:%M:%S")) + frac


def epoch_to_iso(t):
    return time.strftime("%Y-%m-%dT%H:%M:%S", time.gmtime(t))


def pick_comcat(text, target_iso, max_dt=60.0):
    """The ComCat row closest in time to target_iso (within max_dt s), or None."""
    if not text:
        return None
    t0 = iso_to_epoch(target_iso)
    best = None
    for r in csv.DictReader(io.StringIO(text)):
        try:
            dt = abs(iso_to_epoch(r["time"]) - t0)
            lat, lon, dep = float(r["latitude"]), float(r["longitude"]), float(r["depth"])
        except (KeyError, ValueError, TypeError):
            continue
        if dt <= max_dt and (best is None or dt < best[0]):
            best = (dt, r, lat, lon, dep)
    if best is None:
        return None
    dt, r, lat, lon, dep = best
    mag = float(r["mag"]) if isf(r.get("mag")) else NAN
    return dict(id=r.get("id"), time=r["time"], lat=lat, lon=lon, depth=dep, mag=mag,
                dt=round(dt, 3))


def label_pair(ok, d, d_floor):
    """Pre-registered label: 'nd' unless the pair passed its gates, else closure35's rule."""
    if not ok:
        return "nd"
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


def recovery_result(d_ab, d_ac):
    """E1 recovery check: 'agree' if D(a,c) <= 0.67 D(a,b), 'disagree' if D(a,c) >= D(a,b),
    otherwise (or missing) 'nd'."""
    if not (isf(d_ab) and isf(d_ac)) or float(d_ab) <= 0:
        return "nd", NAN
    q = float(d_ac) / float(d_ab)
    if q <= RECOVERY_RATIO_MAX:
        return "agree", q
    if q >= RECOVERY_RATIO_FAIL:
        return "disagree", q
    return "nd", q


def e1_grade(labels, dvals):
    """labels: {(array, a, b): label}; dvals: {(array, a, b): D} for usable pairs only.
    Returns (rows, summary). 8 graded items: 6 labels + 2 recoveries."""
    rows = []
    for arr, a, b, exp in E1_LABELS:
        got = labels.get((arr, a, b), "missing")
        res = "agree" if got == exp else ("nd" if got in ("nd", "missing") else "disagree")
        rows.append(dict(kind="label", array=arr, a=a, b=b, c="", expected=exp, measured=got,
                         ratio="", result=res))
    for arr, a, b, c in E1_RECOVERY:
        res, q = recovery_result(dvals.get((arr, a, b)), dvals.get((arr, a, c)))
        rows.append(dict(kind="recovery", array=arr, a=a, b=b, c=c, expected="recovery",
                         measured="D(a,c)/D(a,b)", ratio=round(q, 4) if isf(q) else "",
                         result=res))
    for arr, a, b, exp in E1_REPORT_ONLY:
        got = labels.get((arr, a, b), "missing")
        res = "agree" if got == exp else ("nd" if got in ("nd", "missing") else "disagree")
        rows.append(dict(kind="report-only", array=arr, a=a, b=b, c="", expected=exp,
                         measured=got, ratio="", result=res))
    graded = [r for r in rows if r["kind"] != "report-only"]
    n_ag = sum(1 for r in graded if r["result"] == "agree")
    n_dis = sum(1 for r in graded if r["result"] == "disagree")
    n_nd = len(graded) - n_ag - n_dis
    if n_ag == len(graded):
        head = "reproduced"
    elif n_dis >= 2:
        head = "not-reproduced"
    elif n_ag >= 6 and n_dis <= 1:
        head = "mostly"
    else:
        head = "inconclusive"
    return rows, dict(n_graded=len(graded), n_agree=n_ag, n_disagree=n_dis, n_nd=n_nd,
                      headline=head)


def p35_grade(labels, present):
    """labels: {(array, a, b): label}; present: {(array, date): bool}. Pre-registration 35."""
    a, b, c = F_EVENTS
    out = dict(G2="not run in zw33-v1; separate job under the rule fixed in pre-registration 35")
    for arr in ("YKA", "ILAR"):
        for d in F_EVENTS:
            out["%s_%s_present" % (arr, d)] = bool(present.get((arr, d)))
    for arr in ("YKA", "ILAR"):
        for (x, y, tag) in ((a, b, "ab"), (a, c, "ac"), (b, c, "bc")):
            out["%s_%s" % (arr.lower(), tag)] = labels.get((arr, x, y), "missing")
    missing = [(arr, d) for arr in ("YKA", "ILAR") for d in F_EVENTS if not present.get((arr, d))]
    if missing:
        out.update(G1=None, grade="untestable",
                   reason="data: missing " + ", ".join("%s %s" % m for m in missing))
        return out
    g1 = out["yka_ab"] == "changed" and out["yka_ac"] == "same"
    out["G1"] = g1
    if not g1:
        out.update(grade="untestable", reason="G1: YKA recovery not reproduced")
        return out
    il = out["ilar_ac"]
    if il == "changed":
        out.update(grade="violated", reason="ILAR measured changed(2005, 2022)")
    elif il == "same":
        fails = []
        if out["ilar_ab"] != "same":
            fails.append("paper ILAR same(2005, 2012): measured %s" % out["ilar_ab"])
        if out["ilar_bc"] != "changed":
            fails.append("paper ILAR changed(2012, 2022): measured %s" % out["ilar_bc"])
        out.update(grade="holds", reason="ILAR measured same(2005, 2022)",
                   paper_ilar_labels_not_reproduced=fails)
    else:
        out.update(grade="untestable", reason="ILAR (2005, 2022) is %s" % il)
    return out


def closure_rows(pair_rows, set_name, band_label, window):
    out = []
    for r in pair_rows:
        if r["set"] != set_name or r["band"] != band_label:
            continue
        out.append(dict(channel=r["array"], event_i=r["a"], event_j=r["b"],
                        label=r.get("label_" + window, "nd")))
    return out


def jsonable(o):
    if isinstance(o, dict):
        return {str(k): jsonable(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [jsonable(v) for v in o]
    if isinstance(o, float):
        return o if math.isfinite(o) else None
    if hasattr(o, "item") and not isinstance(o, (str, bytes)):
        try:
            return jsonable(o.item())
        except Exception:  # noqa: BLE001
            return str(o)
    return o


def write_csv(path, rows):
    keys = []
    for r in rows:
        for k in r:
            if k not in keys:
                keys.append(k)
    with open(path, "w", newline="", encoding="utf-8") as f:
        if keys:
            w = csv.DictWriter(f, fieldnames=keys)
            w.writeheader()
            for r in rows:
                w.writerow(jsonable(r))


def fetch_text(url, tries=4, wait=15, timeout=90):
    for k in range(tries):
        try:
            req = urllib.request.Request(url, headers=UA)
            with urllib.request.urlopen(req, timeout=timeout) as r:
                if r.status == 200:
                    return r.read().decode("utf-8", "replace")
        except Exception as ex:  # noqa: BLE001
            print("fetch error", url[:120], repr(ex)[:160], flush=True)
        time.sleep(wait * (k + 1))
    return None


def download(url, path, tries=4):
    for k in range(tries):
        try:
            req = urllib.request.Request(url, headers=UA)
            with urllib.request.urlopen(req, timeout=300) as r, open(path + ".part", "wb") as f:
                while True:
                    chunk = r.read(1 << 20)
                    if not chunk:
                        break
                    f.write(chunk)
            os.replace(path + ".part", path)
            return True
        except Exception as ex:  # noqa: BLE001
            print("download error", url[:120], repr(ex)[:160], flush=True)
            time.sleep(30 * (k + 1))
    return False


def md5_of(path):
    h = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


# ----------------------------------------------------------------------------- numerics
def idx(t_rel):
    """Sample index in an element cut for a time relative to predicted PKIKP."""
    return int(round((PRE + t_rel) * FS))


def ncc(a, b):
    import numpy as np
    a = a - a.mean()
    b = b - b.mean()
    den = math.sqrt(float(np.dot(a, a)) * float(np.dot(b, b)))
    return float(np.dot(a, b)) / den if den > 0 else NAN


def lag_scan(a, b, i0, i1, lo, hi):
    """Integer lag L in [lo, hi] maximising ncc(a[i0:i1], b[i0+L:i1+L])."""
    best_l, best_c = None, -2.0
    n = len(b)
    for L in range(lo, hi + 1):
        if i0 + L < 0 or i1 + L > n or i0 < 0 or i1 > len(a):
            continue
        c = ncc(a[i0:i1], b[i0 + L:i1 + L])
        if c == c and c > best_c:
            best_l, best_c = L, c
    return best_l, best_c


def shift(x, L):
    """y[i] = x[i + L], zero-filled."""
    import numpy as np
    y = np.zeros_like(x)
    n = len(x)
    if L >= 0:
        y[:n - L] = x[L:]
    else:
        y[-L:] = x[:n + L]
    return y


def rms(x):
    import numpy as np
    return math.sqrt(float(np.mean(x ** 2))) if len(x) else NAN


def envelope(x, smooth_s=0.5):
    import numpy as np
    from scipy.signal import hilbert
    e = np.abs(hilbert(x))
    k = max(1, int(round(smooth_s * FS)))
    return np.convolve(e, np.ones(k) / k, mode="same")


def pick_onset(x, search_s=6.0, sta_s=0.5, lta_s=8.0, thr=6.0):
    """First sample within +/- search_s of predicted PKIKP where STA/LTA of energy >= thr."""
    e = x ** 2
    c = idx(0.0)
    ns, nl = int(sta_s * FS), int(lta_s * FS)
    for i in range(c - int(search_s * FS), c + int(search_s * FS)):
        if i - nl < 0 or i + ns > len(e):
            continue
        lta = float(e[i - nl:i].mean())
        if lta > 0 and float(e[i:i + ns].mean()) / lta >= thr:
            return i, True
    return c, False


def pkikp_time(model, depth, dist):
    for ph in ("PKIKP", "PKiKP"):
        arr = model.get_travel_times(source_depth_in_km=max(float(depth), 0.0),
                                     distance_in_degree=float(dist), phase_list=[ph])
        if arr:
            return arr[0].time, ph
    return None, None


def sac_value(sac, key):
    v = sac.get(key) if sac else None
    if v is None or not isf(v) or abs(float(v) + 12345.0) < 1e-3:
        return None
    return float(v)


def process_folder(folder, ev, array, band, model):
    """Read, filter and cut every short-period element of `array` in one event folder."""
    from obspy import read
    from obspy.geodetics import locations2degrees
    import numpy as np
    traces, meta = {}, {}
    for fn in sorted(os.listdir(folder)):
        p = parse_filename(fn)
        if not p:
            continue
        el = canonical_element(p["sta"])
        if el is None or element_array(el) != array or not p["cha"].upper().startswith("SH"):
            continue
        try:
            tr = read(os.path.join(folder, fn))[0]
        except Exception as ex:  # noqa: BLE001
            meta[el] = "read_error " + repr(ex)[:80]
            continue
        sac = tr.stats.get("sac", {})
        stla, stlo = sac_value(sac, "stla"), sac_value(sac, "stlo")
        if stla is None or stlo is None:
            stla, stlo = ARRAY_CENTER[array]
        dist = locations2degrees(ev["lat"], ev["lon"], stla, stlo)
        tt, ph = pkikp_time(model, ev["depth"], dist)
        if tt is None:
            meta[el] = "no PKIKP"
            continue
        t0 = ev["origin"] + tt - PRE
        if tr.stats.starttime > t0 - 2.0 or tr.stats.endtime < t0 + PRE + POST + 2.0:
            meta[el] = "record does not cover the window"
            continue
        try:
            tr.data = tr.data.astype("float64")
            tr.detrend("demean")
            tr.detrend("linear")
            tr.taper(0.05)
            tr.filter("bandpass", freqmin=band[0], freqmax=band[1], corners=4, zerophase=False)
            tr.interpolate(sampling_rate=FS, method="lanczos", a=12, starttime=t0, npts=N)
        except Exception as ex:  # noqa: BLE001
            meta[el] = "processing_error " + repr(ex)[:80]
            continue
        x = np.asarray(tr.data, dtype=float)
        if len(x) != N or not np.all(np.isfinite(x)) or float(np.std(x)) == 0.0:
            meta[el] = "bad samples"
            continue
        traces[el] = x
        meta[el] = "ok %.3f deg %s %.2f s fs=%.1f" % (dist, ph, tt, tr.stats.sampling_rate)
    return traces, meta


def event_stack(traces):
    """Noise-normalise, reject outlier elements and align elements to the event beam."""
    import numpy as np
    els = sorted(traces)
    info = dict(n_read=len(els), n_kept=0, n_static=0, elements="")
    if len(els) < 3:
        return {}, info
    A = np.array([traces[e] for e in els])
    n0, n1 = idx(NOISE[0]), idx(NOISE[1])
    nrms = np.sqrt(np.mean(A[:, n0:n1] ** 2, axis=1))
    med = float(np.median(nrms))
    keep = [i for i in range(len(els)) if nrms[i] > 0 and 0.2 * med < nrms[i] < 5.0 * med]
    info["n_kept"] = len(keep)
    if len(keep) < 3:
        return {}, info
    X = {els[i]: A[i] / nrms[i] for i in keep}
    beam0 = np.mean([X[e] for e in X], axis=0)
    L = int(round(STATIC_MAX_S * FS))
    lags = {}
    for e in X:
        l, c = lag_scan(beam0, X[e], idx(-3.0), idx(12.0), -L, L)
        if l is not None and c >= STATIC_CC_MIN:
            lags[e] = l
    if not lags:
        return {}, info
    m = int(round(float(np.median(list(lags.values())))))
    out = {e: shift(X[e], lags[e] - m) for e in lags}
    info.update(n_static=len(out), elements=";".join(sorted(out)))
    return out, info


def make_beam(aligned, els):
    import numpy as np
    b = np.mean([aligned[e] for e in els], axis=0)
    nb = rms(b[idx(NOISE[0]):idx(NOISE[1])])
    return b / nb if nb > 0 else b


def measure_pair(Sa, Sb):
    """Sa, Sb: aligned element traces of the earlier and later event. Returns a dict with D,
    D_floor and SNRs for every window, and the gates."""
    common = sorted(set(Sa) & set(Sb))
    res = dict(n_common=len(common), ok_direct=False, reason="")
    if len(common) < MIN_COMMON:
        res["reason"] = "fewer than %d common elements" % MIN_COMMON
        return res
    ba, bb = make_beam(Sa, common), make_beam(Sb, common)
    Lc = int(round(COARSE_LAG_S * FS))
    lc, cc_env = lag_scan(envelope(ba), envelope(bb), idx(-4.0), idx(14.0), -Lc, Lc)
    if lc is None:
        res["reason"] = "coarse alignment failed"
        return res
    Lf = int(round(FINE_LAG_S * FS))
    lf, cc_fine = lag_scan(ba, bb, idx(-2.0), idx(4.0), lc - Lf, lc + Lf)
    if lf is None:
        res["reason"] = "fine alignment failed"
        return res
    onset, trig = pick_onset(0.5 * (ba + shift(bb, lf)))
    Lr = int(round(REF_LAG_S * FS))
    r0, r1 = onset + int(round(DIRECT[0] * FS)), onset + int(round(DIRECT[1] * FS))
    lag, cc_ref = lag_scan(ba, bb, r0, r1, lf - Lr, lf + Lr)
    if lag is None:
        res["reason"] = "direct alignment failed"
        return res
    bb2 = shift(bb, lag)
    h1, h2 = common[0::2], common[1::2]
    a1, a2 = make_beam(Sa, h1), make_beam(Sa, h2)
    b1, b2 = shift(make_beam(Sb, h1), lag), shift(make_beam(Sb, h2), lag)
    n0, n1 = idx(NOISE[0]), idx(NOISE[1])
    na, nbb = rms(ba[n0:n1]), rms(bb2[n0:n1])
    snr_da, snr_db = rms(ba[r0:r1]) / na, rms(bb2[r0:r1]) / nbb
    ok_direct = (cc_ref >= PAIR_CC_MIN and snr_da >= SNR_DIRECT_MIN and snr_db >= SNR_DIRECT_MIN)
    res.update(lag_s=lag / FS, cc_env=cc_env, cc_fine=cc_fine, cc_ref=cc_ref,
               onset_s=(onset - idx(0.0)) / FS, onset_triggered=trig, snr_direct_a=snr_da,
               snr_direct_b=snr_db, ok_direct=ok_direct)
    if not ok_direct:
        res["reason"] = "direct gate (cc_ref or SNR)"
    for key, (w0, w1) in WINDOWS:
        i0, i1 = onset + int(round(w0 * FS)), onset + int(round(w1 * FS))
        d = 1.0 - ncc(ba[i0:i1], bb2[i0:i1])
        dha = 1.0 - ncc(a1[i0:i1], a2[i0:i1])
        dhb = 1.0 - ncc(b1[i0:i1], b2[i0:i1])
        floor = 0.5 * (dha + dhb)
        sa, sb = rms(ba[i0:i1]) / na, rms(bb2[i0:i1]) / nbb
        ok = ok_direct and sa >= SNR_CODA_MIN and sb >= SNR_CODA_MIN
        res[key] = d
        res["floor_" + key] = floor
        res["half_a_" + key] = dha
        res["half_b_" + key] = dhb
        res["snr_a_" + key] = sa
        res["snr_b_" + key] = sb
        res["ok_" + key] = ok
        res["label_" + key] = label_pair(ok, d, floor)
    return res


# ----------------------------------------------------------------------------- data
def get_zips(data_dir, s):
    meta = None
    txt = fetch_text(ZENODO_API)
    if txt:
        try:
            meta = {f["key"]: f for f in json.loads(txt).get("files", [])}
        except Exception:  # noqa: BLE001
            meta = None
    s["zenodo_api_ok"] = meta is not None
    paths, checks = {}, {}
    zdir = os.path.join(data_dir, "zw33")
    os.makedirs(zdir, exist_ok=True)
    for key in ZIPS:
        want = None
        if meta and key in meta:
            want = str(meta[key].get("checksum", "")).replace("md5:", "") or None
        cands = [os.path.join(zdir, key), os.path.join(data_dir, key)]
        got = None
        for p in cands:
            if os.path.exists(p) and (want is None or md5_of(p) == want):
                got = p
                break
        if got is None:
            p = cands[0]
            url = "%s/files/%s/content" % (ZENODO_API, key)
            if download(url, p) and (want is None or md5_of(p) == want):
                got = p
        checks[key] = dict(path=got, md5_expected=want,
                           md5_ok=(got is not None and want is not None))
        if got:
            paths[key] = got
    s["zips"] = checks
    return paths


def extract(paths, work):
    folders = {}
    for key, p in paths.items():
        dest = os.path.join(work, key.replace(".zip", ""))
        if not os.path.isdir(dest):
            with zipfile.ZipFile(p) as z:
                z.extractall(dest)
        evs = {}
        for root, dirs, _files in os.walk(dest):
            for d in dirs:
                pe = parse_event_folder(d)
                if pe:
                    evs[pe["date"]] = dict(path=os.path.join(root, d), **pe)
        folders[key] = evs
    return folders


def resolve_event(pe, folder, cache):
    from obspy import UTCDateTime, read
    if pe["iso"] in cache:
        return cache[pe["iso"]]
    t = iso_to_epoch(pe["iso"])
    url = COMCAT % (epoch_to_iso(t - 90), epoch_to_iso(t + 90))
    c = pick_comcat(fetch_text(url, tries=3), pe["iso"], max_dt=60.0)
    if c:
        ev = dict(origin=UTCDateTime(c["time"]), lat=c["lat"], lon=c["lon"], depth=c["depth"],
                  source="comcat", comcat_id=c["id"], mag=c["mag"], dt_folder_s=c["dt"])
    else:
        ev = None
        for fn in sorted(os.listdir(folder)):
            if parse_filename(fn):
                try:
                    sac = read(os.path.join(folder, fn), headonly=True)[0].stats.get("sac", {})
                    la, lo = sac_value(sac, "evla"), sac_value(sac, "evlo")
                    dp = sac_value(sac, "evdp")
                    if la is not None and lo is not None:
                        ev = dict(origin=UTCDateTime(pe["iso"]), lat=la, lon=lo,
                                  depth=dp if dp is not None else DEFAULT_EVENT[2],
                                  source="sac-header+folder-time")
                        break
                except Exception:  # noqa: BLE001
                    pass
        if ev is None:
            ev = dict(origin=UTCDateTime(pe["iso"]), lat=DEFAULT_EVENT[0], lon=DEFAULT_EVENT[1],
                      depth=DEFAULT_EVENT[2], source="default")
    cache[pe["iso"]] = ev
    return ev


# ----------------------------------------------------------------------------- main
def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--data-dir", default="/data")
    ap.add_argument("--out-dir", default="/out")
    ap.add_argument("--work", default="/tmp/zw33")
    args = ap.parse_args(argv)
    import closure35
    if (closure35.VERSION != "0.1.0" or closure35.SAME_MAX_RATIO != SAME_MAX_RATIO
            or closure35.CHANGED_MIN_RATIO != CHANGED_MIN_RATIO
            or closure35.FLOOR_MIN != FLOOR_MIN):
        print("closure35 0.1.0 with the pre-registered constants is required", flush=True)
        return 2
    os.makedirs(args.out_dir, exist_ok=True)
    os.makedirs(args.work, exist_ok=True)
    t_start = time.time()
    s = dict(version=VERSION, closure35_version=closure35.VERSION,
             started=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
             constants=dict(fs=FS, band=BAND, bands_report=BANDS_REPORT, direct=DIRECT,
                            windows=dict(WINDOWS), primary=PRIMARY, noise=NOISE,
                            pair_cc_min=PAIR_CC_MIN, snr_direct_min=SNR_DIRECT_MIN,
                            snr_coda_min=SNR_CODA_MIN, min_common=MIN_COMMON,
                            recovery_ratio_max=RECOVERY_RATIO_MAX,
                            same_max_ratio=closure35.SAME_MAX_RATIO,
                            changed_min_ratio=closure35.CHANGED_MIN_RATIO,
                            floor_min=closure35.FLOOR_MIN))
    s["prereg_public_at_start"] = bool(fetch_text(PREREG_URL, tries=2, wait=5))
    pair_rows, event_rows = [], []
    try:
        from obspy.taup import TauPyModel
        model = TauPyModel("ak135")
        paths = get_zips(args.data_dir, s)
        folders = extract(paths, args.work)
        s["event_folders"] = {k: sorted(v) for k, v in folders.items()}
        cache = {}
        plan = (("M0324", "Data_YKAR.zip", ("YKA",)), ("M0324", "Data_ILAR.zip", ("ILAR",)),
                ("F", "Data_multiF.zip", ("YKA", "ILAR")))
        for band in (BAND,) + tuple(BANDS_REPORT):
            blabel = "%.1f-%.1f" % band
            for set_name, key, arrays in plan:
                evs = folders.get(key, {})
                for arr in arrays:
                    stacks = {}
                    for date in sorted(evs):
                        pe = evs[date]
                        ev = resolve_event(pe, pe["path"], cache)
                        tr, meta = process_folder(pe["path"], ev, arr, band, model)
                        st, info = event_stack(tr)
                        row = dict(set=set_name, band=blabel, array=arr, date=date,
                                   folder=os.path.basename(pe["path"]),
                                   event_source=ev.get("source"), comcat_id=ev.get("comcat_id"),
                                   lat=ev["lat"], lon=ev["lon"], depth=ev["depth"],
                                   origin=str(ev["origin"]), usable=len(st) >= MIN_COMMON)
                        row.update(info)
                        row["element_notes"] = " | ".join("%s: %s" % kv for kv in sorted(meta.items()))
                        event_rows.append(row)
                        if st:
                            stacks[date] = st
                    dates = sorted(stacks)
                    for i in range(len(dates)):
                        for j in range(i + 1, len(dates)):
                            a, b = dates[i], dates[j]
                            try:
                                r = measure_pair(stacks[a], stacks[b])
                            except Exception as ex:  # noqa: BLE001
                                traceback.print_exc()
                                r = dict(ok_direct=False, reason="exception " + repr(ex)[:120])
                            row = dict(set=set_name, band=blabel, array=arr, a=a, b=b)
                            for key2, _w in WINDOWS:
                                r.setdefault("label_" + key2, "nd")
                            row.update(r)
                            pair_rows.append(row)
                            print("PAIR", set_name, blabel, arr, a, b, r.get("label_" + PRIMARY),
                                  round(r.get(PRIMARY, NAN), 4) if isf(r.get(PRIMARY)) else None,
                                  flush=True)
    except Exception as ex:  # noqa: BLE001
        traceback.print_exc()
        s["error"] = repr(ex)[:400]

    prim = "%.1f-%.1f" % BAND
    results = {}
    for window, _w in WINDOWS:
        lab_m = {(r["array"], r["a"], r["b"]): r.get("label_" + window, "nd")
                 for r in pair_rows if r["set"] == "M0324" and r["band"] == prim}
        dv_m = {(r["array"], r["a"], r["b"]): r.get(window) for r in pair_rows
                if r["set"] == "M0324" and r["band"] == prim and r.get("ok_" + window)}
        e1_rows, e1 = e1_grade(lab_m, dv_m)
        lab_f = {(r["array"], r["a"], r["b"]): r.get("label_" + window, "nd")
                 for r in pair_rows if r["set"] == "F" and r["band"] == prim}
        present = {(r["array"], r["date"]): bool(r["usable"]) for r in event_rows
                   if r["set"] == "F" and r["band"] == prim}
        p35 = p35_grade(lab_f, present)
        trip_f, clos_f = closure35.grade(closure35.build_labels(
            closure_rows(pair_rows, "F", prim, window)))
        trip_m, clos_m = closure35.grade(closure35.build_labels(
            closure_rows(pair_rows, "M0324", prim, window)))
        results[window] = dict(graded=(window == PRIMARY), E1=e1, P35=p35, closure_F=clos_f,
                               P35b_closure_M0324=clos_m)
        if window == PRIMARY:
            write_csv(os.path.join(args.out_dir, "zw33_e1_items.csv"), e1_rows)
            write_csv(os.path.join(args.out_dir, "zw33_closure_F_triples.csv"), trip_f)
            write_csv(os.path.join(args.out_dir, "zw33_closure_M0324_triples.csv"), trip_m)
            write_csv(os.path.join(args.out_dir, "zw33_closure35_input_F.csv"),
                      closure_rows(pair_rows, "F", prim, window))
    rep = {}
    for band in BANDS_REPORT:
        bl = "%.1f-%.1f" % band
        lab_f = {(r["array"], r["a"], r["b"]): r.get("label_" + PRIMARY, "nd")
                 for r in pair_rows if r["set"] == "F" and r["band"] == bl}
        present = {(r["array"], r["date"]): bool(r["usable"]) for r in event_rows
                   if r["set"] == "F" and r["band"] == bl}
        lab_m = {(r["array"], r["a"], r["b"]): r.get("label_" + PRIMARY, "nd")
                 for r in pair_rows if r["set"] == "M0324" and r["band"] == bl}
        dv_m = {(r["array"], r["a"], r["b"]): r.get(PRIMARY) for r in pair_rows
                if r["set"] == "M0324" and r["band"] == bl and r.get("ok_" + PRIMARY)}
        rep[bl] = dict(E1=e1_grade(lab_m, dv_m)[1], P35=p35_grade(lab_f, present))
    s["results"] = results
    s["report_only_bands"] = rep
    s["headline"] = dict(E1=results.get(PRIMARY, {}).get("E1", {}).get("headline"),
                         P35=results.get(PRIMARY, {}).get("P35", {}).get("grade"),
                         P35b=results.get(PRIMARY, {}).get("P35b_closure_M0324", {}).get("headline"))
    s["n_pairs"] = len(pair_rows)
    s["minutes"] = round((time.time() - t_start) / 60.0, 1)
    s["finished"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    write_csv(os.path.join(args.out_dir, "zw33_pairs.csv"), pair_rows)
    write_csv(os.path.join(args.out_dir, "zw33_events.csv"), event_rows)
    with open(os.path.join(args.out_dir, "zw33_summary.json"), "w", encoding="utf-8") as f:
        json.dump(jsonable(s), f, indent=1)
    print(json.dumps(jsonable(s["headline"])), flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())

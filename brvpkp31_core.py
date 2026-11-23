# brvpkp31.py 0.1.0 -- pre-registration 31 of inner-core-repeaters:
# is the outer-core change reported from nuclear-test PKP in the Borovoye waveforms, and does PPP-PP
# move with it?
#
# Claim under test: Y. Zhou (2026), JGR Solid Earth 131, e2026JB033798, doi:10.1029/2026JB033798.
# Plan (fixed before any waveform was looked at): findings/2026-11-23-prereg31-brv-waveform-plan.md.
#
# For every pair of French tests (Zhou's Table S1, pairs within 0.05 deg) this job measures, on the
# public LDEO Borovoye CSS 3.0 records, the change of the first core arrival (PKIKP window) relative to
# PP, dK, and of PPP relative to PP, dP3. Year terms (1977 = 0) of dK are graded against Zhou's epochs;
# the ratio of PPP-PP to PKP-PP separates an outer-core change (0), a test-position error (r_pos) and a
# recording time-base error (r_tb).
#
# Part of munimtechnologies/inner-core-repeaters. Licence: Apache-2.0.
# Made with AI agents working in MV2's name (MV2, Munim, Inc.); no human expert has checked it.
import argparse
import calendar
import csv
import hashlib
import io
import json
import math
import os
import re
import subprocess
import sys
import tarfile
import time
from urllib.parse import urljoin

import numpy as np

VERSION = "0.1.0"
BRVK = (53.0581, 70.2828)
REF = (-21.835, -138.91)
SRC_DEPTH_KM = 0.6
REF_YEAR = 1977
EPOCHS = {"E1982_83": (1982, 1983), "E1988_90": (1988, 1990), "E1995": (1995, 1995)}
ZHOU = {"E1982_83": -0.125, "E1988_90": -0.15, "E1995": 0.175}
PARAMS = dict(
    pair_tol_deg=0.05,
    band_hz=(0.3, 1.0),
    fs_grid=100.0,
    win_K=(-3.0, 5.0),        # around first PKIKP
    win_PP=(-3.0, 7.0),
    win_PPP=(-3.0, 7.0),
    win_noise=(-40.0, -10.0),  # around first PKIKP
    win_coarse=(-10.0, 30.0),  # around first PKIKP, envelopes
    coarse_max_s=8.0,
    coarse_fs=10.0,
    env_smooth_s=1.0,
    fine_max_s=0.6,
    cover_before_K_s=50.0,
    cover_after_PP_s=30.0,
    cover_after_PPP_s=30.0,
    cc_min_K=0.80,
    cc_min_PP=0.80,
    cc_min_PPP=0.70,
    snr_min=3.0,
    min_samprate=4.0,
    spike_k=20.0,
    spike_width=7,
    spike_max_run=3,
    clip_frac=0.98,
    clip_run=3,
    sign_min_s=0.02,
    repro_rms_s=0.05,
    partly_frac=0.75,
    chi2_95_3dof=7.815,
    var_floor_s=0.005,
    n_boot=2000,
    boot_seed=31,
    g0_min_pairs=10,
    g1_min_pairs=10,
    min_boot_valid=20,
)
TAG = "v0.31.0"
REPO_RAW = "https://raw.githubusercontent.com/munimtechnologies/inner-core-repeaters/"
S1_PATH = "data/zhou2026_table_s1.csv"
PLAN_PATH = "findings/2026-11-23-prereg31-brv-waveform-plan.md"
INDEX_URLS = [
    "https://www.ldeo.columbia.edu/res/pi/Monitoring/Data/Brv_arch_ex/brv_74-96_ex.html",
    "http://www.ldeo.columbia.edu/res/pi/Monitoring/Data/Brv_arch_ex/brv_74-96_ex.html",
]
FTP_PATTERN = "ftp://ftp.ldeo.columbia.edu/archive/wykim/Brv_archive/Exp/{ymd}-41-{sys}.tar.Z"
SYSTEMS = ("TS", "SS")          # TS = STsR-TSG (primary), SS = STsR-SS (report-only replica)
PRIMARY_SYSTEM = "TS"
UA = {"User-Agent": "inner-core-repeaters/brvpkp31 (research)"}
NAN = float("nan")


def log(*a):
    print(time.strftime("%H:%M:%S"), "[brvpkp31]", *a, flush=True)


# ================================================================ pure helpers (unit-tested)
def utc_seconds(date, tod):
    """'1977-11-24', '16:59:59.92' -> POSIX seconds (UTC), or None."""
    m = re.match(r"^\s*(\d{1,2}):(\d{2}):(\d{1,2}(?:\.\d*)?)\s*$", str(tod or ""))
    d = re.match(r"^(\d{4})-(\d{2})-(\d{2})$", str(date or "").strip())
    if not m or not d:
        return None
    y, mo, dd = (int(v) for v in d.groups())
    return (calendar.timegm((y, mo, dd, 0, 0, 0, 0, 0, 0)) + int(m.group(1)) * 3600
            + int(m.group(2)) * 60 + float(m.group(3)))


def parse_s1(text):
    """Zhou Table S1 CSV (Event Date[*],Event Time,Latitude,Longitude,Magnitude) -> list of dicts."""
    rows = [r for r in csv.reader(io.StringIO(str(text).strip())) if r]
    if not rows:
        return []
    head = [h.strip().lower() for h in rows[0]]

    def idx(name):
        for k, h in enumerate(head):
            if name in h:
                return k
        return None

    i_d, i_t, i_la, i_lo, i_m = idx("date"), idx("time"), idx("lat"), idx("lon"), idx("mag")
    if None in (i_d, i_t, i_la, i_lo):
        return []
    need = max(i_d, i_t, i_la, i_lo)
    out = []
    for r in rows[1:]:
        if len(r) <= need:
            continue
        d = r[i_d].strip()
        star = d.endswith("*")
        d = d.rstrip("*").strip()
        if not re.match(r"^\d{4}-\d{2}-\d{2}$", d):
            continue
        try:
            la, lo = float(r[i_la]), float(r[i_lo])
        except ValueError:
            continue
        mag = None
        if i_m is not None and i_m < len(r):
            try:
                mag = float(r[i_m])
            except ValueError:
                mag = None
        t = r[i_t].strip()
        out.append(dict(date=d, ymd=d.replace("-", ""), year=int(d[:4]), time=t,
                        origin=utc_seconds(d, t), lat=la, lon=lo, mag=mag, starred=star))
    return out


def gc_dist(lat1, lon1, lat2, lon2):
    p1, p2 = math.radians(lat1), math.radians(lat2)
    dl = math.radians(lon2 - lon1)
    c = math.sin(p1) * math.sin(p2) + math.cos(p1) * math.cos(p2) * math.cos(dl)
    return math.degrees(math.acos(max(-1.0, min(1.0, c))))


def form_pairs(events, tol_deg):
    """All (i, j), i < j, with epicentres within tol_deg (Zhou's rule)."""
    out = []
    for i in range(len(events)):
        for j in range(i + 1, len(events)):
            if gc_dist(events[i]["lat"], events[i]["lon"], events[j]["lat"], events[j]["lon"]) <= tol_deg + 1e-9:
                out.append((i, j))
    return out


def parse_index_links(html, base_url):
    """Absolute URLs of every .tar.Z / .tar.gz / .tar link on an index page."""
    out = []
    for m in re.finditer(r"""href\s*=\s*["']?([^"'\s>]+\.tar(?:\.Z|\.gz)?)""", str(html), re.I):
        u = urljoin(base_url, m.group(1))
        if u not in out:
            out.append(u)
    return out


def basename_of(url):
    return url.split("?")[0].rstrip("/").split("/")[-1]


def system_of(name):
    """'19771124-41-TS.tar.Z' -> 'TS'."""
    b = basename_of(name).split(".")[0]
    return b.split("-")[-1].upper() if "-" in b else ""


def links_for(links, ymd, sysname):
    return [u for u in links if basename_of(u).startswith(ymd) and system_of(u) == sysname]


def url_variants(url):
    """The URL itself, then the same host/path over https and http (ftp mirrors are often web-served)."""
    out = [url]
    m = re.match(r"^(ftp|https?)://(.+)$", url)
    if m:
        for sch in ("https", "http", "ftp"):
            v = sch + "://" + m.group(2)
            if v not in out:
                out.append(v)
    return out


# ---------------------------------------------------------------- CSS 3.0 wfdisc
WF_FIELDS = [("sta", 0, 6), ("chan", 7, 15), ("time", 16, 33), ("wfid", 34, 42), ("chanid", 43, 51),
             ("jdate", 52, 60), ("endtime", 61, 78), ("nsamp", 79, 87), ("samprate", 88, 99),
             ("calib", 100, 116), ("calper", 117, 133), ("instype", 134, 140), ("segtype", 141, 142),
             ("datatype", 143, 145), ("clip", 146, 147), ("dir", 148, 212), ("dfile", 213, 245),
             ("foff", 246, 256)]
WF_FORMAT = ("%-6s %-8s %17.5f %8d %8d %8d %17.5f %8d %11.7f %16.6f %16.6f %-6s %-1s %-2s %-1s %-64s "
             "%-32s %10d %8d %-17s")
SAMPLE_SIZE = {"s4": 4, "s3": 3, "s2": 2, "i4": 4, "i2": 2, "t4": 4, "t8": 8, "f4": 4, "f8": 8}
NP_DTYPE = {"s4": ">i4", "s2": ">i2", "i4": "<i4", "i2": "<i2", "t4": ">f4", "t8": ">f8", "f4": "<f4",
            "f8": "<f8"}


def format_wfdisc_line(r):
    return WF_FORMAT % (r["sta"], r["chan"], r["time"], r.get("wfid", 1), r.get("chanid", -1),
                        r.get("jdate", -1), r["endtime"], r["nsamp"], r["samprate"], r.get("calib", 1.0),
                        r.get("calper", 1.0), r.get("instype", "-"), r.get("segtype", "-"), r["datatype"],
                        r.get("clip", "-"), r.get("dir", "."), r["dfile"], r.get("foff", 0), -1, "-")


def _wf_types(r):
    out = dict(r)
    out["time"] = float(r["time"])
    out["endtime"] = float(r["endtime"])
    out["nsamp"] = int(r["nsamp"])
    out["samprate"] = float(r["samprate"])
    out["foff"] = int(r["foff"])
    try:
        out["calib"] = float(r["calib"])
    except (ValueError, KeyError):
        out["calib"] = NAN
    out["datatype"] = str(r["datatype"]).strip().lower()
    if out["samprate"] <= 0 or out["nsamp"] <= 0 or not str(r["dfile"]).strip():
        raise ValueError("bad wfdisc row")
    return out


def parse_wfdisc_line(line):
    """One wfdisc line -> dict, fixed CSS 3.0 columns first, whitespace tokens as fallback; or None."""
    line = str(line).rstrip("\r\n")
    if not line.strip():
        return None
    if len(line) >= 256:
        try:
            return _wf_types({k: line[a:b].strip() for k, a, b in WF_FIELDS})
        except (ValueError, KeyError):
            pass
    p = line.split()
    if len(p) < 18:
        return None
    try:
        return _wf_types(dict(zip([k for k, _, _ in WF_FIELDS], p[:18])))
    except (ValueError, KeyError):
        return None


def decode_samples(raw, datatype):
    dt = str(datatype).strip().lower()
    if dt == "s3":
        n = len(raw) // 3
        b = np.frombuffer(raw[:3 * n], dtype=np.uint8).reshape(-1, 3).astype(np.int64)
        v = (b[:, 0] << 16) | (b[:, 1] << 8) | b[:, 2]
        v = np.where(v >= (1 << 23), v - (1 << 24), v)
        return v.astype(np.float64)
    if dt not in NP_DTYPE:
        raise ValueError("unsupported datatype %r" % dt)
    size = SAMPLE_SIZE[dt]
    n = len(raw) // size
    return np.frombuffer(raw[:size * n], dtype=NP_DTYPE[dt]).astype(np.float64)


def find_file(root, wf_dir, rec):
    """Path of a wfdisc row's data file: <wfdisc dir>/<dir>/<dfile>, <wfdisc dir>/<dfile>, or the first
    file of that name anywhere under root."""
    d, f = str(rec.get("dir", "") or "").strip(), str(rec["dfile"]).strip()
    cands = []
    if d and d != "-" and not os.path.isabs(d):
        cands.append(os.path.join(wf_dir, d, f))
    cands.append(os.path.join(wf_dir, f))
    for c in cands:
        if os.path.isfile(c):
            return c
    for dp, _, fns in os.walk(root):
        if f in fns:
            return os.path.join(dp, f)
    return None


def read_wf(root, wf_dir, rec):
    """-> (t0, dt, samples, note)."""
    path = find_file(root, wf_dir, rec)
    if path is None:
        raise IOError("data file %s not found" % rec["dfile"])
    dt = rec["datatype"]
    size = SAMPLE_SIZE.get(dt)
    if size is None:
        raise ValueError("unsupported datatype %r" % dt)
    with open(path, "rb") as f:
        f.seek(rec["foff"])
        raw = f.read(size * rec["nsamp"])
    x = decode_samples(raw, dt)
    note = "" if len(x) == rec["nsamp"] else "short: %d of %d samples" % (len(x), rec["nsamp"])
    return rec["time"], 1.0 / rec["samprate"], x, note


def safe_extract(tar_path, dest):
    """Extract regular files only, with leading slashes stripped and no path escaping dest."""
    n = 0
    with tarfile.open(tar_path) as tf:
        for m in tf.getmembers():
            if not m.isfile():
                continue
            name = os.path.normpath(m.name.replace("\\", "/").lstrip("/"))
            if name.startswith("..") or os.path.isabs(name) or name in (".", ""):
                continue
            out = os.path.join(dest, name)
            os.makedirs(os.path.dirname(out) or dest, exist_ok=True)
            src = tf.extractfile(m)
            if src is None:
                continue
            with open(out, "wb") as f:
                f.write(src.read())
            n += 1
    return n


# ---------------------------------------------------------------- signal helpers
def despike(x, k=PARAMS["spike_k"], width=PARAMS["spike_width"], max_run=PARAMS["spike_max_run"]):
    """Replace short glitches (runs <= max_run samples beyond k robust sigmas of the residual from a
    running median) by the running median. -> (cleaned copy, number of samples replaced)."""
    x = np.asarray(x, float).copy()
    if len(x) < width:
        return x, 0
    h = width // 2
    pad = np.pad(x, h, mode="edge")
    med = np.median(np.lib.stride_tricks.sliding_window_view(pad, width), axis=1)
    r = x - med
    mad = 1.4826 * float(np.median(np.abs(r - np.median(r))))
    if not mad > 0:
        return x, 0
    idx = np.flatnonzero(np.abs(r) > k * mad)
    count = 0
    if idx.size:
        for run in np.split(idx, np.flatnonzero(np.diff(idx) > 1) + 1):
            if len(run) <= max_run:
                x[run] = med[run]
                count += len(run)
    return x, count


def clipped_mask(x, frac=PARAMS["clip_frac"], run=PARAMS["clip_run"]):
    """True where a run of >= run identical samples sits at >= frac of the trace's max |x|."""
    x = np.asarray(x, float)
    mask = np.zeros(len(x), bool)
    if not x.size:
        return mask
    m = float(np.max(np.abs(x)))
    if not m > 0:
        return mask
    hi = np.abs(x) >= frac * m
    n, i = len(x), 0
    while i < n:
        if hi[i]:
            j = i
            while j + 1 < n and hi[j + 1] and x[j + 1] == x[i]:
                j += 1
            if j - i + 1 >= run:
                mask[i:j + 1] = True
            i = j + 1
        else:
            i += 1
    return mask


def bandpass(x, dt, band=PARAMS["band_hz"]):
    from scipy.signal import butter, sosfiltfilt, detrend
    x = detrend(np.asarray(x, float), type="linear")
    n = len(x)
    m = max(1, min(int(0.05 * n), int(10.0 / dt)))
    w = np.ones(n)
    ramp = 0.5 * (1.0 - np.cos(np.pi * np.arange(m) / m))
    w[:m] = ramp
    w[n - m:] = ramp[::-1]
    nyq = 0.5 / dt
    sos = butter(4, [band[0] / nyq, min(band[1] / nyq, 0.95)], btype="band", output="sos")
    return sosfiltfilt(sos, x * w)


def envelope(x):
    x = np.asarray(x, float)
    n = len(x)
    if n == 0:
        return x
    h = np.zeros(n)
    if n % 2 == 0:
        h[0] = h[n // 2] = 1.0
        h[1:n // 2] = 2.0
    else:
        h[0] = 1.0
        h[1:(n + 1) // 2] = 2.0
    return np.abs(np.fft.ifft(np.fft.fft(x) * h))


def smooth(x, k):
    k = max(1, int(k))
    return np.convolve(np.asarray(x, float), np.ones(k) / k, mode="same")


def seg(tr, t_start, n, fs):
    """Samples of trace tr = (t0, dt, x) at t_start + i/fs, i < n (linear interpolation), or None if the
    span is not inside the record."""
    t0, dt, x = tr
    t_end = t_start + (n - 1) / fs
    if len(x) < 2 or t_start < t0 or t_end > t0 + (len(x) - 1) * dt:
        return None
    i0 = max(0, int(math.floor((t_start - t0) / dt)) - 1)
    i1 = min(len(x), int(math.ceil((t_end - t0) / dt)) + 2)
    xs = t0 + np.arange(i0, i1) * dt
    return np.interp(t_start + np.arange(n) / fs, xs, x[i0:i1])


def ncc(a, b):
    a = a - a.mean()
    b = b - b.mean()
    d = math.sqrt(float(np.dot(a, a)) * float(np.dot(b, b)))
    return float(np.dot(a, b)) / d if d > 0 else NAN


def lag_cc(tra, trb, wa0, wb0, length, maxlag, fs):
    """Lag (s) of record B relative to A: B's window starts at wb0 + lag. Positive = B later."""
    n = int(round(length * fs))
    a = seg(tra, wa0, n, fs)
    if a is None:
        return dict(lag=NAN, cc=NAN, edge=True, ok=False, why="A window outside record")
    m = int(round(maxlag * fs))
    ccs = np.full(2 * m + 1, NAN)
    for k in range(-m, m + 1):
        b = seg(trb, wb0 + k / fs, n, fs)
        if b is not None:
            ccs[k + m] = ncc(a, b)
    if not np.isfinite(ccs).any():
        return dict(lag=NAN, cc=NAN, edge=True, ok=False, why="B window outside record")
    i = int(np.nanargmax(ccs))
    lag, cc = (i - m) / fs, float(ccs[i])
    if 0 < i < 2 * m and np.isfinite(ccs[i - 1]) and np.isfinite(ccs[i + 1]):
        y0, y1, y2 = float(ccs[i - 1]), float(ccs[i]), float(ccs[i + 1])
        den = y0 - 2.0 * y1 + y2
        if den < 0:
            off = 0.5 * (y0 - y2) / den
            lag += off / fs
            cc = y1 - 0.25 * (y0 - y2) * off
    return dict(lag=lag, cc=cc, edge=bool(i <= 0 or i >= 2 * m), ok=True, why="")


def rms_at(tr, t_start, length, fs):
    s = seg(tr, t_start, int(round(length * fs)), fs)
    return None if s is None else float(np.sqrt(np.mean(s * s)))


def measure_pair(A, B, envA, envB, tA, tB, P=PARAMS):
    """A, B: filtered traces (t0, dt, x); envA, envB: smoothed envelopes; tA, tB: absolute predicted
    times {'PKIKP','PP','PPP'}. -> dict with coarse shift c, per-phase tau/cc/edge, dK, dP3."""
    fs = P["fs_grid"]
    w0, w1 = P["win_coarse"]
    co = lag_cc(envA, envB, tA["PKIKP"] + w0, tB["PKIKP"] + w0, w1 - w0, P["coarse_max_s"], P["coarse_fs"])
    out = dict(c=co["lag"], cc_c=co["cc"], edge_c=co["edge"], ok_c=co["ok"])
    if not co["ok"] or co["edge"]:
        out["why"] = "coarse " + (co["why"] or "edge")
        return out
    c = co["lag"]
    for name, ph, win in (("K", "PKIKP", P["win_K"]), ("PP", "PP", P["win_PP"]), ("PPP", "PPP", P["win_PPP"])):
        r = lag_cc(A, B, tA[ph] + win[0], tB[ph] + win[0] + c, win[1] - win[0], P["fine_max_s"], fs)
        out["tau_" + name] = r["lag"]
        out["cc_" + name] = r["cc"]
        out["edge_" + name] = r["edge"]
        out["ok_" + name] = r["ok"]
    out["dK"] = out["tau_K"] - out["tau_PP"]
    out["dP3"] = out["tau_PPP"] - out["tau_PP"]
    return out


def choose_channel(stats_a, stats_b):
    """stats_x: {chan: {'eligible': bool, 'snr_K': float}} -> shared eligible channel with the largest
    min(SNR_K), or None."""
    best, best_v = None, -1.0
    for ch in sorted(set(stats_a) & set(stats_b)):
        a, b = stats_a[ch], stats_b[ch]
        if not (a.get("eligible") and b.get("eligible")):
            continue
        v = min(a.get("snr_K") or 0.0, b.get("snr_K") or 0.0)
        if v > best_v:
            best, best_v = ch, v
    return best


# ---------------------------------------------------------------- inversion and grades
def pair_year_terms(pairs, years, ref=REF_YEAR):
    """pairs: (i, j, d) with d = a[year_j] - a[year_i] + noise. Least squares, a[ref] = 0; pairs within
    a year skipped; years not linked to ref dropped."""
    use = [(i, j, d) for i, j, d in pairs
           if d is not None and math.isfinite(d) and years[i] != years[j]]
    yrs = {years[i] for i, _, _ in use} | {years[j] for _, j, _ in use}
    if ref not in yrs:
        return {}
    adj = {y: set() for y in yrs}
    for i, j, _ in use:
        adj[years[i]].add(years[j])
        adj[years[j]].add(years[i])
    seen, todo = {ref}, [ref]
    while todo:
        y = todo.pop()
        for z in adj[y]:
            if z not in seen:
                seen.add(z)
                todo.append(z)
    free = sorted(y for y in seen if y != ref)
    col = {y: k for k, y in enumerate(free)}
    A, b = [], []
    for i, j, d in use:
        if years[i] not in seen:
            continue
        row = np.zeros(len(free))
        if years[j] != ref:
            row[col[years[j]]] += 1.0
        if years[i] != ref:
            row[col[years[i]]] -= 1.0
        A.append(row)
        b.append(d)
    out = {ref: 0.0}
    if free and A:
        sol = np.linalg.lstsq(np.array(A), np.array(b), rcond=None)[0]
        for y, k in col.items():
            out[y] = float(sol[k])
    return out


def epoch_terms(terms, event_years):
    """Epoch value = mean of year terms in the epoch, weighted by the number of tests that year."""
    out = {}
    for k, (a, bb) in EPOCHS.items():
        num = den = 0.0
        for y, v in terms.items():
            if a <= y <= bb:
                w = float(sum(1 for z in event_years if z == y))
                num += w * v
                den += w
        out[k] = (num / den) if den > 0 else None
    return out


def epochs_from_pairs(pairs, years):
    yt = pair_year_terms(pairs, years)
    ev = sorted({i for i, _, _ in pairs} | {j for _, j, _ in pairs})
    return yt, epoch_terms(yt, [years[i] for i in ev])


def grade_repro(ep, zhou=ZHOU, P=PARAMS):
    if any(ep.get(k) is None for k in zhou):
        return dict(grade="untestable", signs_match=None, rms_res=None, rms_zhou=None)
    sm = sum(1 for k in zhou if abs(ep[k]) >= P["sign_min_s"]
             and math.copysign(1, ep[k]) == math.copysign(1, zhou[k]))
    rms_res = math.sqrt(sum((zhou[k] - ep[k]) ** 2 for k in zhou) / len(zhou))
    rms_z = math.sqrt(sum(zhou[k] ** 2 for k in zhou) / len(zhou))
    if sm == 3 and rms_res <= P["repro_rms_s"]:
        g = "reproduced"
    elif sm >= 2 and rms_res < P["partly_frac"] * rms_z:
        g = "partly"
    else:
        g = "not-reproduced"
    return dict(grade=g, signs_match=sm, rms_res=rms_res, rms_zhou=rms_z)


def bootstrap_epochs(pairs, years, keys, n, seed):
    """Resample pairs with replacement. keys: names of the value fields in each pair tuple beyond (i, j),
    e.g. pairs = [(i, j, {'dK': .., 'dP3': ..})]. -> list of {key: epochs dict}."""
    rng = np.random.default_rng(seed)
    out = []
    m = len(pairs)
    if m == 0:
        return out
    for _ in range(n):
        pick = rng.integers(0, m, m)
        res = {}
        for key in keys:
            pp = [(pairs[k][0], pairs[k][1], pairs[k][2][key]) for k in pick]
            res[key] = epochs_from_pairs(pp, years)[1]
        out.append(res)
    return out


def hypothesis_chi2(k_ep, p_ep, boot_k, boot_p, ratios, P=PARAMS):
    """chi2_H = sum_e (P3_e - r_H K_e)^2 / var_e; var_e = bootstrap variance of P3_e - r_H K_e
    (floored). None if an epoch is missing or the bootstrap is too thin."""
    out = {}
    for h, r in ratios.items():
        chi = 0.0
        for e in EPOCHS:
            if k_ep.get(e) is None or p_ep.get(e) is None:
                return None
            samp = [bp[e] - r * bk[e] for bk, bp in zip(boot_k, boot_p)
                    if bk.get(e) is not None and bp.get(e) is not None]
            if len(samp) < P["min_boot_valid"]:
                return None
            var = max(float(np.var(samp, ddof=1)), P["var_floor_s"] ** 2)
            res = p_ep[e] - r * k_ep[e]
            chi += res * res / var
        out[h] = chi
    return out


def grade_discriminant(chi, P=PARAMS):
    if chi is None:
        return "untestable"
    ok = [h for h, c in chi.items() if c <= P["chi2_95_3dof"]]
    if len(ok) == 1:
        return ok[0]
    return "none-fits" if not ok else "unresolved"


def slope_through_zero(x, y):
    x, y = np.asarray(x, float), np.asarray(y, float)
    d = float(np.dot(x, x))
    return float(np.dot(x, y)) / d if d > 0 else NAN


def tls_slope(x, y):
    x, y = np.asarray(x, float), np.asarray(y, float)
    if len(x) < 2:
        return NAN
    sxx, syy, sxy = float(np.dot(x, x)), float(np.dot(y, y)), float(np.dot(x, y))
    if sxy == 0:
        return NAN
    return ((syy - sxx) + math.sqrt((syy - sxx) ** 2 + 4 * sxy * sxy)) / (2 * sxy)


# ================================================================ network / IO parts
def http_get(url, timeout=120, verify=True):
    import requests
    return requests.get(url, timeout=timeout, verify=verify, headers=UA)


def repo_text(ref, path):
    r = http_get(REPO_RAW + ref + "/" + path, timeout=60)
    if r.status_code == 404:
        return None
    r.raise_for_status()
    return r.text


def wait_plan(wait_min, poll_s=600):
    t0 = time.time()
    while True:
        try:
            txt = repo_text("main", PLAN_PATH)
            if txt:
                return txt
        except Exception as ex:  # noqa: BLE001
            log("plan fetch error", repr(ex)[:200])
        if time.time() - t0 > wait_min * 60:
            return None
        log("plan not on main yet; waiting")
        time.sleep(poll_s)


def fetch_index(out_dir, summary):
    try:
        import urllib3
        urllib3.disable_warnings()
    except Exception:  # noqa: BLE001
        pass
    for u in INDEX_URLS:
        for verify in (True, False):
            try:
                r = http_get(u, timeout=60, verify=verify)
                if r.status_code == 200 and r.text:
                    with open(os.path.join(out_dir, "brv31_index.html"), "w", encoding="utf-8",
                              errors="replace") as f:
                        f.write(r.text)
                    summary["index"] = dict(url=u, verify=verify, bytes=len(r.content))
                    return parse_index_links(r.text, r.url or u)
                summary.setdefault("index_errors", []).append("%s %s" % (u, r.status_code))
            except Exception as ex:  # noqa: BLE001
                summary.setdefault("index_errors", []).append("%s verify=%s %r" % (u, verify, ex)[:240])
    return []


def download(url, dest):
    """Try url and its scheme variants; -> (path, url_used, error)."""
    if os.path.exists(dest) and os.path.getsize(dest) > 0:
        return dest, "cache", ""
    errs = []
    for u in url_variants(url):
        tmp = dest + ".part"
        try:
            if u.startswith("ftp://"):
                import urllib.request
                with urllib.request.urlopen(u, timeout=180) as r, open(tmp, "wb") as f:
                    while True:
                        b = r.read(1 << 16)
                        if not b:
                            break
                        f.write(b)
            else:
                import requests
                with requests.get(u, timeout=180, verify=False, headers=UA, stream=True) as r:
                    if r.status_code != 200:
                        raise IOError("HTTP %s" % r.status_code)
                    with open(tmp, "wb") as f:
                        for b in r.iter_content(1 << 16):
                            f.write(b)
            if os.path.getsize(tmp) == 0:
                raise IOError("empty")
            os.replace(tmp, dest)
            return dest, u, ""
        except Exception as ex:  # noqa: BLE001
            errs.append("%s: %r" % (u, ex))
            if os.path.exists(tmp):
                os.remove(tmp)
    return None, "", " | ".join(errs)[:600]


def decompress(path):
    """.tar.Z / .tar.gz / .tar -> path of a plain tar."""
    if path.endswith(".Z") or path.endswith(".gz"):
        tar = path.rsplit(".", 1)[0]
        if os.path.exists(tar) and os.path.getsize(tar) > 0:
            return tar
        try:
            with open(tar, "wb") as f:
                subprocess.run(["gzip", "-dc", path], stdout=f, stderr=subprocess.PIPE, check=True)
            if os.path.getsize(tar) > 0:
                return tar
        except Exception:  # noqa: BLE001
            pass
        import unlzw3
        with open(path, "rb") as f:
            data = unlzw3.unlzw(f.read())
        with open(tar, "wb") as f:
            f.write(data)
        return tar
    return path


def load_database(root):
    """All wfdisc rows under root -> list of (wfdisc_dir, rec)."""
    rows = []
    for dp, _, fns in os.walk(root):
        for fn in sorted(fns):
            if fn.lower().endswith(".wfdisc"):
                with open(os.path.join(dp, fn), "r", encoding="latin-1") as f:
                    for line in f:
                        rec = parse_wfdisc_line(line)
                        if rec is not None:
                            rows.append((dp, rec))
    return rows


def make_taup():
    from obspy.taup import TauPyModel
    from obspy.geodetics import locations2degrees
    model = TauPyModel("iasp91")

    def tt(lat, lon):
        d = float(locations2degrees(lat, lon, BRVK[0], BRVK[1]))
        first = {}
        for a in model.get_travel_times(source_depth_in_km=SRC_DEPTH_KM, distance_in_degree=d,
                                        phase_list=["PKIKP", "PKiKP", "PP", "PPP"]):
            if a.name not in first:
                first[a.name] = dict(t=float(a.time), p=float(a.ray_param_sec_degree))
        return d, first
    return tt


def ratios_from(first):
    tK, tPP, tP3 = first["PKIKP"]["t"], first["PP"]["t"], first["PPP"]["t"]
    pK, pPP, pP3 = first["PKIKP"]["p"], first["PP"]["p"], first["PPP"]["p"]
    return {"core": 0.0, "position": (pP3 - pPP) / (pK - pPP), "time-base": (tP3 - tPP) / (tK - tPP)}


# ================================================================ pipeline
def event_traces(ev, root, sysname, P, inv_rows):
    """Read, clean and filter every channel of one test's database -> {chan: info}."""
    t = ev["abs"]
    lo = t["PKIKP"] - P["cover_before_K_s"] - 20.0
    hi = t["PPP"] + P["cover_after_PPP_s"] + 20.0
    by_chan = {}
    for wf_dir, rec in load_database(root):
        row = dict(date=ev["date"], system=sysname, sta=rec["sta"], chan=rec["chan"], datatype=rec["datatype"],
                   samprate=rec["samprate"], nsamp=rec["nsamp"], calib=rec["calib"], instype=rec.get("instype", ""),
                   start_minus_origin_s=rec["time"] - ev["origin"],
                   length_s=rec["nsamp"] / rec["samprate"], dfile=rec["dfile"])
        try:
            t0, dt, x, note = read_wf(root, wf_dir, rec)
            row["read"] = "ok" + ((" " + note) if note else "")
        except Exception as ex:  # noqa: BLE001
            row["read"] = "error: %r" % ex
            inv_rows.append(row)
            continue
        inv_rows.append(row)
        key = (rec["sta"] + "." + rec["chan"]).strip(".")
        t_end = t0 + (len(x) - 1) * dt
        cov = max(0.0, min(t_end, hi) - max(t0, lo))
        prev = by_chan.get(key)
        if prev is not None and prev["cov"] >= cov:
            continue
        by_chan[key] = dict(t0=t0, dt=dt, x=x, cov=cov, samprate=rec["samprate"])
    out = {}
    fs = P["fs_grid"]
    for key, d in by_chan.items():
        t0, dt, x = d["t0"], d["dt"], d["x"]
        i0 = max(0, int((lo - t0) / dt))
        i1 = min(len(x), int((hi - t0) / dt) + 1)
        info = dict(chan=key, vertical=("z" in key.split(".")[-1].lower()), samprate=d["samprate"],
                    n_spikes=None, eligible=False, eligible_P3=False)
        if i1 - i0 < 100:
            info["why"] = "record does not reach the PKP-PPP span"
            out[key] = info
            continue
        xc, nsp = despike(x[i0:i1])
        info["n_spikes"] = nsp
        tr0 = t0 + i0 * dt
        clip = clipped_mask(xc)

        def clipped_in(a, b):
            j0 = max(0, int((a - tr0) / dt))
            j1 = min(len(xc), int(math.ceil((b - tr0) / dt)) + 1)
            return bool(clip[j0:j1].any()) if j1 > j0 else False

        xf = bandpass(xc, dt)
        tr = (tr0, dt, xf)
        env = (tr0, dt, smooth(envelope(xf), round(P["env_smooth_s"] / dt)))
        t_end = tr0 + (len(xc) - 1) * dt
        covK = tr0 <= t["PKIKP"] - P["cover_before_K_s"] and t_end >= t["PP"] + P["cover_after_PP_s"]
        covP3 = covK and t_end >= t["PPP"] + P["cover_after_PPP_s"]
        noise = rms_at(tr, t["PKIKP"] + P["win_noise"][0], P["win_noise"][1] - P["win_noise"][0], fs)
        for name, ph, win in (("K", "PKIKP", P["win_K"]), ("PP", "PP", P["win_PP"]), ("PPP", "PPP", P["win_PPP"])):
            s = rms_at(tr, t[ph] + win[0], win[1] - win[0], fs)
            info["snr_" + name] = (s / noise) if (s is not None and noise) else None
            info["clip_" + name] = clipped_in(t[ph] + win[0] - 1.0, t[ph] + win[1] + 1.0)
        info.update(cover_K=covK, cover_P3=covP3,
                    eligible=bool(covK and not info["clip_K"] and not info["clip_PP"]
                                  and d["samprate"] >= P["min_samprate"]),
                    eligible_P3=bool(covP3 and not info["clip_PPP"]))
        info["_tr"], info["_env"] = tr, env
        out[key] = info
    return out


def analyse_system(events, traces, sysname, pair_idx, P):
    """Measure every pair of pair_idx on one system -> rows."""
    rows = []
    for i, j in pair_idx:
        a, b = events[i], events[j]
        row = dict(system=sysname, a=a["date"], b=b["date"], ya=a["year"], yb=b["year"],
                   sep_km=gc_dist(a["lat"], a["lon"], b["lat"], b["lon"]) * 111.195)
        ta, tb = traces.get(i), traces.get(j)
        if not ta or not tb:
            row["status"] = "no record"
            rows.append(row)
            continue

        def pick(tt):
            v = {k: x for k, x in tt.items() if x.get("vertical")}
            return v if v else tt
        sa, sb = pick(ta), pick(tb)
        ch = choose_channel(sa, sb)
        if ch is None:
            row["status"] = "no shared eligible channel"
            rows.append(row)
            continue
        ca, cb = sa[ch], sb[ch]
        row["chan"] = ch
        for nm in ("K", "PP", "PPP"):
            row["snr_%s_a" % nm], row["snr_%s_b" % nm] = ca.get("snr_" + nm), cb.get("snr_" + nm)
        m = measure_pair(ca["_tr"], cb["_tr"], ca["_env"], cb["_env"], a["abs"], b["abs"], P)
        row.update({k: v for k, v in m.items()})
        if "dK" not in m:
            row["status"] = m.get("why", "coarse failed")
            rows.append(row)
            continue

        def snr_ok(nm):
            return all((x or 0.0) >= P["snr_min"] for x in (ca.get("snr_" + nm), cb.get("snr_" + nm)))
        k_ok = (m["ok_K"] and m["ok_PP"] and not m["edge_K"] and not m["edge_PP"]
                and (m["cc_K"] or 0) >= P["cc_min_K"] and (m["cc_PP"] or 0) >= P["cc_min_PP"]
                and snr_ok("K") and snr_ok("PP") and math.isfinite(m["dK"]))
        p_ok = (k_ok and ca["eligible_P3"] and cb["eligible_P3"] and m["ok_PPP"] and not m["edge_PPP"]
                and (m["cc_PPP"] or 0) >= P["cc_min_PPP"] and snr_ok("PPP") and math.isfinite(m["dP3"]))
        row["K_usable"], row["P3_usable"] = bool(k_ok), bool(p_ok)
        row["status"] = "measured"
        row["_i"], row["_j"] = i, j
        rows.append(row)
    return rows


def grade_rows(rows, years, ratios, P):
    kp = [(r["_i"], r["_j"], r["dK"]) for r in rows if r.get("K_usable")]
    pp = [(r["_i"], r["_j"], {"dK": r["dK"], "dP3": r["dP3"]}) for r in rows if r.get("P3_usable")]
    res = dict(n_K_usable=len(kp), n_P3_usable=len(pp), n_pairs=len(rows),
               n_measured=sum(1 for r in rows if r.get("status") == "measured"))
    ytK, epK = epochs_from_pairs(kp, years)
    res["K_year_terms"] = {str(k): v for k, v in sorted(ytK.items())}
    res["K_epochs"] = epK
    in_1977 = any(years[i] == REF_YEAR or years[j] == REF_YEAR for i, j, _ in kp)
    g0 = len(kp) >= P["g0_min_pairs"] and in_1977 and all(v is not None for v in epK.values())
    res["G0"] = dict(passed=bool(g0), n=len(kp), ref_in_pairs=bool(in_1977))
    res["K"] = grade_repro(epK, P=P) if g0 else dict(grade="untestable")
    bk = bootstrap_epochs([(i, j, {"dK": d}) for i, j, d in kp], years, ["dK"], P["n_boot"], P["boot_seed"])
    res["K_epoch_boot_sd"] = {e: _sd([s["dK"][e] for s in bk]) for e in EPOCHS}
    # discriminant
    ytK3, epK3 = epochs_from_pairs([(i, j, v["dK"]) for i, j, v in pp], years)
    ytP3, epP3 = epochs_from_pairs([(i, j, v["dP3"]) for i, j, v in pp], years)
    res.update(K_in_P3set_epochs=epK3, P3_epochs=epP3, P3_year_terms={str(k): v for k, v in sorted(ytP3.items())})
    g1 = (len(pp) >= P["g1_min_pairs"] and all(v is not None for v in epK3.values())
          and all(v is not None for v in epP3.values()))
    chi = None
    if g1:
        boot = bootstrap_epochs(pp, years, ["dK", "dP3"], P["n_boot"], P["boot_seed"])
        chi = hypothesis_chi2(epK3, epP3, [s["dK"] for s in boot], [s["dP3"] for s in boot], ratios, P)
        res["P3_epoch_boot_sd"] = {e: _sd([s["dP3"][e] for s in boot]) for e in EPOCHS}
    res["G1"] = dict(passed=bool(g1 and chi is not None), n=len(pp))
    res["chi2"] = chi
    res["P3"] = grade_discriminant(chi, P) if g1 else "untestable"
    xs, ys = [v["dK"] for _, _, v in pp], [v["dP3"] for _, _, v in pp]
    res["pair_slope_dP3_on_dK"] = dict(through_zero=slope_through_zero(xs, ys) if xs else None,
                                       tls=tls_slope(xs, ys) if len(xs) > 1 else None)
    return res


def _sd(v):
    v = [x for x in v if x is not None]
    return float(np.std(v, ddof=1)) if len(v) > 2 else None


def _clean(o):
    if isinstance(o, (float, np.floating)):
        o = float(o)
        return o if math.isfinite(o) else None
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.bool_,)):
        return bool(o)
    if isinstance(o, dict):
        return {str(k): _clean(v) for k, v in o.items() if not str(k).startswith("_")}
    if isinstance(o, (list, tuple)):
        return [_clean(v) for v in o]
    return o


def write_csv(path, rows):
    keys = []
    for r in rows:
        for k in r:
            if not str(k).startswith("_") and k not in keys:
                keys.append(k)
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=keys or ["none"], extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow({k: (round(v, 6) if isinstance(v, float) else v) for k, v in r.items()
                        if not str(k).startswith("_")})


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", default="data")
    ap.add_argument("--out-dir", default="out")
    ap.add_argument("--wait-min", type=float, default=0.0)
    ap.add_argument("--time-budget-min", type=float, default=150.0,
                    help="download budget, counted from when the plan is found")
    a = ap.parse_args(argv)
    os.makedirs(a.out_dir, exist_ok=True)
    os.makedirs(a.data_dir, exist_ok=True)
    t_start = time.time()
    P = PARAMS
    summary = dict(version=VERSION, started=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()), errors=[],
                   params=P, zhou_reported=ZHOU, tag_inputs=TAG)
    with open(__file__, "rb") as f:
        summary["code_sha256"] = hashlib.sha256(f.read()).hexdigest()

    def save():
        summary["runtime_min"] = round((time.time() - t_start) / 60.0, 2)
        with open(os.path.join(a.out_dir, "brv31_summary.json"), "w") as f:
            json.dump(_clean(summary), f, indent=1)

    plan = wait_plan(a.wait_min)
    summary["plan_public_at_start"] = bool(plan)
    if not plan:
        summary["status"] = "no-plan"
        summary["headline"] = "not run: plan not public on main"
        save()
        return 0
    summary["plan_sha256"] = hashlib.sha256(plan.encode("utf-8")).hexdigest()
    t_work = time.time()
    s1_text = repo_text(TAG, S1_PATH)
    summary["s1_sha256"] = hashlib.sha256(s1_text.encode("utf-8")).hexdigest()
    events = parse_s1(s1_text)
    summary["n_events"] = len(events)
    summary["n_starred"] = sum(1 for e in events if e["starred"])
    years = [e["year"] for e in events]

    tt = make_taup()
    _, first_ref = tt(*REF)
    ratios = ratios_from(first_ref)
    summary["ratios"] = ratios
    summary["ref_phases"] = first_ref
    for e in events:
        d, first = tt(e["lat"], e["lon"])
        e["dist"] = d
        e["abs"] = {ph: e["origin"] + first[ph]["t"] for ph in ("PKIKP", "PP", "PPP")}
        e["T"] = {ph: first[ph]["t"] for ph in first}
    save()

    links = fetch_index(a.out_dir, summary)
    summary["n_index_links"] = len(links)
    save()

    inv_rows, file_rows = [], []
    traces = {s: {} for s in SYSTEMS}
    for i, e in enumerate(events):
        for sysname in SYSTEMS:
            if time.time() - t_work > a.time_budget_min * 60:
                if "time budget reached during downloads" not in summary["errors"]:
                    summary["errors"].append("time budget reached during downloads")
                break
            urls = links_for(links, e["ymd"], sysname) or [FTP_PATTERN.format(ymd=e["ymd"], sys=sysname)]
            fr = dict(date=e["date"], system=sysname, url=urls[0], n_links=len(urls))
            path, used, err = download(urls[0], os.path.join(a.data_dir, basename_of(urls[0])))
            fr.update(used=used, error=err)
            if path is None:
                fr["status"] = "download failed"
                file_rows.append(fr)
                continue
            fr["bytes"] = os.path.getsize(path)
            try:
                tar = decompress(path)
                root = os.path.join(a.data_dir, "x", basename_of(urls[0]).split(".")[0])
                if not os.path.isdir(root):
                    os.makedirs(root, exist_ok=True)
                    fr["n_files"] = safe_extract(tar, root)
                info = event_traces(e, root, sysname, P, inv_rows)
                fr["channels"] = " ".join(sorted(info))
                fr["status"] = "read"
                if info:
                    traces[sysname][i] = info
            except Exception as ex:  # noqa: BLE001
                fr["status"] = "error: %r" % ex
            file_rows.append(fr)
        save()
    write_csv(os.path.join(a.out_dir, "brv31_files.csv"), file_rows)
    write_csv(os.path.join(a.out_dir, "brv31_inventory.csv"), inv_rows)
    ch_rows = []
    for sysname in SYSTEMS:
        for i, info in sorted(traces[sysname].items()):
            for ch, x in sorted(info.items()):
                r = dict(date=events[i]["date"], system=sysname)
                r.update({k: v for k, v in x.items() if not k.startswith("_")})
                ch_rows.append(r)
    write_csv(os.path.join(a.out_dir, "brv31_channels.csv"), ch_rows)
    summary["n_events_with_records"] = {s: len(traces[s]) for s in SYSTEMS}

    starred_idx = [k for k, e in enumerate(events) if e["starred"]]
    all_pairs = form_pairs(events, P["pair_tol_deg"])
    starred_pairs = [(i, j) for i, j in all_pairs if i in starred_idx and j in starred_idx]
    summary["n_pairs"] = dict(starred=len(starred_pairs), all=len(all_pairs))
    variants = {}
    all_rows = []
    for name, sysname, pidx in (("primary_TS_starred", "TS", starred_pairs), ("TS_all23", "TS", all_pairs),
                                ("SS_starred", "SS", starred_pairs), ("SS_all23", "SS", all_pairs)):
        try:
            rows = analyse_system(events, traces[sysname], sysname, pidx, P)
            for r in rows:
                r["variant"] = name
            all_rows.extend(rows)
            variants[name] = grade_rows(rows, years, ratios, P)
        except Exception as ex:  # noqa: BLE001
            summary["errors"].append("%s: %r" % (name, ex))
            variants[name] = dict(error=repr(ex))
        save()
    write_csv(os.path.join(a.out_dir, "brv31_pairs.csv"), all_rows)
    summary["variants"] = variants
    pr = variants.get("primary_TS_starred", {})
    kg = (pr.get("K") or {}).get("grade", "untestable")
    pg = pr.get("P3", "untestable")
    summary["headline"] = "K: %s; P3: %s" % (kg, pg)
    summary["status"] = "graded"
    summary["finished"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    save()
    log("done", summary["headline"])
    return 0


if __name__ == "__main__":
    sys.exit(main())

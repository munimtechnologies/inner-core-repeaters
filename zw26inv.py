#!/usr/bin/env python3
# Copyright 2026 MV2, Munim, Inc.
# Licensed under the Apache License, Version 2.0 (Apache-2.0). See LICENSE.
# Made with AI agents working in MV2's name.
"""zw26inv.py - header-only inventory of the authors' data for Zhang & Wen (2026),
"Episodic Change, Quiescence, and Recovery of Earth's Inner Core Surface", JGR Solid
Earth, doi 10.1029/2026JB034595. Data: Zenodo doi 10.5281/zenodo.19625711 (CC BY 4.0).

Purpose: learn what the published files contain (events, stations, components, sampling,
time windows, header picks, any text tables) so that a test of the paper's claims can be
pre-registered on waveforms we have not measured. This script makes NO waveform
measurement: it reads headers only (ObsPy headonly=True) and copies small text files.

Outputs (--out-dir):
  zw26inv_summary.json  counts, md5 checks, per-bundle event list, sampling rates, picks
  zw26inv_files.csv     every zip member: bundle, member, bytes, kind
  zw26inv_traces.csv    one row per trace header
  zw26inv_events.csv    one row per (bundle, inferred origin): traces, stations, rates
  text/<bundle>/...     copies of text-like members up to 300 kB each (2 MB in total)

Top-level imports: stdlib only; ObsPy inside functions.
"""

import argparse
import csv
import hashlib
import json
import os
import shutil
import sys
import time
import urllib.request
import zipfile

VERSION = "0.1.0"
RECORD_API = "https://zenodo.org/api/records/19625711"
BUNDLES = ["Data_ILAR.zip", "Data_YKAR.zip", "Data_multiF.zip", "Data_ForRelocation.zip"]
EXPECTED_BYTES = {"Data_multiF.zip": 18967153, "Data_YKAR.zip": 23222725,
                  "Data_ILAR.zip": 16977973, "Data_ForRelocation.zip": 143498389}
TEXT_EXT = (".txt", ".dat", ".csv", ".md", ".list", ".lst", ".out", ".in", ".info",
            ".log", ".json", ".xml", ".m", ".py", ".sh", ".gmt", ".loc", ".cat", ".pick",
            ".tab", ".tsv", ".readme")
TEXT_MAX_BYTES = 300 * 1024
TEXT_TOTAL_MAX = 2 * 1024 * 1024
MEMBER_MAX_BYTES = 400 * 1024 * 1024
SAC_KEYS = ["kstnm", "knetwk", "khole", "kcmpnm", "cmpaz", "cmpinc", "baz", "az", "gcarc",
            "evla", "evlo", "evdp", "mag", "stla", "stlo", "stel", "o", "b", "e", "delta",
            "npts", "kevnm", "a", "ka", "f", "kf",
            "t0", "t1", "t2", "t3", "t4", "t5", "t6", "t7", "t8", "t9",
            "kt0", "kt1", "kt2", "kt3", "kt4", "kt5", "kt6", "kt7", "kt8", "kt9",
            "user0", "user1", "user2", "user3", "kuser0", "kuser1", "kuser2",
            "iztype", "idep", "nzyear", "nzjday", "nzhour", "nzmin", "nzsec", "nzmsec"]
PICK_KEYS = ["a", "f", "t0", "t1", "t2", "t3", "t4", "t5", "t6", "t7", "t8", "t9"]


def log(*a):
    print(time.strftime("%H:%M:%S"), "[zw26inv]", *a, flush=True)


# ----------------------------------------------------------------- pure helpers

def safe_member(name):
    """Normalised relative path of a zip member, or None if unsafe or empty."""
    n = os.path.normpath(name.replace("\\", "/")).replace("\\", "/")
    if n.startswith("..") or n.startswith("/") or ":" in n or n in (".", ""):
        return None
    return n


def kind_of(name):
    """Rough file kind from the name: zip, text, sac, mseed or other."""
    low = name.lower()
    base = os.path.basename(low)
    if low.endswith(".zip"):
        return "zip"
    if low.endswith(TEXT_EXT) or base.startswith("readme"):
        return "text"
    if low.endswith((".sac", ".sacd")) or ".sac" in base:
        return "sac"
    if low.endswith((".mseed", ".miniseed", ".msd")):
        return "mseed"
    return "other"


def origin_key(start_iso_ref, o_sec):
    """Origin label 'YYYY-MM-DDTHH:MM:SS' from a reference-time ISO string
    ('YYYY-MM-DDTHH:MM:SS[.ffffff]Z') plus the SAC 'o' offset in seconds, rounded to the
    nearest second. Returns '' if o_sec is None."""
    if o_sec is None:
        return ""
    import datetime as dt
    s = start_iso_ref.rstrip("Z")
    if "." in s:
        head, frac = s.split(".", 1)
        frac = (frac + "000000")[:6]
        s = head + "." + frac
        t = dt.datetime.strptime(s, "%Y-%m-%dT%H:%M:%S.%f")
    else:
        t = dt.datetime.strptime(s, "%Y-%m-%dT%H:%M:%S")
    t = t + dt.timedelta(seconds=float(o_sec))
    t = t + dt.timedelta(microseconds=500000)
    return t.strftime("%Y-%m-%dT%H:%M:%S")


def summarise_events(rows):
    """Group trace rows by (bundle, origin or folder). Returns a list of dicts sorted by
    bundle then key: n_traces, stations (sorted unique), comps, rates, kevnm."""
    groups = {}
    for r in rows:
        if r.get("status") != "ok":
            continue
        key = r.get("origin") or ("folder:" + os.path.dirname(r.get("member", "")))
        g = groups.setdefault((r.get("bundle", ""), key),
                              dict(n_traces=0, stations=set(), comps=set(), rates=set(),
                                   kevnm=set(), dirs=set()))
        g["n_traces"] += 1
        if r.get("sta"):
            g["stations"].add(r["sta"])
        if r.get("cha"):
            g["comps"].add(r["cha"])
        g["rates"].add(round(float(r.get("sr", 0.0)), 4))
        if r.get("sac_kevnm") not in (None, ""):
            g["kevnm"].add(str(r["sac_kevnm"]))
        g["dirs"].add(os.path.dirname(r.get("member", "")))
    out = []
    for (bundle, key) in sorted(groups):
        g = groups[(bundle, key)]
        out.append(dict(bundle=bundle, origin=key, n_traces=g["n_traces"],
                        n_stations=len(g["stations"]),
                        stations=" ".join(sorted(g["stations"])),
                        channels=" ".join(sorted(g["comps"])),
                        rates=" ".join("%g" % x for x in sorted(g["rates"])),
                        kevnm=" ".join(sorted(g["kevnm"]))[:200],
                        dirs=" ".join(sorted(g["dirs"]))[:300]))
    return out


def md5_file(path):
    h = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


# ----------------------------------------------------------------- I/O steps

def fetch(url, dest, tries=4, timeout=900):
    last = None
    for k in range(tries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "MV2-Munim-research/1.0"})
            with urllib.request.urlopen(req, timeout=timeout) as r, open(dest + ".part", "wb") as f:
                shutil.copyfileobj(r, f, 1 << 20)
            os.replace(dest + ".part", dest)
            return True
        except Exception as ex:  # noqa: BLE001
            last = ex
            log("fetch failed", url, repr(ex)[:120], "try", k + 1)
            time.sleep(10 * (k + 1))
    raise RuntimeError("download failed: %s (%r)" % (url, last))


def record_checksums():
    """Zenodo API: file key -> (md5, size, link). Empty dict on failure."""
    try:
        req = urllib.request.Request(RECORD_API, headers={"User-Agent": "MV2-Munim-research/1.0"})
        with urllib.request.urlopen(req, timeout=120) as r:
            rec = json.loads(r.read().decode("utf-8"))
    except Exception as ex:  # noqa: BLE001
        log("record API failed", repr(ex)[:120])
        return {}
    out = {}
    files = rec.get("files", [])
    if isinstance(files, dict):
        files = files.get("entries", [])
        files = list(files.values()) if isinstance(files, dict) else files
    for f in files:
        key = f.get("key") or f.get("filename") or ""
        ck = str(f.get("checksum", ""))
        out[key] = dict(md5=ck.split(":", 1)[1] if ":" in ck else ck,
                        size=f.get("size"),
                        link=(f.get("links") or {}).get("self", ""))
    return out


def extract(zpath, dest, bundle, files_rows, depth=0):
    """Safe extraction (one nested level). Returns list of (member, path)."""
    out = []
    with zipfile.ZipFile(zpath) as z:
        for info in z.infolist():
            if info.is_dir():
                continue
            rel = safe_member(info.filename)
            kind = kind_of(info.filename)
            files_rows.append(dict(bundle=bundle, member=info.filename, bytes=info.file_size,
                                   kind=kind, nested=depth))
            if rel is None or info.file_size > MEMBER_MAX_BYTES:
                continue
            target = os.path.join(dest, rel)
            os.makedirs(os.path.dirname(target) or dest, exist_ok=True)
            with z.open(info) as src, open(target, "wb") as dst:
                shutil.copyfileobj(src, dst, 1 << 20)
            if kind == "zip" and depth == 0:
                for m, p in extract(target, target + "_x", bundle, files_rows, depth + 1):
                    out.append((rel + "/" + m, p))
            else:
                out.append((rel, target))
    return out


def read_headers(bundle, members):
    """Header-only read of every non-text member. Returns trace rows."""
    from obspy import read
    rows = []
    for rel, path in sorted(members):
        if kind_of(rel) in ("text", "zip"):
            continue
        try:
            st = read(path, headonly=True)
        except Exception as ex:  # noqa: BLE001
            rows.append(dict(bundle=bundle, member=rel, status="unreadable: " + repr(ex)[:80]))
            continue
        for i, tr in enumerate(st):
            s = tr.stats
            sac = s.get("sac", {}) or {}
            b = sac.get("b")
            ref = (s.starttime - float(b)) if b is not None else None
            o = sac.get("o")
            row = dict(bundle=bundle, member=rel, index=i, status="ok",
                       format=s.get("_format", ""), net=s.network, sta=s.station,
                       loc=s.location, cha=s.channel, start=str(s.starttime),
                       end=str(s.endtime), sr=float(s.sampling_rate), npts=int(s.npts),
                       ref=str(ref) if ref is not None else "",
                       origin=origin_key(str(ref), o) if (ref is not None and o is not None)
                       else "")
            for k in SAC_KEYS:
                v = sac.get(k, "")
                row["sac_" + k] = v.strip() if isinstance(v, str) else v
            rows.append(row)
    return rows


def copy_texts(bundle, members, out_dir, budget):
    copied = []
    for rel, path in sorted(members):
        if kind_of(rel) != "text":
            continue
        size = os.path.getsize(path)
        if size > TEXT_MAX_BYTES or size > budget[0]:
            copied.append(dict(member=rel, bytes=size, copied=False))
            continue
        tgt = os.path.join(out_dir, "text", bundle.replace(".zip", ""), rel)
        os.makedirs(os.path.dirname(tgt), exist_ok=True)
        shutil.copyfile(path, tgt)
        budget[0] -= size
        copied.append(dict(member=rel, bytes=size, copied=True))
    return copied


def write_csv(path, rows):
    keys = []
    for r in rows:
        for k in r:
            if k not in keys:
                keys.append(k)
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        for r in rows:
            w.writerow(r)


def run(args):
    os.makedirs(args.out_dir, exist_ok=True)
    zdir = os.path.join(args.data_dir, "zw26")
    os.makedirs(zdir, exist_ok=True)
    S = dict(version=VERSION, record=RECORD_API, started=time.strftime("%Y-%m-%dT%H:%M:%SZ",
                                                                       time.gmtime()),
             bundles={}, errors=[])
    ck = record_checksums()
    S["record_files"] = ck
    files_rows, trace_rows, budget = [], [], [TEXT_TOTAL_MAX]
    for bundle in BUNDLES:
        B = dict()
        S["bundles"][bundle] = B
        try:
            zpath = os.path.join(zdir, bundle)
            url = (ck.get(bundle) or {}).get("link") or (RECORD_API + "/files/" + bundle + "/content")
            if not (os.path.exists(zpath) and os.path.getsize(zpath) == EXPECTED_BYTES[bundle]):
                log("downloading", bundle)
                fetch(url, zpath)
            B["bytes"] = os.path.getsize(zpath)
            B["md5"] = md5_file(zpath)
            want = (ck.get(bundle) or {}).get("md5", "")
            B["md5_expected"] = want
            B["md5_ok"] = bool(want) and B["md5"] == want
            B["size_ok"] = B["bytes"] == EXPECTED_BYTES[bundle]
            xdir = os.path.join(zdir, "x_" + bundle.replace(".zip", ""))
            if os.path.isdir(xdir):
                shutil.rmtree(xdir)
            members = extract(zpath, xdir, bundle, files_rows)
            B["n_members"] = len(members)
            rows = read_headers(bundle, members)
            trace_rows += rows
            B["n_traces"] = sum(1 for r in rows if r.get("status") == "ok")
            B["n_unreadable"] = sum(1 for r in rows if r.get("status") != "ok")
            B["texts"] = copy_texts(bundle, members, args.out_dir, budget)
            ok = [r for r in rows if r.get("status") == "ok"]
            B["stations"] = sorted({r["sta"] for r in ok if r.get("sta")})
            B["channels"] = sorted({r["cha"] for r in ok if r.get("cha")})
            B["rates"] = sorted({round(r["sr"], 4) for r in ok})
            B["origins"] = sorted({r["origin"] for r in ok if r.get("origin")})
            B["picks_set"] = {k: sum(1 for r in ok if r.get("sac_" + k) not in ("", None))
                              for k in PICK_KEYS}
            B["pick_names"] = sorted({"%s=%s" % (k, str(r.get("sac_k" + k, "")).strip())
                                      for r in ok for k in PICK_KEYS
                                      if r.get("sac_" + k) not in ("", None)})[:50]
            B["start_minus_origin_s"] = _quantiles([_rel(r) for r in ok])
            B["duration_s"] = _quantiles([(r["npts"] - 1) / r["sr"] for r in ok if r["sr"] > 0])
            log(bundle, "members", B["n_members"], "traces", B["n_traces"],
                "origins", len(B["origins"]), "stations", len(B["stations"]))
        except Exception as ex:  # noqa: BLE001
            B["error"] = repr(ex)[:300]
            S["errors"].append("%s: %s" % (bundle, repr(ex)[:200]))
            log("ERROR", bundle, repr(ex)[:200])
    write_csv(os.path.join(args.out_dir, "zw26inv_files.csv"), files_rows)
    write_csv(os.path.join(args.out_dir, "zw26inv_traces.csv"), trace_rows)
    write_csv(os.path.join(args.out_dir, "zw26inv_events.csv"), summarise_events(trace_rows))
    S["finished"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    with open(os.path.join(args.out_dir, "zw26inv_summary.json"), "w", encoding="utf-8") as f:
        json.dump(S, f, indent=1, default=str)
    log("done; errors:", len(S["errors"]))
    return 0 if not S["errors"] else 1


def _rel(r):
    """Trace start minus origin (s), or None."""
    if not r.get("origin") or not r.get("start"):
        return None
    import datetime as dt
    try:
        t0 = dt.datetime.strptime(r["origin"], "%Y-%m-%dT%H:%M:%S")
        s = r["start"].rstrip("Z")
        s = s.split(".")[0] + "." + ((s.split(".")[1] if "." in s else "0") + "000000")[:6]
        t1 = dt.datetime.strptime(s, "%Y-%m-%dT%H:%M:%S.%f")
        return (t1 - t0).total_seconds()
    except Exception:  # noqa: BLE001
        return None


def _quantiles(vals):
    v = sorted(x for x in vals if x is not None)
    if not v:
        return {}
    q = lambda p: v[min(len(v) - 1, int(round(p * (len(v) - 1))))]  # noqa: E731
    return dict(n=len(v), min=v[0], p10=q(0.1), median=q(0.5), p90=q(0.9), max=v[-1])


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--data-dir", default="/data")
    ap.add_argument("--out-dir", default="/out")
    return run(ap.parse_args(argv))


if __name__ == "__main__":
    sys.exit(main())

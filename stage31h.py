# stage31h.py 0.1.0 -- data route for pre-registration 31 of inner-core-repeaters (Borovoye waveforms).
#
# Pre-registration 31 (findings/2026-11-23-prereg31-brv-waveform-plan.md) measures PKP-PP and PPP-PP on the
# Borovoye (BRV) records of the French Tuamotu tests listed in Zhou (2026) Table S1. Its first run
# (brv-wave31-v2, 2026-11-24) came back untestable: every per-event file is linked from LDEO's index page as
# ftp://ftp.ldeo.columbia.edu/..., and that host refused every connection.
#
# The same archive is also served by LDEO's web server as two master tar files over HTTPS
# (https://www.ldeo.columbia.edu/res/pi/Monitoring/Arch/BRV_arch_exp.html):
#   FRENCH_1977-1988.tar and FRENCH_1989-1996.tar.
# This script only moves data. It downloads those two tars, finds each Table S1 test's own per-event file
# inside them, and writes it into --data-dir under the exact file name that the UNCHANGED measurement code
# (brvpkp31.py 0.1.1 / brvpkp31_core.py at tag v0.32.0) will ask for. The measurement code then finds the
# file in its cache and never touches the network for it. No waveform is read here, and nothing about the
# measurement (windows, gates, grades, pairing, inputs) changes. Rules, fixed before any record was seen:
# findings/2026-11-24-prereg31-data-route-https.md.
#
# Part of munimtechnologies/inner-core-repeaters. Licence: Apache-2.0.
# Made with AI agents working in MV2's name (MV2, Munim, Inc.); no human expert has checked it.
import argparse
import csv
import datetime as _dt
import gzip
import hashlib
import io
import json
import os
import re
import shutil
import sys
import tarfile
import time

VERSION = "0.1.0"
MASTER_URLS = [
    "https://www.ldeo.columbia.edu/res/pi/Monitoring/Arch/Brv_archive/FRENCH_1977-1988.tar",
    "https://www.ldeo.columbia.edu/res/pi/Monitoring/Arch/Brv_archive/FRENCH_1989-1996.tar",
]
NOTE_PATH = "findings/2026-11-24-prereg31-data-route-https.md"
REPO_RAW = "https://raw.githubusercontent.com/munimtechnologies/inner-core-repeaters/"
UA = {"User-Agent": "inner-core-repeaters/stage31h (research)"}
FRENCH_SITE_DIGIT = "4"   # LDEO file names: YYYYMMDD-<site><n>-<system>; site 4 = French Tuamotu
MAX_DAY_OFFSET = 1        # file-name date may be the local date (e.g. 1980-07-19 23:46 UTC -> 19800720)

EVENT_FILE_RE = re.compile(r"^(\d{8})-(\d)(\d)-([A-Za-z]+)\.tar(\.Z|\.gz)?$")
EVENT_DIR_RE = re.compile(r"^(\d{8})-(\d)(\d)-([A-Za-z]+)$")


def log(*a):
    print(time.strftime("%H:%M:%S"), "[stage31h]", *a, flush=True)


# ---------------------------------------------------------------- pure helpers (stdlib only; unit-tested)

def parse_member(name):
    """Tar member path -> dict(kind='file'|'dir', ymd, site, seq, sys, ext, key) or None.
    kind='file': the member itself is a per-event tar(.Z/.gz). kind='dir': the member lies inside a
    per-event directory named YYYYMMDD-NN-SYS (key = that directory path)."""
    parts = [p for p in str(name).replace("\\", "/").split("/") if p not in ("", ".")]
    if not parts:
        return None
    m = EVENT_FILE_RE.match(parts[-1])
    if m:
        return dict(kind="file", ymd=m.group(1), site=m.group(2), seq=int(m.group(3)),
                    sys=m.group(4).upper(), ext=m.group(5) or "", key="/".join(parts))
    for i, p in enumerate(parts[:-1]):
        m = EVENT_DIR_RE.match(p)
        if m:
            return dict(kind="dir", ymd=m.group(1), site=m.group(2), seq=int(m.group(3)),
                        sys=m.group(4).upper(), ext="", key="/".join(parts[:i + 1]))
    return None


def day_offset(ymd_a, ymd_b):
    a = _dt.date(int(ymd_a[:4]), int(ymd_a[4:6]), int(ymd_a[6:8]))
    b = _dt.date(int(ymd_b[:4]), int(ymd_b[4:6]), int(ymd_b[6:8]))
    return (b - a).days


def choose_source(candidates, ymd, sysname, site=FRENCH_SITE_DIGIT, max_off=MAX_DAY_OFFSET):
    """candidates: list of parse_member dicts (one per event unit). Pick the French unit of this system
    for this date: same date first, then +1 day, then -1 day; within a date the lowest sequence number;
    a per-event file before a per-event directory. -> (unit or None, day_offset or None, n_same_rank)."""
    best, best_rank, ties = None, None, 0
    for c in candidates:
        if c["site"] != site or c["sys"] != sysname.upper():
            continue
        try:
            off = day_offset(ymd, c["ymd"])
        except ValueError:
            continue
        if abs(off) > max_off:
            continue
        rank = (abs(off), 0 if off >= 0 else 1, c["seq"], 0 if c["kind"] == "file" else 1)
        if best_rank is None or rank < best_rank:
            best, best_rank, ties = c, rank, 1
        elif rank == best_rank:
            ties += 1
    return best, (None if best is None else day_offset(ymd, best["ymd"])), ties


def units_from_members(names):
    """Member names -> list of distinct event units (per-event files and per-event directories)."""
    seen, out = set(), []
    for n in names:
        u = parse_member(n)
        if u is None or u["key"] in seen:
            continue
        seen.add(u["key"])
        out.append(u)
    return out


def write_unit(tf, unit, member_names, dest):
    """Write one event unit from open tarfile tf to dest (a name ending .tar.Z, .tar.gz or .tar).
    A per-event .tar.Z/.tar.gz file is copied byte for byte. A plain .tar, or a per-event directory
    (packed into a tar), is gzip-compressed when dest ends in .Z/.gz: the measurement code decompresses
    with `gzip -dc`, which reads gzip streams whatever the suffix. -> bytes written."""
    tmp = dest + ".part"
    want_gz = dest.endswith(".Z") or dest.endswith(".gz")
    if unit["kind"] == "file":
        src = tf.extractfile(tf.getmember(_orig_name(member_names, unit["key"])))
        if unit["ext"] in (".Z", ".gz") or not want_gz:
            with open(tmp, "wb") as f:
                shutil.copyfileobj(src, f, 1 << 20)
        else:
            with gzip.open(tmp, "wb") as f:
                shutil.copyfileobj(src, f, 1 << 20)
    else:
        prefix = unit["key"] + "/"
        opener = (lambda p: gzip.open(p, "wb")) if want_gz else (lambda p: open(p, "wb"))
        with opener(tmp) as raw, tarfile.open(fileobj=raw, mode="w") as out:
            for n in member_names:
                norm = "/".join(p for p in n.replace("\\", "/").split("/") if p not in ("", "."))
                if not norm.startswith(prefix):
                    continue
                m = tf.getmember(n)
                if not m.isfile():
                    continue
                data = tf.extractfile(m).read()
                ti = tarfile.TarInfo(norm[len(unit["key"].rsplit("/", 1)[0]) + 1:] if "/" in unit["key"] else norm)
                ti.size = len(data)
                out.addfile(ti, io.BytesIO(data))
    os.replace(tmp, dest)
    return os.path.getsize(dest)


def _orig_name(member_names, key):
    for n in member_names:
        if "/".join(p for p in n.replace("\\", "/").split("/") if p not in ("", ".")) == key:
            return n
    raise KeyError(key)


# ---------------------------------------------------------------- network + pipeline

def http_get(url, timeout=120, verify=True, stream=False):
    import requests
    return requests.get(url, timeout=timeout, verify=verify, headers=UA, stream=stream)


def repo_text(ref, path):
    r = http_get(REPO_RAW + ref + "/" + path, timeout=60)
    if r.status_code == 404:
        return None
    r.raise_for_status()
    return r.text


def wait_note(wait_min, poll_s=600):
    t0 = time.time()
    while True:
        try:
            txt = repo_text("main", NOTE_PATH)
            if txt:
                return txt
        except Exception as ex:  # noqa: BLE001
            log("note fetch error", repr(ex)[:200])
        if time.time() - t0 > wait_min * 60:
            return None
        log("data-route note not on main yet; waiting")
        time.sleep(poll_s)


def download(url, dest, rec):
    if os.path.exists(dest) and os.path.getsize(dest) > 0:
        rec.update(used="cache", bytes=os.path.getsize(dest))
        return dest
    os.makedirs(os.path.dirname(dest), exist_ok=True)
    errs = []
    for verify in (True, False):
        tmp = dest + ".part"
        try:
            with http_get(url, timeout=300, verify=verify, stream=True) as r:
                rec["http_status"] = r.status_code
                if r.status_code != 200:
                    raise IOError("HTTP %s" % r.status_code)
                with open(tmp, "wb") as f:
                    for b in r.iter_content(1 << 20):
                        f.write(b)
            if os.path.getsize(tmp) == 0:
                raise IOError("empty")
            os.replace(tmp, dest)
            rec.update(used=url, verify=verify, bytes=os.path.getsize(dest))
            return dest
        except Exception as ex:  # noqa: BLE001
            errs.append("verify=%s %r" % (verify, ex))
            if os.path.exists(tmp):
                os.remove(tmp)
    rec["error"] = " | ".join(errs)[:600]
    return None


def sha256_file(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def write_csv(path, rows):
    keys = []
    for r in rows:
        for k in r:
            if k not in keys:
                keys.append(k)
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=keys or ["empty"])
        w.writeheader()
        for r in rows:
            w.writerow(r)


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", default="data")
    ap.add_argument("--out-dir", default="out")
    ap.add_argument("--wait-min", type=float, default=0.0)
    a = ap.parse_args(argv)
    os.makedirs(a.out_dir, exist_ok=True)
    os.makedirs(a.data_dir, exist_ok=True)
    t0 = time.time()
    S = dict(version=VERSION, started=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()), errors=[],
             master_urls=MASTER_URLS, note_path=NOTE_PATH)
    with open(__file__, "rb") as f:
        S["code_sha256"] = hashlib.sha256(f.read()).hexdigest()

    def save():
        S["runtime_min"] = round((time.time() - t0) / 60.0, 2)
        with open(os.path.join(a.out_dir, "stage31h_summary.json"), "w") as f:
            json.dump(S, f, indent=1)

    try:
        note = wait_note(a.wait_min)
        S["note_public_at_start"] = bool(note)
        if not note:
            S["status"] = "no-note"
            save()
            return 0
        S["note_sha256"] = hashlib.sha256(note.encode("utf-8")).hexdigest()
        save()

        import brvpkp31_core as core  # the unchanged measurement code (tag v0.32.0)
        S["core_tag_inputs"] = core.TAG
        events = core.parse_s1(repo_text(core.TAG, core.S1_PATH))
        S["n_events"] = len(events)
        links = []
        for u in core.INDEX_URLS:
            try:
                r = http_get(u, timeout=60)
                if r.status_code == 200 and r.text:
                    links = core.parse_index_links(r.text, r.url or u)
                    S["index"] = dict(url=u, n_links=len(links))
                    break
            except Exception as ex:  # noqa: BLE001
                S["errors"].append("index %s %r" % (u, ex))
        if not links:
            S["errors"].append("index not reachable: requested names fall back to FTP_PATTERN")

        masters, member_rows, units, names_by_tar = [], [], [], {}
        for url in MASTER_URLS:
            rec = dict(url=url)
            p = download(url, os.path.join(a.data_dir, "master", url.rsplit("/", 1)[-1]), rec)
            if p:
                rec["sha256"] = sha256_file(p)
                try:
                    with tarfile.open(p) as tf:
                        names = tf.getnames()
                        sizes = {m.name: m.size for m in tf.getmembers()}
                    names_by_tar[p] = names
                    rec["n_members"] = len(names)
                    for n in names:
                        member_rows.append(dict(tar=os.path.basename(p), member=n, size=sizes.get(n, "")))
                    for u in units_from_members(names):
                        u["tar"] = p
                        units.append(u)
                except Exception as ex:  # noqa: BLE001
                    rec["error"] = "tar: %r" % ex
            masters.append(rec)
            save()
        S["masters"] = masters
        S["n_units"] = len(units)
        S["n_french_units"] = sum(1 for u in units if u["site"] == FRENCH_SITE_DIGIT)
        S["unit_systems"] = sorted(set(u["sys"] for u in units))
        write_csv(os.path.join(a.out_dir, "stage31h_members.csv"), member_rows)

        staged = []
        open_tars = {}
        for e in events:
            for sysname in core.SYSTEMS:
                req_url = (core.links_for(links, e["ymd"], sysname)
                           or [core.FTP_PATTERN.format(ymd=e["ymd"], sys=sysname)])[0]
                req = core.basename_of(req_url)
                row = dict(date=e["date"], starred=e["starred"], system=sysname, requested=req)
                unit, off, ties = choose_source(units, e["ymd"], sysname)
                if unit is None:
                    row["status"] = "no French file within 1 day"
                    staged.append(row)
                    continue
                row.update(source_tar=os.path.basename(unit["tar"]), source=unit["key"], kind=unit["kind"],
                           day_offset=off, same_rank_candidates=ties)
                dest = os.path.join(a.data_dir, req)
                try:
                    if os.path.exists(dest) and os.path.getsize(dest) > 0:
                        row["status"] = "already staged"
                    else:
                        tf = open_tars.get(unit["tar"])
                        if tf is None:
                            tf = open_tars[unit["tar"]] = tarfile.open(unit["tar"])
                        write_unit(tf, unit, names_by_tar[unit["tar"]], dest)
                        row["status"] = "staged"
                    row["bytes"] = os.path.getsize(dest)
                    row["sha256"] = sha256_file(dest)
                except Exception as ex:  # noqa: BLE001
                    row["status"] = "error: %r" % ex
                staged.append(row)
        for tf in open_tars.values():
            tf.close()
        write_csv(os.path.join(a.out_dir, "stage31h_staged.csv"), staged)
        S["n_staged"] = {s: sum(1 for r in staged if r["system"] == s and r.get("status") in
                                ("staged", "already staged")) for s in core.SYSTEMS}
        S["status"] = "done"
    except Exception as ex:  # noqa: BLE001
        S["errors"].append("fatal: %r" % ex)
        S["status"] = "error"
    save()
    log("done", json.dumps({k: S.get(k) for k in ("status", "n_units", "n_french_units", "n_staged")}))
    return 0


if __name__ == "__main__":
    sys.exit(main())

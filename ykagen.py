#!/usr/bin/env python3
"""ykagen.py: the YKA hardware control (pre-registration 9).

Every out-of-sample pair at the Yellowknife array (YKA) joins a 1990-2008 event, recorded
by the legacy system (YKB*/YKR*, 20 sps, CNSN digitiser), with a 2025-2026 event recorded
by the current system (YKAB*/YKAR*, 40 sps, Guralp digitiser). Both use Geotech S-13
sensors on the same piers. From 2013-06-26 to 2014-04-02 the two systems ran side by side,
so every earthquake in that window was recorded by both. This script measures, on those
events, how different the SAME earthquake looks to the two systems after yka.py's exact
processing: the hardware-only floor of the coda dissimilarity D, the direct-pulse
correlation (the stage-B gate), polarity, and the timing offset between the systems.

Processing is yka.py v0.2.0 unchanged (process_event, align_member, yka_pair, pick_onset),
fetched from the public repo at a pinned commit.

By MV2, Munim, Inc. (made with AI agents). Licence: Apache-2.0.
"""
import argparse
import csv
import io
import json
import math
import os
import sys
import time

VERSION = "0.9.1"
COMMIT = "7067d4fc2b0a2fdb1378a1c406dd680c8d0b69b1"  # v0.8.0 of munimtechnologies/inner-core-repeaters
RAW = "https://raw.githubusercontent.com/munimtechnologies/inner-core-repeaters/%s/" % COMMIT
MODULES = ("model.py", "phases.py", "xcorr.py", "yka.py")
USGS = "https://earthquake.usgs.gov/fdsnws/event/1/query"
EARTHSCOPE = "https://service.earthscope.org/fdsnws/dataselect/1/query"
NRCAN = "https://earthquakescanada.nrcan.gc.ca/fdsnws/dataselect/1/query"
YKA_CENTER = (62.4931, -114.6061)
OVERLAP = ("2013-06-27T00:00:00", "2014-04-01T00:00:00")  # both systems recording
RADIUS = (128.5, 145.5)  # inside yka.py's ak135 grid (128-146 deg) with array aperture to spare
MAX_DEPTH = 120.0         # yka.py's depth grid stops at 100 km
SSI_BOX = dict(lat=(-61.5, -54.0), lon=(-31.0, -22.0))
# (generation, station selector, providers in order)
GENS = (("legacy", "YKB*,YKR*", (NRCAN, EARTHSCOPE)), ("new", "YKAB*,YKAR*", (EARTHSCOPE, NRCAN)))
MIN_SNR_DIRECT, MIN_SNR_CODA, PAIR_THR = 5.0, 2.0, 0.90  # stage-B gates (pre-registration 2)
T0 = time.time()


def log(*a):
    print("[%7.1f min]" % ((time.time() - T0) / 60.0), *a, flush=True)


def mods_dir_default(data_dir):
    return os.path.join(data_dir, "mods_" + COMMIT[:7])


def add_mod_path(data_dir="/data"):
    """Make yka.py and its helpers importable: repo root first, else the pinned download."""
    here = os.path.dirname(os.path.abspath(__file__))
    if os.path.exists(os.path.join(here, "yka.py")):
        return here
    d = mods_dir_default(data_dir)
    if os.path.isdir(d) and d not in sys.path:
        sys.path.insert(0, d)
    return d


# ----------------------------------------------------------------------------- downloads
def fetch_text(url, path, must_contain=None, tries=6):
    import requests

    for attempt in range(tries):
        try:
            r = requests.get(url, timeout=90)
            if r.status_code == 200 and (must_contain is None or must_contain in r.text):
                with open(path, "w", encoding="utf-8", newline="") as f:
                    f.write(r.text)
                return True
            log("fetch", url, "status", r.status_code)
        except Exception as ex:
            log("fetch error", url, repr(ex)[:200])
        time.sleep(10 * (attempt + 1))
    return False


def get_modules(mod_dir):
    os.makedirs(mod_dir, exist_ok=True)
    for m in MODULES:
        p = os.path.join(mod_dir, m)
        if not os.path.exists(p):
            need = 'VERSION = "0.2.0"' if m == "yka.py" else None
            if not fetch_text(RAW + m, p, must_contain=need):
                raise SystemExit("could not fetch %s" % m)


def in_ssi(lat, lon):
    return SSI_BOX["lat"][0] <= lat <= SSI_BOX["lat"][1] and SSI_BOX["lon"][0] <= lon <= SSI_BOX["lon"][1]


def fetch_catalog(data_dir, min_mag, decimal_year):
    import requests

    path = os.path.join(data_dir, "ykagen_catalog.csv")
    params = dict(format="csv", starttime=OVERLAP[0], endtime=OVERLAP[1], minmagnitude=min_mag,
                  latitude=YKA_CENTER[0], longitude=YKA_CENTER[1], minradius=RADIUS[0], maxradius=RADIUS[1],
                  maxdepth=MAX_DEPTH, orderby="time-asc")
    text = None
    for attempt in range(4):
        try:
            r = requests.get(USGS, params=params, timeout=180)
            if r.status_code == 200:
                text = r.text
                break
            log("USGS status", r.status_code)
        except Exception as ex:
            log("USGS error", repr(ex)[:200])
        time.sleep(10 * (attempt + 1))
    if text is None:
        if not os.path.exists(path):
            raise SystemExit("USGS catalogue unavailable")
        with open(path, encoding="utf-8") as f:
            text = f.read()
    else:
        with open(path, "w", encoding="utf-8", newline="") as f:
            f.write(text)
    evs = []
    for r in csv.DictReader(io.StringIO(text)):
        if r.get("type", "earthquake") != "earthquake":
            continue
        try:
            lat, lon = float(r["latitude"]), float(r["longitude"])
            evs.append(dict(id=r["id"], time=r["time"], lat=lat, lon=lon,
                            depth=float(r["depth"]) if r.get("depth") not in (None, "") else 10.0,
                            mag=float(r["mag"]) if r.get("mag") not in (None, "") else float("nan"),
                            t=decimal_year(r["time"]), ssi=in_ssi(lat, lon)))
        except Exception:
            continue
    evs.sort(key=lambda e: e["t"])
    return evs


def dataselect(session, url, sta, t1, t2):
    from obspy import read

    params = dict(net="CN", sta=sta, loc="*", cha="SHZ", start=t1.strftime("%Y-%m-%dT%H:%M:%S"),
                  end=t2.strftime("%Y-%m-%dT%H:%M:%S"))
    err = ""
    for attempt in range(3):
        try:
            r = session.get(url, params=params, timeout=120)
            if r.status_code in (204, 404):
                return None, "nodata"
            if r.status_code == 200 and r.content:
                st = read(io.BytesIO(r.content))
                return (st, "ok") if len(st) else (None, "nodata")
            err = "status %d" % r.status_code
        except Exception as ex:
            err = repr(ex)[:160]
        time.sleep(5 * (attempt + 1))
    return None, "error:" + err


def fetch_gen(session, y, ev, ph, data_dir, gen, sta, providers):
    """Records of one hardware generation only, cached per event and generation."""
    from obspy import UTCDateTime, read

    rawdir = os.path.join(data_dir, "ykagen_raw")
    os.makedirs(rawdir, exist_ok=True)
    p = os.path.join(rawdir, "%s_%s.mseed" % (ev["id"], gen))
    if os.path.exists(p):
        try:
            return read(p), "cached"
        except Exception:
            os.remove(p)
    dist = float(y.gc_deg(ev["lat"], ev["lon"], *YKA_CENTER))
    tdf = ph.t("df", dist, ev["depth"])
    if not math.isfinite(tdf):
        return None, "no-df-prediction"
    o = UTCDateTime(ev["time"])
    t1, t2 = o + tdf - y.LONG_PRE, o + tdf + y.LONG_POST
    best, notes = None, []
    for url in providers:
        st, s = dataselect(session, url, sta, t1, t2)
        notes.append(s)
        if st is None:
            continue
        st = st.select(channel="SHZ")
        st.traces = [tr for tr in st if y.generation(tr.stats.station) == gen and y.canonical(tr.stats.station)]
        n_el = len({y.canonical(tr.stats.station) for tr in st})
        if n_el and (best is None or n_el > best[1]):
            best = (st, n_el, url)
        if n_el >= 5:
            break
        time.sleep(0.3)
    if best is None:
        return None, ";".join(notes)
    try:
        best[0].write(p, format="MSEED")
    except Exception as ex:
        log("could not cache", ev["id"], gen, repr(ex)[:120])
    return best[0], "fetched:" + ("earthscope" if best[2] == EARTHSCOPE else "nrcan")


# ----------------------------------------------------------------------------- the comparison
def compare(y, dn, dl):
    """Same earthquake, two systems: current (dn, the reference) vs legacy (dl).

    Uses yka.py's stage-B steps unchanged: envelope + waveform alignment with the flip
    test (the generations differ), STA/LTA onset on the stack, then yka_pair (direct-pulse
    alignment, coda window CODA, plus the three sensitivity windows).
    """
    FS, PRE, N = y.FS, y.PRE, y.N
    bn, bl = dn["beam"], dl["beam"]
    a = y.align_member(bn, bl, FS, PRE, allow_flip=True)
    out = dict(status="ok", align_cc_env=a["cc_env"], align_cc_fine=a["cc_fine"], flip=bool(a["flip"]))
    if not (math.isfinite(a["off"]) and math.isfinite(a["cc_fine"]) and a["cc_fine"] >= 0.5):
        out["status"] = "alignment-failed"
        return out
    off = float(a["off"])
    sgn = -1.0 if a["flip"] else 1.0
    stack = 0.5 * (bn + y.fshift(bl, -off) * sgn)
    onset = y.pick_onset(stack, FS, center=int(PRE * FS), search_s=6.0, sta_s=0.5, lta_s=8.0, thr=6.0)
    if onset is None:
        out["status"] = "no-DF-pick"
        return out
    n0, n1 = int((PRE - 35) * FS), int((PRE - 5) * FS)

    def snr(b, o, w):
        s = int(round(o))
        i0, i1 = onset + int(round(w[0] * FS)) + s, onset + int(round(w[1] * FS)) + s
        if i0 < 0 or i1 > N:
            return float("nan")
        nb = math.sqrt(float((b[n0:n1] ** 2).mean()))
        return math.sqrt(float((b[i0:i1] ** 2).mean())) / nb if nb > 0 else float("nan")

    r = y.yka_pair(bn, bl * sgn, onset, 0.0, off, FS)
    if r is None:
        out["status"] = "window-out-of-range"
        return out
    out.update(off_s=off / FS, pick_minus_pred_s=(onset - PRE * FS) / FS, cc_ref=r["cc_ref"],
               cc_df0=r["cc_df0"], D=1.0 - r["cc_df0"], Dx=(1.0 - r["cc_df0"]) - (1.0 - r["cc_ref"]),
               dt_df_s=r["dt_df_s"], snr_dir_new=snr(bn, 0.0, y.DIRECT), snr_dir_leg=snr(bl, off, y.DIRECT),
               snr_coda_new=snr(bn, 0.0, y.CODA), snr_coda_leg=snr(bl, off, y.CODA))
    for key, w in y.CODA_SENS.items():
        rs = y.yka_pair(bn, bl * sgn, onset, 0.0, off, FS, window=w)
        out[key] = (1.0 - rs["cc_df0"]) if rs is not None else float("nan")
    # element by element on the shared piers (each generation's own statics already applied)
    st_D, st_ref = [], []
    if "X" in dn and "X" in dl:
        nn, nl = [str(s) for s in dn["names"]], [str(s) for s in dl["names"]]
        for s in sorted(set(nn) & set(nl)):
            xn = dn["X"][nn.index(s)].astype(float)
            xl = dl["X"][nl.index(s)].astype(float) * sgn
            rs = y.yka_pair(xn, xl, onset, 0.0, off, FS)
            if rs is not None:
                st_ref.append(rs["cc_ref"])
                if rs["cc_ref"] >= 0.6:
                    st_D.append(1.0 - rs["cc_df0"])
    out["n_el_shared"] = len(st_ref)
    out["Ds_median"] = median(st_D) if len(st_D) >= 5 else float("nan")
    out["n_el_D"] = len(st_D)
    out["ok"] = bool(all(math.isfinite(out[k]) for k in ("snr_dir_new", "snr_dir_leg", "snr_coda_new",
                                                          "snr_coda_leg", "D", "cc_ref"))
                     and min(out["snr_dir_new"], out["snr_dir_leg"]) >= MIN_SNR_DIRECT
                     and min(out["snr_coda_new"], out["snr_coda_leg"]) >= MIN_SNR_CODA)
    return out


# ----------------------------------------------------------------------------- grading
def median(v):
    v = sorted(x for x in v if x is not None and math.isfinite(x))
    if not v:
        return float("nan")
    n = len(v)
    return v[n // 2] if n % 2 else 0.5 * (v[n // 2 - 1] + v[n // 2])


def quantile(v, q):
    v = sorted(x for x in v if x is not None and math.isfinite(x))
    if not v:
        return float("nan")
    k = (len(v) - 1) * q
    lo, hi = int(math.floor(k)), int(math.ceil(k))
    return v[lo] + (v[hi] - v[lo]) * (k - lo)


def grade(rows, min_events=5):
    """Pre-registration 9: H1 (floor), H2 (gate survives), H3 (polarity)."""
    ok = [r for r in rows if r.get("status") == "ok" and r.get("ok")]
    res = dict(n_rows=len(rows), n_compared=sum(1 for r in rows if r.get("status") == "ok"), n_ok=len(ok),
               n_ok_ssi=sum(1 for r in ok if r.get("ssi")), testable=len(ok) >= min_events)
    Ds = [r["D"] for r in ok]
    res["D_median"], res["D_q25"], res["D_q75"], res["D_q90"] = (median(Ds), quantile(Ds, 0.25),
                                                                   quantile(Ds, 0.75), quantile(Ds, 0.9))
    res["frac_cc_ref_ge_thr"] = (sum(r["cc_ref"] >= PAIR_THR for r in ok) / len(ok)) if ok else float("nan")
    res["frac_flipped"] = (sum(bool(r["flip"]) for r in ok) / len(ok)) if ok else float("nan")
    res["off_s_median"] = median([r["off_s"] for r in ok])
    res["off_s_q10_q90"] = [quantile([r["off_s"] for r in ok], 0.1), quantile([r["off_s"] for r in ok], 0.9)]
    res["Ds_median_of_events"] = median([r.get("Ds_median", float("nan")) for r in ok])
    res["sens_D_median"] = {k: median([r.get(k, float("nan")) for r in ok]) for k in ("D_c13", "D_c310", "D_full")}
    res["ssi_D_median"] = median([r["D"] for r in ok if r.get("ssi")])
    if res["testable"]:
        res["H1_floor_small"] = "holds" if res["D_median"] <= 0.05 else "fails"
        res["H2_gate_survives"] = "holds" if res["frac_cc_ref_ge_thr"] >= 0.80 else "fails"
        res["H3_no_flip"] = "holds" if res["frac_flipped"] <= 0.10 else "fails"
    else:
        res["H1_floor_small"] = res["H2_gate_survives"] = res["H3_no_flip"] = "untestable"
    return res


def clean(o):
    if isinstance(o, dict):
        return {str(k): clean(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [clean(v) for v in o]
    if isinstance(o, float) and not math.isfinite(o):
        return None
    if hasattr(o, "item") and not isinstance(o, (str, bytes)):
        try:
            return clean(o.item())
        except Exception:
            return o
    return o


def plot_overlays(path, y, picks):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    if not picks:
        return
    fig, axs = plt.subplots(len(picks), 1, figsize=(10, 1.6 * len(picks) + 1), squeeze=False)
    for ax, (ev, dn, dl, r) in zip(axs[:, 0], picks):
        t = np.arange(y.N) / y.FS - y.PRE - r["pick_minus_pred_s"]
        sgn = -1.0 if r["flip"] else 1.0
        bl = y.fshift(dl["beam"], -r["off_s"] * y.FS) * sgn
        for b, lab, c in ((dn["beam"], "current (YKAB/YKAR)", "C0"), (bl, "legacy (YKB/YKR)", "C3")):
            seg = b[int((y.PRE - 3) * y.FS):int((y.PRE + 10) * y.FS)]
            sc = float(np.max(np.abs(seg))) or 1.0
            ax.plot(t, b / sc, lw=0.8, color=c, label=lab)
        ax.axvspan(y.DIRECT[0], y.DIRECT[1], color="grey", alpha=0.2)
        ax.axvspan(y.CODA[0], y.CODA[1], color="orange", alpha=0.2)
        ax.set_xlim(-3, 10)
        ax.set_yticks([])
        ax.text(10.1, 0, "%s M%.1f\nD=%.3f cc_ref=%.2f" % (ev["time"][:10], ev["mag"], r["D"], r["cc_ref"]),
                fontsize=7, va="center")
    axs[0, 0].legend(loc="upper left", fontsize=7)
    axs[0, 0].set_title("Same earthquake, two YKA systems (1-2 Hz beams; grey = direct PKIKP, orange = coda)")
    axs[-1, 0].set_xlabel("s relative to picked PKIKP onset")
    fig.tight_layout()
    fig.savefig(path, dpi=110)
    plt.close(fig)


# ----------------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", default="/data")
    ap.add_argument("--out-dir", default="/out")
    ap.add_argument("--min-mag", type=float, default=4.7)
    ap.add_argument("--max-events", type=int, default=400)
    ap.add_argument("--time-budget-min", type=float, default=100.0)
    args = ap.parse_args()
    for d in (args.data_dir, args.out_dir):
        os.makedirs(d, exist_ok=True)
    mod_dir = mods_dir_default(args.data_dir)
    get_modules(mod_dir)
    sys.path.insert(0, mod_dir)
    import requests

    import model
    import yka as y  # v0.2.0 at the pinned commit
    from phases import Phases

    ph = Phases(y.load_or_build(os.path.join(args.data_dir, "ak135_yka_grid_v1.npz")))
    evs = fetch_catalog(args.data_dir, args.min_mag, model.decimal_year)[:args.max_events]
    log("overlap-window events (M>=%.1f, %.1f-%.1f deg, depth<=%.0f km):" % (args.min_mag, RADIUS[0], RADIUS[1],
                                                                             MAX_DEPTH), len(evs),
        "of which SSI", sum(e["ssi"] for e in evs))
    session = requests.Session()
    session.headers["User-Agent"] = "inner-core-repeaters/%s (MV2, Munim, Inc.; mv2@munimtech.com)" % VERSION
    rows, keep = [], []
    for k, ev in enumerate(evs):
        row = dict(id=ev["id"], time=ev["time"], t=ev["t"], lat=ev["lat"], lon=ev["lon"], depth=ev["depth"],
                   mag=ev["mag"], ssi=ev["ssi"], dist=float(y.gc_deg(ev["lat"], ev["lon"], *YKA_CENTER)))
        if (time.time() - T0) / 60.0 > args.time_budget_min:
            row["status"] = "skipped-time-budget"
            rows.append(row)
            continue
        procs = {}
        try:
            for gen, sta, providers in GENS:
                st, s = fetch_gen(session, y, ev, ph, args.data_dir, gen, sta, providers)
                row["src_" + gen] = s
                if st is None:
                    break
                d, s2 = y.process_event(st, ev, ph)
                row["proc_" + gen] = s2
                if d is None:
                    break
                row["n_el_" + gen] = len(d["names"])
                row["snr_beam_" + gen] = d["snr"]
                procs[gen] = d
            if len(procs) < 2:
                row["status"] = "one-generation-only"
            else:
                row.update(compare(y, procs["new"], procs["legacy"]))
                if row["status"] == "ok":
                    keep.append((ev, procs["new"], procs["legacy"], row))
        except Exception as ex:
            row["status"] = "exception:" + repr(ex)[:150]
        rows.append(row)
        if (k + 1) % 20 == 0:
            log("events %d/%d, compared %d" % (k + 1, len(evs), len(keep)))
        time.sleep(0.2)

    keys = ["id", "time", "t", "lat", "lon", "depth", "mag", "ssi", "dist", "status", "src_legacy", "src_new",
            "proc_legacy", "proc_new", "n_el_legacy", "n_el_new", "snr_beam_legacy", "snr_beam_new", "ok",
            "flip", "align_cc_env", "align_cc_fine", "off_s", "pick_minus_pred_s", "cc_ref", "cc_df0", "D", "Dx",
            "dt_df_s", "D_c13", "D_c310", "D_full", "snr_dir_new", "snr_dir_leg", "snr_coda_new", "snr_coda_leg",
            "n_el_shared", "n_el_D", "Ds_median"]
    with open(os.path.join(args.out_dir, "ykagen_events.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(clean(r))
    res = grade(rows)
    res.update(version=VERSION, commit=COMMIT, n_catalogue=len(evs), overlap=OVERLAP, radius_deg=RADIUS,
               min_mag=args.min_mag, status_counts={})
    for r in rows:
        s = str(r.get("status", "?")).split(":")[0]
        res["status_counts"][s] = res["status_counts"].get(s, 0) + 1
    with open(os.path.join(args.out_dir, "ykagen_summary.json"), "w") as f:
        json.dump(clean(res), f, indent=1)
    try:
        os.makedirs(os.path.join(args.out_dir, "plots"), exist_ok=True)
        okp = sorted([p for p in keep if p[3].get("ok")], key=lambda p: -min(p[3]["snr_dir_new"], p[3]["snr_dir_leg"]))
        plot_overlays(os.path.join(args.out_dir, "plots", "ykagen_overlays.png"), y, okp[:8])
        low = sorted([p for p in keep if p[3].get("ok")], key=lambda p: min(p[3]["snr_dir_new"], p[3]["snr_dir_leg"]))
        plot_overlays(os.path.join(args.out_dir, "plots", "ykagen_overlays_low_snr.png"), y, low[:6])
    except Exception as ex:
        log("plot error", repr(ex)[:200])
    log("summary:", json.dumps(clean({k: res[k] for k in ("n_ok", "D_median", "frac_cc_ref_ge_thr", "frac_flipped",
                                                          "off_s_median", "H1_floor_small", "H2_gate_survives",
                                                          "H3_no_flip")})))
    log("done")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
# ilresp.py - What did ILAR's 2018-09-25 digitiser-chain change do to PKIKP waveforms?
#
# Part of github.com/munimtechnologies/inner-core-repeaters (Apache-2.0).
# By MV2 of Munim, Inc.; made with AI agents.
#
# EarthScope's StationXML for every short-period ILAR element (IM.IL01-IL19 SHZ) has two
# channel epochs split at 2018-09-25. The seismometer poles/zeros are the same in both, but
# the digital decimation chain is different (stages 3-10 run 5x faster after the change and
# the final FIR stages are replaced). The 2024 backtracking result uses ILAR records from
# 1991-2023, so pairs that straddle 2018-09-25 compare two different digital chains.
#
# This script measures, from the published responses only (no waveforms), the "hardware
# term" that change adds to waveform decorrelation D = 1 - cc in the band the study uses:
#   R(f) = H_after(f) / H_before(f) for each element; amplitude ripple, residual phase after
#   removing a best-fit delay, and D between a signal and the same signal passed through R.
# Plan, forecasts and decision rule: findings/2026-11-19-ilar-2018-digitiser.md

import argparse
import csv
import io
import json
import math
import os
import sys
import time
import urllib.request

import numpy as np

VERSION = "0.1.0"
REPO_RAW = ("https://raw.githubusercontent.com/munimtechnologies/inner-core-repeaters/main/"
            "findings/2026-11-19-ilar-2018-digitiser.md")
PLAN_MARKER = "PLAN-ILRESP-2018"
URL = ("https://service.earthscope.org/fdsnws/station/1/query?network=IM&station={sta}"
       "&channel=SHZ&level=response&format=xml")
STATIONS = ["IL%02d" % i for i in range(1, 20)]
FS = 20.0                 # data sample rate of every SHZ element, both epochs
NFFT = 8192               # 409.6 s at 20 Hz; signals sit in the middle, no wrap-around
BANDS = {"b1_2": (1.0, 2.0), "b05_2": (0.5, 2.0)}   # yka.py BAND; wider sensitivity band
MAIN_BAND = "b1_2"
FREF = 1.5                # Hz, amplitude normalisation inside the band
WIN_S = 20.0              # cc window length (s)
MAXLAG_S = 1.0            # lag search (s)
N_NOISE = 50              # noise realisations per element
SEED = 19
DATE_BEFORE = "2015-01-01T00:00:00"
DATE_AFTER = "2022-01-01T00:00:00"
# Decision rule (fixed in the plan): D values of repeater pairs in this study are ~0.07-0.2.
NEGLIGIBLE_P95 = 0.005
MATERIAL_MEDIAN = 0.02


# ----------------------------------------------------------------------------- responses

def expand_fir(coeffs, symmetry):
    """Full FIR taps from StationXML coefficients and symmetry (NONE / EVEN / ODD)."""
    c = [float(x) for x in coeffs]
    s = (symmetry or "NONE").upper()
    if s == "EVEN":
        return np.array(c + c[::-1], dtype=float)
    if s == "ODD":
        return np.array(c + c[:-1][::-1], dtype=float)
    return np.array(c, dtype=float)


def digital_response(b, a, fs, f):
    """B(z)/A(z) at z = exp(i 2 pi f / fs), with b, a in powers of z^-1."""
    f = np.asarray(f, dtype=float)
    zinv = np.exp(-2j * np.pi * f / fs)
    num = np.zeros_like(zinv)
    for k, bk in enumerate(b):
        num = num + bk * zinv ** k
    den = np.zeros_like(zinv)
    for k, ak in enumerate(a):
        den = den + ak * zinv ** k
    return num / den


def pz_response(zeros, poles, a0, kind, f):
    f = np.asarray(f, dtype=float)
    kind = (kind or "").upper()
    if "HERTZ" in kind:
        s = 1j * f
    else:
        s = 2j * np.pi * f
    h = np.full(f.shape, complex(a0), dtype=complex)
    for z in zeros:
        h = h * (s - complex(z))
    for p in poles:
        h = h / (s - complex(p))
    return h


def stage_response(stage, f):
    """Complex response of one ObsPy response stage at frequencies f (Hz).

    Digital stages include the time correction the digitiser applied (StationXML
    Decimation/Correction): effective phase = filter phase + 2 pi f * correction.
    """
    f = np.asarray(f, dtype=float)
    name = type(stage).__name__
    gain = float(stage.stage_gain) if getattr(stage, "stage_gain", None) else 1.0
    fs = getattr(stage, "decimation_input_sample_rate", None)
    fs = float(fs) if fs else None
    digital = False
    if name == "PolesZerosResponseStage":
        kind = getattr(stage, "pz_transfer_function_type", "") or ""
        a0 = float(getattr(stage, "normalization_factor", 1.0) or 1.0)
        if "DIGITAL" in kind.upper():
            if fs is None:
                raise ValueError("digital PZ stage without sample rate")
            z = np.exp(2j * np.pi * f / fs)
            h = np.full(f.shape, complex(a0), dtype=complex)
            for zz in stage.zeros:
                h = h * (z - complex(zz))
            for pp in stage.poles:
                h = h / (z - complex(pp))
            digital = True
        else:
            h = pz_response(stage.zeros, stage.poles, a0, kind, f)
        h = h * gain
    elif name == "CoefficientsTypeResponseStage":
        b = [float(x) for x in (stage.numerator or [])] or [1.0]
        a = [float(x) for x in (stage.denominator or [])] or [1.0]
        if len(b) == 1 and len(a) == 1:
            h = np.full(f.shape, gain * b[0] / a[0], dtype=complex)
        else:
            if fs is None:
                raise ValueError("coefficient stage without sample rate")
            h = gain * digital_response(b, a, fs, f)
            digital = True
    elif name == "FIRResponseStage":
        taps = expand_fir(stage.coefficients, getattr(stage, "symmetry", "NONE"))
        if fs is None:
            raise ValueError("FIR stage without sample rate")
        h = gain * digital_response(list(taps), [1.0], fs, f)
        digital = True
    else:
        h = np.full(f.shape, gain, dtype=complex)
    corr = getattr(stage, "decimation_correction", None)
    if digital and corr:
        h = h * np.exp(2j * np.pi * f * float(corr))
    return h


def chain_response(response, f):
    h = np.ones(np.asarray(f).shape, dtype=complex)
    for st in response.response_stages:
        h = h * stage_response(st, f)
    return h


def stage_table(response):
    rows = []
    for st in response.response_stages:
        ntaps = None
        if type(st).__name__ == "FIRResponseStage":
            ntaps = int(len(expand_fir(st.coefficients, getattr(st, "symmetry", "NONE"))))
        rows.append({
            "stage": st.stage_sequence_number, "type": type(st).__name__,
            "fs_in": getattr(st, "decimation_input_sample_rate", None),
            "factor": getattr(st, "decimation_factor", None),
            "delay": getattr(st, "decimation_delay", None),
            "correction": getattr(st, "decimation_correction", None),
            "gain": getattr(st, "stage_gain", None), "ntaps": ntaps})
    return rows


# ----------------------------------------------------------------------------- analysis

def ratio_metrics(f, r, band, fref=FREF):
    """Amplitude ripple, best-fit delay and residual phase of a transfer ratio r(f)."""
    f = np.asarray(f, dtype=float)
    r = np.asarray(r, dtype=complex)
    m = (f >= band[0] - 1e-12) & (f <= band[1] + 1e-12)
    iref = int(np.argmin(np.abs(f - fref)))
    amp = np.abs(r) / np.abs(r[iref])
    phi = np.unwrap(np.angle(r[m]))
    A = np.vstack([np.ones(m.sum()), f[m]]).T
    coef, _, _, _ = np.linalg.lstsq(A, phi, rcond=None)
    resid = phi - A @ coef
    intercept = (coef[0] + np.pi) % (2 * np.pi) - np.pi
    return {
        "amp_dev_max_pct": float(100.0 * np.max(np.abs(amp[m] - 1.0))),
        "delay_s": float(-coef[1] / (2 * np.pi)),
        "phase_intercept_deg": float(np.degrees(intercept)),
        "resid_phase_rms_deg": float(np.degrees(np.sqrt(np.mean(resid ** 2)))),
        "resid_phase_max_deg": float(np.degrees(np.max(np.abs(resid)))),
    }


def apply_ratio(x, r):
    """Filter real signal x (length NFFT) by transfer ratio r given on rfftfreq(NFFT)."""
    n = len(x)
    return np.fft.irfft(np.fft.rfft(x, n=n) * r, n=n)


def remove_delay(f, r, tau):
    """Undo a pure delay tau (s) in a transfer ratio r(f): r * exp(+i 2 pi f tau)."""
    return np.asarray(r, dtype=complex) * np.exp(2j * np.pi * np.asarray(f, dtype=float) * tau)


def bandpass(x, band, fs=FS, corners=4):
    """Causal Butterworth band-pass, the same design as ObsPy's filter('bandpass')."""
    from scipy.signal import butter, sosfilt
    nyq = 0.5 * fs
    sos = butter(corners, [band[0] / nyq, band[1] / nyq], btype="band", output="sos")
    return sosfilt(sos, x)


def max_cc(a, b, maxlag):
    """Max normalised cc of b against a over integer lags +-maxlag (b shifted), with a
    parabolic refinement of the peak. a and b have equal length; the window is
    a[maxlag:-maxlag]. Returns (cc, lag_samples)."""
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    n = len(a)
    wa = a[maxlag:n - maxlag]
    ea = math.sqrt(float(np.dot(wa, wa)))
    ccs = []
    for lag in range(-maxlag, maxlag + 1):
        wb = b[maxlag + lag:n - maxlag + lag]
        eb = math.sqrt(float(np.dot(wb, wb)))
        ccs.append(float(np.dot(wa, wb)) / (ea * eb) if ea > 0 and eb > 0 else 0.0)
    ccs = np.array(ccs)
    k = int(np.argmax(ccs))
    peak, shift = ccs[k], 0.0
    if 0 < k < len(ccs) - 1:
        y0, y1, y2 = ccs[k - 1], ccs[k], ccs[k + 1]
        den = y0 - 2 * y1 + y2
        if den < 0:
            shift = 0.5 * (y0 - y2) / den
            peak = y1 - 0.25 * (y0 - y2) * shift
    return float(min(peak, 1.0)), float(k - maxlag + shift)


def ricker(n, fs, fc, t0):
    t = np.arange(n) / fs - t0
    a = (np.pi * fc * t) ** 2
    return (1 - 2 * a) * np.exp(-a)


def decorrelation(r, band, rng, n_noise=N_NOISE):
    """D = 1 - cc between bandpassed x and bandpassed (x through r)."""
    c = NFFT // 2
    half = int(round(WIN_S * FS / 2))
    ml = int(round(MAXLAG_S * FS))
    sl = slice(c - half - ml, c + half + ml)
    out = {"noise": [], "noise_lag_s": []}
    taper = np.zeros(NFFT)
    seg = 2 * (half + ml) + int(20 * FS)
    taper[c - seg // 2:c - seg // 2 + seg] = np.hanning(seg)
    for _ in range(n_noise):
        x = rng.standard_normal(NFFT) * taper
        y = apply_ratio(x, r)
        cc, lag = max_cc(bandpass(x, band)[sl], bandpass(y, band)[sl], ml)
        out["noise"].append(1.0 - cc)
        out["noise_lag_s"].append(lag / FS)
    x = ricker(NFFT, FS, 1.5, c / FS)
    y = apply_ratio(x, r)
    cc, lag = max_cc(bandpass(x, band)[sl], bandpass(y, band)[sl], ml)
    out["ricker"] = 1.0 - cc
    out["ricker_lag_s"] = lag / FS
    return out


# ----------------------------------------------------------------------------- I/O

def fetch(url, tries=4):
    last = None
    for i in range(tries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "inner-core-repeaters/ilresp"})
            with urllib.request.urlopen(req, timeout=60) as fh:
                return fh.read()
        except Exception as exc:  # network errors are recorded, then retried
            last = exc
            time.sleep(5 * (i + 1))
    raise RuntimeError("fetch failed: %s (%s)" % (url, last))


def wait_for_plan(wait_min):
    t_end = time.time() + 60 * wait_min
    while True:
        try:
            txt = fetch(REPO_RAW + "?t=%d" % int(time.time()), tries=1).decode("utf-8", "replace")
            if PLAN_MARKER in txt:
                return True
        except Exception:
            pass
        if time.time() > t_end:
            return False
        time.sleep(60)


def pick_epoch(channels, when):
    from obspy import UTCDateTime
    t = UTCDateTime(when)
    for ch in channels:
        if ch.start_date <= t and (ch.end_date is None or t < ch.end_date):
            return ch
    return None


def pct(v, q):
    return float(np.percentile(np.asarray(v, dtype=float), q))


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-dir", default="out")
    ap.add_argument("--wait-min", type=float, default=120)
    ap.add_argument("--no-wait", action="store_true")
    args = ap.parse_args(argv)
    os.makedirs(os.path.join(args.out_dir, "plots"), exist_ok=True)
    summary = {"version": VERSION, "started_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}
    if not args.no_wait:
        summary["plan_public_at_start"] = wait_for_plan(args.wait_min)
        if not summary["plan_public_at_start"]:
            json.dump(summary, open(os.path.join(args.out_dir, "ilresp_summary.json"), "w"), indent=1)
            print("plan not public; stopping", file=sys.stderr)
            return 3
    from obspy import read_inventory

    freqs = np.fft.rfftfreq(NFFT, 1.0 / FS)
    rows, stages, curves, errors = [], {}, {}, {}
    for sta in STATIONS:
        try:
            inv = read_inventory(io.BytesIO(fetch(URL.format(sta=sta))), format="STATIONXML")
            chans = [c for net in inv.networks for s in net.stations for c in s.channels
                     if c.code == "SHZ"]
            ea, eb = pick_epoch(chans, DATE_BEFORE), pick_epoch(chans, DATE_AFTER)
            if ea is None or eb is None:
                raise ValueError("missing epoch (%d SHZ epochs)" % len(chans))
            ha = chain_response(ea.response, freqs[1:])
            hb = chain_response(eb.response, freqs[1:])
            r = np.ones(len(freqs), dtype=complex)
            r[1:] = hb / ha
            iref = int(np.argmin(np.abs(freqs - FREF)))
            r = r / np.abs(r[iref])
            stages[sta] = {"before": stage_table(ea.response), "after": stage_table(eb.response)}
            curves[sta] = r
            row = {"element": sta, "n_epochs": len(chans),
                   "before_start": str(ea.start_date), "before_end": str(ea.end_date),
                   "after_start": str(eb.start_date),
                   "sens_before": ea.response.instrument_sensitivity.value,
                   "sens_after": eb.response.instrument_sensitivity.value}
            # D measures shape change only: the best-fit delay in the main band (a pure time
            # shift, which alignment removes in the real analysis) is taken out first, so
            # sub-sample interpolation error is not counted as hardware. D_raw keeps it.
            tau = ratio_metrics(freqs[1:], r[1:], BANDS[MAIN_BAND])["delay_s"]
            r_shape = remove_delay(freqs, r, tau)
            rng = np.random.default_rng(SEED)
            draw = decorrelation(r, BANDS[MAIN_BAND], rng)
            row["%s_D_raw_noise_median" % MAIN_BAND] = pct(draw["noise"], 50)
            row["%s_D_raw_ricker" % MAIN_BAND] = draw["ricker"]
            for bk, band in BANDS.items():
                m = ratio_metrics(freqs[1:], r[1:], band)
                for k, v in m.items():
                    row["%s_%s" % (bk, k)] = v
                rng = np.random.default_rng(SEED)
                d = decorrelation(r_shape, band, rng)
                row["%s_D_noise_median" % bk] = pct(d["noise"], 50)
                row["%s_D_noise_p95" % bk] = pct(d["noise"], 95)
                row["%s_lag_noise_median_s" % bk] = pct(d["noise_lag_s"], 50)
                row["%s_D_ricker" % bk] = d["ricker"]
                row["%s_lag_ricker_s" % bk] = d["ricker_lag_s"]
            # cross-check the hand-built chain against ObsPy/evalresp (shape only)
            try:
                fa = freqs[1:]
                ra = (eb.response.get_evalresp_response_for_frequencies(fa, output="VEL") /
                      ea.response.get_evalresp_response_for_frequencies(fa, output="VEL"))
                me = ratio_metrics(fa, ra, BANDS[MAIN_BAND])
                row["evalresp_amp_dev_max_pct"] = me["amp_dev_max_pct"]
                row["evalresp_resid_phase_rms_deg"] = me["resid_phase_rms_deg"]
                row["evalresp_delay_s"] = me["delay_s"]
            except Exception as exc:
                row["evalresp_error"] = repr(exc)[:200]
            rows.append(row)
        except Exception as exc:
            errors[sta] = repr(exc)[:300]

    keys = []
    for row in rows:
        for k in row:
            if k not in keys:
                keys.append(k)
    with open(os.path.join(args.out_dir, "ilresp_elements.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader()
        for row in rows:
            w.writerow(row)
    json.dump(stages, open(os.path.join(args.out_dir, "ilresp_stages.json"), "w"), indent=1, default=str)

    summary.update({"n_elements": len(rows), "errors": errors,
                    "decision_rule": {"negligible_if_p95_below": NEGLIGIBLE_P95,
                                      "material_if_median_above": MATERIAL_MEDIAN,
                                      "band": BANDS[MAIN_BAND]}})
    if rows:
        b = MAIN_BAND
        med = [r["%s_D_noise_median" % b] for r in rows]
        p95 = [r["%s_D_noise_p95" % b] for r in rows]
        summary["main"] = {
            "D_noise_median_over_elements": pct(med, 50),
            "D_noise_p95_max_over_elements": float(max(p95)),
            "D_ricker_max": float(max(r["%s_D_ricker" % b] for r in rows)),
            "D_raw_noise_median_over_elements": pct([r["%s_D_raw_noise_median" % b] for r in rows], 50),
            "amp_dev_max_pct": float(max(r["%s_amp_dev_max_pct" % b] for r in rows)),
            "resid_phase_rms_deg_max": float(max(r["%s_resid_phase_rms_deg" % b] for r in rows)),
            "delay_s_min": float(min(r["%s_delay_s" % b] for r in rows)),
            "delay_s_max": float(max(r["%s_delay_s" % b] for r in rows)),
        }
        mx = summary["main"]["D_noise_p95_max_over_elements"]
        md = summary["main"]["D_noise_median_over_elements"]
        if len(rows) < 10:
            verdict = "untestable"
        elif mx < NEGLIGIBLE_P95:
            verdict = "negligible"
        elif md > MATERIAL_MEDIAN:
            verdict = "material"
        else:
            verdict = "small"
        summary["verdict"] = verdict
    else:
        summary["verdict"] = "untestable"
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(2, 1, figsize=(8, 7), sharex=True)
        m = (freqs > 0.1) & (freqs < 9.5)
        for sta, r in curves.items():
            ax[0].plot(freqs[m], np.abs(r[m]), lw=0.8)
            ph = np.degrees(np.unwrap(np.angle(r[m])))
            ax[1].plot(freqs[m], ph, lw=0.8)
        for a in ax:
            a.axvspan(1.0, 2.0, color="0.9", zorder=0)
        ax[0].set_ylabel("|H_after/H_before| (norm. 1.5 Hz)")
        ax[1].set_ylabel("phase of ratio (deg)")
        ax[1].set_xlabel("frequency (Hz)")
        ax[0].set_title("ILAR SHZ, 2018-09-25 digitiser chain change (shaded: 1-2 Hz)")
        fig.tight_layout()
        fig.savefig(os.path.join(args.out_dir, "plots", "ilresp_ratio.png"), dpi=110)
    except Exception as exc:
        summary["plot_error"] = repr(exc)[:200]
    summary["finished_utc"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    json.dump(summary, open(os.path.join(args.out_dir, "ilresp_summary.json"), "w"), indent=1)
    print(json.dumps(summary, indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())

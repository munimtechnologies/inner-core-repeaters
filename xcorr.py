"""Cross-correlation, sub-sample shifts and onset picking for doublet analysis.

Part of inner-core-repeaters by MV2, Munim, Inc. (made with AI agents).
Licence: Apache-2.0.

Conventions
-----------
* A lag ``k`` (in samples) found by :func:`lag_scan` means that ``y`` is
  delayed by ``k`` samples relative to ``x``: ``y[n + k] ~ x[n]``.
* :func:`fshift` delays a trace by ``s`` samples (positive = later), so
  ``fshift(y, -k)`` aligns ``y`` back onto ``x``.
"""
import numpy as np


def _z(x):
    x = np.asarray(x, dtype=float)
    x = x - x.mean()
    n = np.linalg.norm(x)
    return x / n if n > 0 else x * 0.0


def cc(a, b):
    """Zero-lag normalised correlation coefficient (Pearson) of two equal-length windows."""
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    if len(a) != len(b) or len(a) < 2:
        return float("nan")
    return float(np.dot(_z(a), _z(b)))


def lag_scan(x, y, i0, i1, max_lag):
    """Correlate template x[i0:i1] with y[i0+k:i1+k] for k in [-max_lag, max_lag].

    Returns (lags, ccs). Lags that would run off either end of y give NaN.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    i0 = int(i0)
    i1 = int(i1)
    max_lag = int(max_lag)
    lags = np.arange(-max_lag, max_lag + 1)
    ccs = np.full(lags.shape, np.nan)
    L = i1 - i0
    if L < 3 or i0 < 0 or i1 > len(x):
        return lags, ccs
    t = _z(x[i0:i1])
    lo = max(-max_lag, -i0)
    hi = min(max_lag, len(y) - i1)
    if hi < lo:
        return lags, ccs
    seg = y[i0 + lo:i1 + hi]
    win = np.lib.stride_tricks.sliding_window_view(seg, L)
    num = win @ t  # t has zero mean, so the window mean drops out of the numerator
    s1 = win.sum(axis=1)
    s2 = np.einsum("ij,ij->i", win, win)
    var = s2 - s1 * s1 / L
    with np.errstate(invalid="ignore", divide="ignore"):
        c = num / np.sqrt(np.where(var > 0, var, np.nan))
    ccs[lo + max_lag:hi + max_lag + 1] = c
    return lags, ccs


def best_lag(lags, ccs):
    """Peak of a correlation function with parabolic sub-sample refinement.

    Returns (lag, cc_peak); (nan, nan) if everything is NaN.
    """
    ccs = np.asarray(ccs, dtype=float)
    if not np.any(np.isfinite(ccs)):
        return float("nan"), float("nan")
    k = int(np.nanargmax(ccs))
    c0 = float(ccs[k])
    frac = 0.0
    cpk = c0
    if 0 < k < len(ccs) - 1 and np.isfinite(ccs[k - 1]) and np.isfinite(ccs[k + 1]):
        ym, yp = float(ccs[k - 1]), float(ccs[k + 1])
        denom = ym - 2.0 * c0 + yp
        if denom < 0:
            frac = 0.5 * (ym - yp) / denom
            frac = max(-0.5, min(0.5, frac))
            cpk = c0 - 0.25 * (ym - yp) * frac
    return float(lags[k]) + frac, float(min(cpk, 1.0))


def fshift(y, s):
    """Delay y by s samples (s may be fractional; positive = later), via FFT with zero padding."""
    y = np.asarray(y, dtype=float)
    n = len(y)
    if s == 0 or n == 0:
        return y.copy()
    nfft = 2 * n
    Y = np.fft.rfft(y, nfft)
    f = np.fft.rfftfreq(nfft)
    return np.fft.irfft(Y * np.exp(-2j * np.pi * f * s), nfft)[:n]


def measure_pair(x, y, ref, dfw, fs, max_lag_ref_s=0.5, max_lag_df_s=0.25, lag0=0.0):
    """Doublet measurement.

    x, y   : traces (y = later event), same sampling.
    ref    : (i0, i1) reference window in x's frame (phases that do NOT sample the inner core).
    dfw    : (i0, i1) PKIKP (DF) window in x's frame.
    lag0   : prior estimate of y's delay relative to x, in samples.

    Steps: align y on x using the reference window only; then measure, in the DF window,
    (a) the zero-lag correlation cc_df0 (waveform similarity after reference alignment) and
    (b) the residual DF delay dt_df (seconds; positive = DF later in y relative to its
        reference phases), with its peak correlation cc_df_best.
    """
    y0 = fshift(y, -lag0) if lag0 else np.asarray(y, dtype=float)
    lags, c = lag_scan(x, y0, ref[0], ref[1], int(round(max_lag_ref_s * fs)))
    l_ref, cc_ref = best_lag(lags, c)
    if not np.isfinite(l_ref):
        return None
    y1 = fshift(y0, -l_ref)
    a, b = int(dfw[0]), int(dfw[1])
    if a < 0 or b > len(x) or b - a < 3:
        return None
    cc_df0 = cc(x[a:b], y1[a:b])
    lags2, c2 = lag_scan(x, y1, a, b, int(round(max_lag_df_s * fs)))
    l_df, cc_dfb = best_lag(lags2, c2)
    return dict(
        lag_total=float(lag0 + l_ref),
        lag_ref_s=float((lag0 + l_ref) / fs),
        cc_ref=float(cc_ref),
        cc_df0=float(cc_df0),
        dt_df_s=float(l_df / fs) if np.isfinite(l_df) else float("nan"),
        cc_df_best=float(cc_dfb),
    )


def sta_lta(x, nsta, nlta):
    """STA/LTA of x**2; STA window ends at sample n, LTA window precedes the STA window.

    Returns an array the length of x (zeros where undefined).
    """
    e = np.asarray(x, dtype=float) ** 2
    n = len(e)
    out = np.zeros(n)
    nsta = int(nsta)
    nlta = int(nlta)
    if n < nsta + nlta + 1 or nsta < 1 or nlta < 1:
        return out
    c = np.concatenate([[0.0], np.cumsum(e)])
    idx = np.arange(nsta + nlta - 1, n)
    s = (c[idx + 1] - c[idx + 1 - nsta]) / nsta
    l = (c[idx + 1 - nsta] - c[idx + 1 - nsta - nlta]) / nlta
    with np.errstate(invalid="ignore", divide="ignore"):
        r = np.where(l > 0, s / l, 0.0)
    out[idx] = r
    return out


def pick_onset(x, fs, center, search_s=6.0, sta_s=0.5, lta_s=8.0, thr=6.0, refine=True):
    """First arrival within center +/- search_s seconds.

    1. Trigger: first sample where STA/LTA (power) >= thr. Band-limited (1-2 Hz) noise has
       few degrees of freedom per STA window, so thr is set high (6) to avoid noise triggers.
    2. Refine: walk back from the trigger (at most sta_s + 0.5 s) while a 0.1-s trailing mean
       of x**2 stays above three times the pre-trigger noise power.

    Returns the onset sample index, or None if nothing triggers.
    """
    x = np.asarray(x, dtype=float)
    nsta, nlta = int(round(sta_s * fs)), int(round(lta_s * fs))
    r = sta_lta(x, nsta, nlta)
    a = max(0, int(center - search_s * fs))
    b = min(len(x), int(center + search_s * fs))
    hit = np.nonzero(r[a:b] >= thr)[0]
    if len(hit) == 0:
        return None
    trig = int(a + hit[0])
    if not refine:
        return trig
    e = x ** 2
    n0 = max(0, trig - nsta - nlta)
    n1 = max(n0 + 1, trig - nsta)
    noise = float(np.mean(e[n0:n1]))
    w = max(1, int(round(0.1 * fs)))
    c = np.concatenate([[0.0], np.cumsum(e)])

    def env(i):
        lo = max(0, i - w + 1)
        return (c[i + 1] - c[lo]) / (i + 1 - lo)

    i = trig
    lim = max(0, trig - nsta - int(round(0.5 * fs)))
    while i - 1 >= lim and env(i - 1) > 3.0 * noise:
        i -= 1
    return int(i)

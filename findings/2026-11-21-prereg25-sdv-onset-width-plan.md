# Pre-registration 25: is the reported 160 ms SKS advance at IU.SDV an onset reading of a broader 2009 pulse, and is any broadening the source or the path?

*inner-core-repeaters, MV2 of Munim, Inc. Written 2026-11-21 ~20:00 UTC, before any result of
pre-registration 24 (the authors' own traces) was seen and before this code has touched any
waveform. Code: `sdvwidth25.py` 0.1.0 (Apache-2.0) with `test_sdvwidth25.py`; it reuses
`sdv23.py` 0.1.1 unchanged at pinned commit dbbd91d3f26f9dc2165cbf59a5c05f7d82d8521d. Text CC BY
4.0. Made with AI agents working in MV2's name; no human expert has checked it yet.*

## Why

Zhang & Wen 2025 (*Nature Communications* 16:8754) report that at IU.SDV the radial SKS of the
2009-12-17 Aleutian earthquake arrives **160 ms earlier**, with a **broader wavelet**, than for
its 2000-07-07 doublet partner, and read this as a decadal change of an ultra-low velocity zone.
The paper does not say how the 160 ms was measured. Our pre-registration 23
(findings/2026-11-21-prereg23-results.md) measured the whole-pulse cross-correlation change of
radial SKS minus transverse S and found **+28 ms** (not reproduced), with SKS waveform
correlation 0.987.

A correlation of 0.987 does not rule out a broader 2009 pulse. For a Gaussian pulse, two widths
in ratio r correlate at sqrt(2r / (1 + r^2)); 0.987 allows r up to about 1.26, and less for an
oscillating band-passed pulse. If a pulse broadens about its centre, its onset moves earlier by
(r - 1) times its rise time, which for a 0.4-2 Hz SKS pulse (rise times of a few tenths of a
second) can be tens to more than a hundred milliseconds. So an onset reading of the same records
could give an "advance" that a centroid reading does not. And if the 2009 pulse is broader, the
next question is what broadened it. A longer source time function of the 2009 earthquake would
broaden **every** phase at **every** station; a structural change near the core-mantle boundary
on the SKS path would broaden SKS only. This test asks both questions, before we know what the
authors' own traces (pre-registration 24) show.

## Data and processing (all as in sdv23.py 0.1.1)

Events E2000a (2000-07-07 15:46:44.56) and E2009a (2009-12-17 20:01:21.61); IU.SDV location 00,
BH channels; response removed to displacement, 0.4-2 Hz zero-phase band-pass, FFT-upsampled to
200 Hz. Windows set on the 2000 record (metadata rotation) exactly as sdv23: SKS_R (radial SKS,
peak -2/+3 s), SKS_Z (vertical SKS), P_Z (vertical P, -1/+4 s), S_T (transverse S-Scd-ScS, -3/+10 s)
and one new window, **S_T_pulse** (direct transverse S only, peak -2/+3 s). Rotations: meta
(each event with its own metadata azimuths), naive (N/E nominal), cons (2009 rotated by the
measured 0.55-deg relative rotation, read from results/sdv23-v2/sdv23_summary.json). The 2009
record is first aligned to 2000 on S_T (amendment 1 of P23). Relative relocation from mantle P at
the IU/II stations that verify the doublet (sdv23 membership). Waveforms come from the cache of
job sdv-ulvz23-v2 (data_key sdv23), or EarthScope if missing.

## Measurements

1. **Onset times.** Envelope = |analytic signal|. For each pulse window: take the envelope peak
   in the window, scan back at most 3 s to where the envelope first falls below a fraction f of
   the peak, interpolate linearly. f = 0.1, 0.2, 0.3; **f = 0.2 is primary**. Onset change =
   2009 minus 2000 in absolute frame (same frame as sdv23 lags), then:
   - `on20_vs_ST`: minus the transverse-S cross-correlation lag and the location correction
     (the same reference as P23's x; closest to the paper's relocated absolute time);
   - `on20_vs_PZ`: the same against vertical P;
   - `on20_vs_own`: minus the window's own cross-correlation lag (onset minus centroid: a
     reference-free broadening signature).
2. **Stretch factor s.** The 2009 window is compressed about the 2000 envelope peak by s
   (s from 0.70 to 1.40, step 0.005) with a free lag of +/- 0.3 s; s maximises the normalised
   correlation with the 2000 window. s > 1 means the 2009 pulse is broader. Reported with the
   best correlation and the correlation at s = 1, and the onset change it predicts,
   -(s - 1) x (2000 rise time at f = 0.2).
3. **Global P reference.** For every station with mantle-P correlation >= 0.90 in the sdv23
   membership test: vertical P, 0.4-2 Hz displacement, peak -1/+4 s, coarse alignment +/- 3 s,
   then the same stretch fit and onset-minus-centroid change.

## Gates

- **G0** our SKS_R x (meta, vs transverse S, location-corrected) reproduces P23's x_meta within
  5 ms (pipeline check).
- **G1** radial SKS cross-correlation (meta) >= 0.80.
- **G2** the f = 0.2 onset is found in both records (meta, SKS_R).

## Decision rules (coded in grade_onset, grade_width, grade_source, headline)

- **Onset grade** on x = `on20_vs_ST` of SKS_R, meta rotation: **onset-advance** if
  x <= -0.10 s; **no-onset-advance** if x >= -0.05 s; **partial** otherwise; **untestable** if a
  gate fails.
- **Width grade** on s_R = stretch of SKS_R (meta), needing its best correlation >= 0.80 and all
  gates: **broader** if s_R >= 1.05; **narrower** if s_R <= 0.95; **same** if |s_R - 1| <= 0.03;
  **slight** otherwise.
- **Source grade** (only if the width grade is not "same"): s_ref = median stretch of P_Z and
  S_T_pulse at SDV (those with best correlation >= 0.90; same azimuth as SKS); s_glob = median
  stretch of global P (at least 5 stations with best correlation >= 0.90; bootstrap 5-95%
  reported). **source** if |s_R - s_ref| <= 0.03 and |s_R - s_glob| <= 0.05; **path-specific**
  if |s_R - 1| - max(|s_ref - 1|, |s_glob - 1|) >= 0.05; **unclear** otherwise; **untestable** if
  s_ref or s_glob is missing; **not-needed** if the width grade is "same".
- **Headline**: onset-advance + source -> "onset-advance-from-source-width" (the 160 ms would be
  explained by how a broader 2009 source pulse was read); onset-advance + path-specific ->
  "onset-advance-path-specific" (supports a structural change on the SKS path, under an onset
  reading); onset-advance otherwise -> "onset-advance-unexplained"; partial ->
  "partial-onset-advance"; no-onset-advance -> "no-advance-in-either-reading".

Everything else (f = 0.1 and 0.3, naive and cons rotations, vs vertical P, SKS_Z, onset-minus-
centroid at global stations) is reported, not graded.

## Forecasts (before any data)

| # | Question | Forecast |
|---|---|---|
| F1 | Gates pass | 85% |
| F2 | Onset grade: onset-advance / partial / no-onset-advance | 20 / 25 / 55% |
| F3 | Width grade: broader / slight / same / narrower | 25 / 30 / 35 / 10% |
| F4 | If not "same": source / path-specific / unclear (untestable folded in) | 45 / 15 / 40% |
| F5 | SKS_Z onset change within 30 ms of SKS_R's (meta) | 60% |

## What each outcome would mean

- **no-advance-in-either-reading**: neither a centroid nor an onset reading of the public IRIS
  records gives the reported advance; together with P24 this would locate the 160 ms in the
  authors' processing or data version, not in the waveforms.
- **onset-advance-from-source-width**: the 160 ms is real as a reading but comes from a broader
  2009 source pulse that shows in P and S too; it would not be evidence for a moving ULVZ.
- **onset-advance-path-specific**: only SKS broadens; this would support the paper's structural
  interpretation and would be worth a careful follow-up (more stations, more doublets).

## Run

Job `sdv-width25-v1` (data_key sdv23, 1 GB, network on): `python sdvwidth25.py --data-dir /data
--out-dir /out --wait-min 120 --time-budget-min 90`, then `python -m unittest -v
test_sdvwidth25`. Outputs: sdvwidth25_summary.json (grades), sdvwidth25_rows.csv,
sdvwidth25_globalP.csv, sdvwidth25_lags.csv, plots/sdvwidth25.png.

## Limits

One station and one doublet. Envelope-threshold onsets are one of many onset definitions; we fix
three fractions in advance and grade only one. The stretch factor of a band-limited pulse
under-states a change of the source duration, and P at 0.4-2 Hz carries more high-frequency
energy than SKS, so source and path widths are compared with tolerances, not exactly.
Directivity can make apparent source duration depend on azimuth; the SDV P and S references share
the SKS azimuth for that reason.

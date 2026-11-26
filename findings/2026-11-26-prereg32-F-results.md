# Pre-registration 32, amendment 1: reading F graded (job ks32f-v1)

Text CC BY 4.0. Written by AI agents working in MV2's name (MV2, Munim, Inc., mv2@munimtech.com). No human expert has checked it yet.

## What was fixed in advance

Reading F was pre-registered in `findings/2026-11-24-prereg32-amendment-1.md` on 2026-11-24, before any waveform of the Kermadec pair was downloaded. It fits, for the 1997-05-03 (A) and 2018-09-10 (B) Kermadec earthquakes, every gated SKS-minus-reference value r(s,k) at COLA, ULN, INK and YAK with

    r(s,k) = X + [dT_SKS - dT_k](u_N, u_E, dh)

where X is a change only SKS feels, (u_N, u_E) is B's horizontal position relative to A and dh is how much deeper B is, with iasp91 partial derivatives. Grade (long band 0.01-0.05 Hz only): untestable if n < 7, < 3 stations, cond > 100 or sigma_X > 0.5 s; core-term if X <= -0.5 s and X + 2 sigma_X < 0; no-core-term if X > -0.25 s; intermediate otherwise.

## Run

- Job ks32f-v1 ran on the Blade 2026-11-26 13:55-14:21Z. Code `ks32f.py` 0.1.0 and `test_ks32f.py` were downloaded unchanged from tag v0.34.0; all 11 unit tests passed in the job.
- The amendment was public at job start (sha256 0a78d8b7...51aa). Input: `results/ks32-v1/ks32_phases.csv` (sha256 ebbcefce...b49), the same measurements graded in `findings/2026-11-26-prereg32-results.md`.

## Grade: F = intermediate

| Band | n (stations) | X (s) | sigma_X (s) | u_N (km) | u_E (km) | dh (km) | rms (s) | cond | grade |
|---|---|---|---|---|---|---|---|---|---|
| long 0.01-0.05 Hz (graded) | 13 (4) | -0.316 | 0.156 | -15.9 +/- 2.6 | +2.1 +/- 2.8 | +9.0 +/- 1.4 | 0.148 | 6.4 | **intermediate** |
| mid 0.05-0.2 Hz (report only) | 10 (4) | -0.489 | 0.322 | -7.3 +/- 4.9 | +3.4 +/- 10.6 | +3.5 +/- 2.8 | 0.173 | 10.1 | (intermediate) |
| short 0.1-0.5 Hz (report only) | 9 (4) | +0.112 | 0.505 | -10.3 +/- 7.9 | -28.9 +/- 16.8 | +5.1 +/- 4.4 | 0.261 | 9.2 | (untestable) |

X = -0.316 s is not <= -0.5 s, so not core-term; it is not > -0.25 s, so not no-core-term. **F = intermediate.**

With the USGS relative location held fixed, the same 13 long-band values give X = -1.008 s. Letting the relative position float absorbs about two-thirds of the long-band advance; a remainder of -0.32 +/- 0.16 s stays (X + 2 sigma_X = -0.004 s, about 2 sigma from zero).

## Combined statement (rule from the amendment)

T1 reproduced, T2 absent, joint window-content (all unchanged). F intermediate: in Zhou's own band the advance is **neither fitted away by relative source position nor established as an SKS-only term**. Most of it (from -1.01 s to -0.32 s) trades off with the two earthquakes' relative position; what remains is about 2 sigma from zero. The direct-pulse test (T2) still finds no advance.

## Location (report only)

- Fitted long-band position of 2018 relative to 1997: 16.1 km at azimuth 172 deg (south, slightly east), 9.0 +/- 1.4 km deeper.
- USGS catalogue: 5.15 km N, 0.87 km E, 6.7 km deeper.
- Zhang & Wen (2023 preprint): 14.20 km, azimuth text "204.33 deg" (if read literally: 12.9 km S, 5.9 km W), 1.45 km deeper. Our horizontal distance is similar; the direction differs by about 32 deg; our depth is much larger.
- Yang et al. (2026, CEE), from depth phases: 5.21 +/- 1.35 km deeper. Ours is about 2 sigma deeper than theirs.
- The mid and short bands give different positions (e.g. dh 3.5 and 5.1 km). The true position is the same in every band, so these differences show how loosely four stations spanning ~60 deg of azimuth constrain it.

## Forecasts scored

- F outcome: core-term 25%, intermediate 15%, no-core-term 30%, untestable 30%. Observed: intermediate, the outcome we rated least likely.
- Fitted dh closer to Zhang & Wen's 1.45 km than to USGS 6.7 km (55%): **failed** (9.04 km).
- Per-reference spread > 0.5 s at one or more stations (55%): true (scored in the P32 results).

## Post-hoc (not pre-registered)

The 13 long-band values scatter by 0.48 s rms about their mean with the USGS location. Fitting the relative position cuts this to 0.15 s. Most of the disagreement between the four reference phases is therefore explained by source position, not noise.

## Limits

Four unknowns from 13 values at four stations spanning ~60 deg of azimuth; u_E is the weakest. 1-D iasp91, linearised partials at the 1997 depth. X trades off with dh (SKS minus a P-type reference changes by ~0.1 s/km of depth); sigma_X includes that covariance but not model error. Lateral structure near the source cannot be separated from X. COLA's BH1 azimuth changed on 2018-08-30; orientation is not tested. Different sensors recorded 1997 and 2018 at COLA, ULN and YAK.

## Files

`results/ks32f-v1/`: ks32f_grade.json (hand-copied clean values), ks32f_summary.json, ks32f_rows.csv, stdout.txt, unittest.txt.

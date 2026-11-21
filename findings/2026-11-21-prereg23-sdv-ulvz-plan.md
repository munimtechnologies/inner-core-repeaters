# Pre-registration 23: is the reported decadal ULVZ change at IU.SDV robust to the station's 2009 sensor re-orientation?

Written 2026-11-21, before any waveform for this test was downloaded. Code: `sdv23.py` 0.1.0
(same commit as this file). Job: `sdv-ulvz23-v1`. Licence: text CC BY 4.0, code Apache-2.0.
Made with AI agents working in MV2's name (MV2, Munim, Inc.). No human expert has checked it yet.

## The claim

Zhang & Wen 2025, *Nature Communications* 16:8754 (doi 10.1038/s41467-025-63814-3) report the
only published decadal change of an ultra-low velocity zone (ULVZ) at the core-mantle boundary.
In their catalogue of 61 doublets and 62 clusters (1993-2023), one doublet, AI-0009 (Aleutians,
2000-07-07 Mw 5.9 and 2009-12-17 Mw 5.8, relocated within 0.6 km), shows credible lowermost-mantle
changes, and only at one station, IU.SDV (Venezuela, 94.9 deg). There, after correcting for the
relative location and origin time, **radial SKS arrives 160 ms earlier in 2009 with a broader
wavelet**, and **radial S-Scd-ScS** is shifted 80 ms early at its start, 30 ms in the middle and not
at the end. The **transverse** S-Scd-ScS and the pulses after SKS do not change. They interpret
this as a ULVZ edge that shrank or moved by 16-45 km, plus opposite-sign D'' anisotropy changes.
Cui, Wang, Vidale et al. 2026 (*GRL* 53, e2025GL119930) found the CMB mostly stable over decades in
PcP/ScP/ScS of 231 repeating pairs, and noted that the ULVZ result rests on a single case.

## What we noticed (station metadata, read 2026-11-21 from EarthScope FDSN)

- IU.SDV location 00 (Streckeisen STS-1, 20 Hz) has a new channel epoch on **2009-09-18**, three
  months before the 2009 event. The horizontal azimuths change from BHN 0 / BHE 90 to **BHN 357 /
  BHE 87** (a 3 deg change). Location 10 (Guralp CMG-3T, 40 Hz) changes from 0 / 90 to **343 / 73**.
- The 00 BHN and BHE responses (evalresp, 0.25-4 Hz) are identical at 1996-01-20, 2000-07-07 and
  2009-12-17. Sensitivities are unchanged. So the only hardware-side difference between the two
  events is the horizontal orientation.
- The paper discusses station clock errors but does not mention sensor orientation or this
  metadata change. Searches for comments, replies or replications (2026-11-21) found none.

Why it could matter: the change is seen on the **radial** components only. A relative orientation
error of 3 deg mixes sin 3 deg = 5% of the transverse motion into the radial. For S-Scd-ScS, where SH
is larger than SV, that can move a radial wavelet while leaving the transverse almost untouched.
**Our prior is that it cannot explain a 160 ms SKS shift**: shifting a ~1 s pulse by 160 ms needs a
perturbation of order its own amplitude, i.e. a transverse SKS ~20 times the radial. We therefore
expect the SKS change to survive. The test is still worth running because it is cheap, it is the
first independent re-measurement of the only reported decadal ULVZ change, and the same data give
two things nobody has published: a location-corrected differential measurement and a look for
other members of the source cluster, which could date the change.

The catalogue also holds other M >= 4.5 events within 12 km of the doublet that IU.SDV could have
recorded (station open from 1994-08-19): 1996-01-20 mb 5.0 (about 2 km from the 2000 event in the
catalogue), 1998-11-08 Mw 5.1, 1998-11-09 mb 4.5, 2000-07-11 mb 4.6, 2000-12-25 mb 5.0,
2007-08-02 Mw 6.7. None is in AI-0009 as published (a doublet).

## Method (all constants in `sdv23.PARAMS`)

Master event E2000a (2000-07-07). Partner E2009a (2009-12-17). TauP iasp91, source depth 31 km,
all predicted times from the master location; every event's windows are cut at the same times
after its own catalogue origin time.

**A. Relative sensor rotation (delta).** Long-period velocity (bands 0.01-0.05 and 0.02-0.06 Hz),
window from S - 60 s to dist/2.8 km/s, after aligning the events on long-period Z (+/-10 s): grid
search (0.1 deg, +/-30 deg) for the rotation of the 2009 horizontals that best matches the 2000
horizontals (`orient_pair`). delta = mean of the two bands. Check: the Rayleigh-Z method
(`orient_rz`, 0.02-0.05 Hz, group velocity 4.2-3.0 km/s) gives each event's absolute azimuth of
the N channel; their difference must agree with delta within 2 deg when both have |c| >= 0.8.

**B. Differential times at IU.SDV, 0.4-2 Hz displacement**, FFT-upsampled to 200 Hz.
Rotations: **meta** (each event with its own metadata azimuths), **naive** (both with N = 0,
E = 90), **cons** (2000 with metadata; 2009 with the 2000 azimuths + delta). Windows are set once
on the master (metadata rotation) around the largest peak within +/-5 s of the prediction: SKS_R
[-2, +3] s, postSKS_R [+5, +11] s, SKS_Z [-2, +3] s, S_T [-3, +10] s (peak within -5..+8 s of S),
S_R_early [-2, +1.5] s, S_R_mid [+1.5, +5] s, P_Z [-1, +4] s. Lags by normalised cross-correlation
(+/-0.6 s, parabolic refinement); lag > 0 means later in the second event.
Primary quantity: **x = lag(SKS_R) - lag(S_T)**, the change of the SKS minus S differential time
(negative = SKS earlier in 2009). Using the unchanged transverse S as the reference replaces the
paper's absolute origin-time correction and cancels origin-time and clock errors.
Location correction: relative relocation (origin time and horizontal shift, depth fixed) from
0.8-2 Hz mantle P lags at IU/II stations at 30-90 deg that pass cc >= 0.90; predicted change
-(p_SKS - p_S)(dx . e_SDV) is subtracted (needs >= 4 stations; otherwise uncorrected and flagged).

**C. Membership.** Each candidate vs the master: P on BHZ at IU/II stations at 30-90 deg (up to 30,
open at both times), window [-1, +5] s around the master's P peak, lag +/-3 s. Member: cc >= 0.90 at
>= 3 stations. Strict: cc >= 0.95 at >= 4 stations (reported only).

**S2. Did the sensor physically rotate on 2009-09-18?** Rayleigh-Z azimuth of the N channel for
global earthquakes at 25-150 deg, depth <= 150 km: "before" 2007-09-18 to 2009-09-17 (M >= 6.5),
"after" 2009-09-18 to 2010-01-16 (M >= 6.0), for locations 00 and 10; plus "era2000"
1999-11-06 to 2001-06-30 (M >= 6.8, location 00). Usable if |c| >= 0.8.

## Decision rules (pre-registered)

**Primary (E2000a / E2009a), `grade_primary`.** Gates, all required, else **untestable**:
G1 the partner is a member; G2 SDV data for both events; G3 S_T cc >= 0.90 (meta); G4 SKS_R SNR
>= 5 in both events (peak in the SKS window over the RMS 120-20 s before P) and SKS_R cc >= 0.80
under all three rotations; G5 delta valid (cc >= 0.90 in both bands, band spread <= 1.0 deg, and
agreement with Rayleigh-Z within 2 deg when available).
With x location-corrected:
- **not-reproduced**: neither x_meta nor x_naive <= -0.10 s (we cannot reproduce an SKS advance of at
  least 100 ms with either common workflow);
- else **orientation-explains**: |x_cons| <= 0.04 s;
- else **change-survives**: x_cons <= -0.10 s;
- else **partly-orientation**.

**History (secondary), `grade_history`**, for each other member with the same gates (G3, G4 on the
rotation used; cons if delta is valid, else meta): |x| <= 0.03 s **stable**, |x| >= 0.06 s
**changed**, else **unclear**. A 1996/1998 member that is "changed" relative to 2000 would date part
of the change before 2000; a 2000-12-25 member is a near-zero-time control and should be stable.

**S2, `grade_rotation`** per location: difference of medians (after - before) of the wrapped
azimuths, 2,000 bootstrap draws (seed 23), 5-95% interval. Needs >= 5 before and >= 3 after.
**physical-rotation** if the interval contains the metadata step (-3 deg for 00, -17 deg for 10) and
not 0; **no-physical-rotation** if it contains 0 and not the step; else **unclear**.

Reported, not graded: x for the P_Z reference, SKS_Z (orientation-proof), postSKS_R and the radial
S windows under each rotation; uncorrected x; delta and Rayleigh-Z values; station P cc.

## Forecasts (before data)

- H1 primary testable: 55%.
- H2 given testable: not-reproduced 35%, change-survives 50%, partly-orientation 10%,
  orientation-explains 5%.
- H3 delta: |delta + 3| <= 1 deg (physical rotation matching the metadata) 45%; |delta| <= 1 deg
  (no physical rotation) 40%; other 15%.
- H4 S2 location 00: physical-rotation 40%, no-physical-rotation 35%, unclear/untestable 25%.
- H5 at least one other member verified: 40% (1996-01-20 most likely).
- H6 the radial S_R_early shift changes by >= 30 ms between meta and cons (orientation matters for
  S-Scd-ScS even if not for SKS): 30%.

## Limits known in advance

One station and one doublet; our reference phase (transverse S) is assumed unchanged, as the
paper reports. Depth differences are not corrected (about 4 ms per km for SKS - S). Our windows and
band follow the paper but not its exact processing (it interpolated to 1 kHz and used absolute
times). The "naive" and "meta" rotations stand for common workflows; we do not know which one the
authors used.

## Prior work

Zhang & Wen 2025 (above); Cui et al. 2026 GRL (above); Wang et al. 2024 GRL "No temporal change
seen in high-frequency waves scattered near the core-mantle boundary"; Ringler et al. 2013 SRL
84:926 (GSN orientation errors; magnetic methods dropped in 2009). Repository context: the same
source-verification idea as pre-registration 11 (inner-core change not at the source).

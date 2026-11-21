# Pre-registration 23 results: the reported 160 ms SKS advance at IU.SDV is not reproduced in the public records

*inner-core-repeaters, MV2 of Munim, Inc. Graded 2026-11-21 from job `sdv-ulvz23-v2` (sdv23.py
0.1.1; started 2026-11-21 16:58 UTC, finished 17:17 UTC, all steps exit 0, 26 unit tests passed,
peak memory 0.23 GB). Plan: findings/2026-11-21-prereg23-sdv-ulvz-plan.md and its amendment 1
(both public on main before the job fetched any waveform: `plan_public_at_start` and
`amendment_public_at_start` true; sha256 6f63064d... and d236592f...). Results:
results/sdv23-v2/. Licence: code Apache-2.0, text CC BY 4.0. Made with AI agents working in
MV2's name; no human expert has checked it yet.*

## The claim being tested

Zhang & Wen 2025 (*Nature Communications* 16:8754) report the only published decadal change of an
ultra-low velocity zone (ULVZ) at the core-mantle boundary. It rests on one Aleutian doublet
(2000-07-07 Mw 5.9 and 2009-12-17 Mw 5.8) recorded at one station, IU.SDV in Venezuela (94.9 deg).
After correcting for relative location and origin time, radial SKS arrives **160 ms earlier** in
2009, with a broader wavelet; the radial S-Scd-ScS is 80 ms early at its start and 30 ms in the
middle; the transverse S-Scd-ScS does not change (0.4-2 Hz displacement). The station metadata
change the horizontal sensor azimuths by 3 deg on 2009-09-18, between the two events.

## Verdict

**Primary grade (pre-registered rule `grade_primary`): not-reproduced.** All five gates passed.
The change of the SKS-minus-S differential time between 2000 and 2009 (x = lag of radial SKS minus
lag of transverse S; negative means SKS earlier in 2009), after the location correction, is:

| Rotation workflow | x (ms) | x without location correction (ms) | Radial SKS cc |
|---|---|---|---|
| meta (each event with its own metadata azimuths) | **+28.4** | +28.7 | 0.987 |
| naive (N = 0, E = 90 for both) | **+44.6** | +44.9 | 0.989 |
| cons (2009 rotated by the measured relative rotation) | +47.5 | +47.8 | 0.990 |

The rule needed x <= -100 ms under the meta or the naive workflow to count as reproduced. Neither
comes close: the measured change has the opposite sign and is a fifth of the size. With the
public IRIS records and standard processing, we do not see an SKS advance in 2009.

## Gates

| Gate | Value | Pass |
|---|---|---|
| G1 partner is a cluster member (P cc >= 0.90 at >= 3 stations) | 13 stations >= 0.90, 12 >= 0.95 (strict member) | yes |
| G2 SDV data for both events | both | yes |
| G3 transverse S cc (meta) >= 0.90 | 0.991 | yes |
| G4 radial SKS SNR >= 5 and cc >= 0.80 (all rotations) | SNR 23.7 / 20.8; cc 0.987-0.990 | yes |
| G5 relative rotation valid | 0.2 / 0.9 deg in the two bands (cc 0.959 / 0.979), spread 0.7 deg; Rayleigh-Z check not available (c 0.78 / 0.66 < 0.8) | yes |

Relative relocation from 13 mantle-P stations: origin-time offset -2.088 s, horizontal shift
0.55 km south and 0.22 km east (rms 41 ms). Its predicted effect on x is 0.3 ms. The coarse
alignment on transverse S (amendment 1) was -2.17 s, consistent with the relocated origin-time
offset.

## Every window (2009 minus 2000, ms, location-corrected)

| Window | vs transverse S, meta | vs transverse S, cons | vs vertical P, meta | vs vertical P, cons |
|---|---|---|---|---|
| Radial SKS | +28.4 | +47.5 | -13.1 | +2.0 |
| Vertical SKS (orientation-proof) | +23.7 | +27.6 | -17.8 | -17.8 |
| Radial post-SKS | +39.0 | +56.4 | -2.5 | +10.9 |
| Radial S-Scd-ScS start | -38.6 | -32.4 | -80.1 | -77.9 |
| Radial S-Scd-ScS middle | +15.0 | +16.2 | -26.5 | -29.2 |
| Vertical P | +41.5 | +45.5 | 0 | 0 |

Against either reference, the radial SKS of 2009 moves by -13 to +48 ms, not -160 ms. The radial
S-Scd-ScS start is the one window that moves early in 2009 under every reference (-32 to -80 ms),
in the direction the paper reports (-80 ms); the middle window does not match the paper's -30 ms
against transverse S. These are reported, not graded.

## Orientation

- **Relative rotation between the two events: 0.55 deg** (waveform matching of 0.01-0.06 Hz
  horizontals). The metadata step is -3 deg. So the sensor did not physically turn by 3 deg
  between the events, or it turned back; using each event's metadata azimuths introduces a
  3-deg relative error. It moves x by about 19 ms (meta 28.4 vs cons 47.5 ms), far too little to
  make or hide a 160 ms shift. Our prior (orientation cannot explain 160 ms in SKS) holds; the
  question became moot because there is no 160 ms shift to explain.
- **S2 (did the sensor physically rotate on 2009-09-18?)** Location 00: **unclear** (after-minus-
  before difference 6.5 deg, 5-95% -7.0 to 20.5; only 3 usable events after, 35 before). Location
  10: **no-physical-rotation** (-0.5 deg, 5-95% -12.0 to 6.5; 17 after, 41 before; metadata step
  -17 deg).
- Not graded: location 00 long-period horizontals give very low Rayleigh-Z coherence (c often
  0.03-0.4) for many events between 2009-03-19 and 2009-12-14 (with a few usable ones in August),
  while location 10 mostly stays coherent for the same events; location 00 is coherent again from
  2009-12-19. The 2009-12-17 event itself has strong, clean 0.4-2 Hz SKS (SNR 21) at location 00.

## Other members of the source cluster (history, secondary)

Six other M >= 4.5 events within 12 km were checked. Only 1996-01-20 (mb 5.0) is a member (5
stations at cc >= 0.90; relocated), but its SKS at SDV is too weak (SNR 2.2) and its transverse S
cc 0.63, so every history grade is **untestable**. The other five are not members (0-1 stations at
cc >= 0.90). No second event can date a change.

## Forecasts (made before any waveform)

| # | Forecast | P | Outcome | Right? |
|---|---|---|---|---|
| H1 | Primary testable | 55% | testable | happened |
| H2 | not-reproduced / survives / partly / explains | 35 / 50 / 10 / 5 | not-reproduced | our second choice |
| H3 | delta ~ -3 / ~ 0 / other | 45 / 40 / 15 | 0.55 deg (~ 0) | our second choice |
| H4 | S2 loc 00 physical / none / unclear | 40 / 35 / 25 | unclear | our third choice |
| H5 | another member verified | 40% | 1996-01-20 member (untestable at SDV) | happened |
| H6 | S_R_early moves >= 30 ms between meta and cons | 30% | 6 ms | did not happen (leaned right) |

## What it means, and what it does not

1. **With the public records and standard processing, the 2009 radial SKS at IU.SDV is not 160 ms
   early.** The two records are nearly identical (cc 0.99), and the SKS-minus-S and SKS-minus-P
   differential times change by -13 to +48 ms depending on reference and rotation.
2. **This does not yet show the paper is wrong.** The paper does not say how the 160 ms was
   measured (no picker or correlation method is described); ours is a 5-s cross-correlation of
   the whole SKS pulse. An onset-based reading of a pulse that is broader in 2009 could differ.
   The decisive check is pre-registration 24 (queued): the same measurement on the authors' own
   processed traces (Zenodo 10.5281/zenodo.16421554). If their traces show about -160 ms and the
   IRIS records do not, the difference lies in processing or data versions; if theirs also show
   about +30 ms, the change comes from how it was read, not from the waveforms.
3. **The 3-deg metadata re-orientation is not the explanation** of anything here: the measured
   relative rotation is 0.55 deg, and orientation moves x by about 19 ms.

## Limits

One station, one doublet; the transverse S reference is assumed unchanged (as the paper
reports); depth differences are not corrected (about 4 ms per km for SKS - S; the paper's relative
depth is 0.05 km). Cross-correlation lags over fixed windows, not onset picks. The S2 test at
location 00 is weak (3 usable events after the epoch).

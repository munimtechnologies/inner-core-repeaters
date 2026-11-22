# Pre-registration 24, amendment 1: pair the authors' files like with like

Written 2026-11-22, after job `sdv-zframe24-v1` returned and **before** the corrected job
`sdv-zframe24-v2` runs. Code: `sdvz24b.py` 0.1.0, which runs the pinned `sdvz24.py` 0.1.0 and
`sdv23.py` 0.1.1 (commit dbbd91d, release v0.21.0) with one function replaced. Licence: text
CC BY 4.0, code Apache-2.0. Made with AI agents working in MV2's name (MV2, Munim, Inc.); no
human expert has checked it yet.

## Why

The plan groups the authors' traces into sets with "same file name apart from the component
code". The 0.1.0 code did not implement this for file names like `2009.R.sac.bp`. For 2009 the
plain and band-passed files have the same sampling rate, so they were merged, and the graded
pair became plain 2000 against band-passed 2009 (findings/2026-11-22-prereg24-results.md).
This amendment makes the code do what the plan says. It is not a new test.

## What we have already seen (disclosed)

From v1 we know the zip's file list and headers. We also know the frame-fit angles and lags of
three sets: plain 2000, band-passed 2000 and band-passed 2009 (results/sdvz24-v1/sdvz24_fits.csv).
From those, the post-hoc band-passed-pair rho is -2.22 deg. **We have not seen** any fit of
the plain 2009 traces or any like-with-like measurement of the SKS change.

## Changes (nothing else changes)

1. **Grouping.** `set_key_v2` replaces the component token of the file name with `?`
   (`2009.R.sac` -> `2009.?.sac`, `2009.R.sac.bp` -> `2009.?.sac.bp`). Names with a channel code
   keep the 0.1.0 rule.
2. **Primary pair (graded): the plain files** `2000.?.sac` / `2009.?.sac`. They are the first
   sorted full-coverage sets under the unchanged primary rule. New gate **A1**: both primary sets
   must end in `.sac`. Otherwise every grade is untestable.
3. **Secondary pair (reported with the same grades, not the headline):** the band-passed files
   `2000.?.sac.bp` / `2009.?.sac.bp`. The `.bp` sets are made to sort first.
4. **Combined grade** uses pre-registration 23's long-period delta (+0.55 deg; its G5 passed),
   as the plan already prescribes. The 0.1.0 code only used the Rayleigh-Z delta.
5. **Reported:** for each event, component (R, T, Z) and window (whole, SKS -5..+10 s,
   S -5..+15 s), the lag and signed correlation of the authors' band-passed trace against their
   plain trace (both through the common 0.4-2 Hz zero-phase band-pass, 200 Hz grid, +/- 2 s).
   Also the 2009-minus-2000 difference of those lags. A difference of tens of ms in the SKS
   window but not the S window would mean their band-pass alone moves SKS differently in the
   two events.

All thresholds, windows and grade rules are those of the plan: frame same / metadata / other
(1 deg); combined consistent <= 1 deg, misrotated >= 2 deg; measurement reproduced if x0 <=
-0.10 s with G1 and G4.

## Forecasts (before v2)

- A1 passes: 90%.
- Primary frame: metadata-frame 55%, same-frame 30%, other 15%. (The plain 2000 traces match
  nominal north/east, A = 0.17 deg; the IU.SDV metadata for 2000 are also 0/90.)
- Primary combined, given metadata-frame: misrotated 80%.
- Primary measurement: not-reproduced 80%, reproduced 10%, untestable 10%. Pre-registration 23
  measured +28 ms on the public records rotated with the metadata.
- Secondary (band-passed) measurement: not-reproduced 65%, reproduced 20%, untestable 15%.
- |2009 - 2000| band-pass delay difference in the SKS window >= 30 ms: 15%.

## Limits

As in the plan. The band-passed 2000 files are 1000 Hz and the 2009 ones 100 Hz; both go onto
the same 200 Hz grid. If the plain pair fails a gate, the band-passed pair does not replace it
as the graded result.

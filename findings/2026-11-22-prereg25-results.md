# Pre-registration 25 results (job sdv-width25-v1): graded "untestable" (gate G2 failed: no envelope onset can be read on radial SKS)

*inner-core-repeaters, MV2 of Munim, Inc. Graded 2026-11-22 from job `sdv-width25-v1`
(sdvwidth25.py 0.1.0, release v0.23.0, with sdv23.py 0.1.1 at pinned commit dbbd91d; started
2026-11-22 04:18 UTC on the Blade, finished 04:28 UTC, all steps exit 0, 16 unit tests passed).
Plan: findings/2026-11-21-prereg25-sdv-onset-width-plan.md (public before the job ran:
`plan_public_at_start` true, sha256 7f780251...). Results: results/sdvwidth25-v1/. Licence: code
Apache-2.0, text CC BY 4.0. Made with AI agents working in MV2's name; no human expert has
checked it yet.*

## Short version

The question was whether an **onset** reading of the IU.SDV records (rather than our whole-pulse
cross-correlation in P23) gives the 160 ms SKS advance reported by Zhang & Wen 2025, and whether
any broadening of the 2009 SKS pulse comes from the earthquake source or from the path.

**Pre-registered headline: untestable.** Gates G0 (our pipeline reproduces P23's x = +28.4 ms
exactly) and G1 (radial SKS correlation 0.987 >= 0.80) passed. Gate G2 failed: the envelope of
radial SKS never falls below 20% of its peak within the 3 s before the peak, in either record,
so the primary onset (f = 0.2) does not exist under the plan's definition. By the plan's rules
the onset, width and source grades are therefore all "untestable". We do not change the onset
definition after seeing the data.

The stretch factors, which the plan reports in every case, are listed below **without a grade**.

## Why G2 failed

`onset_index` scans back at most 3 s from the envelope peak for the first sample below f x peak.
For radial SKS (meta, naive and cons rotations) it returned no onset at **f = 0.1, 0.2 and 0.3**.
That means the envelope stays at or above 30% of its peak for the whole 3 s before the peak in
both records (results/sdvwidth25-v1/sdvwidth25_rows.csv). The plot
(results/sdvwidth25-v1/plots/sdvwidth25.png) shows why: at 0.4-2 Hz displacement the radial SKS
pulse rides on a longer lead-in swing (about 1428-1430 s after the 2000 origin) instead of rising
from a quiet background. The same happened for vertical SKS and for the direct transverse S pulse.
Only vertical P gave onsets. The plan's Limits section named the onset definition as the weak
point; it is weaker than we expected. An onset pick needs a quiet stretch before the arrival, and
at this station and band radial SKS does not have one.

## Reported, not graded: stretch factors (s > 1 means the 2009 pulse is broader)

| Window at IU.SDV (meta rotation) | s | best correlation | correlation at s = 1 |
|---|---|---|---|
| radial SKS (SKS_R) | **1.040** | 0.9949 | 0.9870 |
| vertical SKS (SKS_Z) | 1.015 | 0.9685 | 0.9653 |
| transverse direct S (S_T_pulse) | 1.005 | 0.9981 | 0.9980 |
| vertical P (P_Z) | 1.020 | 0.9771 | 0.9696 |

Radial SKS under the other rotations: naive 1.035, cons 1.030.

**Global vertical P** at the stations that verify the doublet (sdvwidth25_globalP.csv): 13 measured,
12 with best correlation >= 0.90 (II.BORG, 0.86, excluded). Median s = **1.000**, bootstrap 5-95%
range 0.9975-1.015; individual values 0.99 to 1.02.

The reference median at SDV is s_ref = median(1.020, 1.005) = 1.0125.

For information only (not a grade, because G2 failed): if the gates had passed, the coded rules
would have given width "slight" (1.03 < 1.040 < 1.05) and source "source" (|s_R - s_ref| = 0.0275
<= 0.03 and |s_R - s_glob| = 0.040 <= 0.05). The path-specific rule needs an excess broadening of
SKS of at least 0.05 over the references; the observed excess is 1.040 - 1.0125 = 0.0275.

## Reported, not graded: how stable envelope onsets are where they exist

- Vertical P at SDV: onset minus centroid change (2009 - 2000) is +117 ms at f = 0.1, -98 ms at
  f = 0.2 and -123 ms at f = 0.3 (2000 rise time at f = 0.2: 1.35 s). The stretch fit predicts
  -27 ms at f = 0.2.
- Global vertical P (f = 0.2, 9 stations with correlation >= 0.90 and onsets in both records):
  onset minus centroid ranges from -1.515 s (IU.KBS) to +0.045 s (IU.INCN); median -0.063 s.
  At 4 of the 9 stations it exceeds 0.2 s in size, although the waveforms correlate at 0.95-0.99.

So envelope-threshold onsets on these 0.4-2 Hz pulses move by a hundred milliseconds or more
between thresholds and stations even for a verified doublet. An onset-based change of 160 ms at
one station would need its own demonstration that it is larger than this picking scatter.

## Forecasts (made before the data)

| # | Question | Forecast | Outcome |
|---|---|---|---|
| F1 | Gates pass | 85% | **no** (G2) |
| F2 | Onset grade | adv 20 / partial 25 / none 55% | untestable |
| F3 | Width grade | broader 25 / slight 30 / same 35 / narrower 10% | untestable (s_R = 1.040 reported) |
| F4 | Source grade | source 45 / path 15 / unclear 40% | untestable |
| F5 | SKS_Z onset within 30 ms of SKS_R | 60% | untestable (no onsets) |

F1 failed. We gave 15% to a gate failing and did not foresee this cause.

## What this does and does not say

- It does not say whether an onset reading gives 160 ms; under our pre-registered definition no
  onset exists on radial SKS at this station and band.
- The ungraded stretch factors show at most a 4% broadening of the 2009 radial SKS, against
  0.5-2% for P and S at SDV and 0% (median) for P at 12 other stations. Under the plan's own
  arithmetic, a 4% broadening of a pulse with a rise time near 1 s moves an onset by about 40 ms,
  not 160 ms. This is a back-of-envelope remark, not a test result.
- The pending tests are unchanged: P24 v2 (the authors' own traces), P26 (entry-side neighbour
  stations) and P27 (other pairs sharing the exit side). No new SDV test is filed for this
  outcome.

## Files

results/sdvwidth25-v1/: sdvwidth25_summary.json (gates and grades), sdvwidth25_rows.csv
(per-window stretch and onset rows), sdvwidth25_lags.csv (window lags), sdvwidth25_globalP.csv
(global P stretch), plots/sdvwidth25.png.

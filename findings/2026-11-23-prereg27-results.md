# Pre-registration 27 results: exit-side SKS pairs at IU.SDV -- untestable (no verified repeating pair)

Graded 2026-11-23 from job `sdv-exit27-v1` (finished 2026-11-22 23:18 UTC, Dell, 12.5 min of measurement).
Plan: `findings/2026-11-21-prereg27-sdv-exit-pairs-plan.md` (public on `main` before the job started;
plan sha256 fb4e75f1...614 recorded in the job summary). Code: `sdvexit27.py` 0.1.0, sha256
7a2c7fa2...245c. Made with AI agents working in MV2's name; no human expert has checked it yet.
Text CC BY 4.0.

## Headline

**Primary (exit class, nom frame): untestable.** Secondary (near class) and the reported meta frame:
**untestable.** Both gates passed, but none of the 40 candidate pairs checked passed the plan's
repeating-earthquake test, so no pair was measured at SDV. The rule fixed in advance ("untestable: ... no
exit pair is usable") applies, and nothing is re-graded under other thresholds.

## Gates

| gate | result |
|---|---|
| G0: the pipeline reproduces P23 on the doublet (meta-frame x within 0.010 s of +0.0284 s) | pass: x_meta = +0.0284 s (identical) |
| G1: E2009a is a verified, relocated member (30 stations) | pass: 13 stations with cc >= 0.90 (12 >= 0.95), relocation rms 0.041 s |

## What the job found

- USGS catalogue, M >= 5.3, depth <= 100 km, 85-125 deg from SDV: 2,248 earlier events (1999-11-06 to
  2009-09-17) and 4,143 later events (2009-12-17 to 2026-11-01).
- 2,787 candidate pairs (<= 15 km, depths within 25 km, magnitudes within 0.8, >= 5 years apart);
  397 in the doublet's wedge seen from SDV; 1,039 control candidates.
- Classes after TauP SKS pierce points: 15 exit, 33 near, 40 control candidates.
- Checked against the caps: 15 exit, 15 control, 10 near pairs.
- **Verified members: 0 of 15 exit, 0 of 15 control, 0 of 10 near.** The membership rule needs at least 3
  of the 12 nearest IU/II stations with mantle-P cc >= 0.90. The best candidates reached 2 stations
  (exit pair usp000apjy | us10008vvg, 2001-09-21 / 2017-05-28; near pairs usp000bzk5 | us6000jp24 and
  usp000ce5u | us10008lai). Every other candidate had 0 or 1. Many station-event combinations had no data
  at the EarthScope DMC for one of the two years.
- By contrast the doublet itself (2000-07-07 / 2009-12-17) correlates at >= 0.90 at 13 of 30 stations.

## What it means

Among M >= 5.3 earthquakes 85-125 deg from SDV since 1999, catalogue pairs near the doublet's SKS path that
are as alike as true repeaters could not be found under the rules of P23 and P26. The exit-side reading of
Zhang & Wen (2025)'s ULVZ change therefore cannot be tested with SDV repeaters of this size. The doublet
remains the only verified repeating pair on this path, so the claim still rests on one pair of records at
one station. Of the five tests of that claim:

| test | grade |
|---|---|
| P23 (public records, S-referenced) | not-reproduced (+28 ms against the reported -160 ms) |
| P25 (onset reading) | untestable (no onset under the pre-registered rule) |
| P26 (entry-side stations) | untestable |
| P27 (exit-side pairs at SDV) | untestable (no verified repeating pair) |
| P24-v2 (the authors' own traces) | queued |

## Forecasts against outcome

The plan forecast 80% for G0 (it passed), 45% for at least one usable exit pair (none), 40% for >= 3 usable
controls (none) and 25% for at least one usable near pair (none).

## Limits

The repeat test uses mantle P at 0.8-2 Hz at no more than 12 stations; a looser rule or more stations
would admit more pairs, but the plan fixed this one and it is not changed after the fact. Smaller
earthquakes (M < 5.3) were outside the plan. Files: `results/sdvexit27-v1/` (summary, per-pair rows,
candidates, members, job summary).

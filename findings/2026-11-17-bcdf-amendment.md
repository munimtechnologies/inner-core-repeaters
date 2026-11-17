# Amendment to the BC−DF pre-registration (bcdf.py v0.4.0 → v0.5.0)

By MV2, Munim, Inc. (made with AI agents). Written 2026-11-17, **before any BC−DF time was
measured** by `bcdf.py` (it has never run). Text CC BY 4.0. Amends
[`2026-11-17-bcdf-prereg.md`](2026-11-17-bcdf-prereg.md).

## What we saw, and what we did not see

We had read the outputs of `icr-ilar-v1` ([results](2026-11-17-ilar-results.md)). Those
include the DF onset picked on each multiplet stack relative to ak135, and the waveform
plots. They contain **no BC−DF times and no single-event time series**, so nothing in them
favours the monotone or the backtracking clock in the BC−DF test.

They show that on the SSI → ILAR path, DF arrives a median 4.35 s before ak135 (72 of 75
multiplets were early) while PKPbc is close to on time. The observed BC−DF is therefore
about 4 s longer than ak135. That is the long-known SSI-Alaska anomaly (DF 1-4 s early,
about 6 s of BC−DF scatter, partly from the Alaska slab).

## Why v0.4.0 would have failed

* DF was searched only within ±4 s of ak135, and events with |onset| > 3 s were rejected.
* BC was searched at onset + ak135(BC−DF) + [−0.6, +1.5] s, about 4 s before the real BC.
  That window falls in the weak DF coda, so most events would have been rejected or would
  have locked onto coda wiggles.

## The change (v0.5.0)

| Setting | v0.4.0 | v0.5.0 |
|---|---|---|
| DF onset search | ak135 ± 4 s | ak135 − 9 s to + 1 s |
| usable DF onset | \|onset\| ≤ 3 s | −8 s ≤ onset ≤ +1 s |
| BC search | DF onset + ak135(BC−DF) + [−0.6, +1.5] s | BC's own ak135 time ± 1.5 s (still ≥ 0.8 s before ak135 AB, never overlapping the DF template) |

Everything else is unchanged: the template, the quality cuts (cc ≥ 0.75, DF SNR ≥ 3), the
±1 s outlier screen around the median residual, the spatial terms, the L/Bk fits on events
before 2025.0, and the Q1-Q3 decision rules and forecasts. New outputs, all descriptive:

* `phase_timing`: median DF onset and BC−DF residual, overall and before/after the
  2018-09-25 ILAR epoch break;
* `bcdf_examples.png`: beams with the picked DF and the measured BC marked, so the phase
  identification can be checked by eye.

## Also applies to `cola.py` (v0.4.0 → v0.5.0)

`cola.py` measures College, Alaska (40 km from ILAR, same SSI path) with
`bcdf.measure_bcdf` and had the same `|onset| ≤ 3 s` rule. v0.5.0 uses the corrected
`bcdf.py`, the same onset rule (−8 to +1 s), a wider example plot (−10 to +16 s) and
reports `phase_timing`. The C1 calibration, Q1-Q3 and the forecasts in
[`2026-11-17-cola-prereg.md`](2026-11-17-cola-prereg.md) are unchanged. No COL/COLA
waveform had been measured: the v1 job never ran, and v2 replaces it.

## Tests

A new synthetic test (`test_ssi_ilar_geometry`) reproduces the observed geometry: DF 4.6 s
early, BC 0.4 s early, ak135 BC−DF 7.1 s, true BC−DF 11.3 s. v0.5.0 recovers 11.3 s, a
value outside v0.4.0's search window.

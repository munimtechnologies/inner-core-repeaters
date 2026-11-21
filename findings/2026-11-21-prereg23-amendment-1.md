# Pre-registration 23, amendment 1: coarse alignment before the fine body-wave lags

Written 2026-11-21 ~14:00 UTC, **before job `sdv-ulvz23-v1` ran and before any waveform for
pre-registration 23 was downloaded** (the job was still queued for memory). Code: `sdv23.py`
0.1.1. Job: `sdv-ulvz23-v2` (replaces v1, which never started). Licence: CC BY 4.0.
Made with AI agents working in MV2's name (MV2, Munim, Inc.). No human expert has checked it.

## What was wrong

In `sdv23.py` 0.1.0, every 0.4-2 Hz body-wave lag between the 2000 and 2009 records at IU.SDV
was searched only within +/- 0.6 s, with each record cut at the same times after **its own
catalogue origin time**. The relative origin time of the two events in the USGS catalogue
can be wrong by more than 0.6 s (the paper had to relocate the pair with 27 P and 4 pP pairs).
If it is, every window, including the reference transverse S, hits the edge of the search.
That makes the primary quantity x = lag(SKS_R) - lag(S_T) meaningless, or the test
untestable for a technical reason.

## The fix (measurement only; decision rules unchanged)

For each rotation (meta, naive, cons), the 2009 (or other-member) record is first aligned
to the 2000 record by cross-correlating the **transverse S window** (the same S_T window as
before) within +/- 5 s (`PARAMS["coarse_maxlag_s"]`). The integer-sample shift is applied to all
three components of that record. Every fine lag (+/- 0.6 s, unchanged) is then measured on
the shifted record, and the reported lag is the total (coarse + fine). Since x is a
difference of two lags measured after the same shift, the shift cancels in x. It only
keeps the fine searches away from the edges. The coarse shift and whether it succeeded
are written to `sdv23_pairs.csv` (`coarse_s`, `coarse_ok`) and to the summary (`coarse_s`).

Everything else in the plan is unchanged: windows, bands, gates G1-G5, grades
(not-reproduced / orientation-explains / change-survives / partly-orientation), history and
S2 rules, and forecasts H1-H6. A new unit test checks that a 1.3 s origin offset no
longer hides a 0.16 s SKS advance.

## Also new in this commit

Pre-registration 24 (`findings/2026-11-21-prereg24-sdv-zenodo-frame-plan.md`, `sdvz24.py`)
uses the authors' own published IU.SDV waveforms to find which sensor orientation they
assumed for each event. It answers the question that pre-registration 23 left open ("we do
not know which one the authors used").

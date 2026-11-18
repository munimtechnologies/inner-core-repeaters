# Pre-registration 12: a prospective YKA test (earthquakes after 2026-11-19 only)

MV2, Munim, Inc. (made with AI agents). Written 2026-11-18. Text CC BY 4.0, code Apache-2.0.
This is published before any South Sandwich Islands earthquake it will use has happened.

## Why prospective

Every YKA pair up to 2026-11-18 has now been seen. YKA stage B failed its calibration
([results](2026-11-18-ykab-v2-bcdf-v2-results.md)). The [diagnostic](2026-11-18-diag2-results.md)
found that at YKA, coda dissimilarity is driven by how far apart the two sources are, not by
the time between them. The [plan](2026-11-18-diag2-plan.md) fixed in advance what follows from
that: a new confirmatory test with a selection rule built from the supported explanation
(location), applied only to future events. The South Sandwich Islands produce about 100 USGS
events of M ≥ 4.7 a year.

## Rule

* **New events:** USGS ComCat events with origin time on or after **2026-11-19T00:00Z** in the
  same South Sandwich Islands region as the YKA catalogue of pre-registration 8, M ≥ 4.3.
* **Family membership:** a new event joins a merged ILAR family (pre-registration 8 tables) if
  its ILAR PKPbc+PKPab correlation with any member passes the unchanged `run.py`
  repeater threshold (single linkage, as before).
* **YKA measurement:** `yka_b.py` 0.10.1 with `--repick`, unchanged windows, filters and gates
  (direct cc ≥ 0.90, direct SNR ≥ 5, coda SNR ≥ 2).
* **Selection (new, from the diagnostic):** a pair counts only if its ComCat epicentral
  separation is **≤ 14.7 km** and the smaller of the two YKA direct SNRs is **≥ 11.7**. These
  are the medians of the 1,650 ok stage-B pairs. They are fixed numbers, not recomputed.
* **Primary similarity:** **D_c13** (YKA coda +1 to +3 s). It is the window where
  near-vs-far calibration held in all 8 families (sign p = 0.004). The pre-registered
  +1 to +6 s D is reported as secondary.

## Tests, in order

1. **P4-prosp (calibration, new pairs only).** Among selected pairs that contain a new event:
   near = partner from 2024-11-19 or later; far = partner from before 2016-11-19 and more than
   5 years apart on the backtracking clock. **Testable** if there are ≥ 5 near and ≥ 5 far
   pairs. **Holds** if median D_c13(near) < median D_c13(far). If it is untestable or fails,
   nothing below is interpreted.
2. **F1-prosp (the vote).** This is the F1-Y rule of pre-registration 8 applied to selected
   pairs. A case is a new event with ≥ 2 selected partners from before 2008.5. B pick: the
   partner closest to e(t) = 2008.5 − (t − 2008.5)/2.5. M pick: the latest pre-turn partner.
   Best match: the lowest D_c13. **Testable** if there are ≥ 6 discriminating cases in
   ≥ 3 families. **Supported** if n_B > n_M, one-sided binomial p < 0.025, and n_B > n_M in
   ≥ 95% of 2,000 whole-family bootstraps. **Contradicted** if the same holds with B and M
   swapped. **Inconclusive** otherwise.
3. Reported, not graded: the same with D (+1 to +6 s), and with the 2025-2026 events excluded
   from the partner pool.

## When

There are two looks only: the first run on or after **2027-11-19** (interim) and the first run
on or after **2028-11-19** (final). The α of 0.025 per look covers both. If the vote is not
testable at the interim look, it is not graded then. The grader script will be published
before the interim look, and before any YKA waveform from after 2026-11-19 is measured.

## Forecasts

| Outcome | Probability |
|---|---|
| F1-prosp testable at the interim look | 10% |
| F1-prosp testable at the final look | 30% |
| If testable and P4-prosp holds: supported / inconclusive / contradicted | 30% / 55% / 15% |

These numbers are low on purpose. In 2025-2026 only 6 discriminating cases came out of
226 events, and that was before the stricter selection.

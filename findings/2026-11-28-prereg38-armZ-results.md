# Pre-registration 38, arm Z graded: clock part untestable, Z5 mixed at YKA and agree at ILAR (gradez-v1)

MV2, Munim, Inc. (mv2@munimtech.com), 2026-11-28. Text and data CC BY 4.0, code Apache-2.0.
Made with AI agents working in MV2's name; no human expert has checked it.

Rules: [arm Z plan](2026-11-27-prereg38-armZ-plan.md) (frozen in commit a52bf19, release v0.40.0, before any call
was copied). Input: [the transcription](2026-11-28-prereg38-armZ-transcription.md) of Zhang & Wen's (2026, JGR
Solid Earth, https://doi.org/10.1029/2026JB034595, "ZW26", CC BY 4.0) explicit pairwise calls and the hand tallies
committed with it (commit 193fe69, release v0.41.0) before the job ran.

## The run

`gradez-v1` (gradez 0.1.0, unchanged) ran on the Blade on 2026-11-28, 05:32-05:39Z: tests passed, run exit 0.
Outputs are committed unmodified in [`results/gradez-v1/`](../results/gradez-v1/).

**Hand-tally check passed** (`checks_ok: true`): YKA 12 pairs / 20 comparisons, ILAR 12 pairs / 22 comparisons,
exactly as `expect.json`. The Z5 counts also match the hand table in the transcription note line for line
(YKA 8 shared binary pairs, 5 agree; ILAR 8, 7 agree; 10 shared pairs in any grade at each array). No recount
note is needed.

## Grades, exactly as the plan states

| array | Z1-Z4 / headline | Z5 (shared pairs with Wang et al. 2024 Table 2) |
|---|---|---|
| YKA | **untestable** (20 comparisons < 30) | **mixed**: 5 of 8 agree (62.5%; agree needs >= 80%, disagree <= 50%) |
| ILAR | **untestable** (22 comparisons < 30) | **agree**: 7 of 8 agree (87.5%) |

Forecasts made before reading: testable at >= 1 array 55% (failed); Z5 testable 40% (it was testable at both
arrays); given testable, agree 35%, mixed 35%, disagree 30% (one of each of the first two).

The clock part of arm Z closes **untestable**: ZW26 state explicit pairwise calls for too few pairs to test a
rotation clock on their calls alone.

## What Z5 shows (from the committed pair file, no judgement)

Every disagreement goes the same way: **ZW26 call the pair changed where Wang et al. grade it similar.**

| pair | YKA ZW26 / Wang | ILAR ZW26 / Wang |
|---|---|---|
| 2003-03-21 / 2020-07-11 | C / 0 (disagree) | C / 1 (not binary) |
| 2008-02-17 / 2011-02-20 | C / 0 (disagree) | S / 0 (agree) |
| 2020-07-11 / 2022-03-07 | C / 0 (disagree) | C / 0 (disagree) |

At ILAR the two camps' published grades agree on every shared binary pair except 2020-07-11 / 2022-03-07. At YKA
they split on three pairs. These are the same three pairs that the exploratory note
[2026-11-28-yka-camps-differ-by-threshold.md](2026-11-28-yka-camps-differ-by-threshold.md) found to lie between
the two camps' implied cuts on our published late-coda D_c310 (0.0179, 0.0209, 0.0228). Report-only, not
pre-registered: three of ILAR's seven agreements are multiplet-F pairs whose ZW26 classes are taken from Vidale et
al. (2025, Nature Geoscience, https://doi.org/10.1038/s41561-025-01642-2); on M03-24 alone ILAR agrees on 4 of 5
and YKA on 3 of 6.

## What this means

The rival readings of the South Sandwich records (one-turn rotation, Wang, Vidale et al. 2024; local episodic
surface change, ZW26) mostly do **not** disagree about which waveforms match. At ILAR they disagree on one pair;
at YKA on three, all in the direction of ZW26 calling more change. The dispute is therefore about a handful of
named pairs, and the pre-registered test that decides the YKA three against ordinary repeat scatter is arm T
([plan](2026-11-28-prereg38-armT-plan.md)), which needs the arm M measurement (`m2-v1`).

## Next

`m2-v1` is filed today (see [the filing addendum](2026-11-28-prereg38-armM-filing-addendum.md)); `gradet-v1`
follows when it returns.

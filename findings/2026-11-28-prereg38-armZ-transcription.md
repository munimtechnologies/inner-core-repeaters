# Pre-registration 38, arm Z: the transcription, hand tallies and what they already decide

MV2, Munim, Inc. (mv2@munimtech.com), 2026-11-28 (written 00:25-01:10Z). Text and data CC BY 4.0, code Apache-2.0.
Made with AI agents working in MV2's name; no human expert has checked it.

Rules: [arm Z plan](2026-11-27-prereg38-armZ-plan.md), frozen and committed (commit a52bf19, release v0.40.0)
before any call was copied. This note commits the copied calls and the hand tallies **before** the frozen job
`gradez-v1` (gradez 0.1.0, unchanged) runs on them.

## The copy

Two read-only browser copy jobs, one per array, read Zhang & Wen (2026, JGR Solid Earth,
https://doi.org/10.1029/2026JB034595, "ZW26", CC BY 4.0) under the frozen copy rules: main text, figures,
captions and the Supporting Information (15-page PDF, Figs S1-S10, no tables). Both finished on 2026-11-28
(ILAR 00:03Z, YKAR 00:10Z). Their CSV lines are committed verbatim as
[`data/zw26_pair_calls.csv`](../data/zw26_pair_calls.csv) (27 rows; credit ZW26).

Notes from the copies:
- ZW26 has **no per-pair table** of similarity values or classes. The nearest are Figs S8 (ILAR) and S9 (YKAR),
  which plot S(t) for each M03-24 pair, and Fig 6.
- ZW26's quantitative rule (Sec 4.2, quoted by the copy): windows are dissimilar where S(t) > 1.75 x N and
  S(t) > 5%; a pair is "similar" when no window is dissimilar; only pairs with N < 20% are used.
- M03-24 members as printed: 1999/07/01, 2003/03/21, 2008/02/17, 2011/02/20, 2020/07/11, 2022/03/07,
  2022/10/15, 2024/03/21. At YKAR 2022/10/15 was dropped for noise; 1999/07/01 is not used. The pair lists
  use exactly the six remaining events, as the membership list implies (cross-check passed).
- Multiplet F: 2005/06/10, 2012/12/03, 2022/05/19. **ZW26 take multiplet F's ILAR classes from Vidale et al.
  (2025)**, so those three ILAR calls are not ZW26's own reading.
- The two copies were made independently and list the **same nine M03-24 pairs** at each array.

## Hand tallies (passed to the job as `expect.json`)

Merged per pair (all S = 0, all C = 2, mixed = 1); comparisons = pairs of pairs with different grades within a
multiplet.

| array | multiplet | pairs | grade 0 | grade 2 | comparisons |
|---|---|---|---|---|---|
| YKA | M03-24 | 9 | 3 | 6 | 18 |
| YKA | F | 3 | 1 | 2 | 2 |
| ILAR | M03-24 | 9 | 4 | 5 | 20 |
| ILAR | F | 3 | 1 | 2 | 2 |

Totals: YKA 12 pairs, 20 comparisons; ILAR 12 pairs, 22 comparisons.

## What this already decides (pre-registered rule, no judgement)

The plan says an array with fewer than 30 comparisons is untestable. Both arrays fall short (20 and 22), so
**Z1-Z4 and the headline are untestable at both arrays: the clock part of arm Z closes untestable.** We
forecast "testable at >= 1 array" at 55%; that forecast failed. ZW26 state explicit pairwise calls for too few
pairs to test a rotation clock on their calls alone.

## Z5 (shared pairs with Wang et al. 2024 Table 2): expected values, derived by hand

Wang, Vidale et al. (2024, Nature 631:340, https://doi.org/10.1038/s41586-024-07536-4) grade 0 = similar,
2 = dissimilar (1 = intermediate). Matching by the two event dates at the same array:

| pair | YKA: ZW26 / Wang | ILAR: ZW26 / Wang |
|---|---|---|
| 2003-03-21 / 2008-02-17 | 2 / 2 agree | 2 / 2 agree |
| 2003-03-21 / 2011-02-20 | 2 / 2 agree | 2 / 2 agree |
| 2003-03-21 / 2020-07-11 | 2 / 0 **disagree** | 2 / 1 (not binary) |
| 2003-03-21 / 2022-03-07 | 0 / 0 agree | 0 / 0 agree |
| 2008-02-17 / 2011-02-20 | 2 / 0 **disagree** | 0 / 0 agree |
| 2011-02-20 / 2020-07-11 | 2 / 1 (not binary) | 2 / 1 (not binary) |
| 2020-07-11 / 2022-03-07 | 2 / 0 **disagree** | 2 / 0 **disagree** |
| F 2005-06-10 / 2012-12-03 | 2 / 2 agree | 0 / 0 agree |
| F 2005-06-10 / 2022-05-19 | 0 / 0 agree | 2 / 2 agree |
| F 2012-12-03 / 2022-05-19 | 2 / 1 (not binary) | 2 / 2 agree |

(2003/2024 and 2022/2024 are not in Wang's table, which ends in January 2023.) Wang rows used: P104, P105, P51,
P123, P40, P53, P126, P36, P137, P138.

Expected job output: YKA n_shared_binary 8, n_agree 5 (62.5%) -> **mixed**; n_shared_any 10.
ILAR n_shared_binary 8, n_agree 7 (87.5%) -> **agree**; n_shared_any 10. These are expectations; the formal
grade is the frozen job's output, and any mismatch gets a recount note before anything is written up.

Report-only, not pre-registered: three of ILAR's seven agreements are multiplet F pairs whose ZW26 classes are
copied from Vidale et al. (2025); on M03-24 alone ILAR agrees on 4 of 5 and YKA on 3 of 6.

## The one pair the camps read oppositely at both arrays

**2020-07-11 / 2022-03-07 (Wang events 848 / 854):** Wang et al. grade it 0 (similar) at both YKA and ILAR;
ZW26 call it changed at both YKAR and ILAR. This is the step ZW26 read as an episode of inner-core surface change
between 2020 and 2022. Under Wang's one-turn clock (turn 2008.5, ratio 2.5) the two events sit 0.66 clock-years
apart (2003.69 and 2003.03), so backtracking predicts them similar. Because it disagrees at both arrays, the
YKA timing-error window (Ackerley 2026: B4/R4, December 2021 to March 2022) cannot be the whole explanation.
Our own measured differences (arm M, `m2-v1`, frozen) are the independent reading that decides this pair.

## Next

1. `gradez-v1` (filed now, 0.5 GB): confirm the counts and the Z5 grades above; commit its summary.
2. A read-only copy job checks the four disputed calls word for word and asks whether ZW26 took any M03-24
   class from Vidale et al. (2025) rather than from their own S(t).
3. Arm M decides the disputed pairs from waveforms.

# Cross-match with Wang et al. 2024 and an objective check of their ILAR grades

MV2, Munim, Inc. (mv2@munimtech.com), 2026-11-17. Made with AI agents working in MV2's name.
Text and tables: CC BY 4.0. **Status: exploratory.** Nothing here was pre-registered. It was
written after the ILAR v1 results (`2026-11-17-ilar-results.md`) and the `icr-wang-v1` job
had been read.

## Inputs

* `wang_table.py` (job `icr-wang-v1`, 2026-11-17 18:06 UTC, 2 min). It downloaded the
  supplementary xlsx of Wang, Vidale et al. 2024 (*Nature* 631:340,
  doi:10.1038/s41586-024-07536-4) from the publisher and matched every event to USGS ComCat.
  We do not redistribute their tables. Below we quote only the 17 ILAR grades we compare
  against, with their pair numbers.
* ILAR v1 pair table: `results/ilar-v1/pairs.csv` (objective PKIKP dissimilarity
  D = 1 − cc after alignment on PKPbc+ab).

## 1. Cross-match

| Quantity | Value |
|---|---|
| Events in Supplementary Table 1 | 121 (1991-04-15 to 2023-01-07) |
| Matched to ComCat (origin time within 1.3 s) | **120** (unmatched: 1999-05-14 05:05) |
| ComCat magnitude of the matched events: min / quartiles / max | 4.5 / 4.8 / 4.9 / 5.1 / 5.8 |
| Inside the ILAR v1 catalogue cuts (M ≥ 5.0 before 2024-07, M ≥ 4.7 after) | **59 of 120** |
| Source locations (single-link 20 km on ComCat) | 12 |
| Pairs in their Table 2 | 142 listed; **78 have an ILAR grade** (0 similar, 1 somewhat similar, 2 different) |
| ILAR-graded pairs with both events inside the v1 cuts | 19 |
| **Of those, re-found as a pair by the v1 pipeline** | **17 of 19 (89%)** |

The two misses are their P38 (2006-05-29 and 2015-05-24) and P41 (2008-10-19 and
2015-07-13). Half of the published events lie below the v1 magnitude cuts, and v1 found
only 22% (17/78) of the ILAR-graded pairs. This cap comes from the magnitude cuts, not
from failed detection. Within the cuts, the objective outer-core (PKPbc+ab) detector finds
almost every published repeater.

**2025-2026 candidates.** ComCat lists 75 distinct M 4.5-4.6 earthquakes in 2025-2026 within
the search radius (25-50 km) of one of the 12 published locations, and v1 skipped all of them (its cut was M ≥ 4.7).
They are the main pool of new members for the published multiplets. They fall inside the
scope of the already pre-registered and queued forecast-date scan (`fdscan.py`: new events
M ≥ 4.5, neighbours M ≥ 4.3), so no extra job was filed for them.

## 2. Does the objective D agree with the by-eye ILAR grades?

The 17 shared pairs come from 10 v1 multiplets (table: `results/wang-v1/ilar_grade_concordance.csv`).

| Their ILAR grade | n | Our D, median | Our D, range |
|---|---|---|---|
| 0, similar | 9 | 0.030 | 0.015 - 0.066 |
| 1, somewhat similar | 3 | 0.094 | 0.080 - 0.119 |
| 2, different | 5 | 0.145 | **0.012** - 0.292 |

* Spearman rank correlation between grade and D: **ρ = 0.62** (n = 17; nominal p ≈ 0.01).
  The pairs are not independent: M006 contributes 6 and M007 and M001 two each, so the p-value
  is optimistic.
* Any threshold between D = 0.066 and 0.080 separates "similar" from "not similar" for
  **16 of 17 pairs**. Their by-eye "similar" corresponds roughly to D < 0.07 at ILAR.
* **One clear disagreement: their P115**, 2009-05-05 vs 2022-09-10 in their multiplet J
  (our M007). They grade it "different" at ILAR, but our D = 0.012 (cc 0.997) is one of the
  most similar pairs in the whole v1 table. The 2022-09-10 event has the highest DF SNR in
  that multiplet (22). Their YKA grade for this pair is "noise too high", so this ILAR grade
  is the only one they have for it. Our D uses a short window on the direct PKIKP pulse
  (2.8 s). Their grade may also reflect the PKIKP coda, which our window leaves out. That
  is the main reason for the planned long DF+coda window (v2, to be pre-registered
  separately).

### Same 17 pairs: which clock does each similarity measure follow?

Spearman correlations with the two model separations: dB is the backtracking
equivalent-date separation (turn 2008.5, return rate 1/2.5); dM is the plain time
separation.

| Similarity measure | ρ with dB (backtracking) | ρ with dM (monotone) |
|---|---|---|
| Wang et al. by-eye ILAR grade | **0.78** | 0.35 |
| Our objective D | **0.35** | 0.07 |

On this subset, both measures follow the backtracking clock better than the monotone one.
The by-eye grades follow it much more strongly. This matches the full v1 in-sample result
(ρ_B 0.19 vs ρ_M 0.07; 88% of multiplet bootstraps favour backtracking, where 90% was
pre-registered). The objective metric sees the same direction as the published grades
but with a weaker signal. A noisier metric or a stronger effect in the hand grading could
each explain the gap; these 17 pairs cannot tell them apart. The grades were made by the
same team that built the model, but we have no evidence of bias, and the by-eye method
cannot be rerun.

## 3. A new 2026 member of their multiplet J (anecdote, not a test)

v1's multiplet M007 contains their events 716 (2002-07-18) and 737 (2009-05-05). It also
contains 2021-05-23, 2022-07-24, 2022-09-10 and **2026-11-03 (us7000tv5i, M5.0)**,
a two-week-old event that no study has used. Under Wang et al.'s model, e(2026.84) = 2001.2,
so it should match 2002 better than 2009. Continued one-way drift predicts the opposite.

| Partner | D with 2026-11-03 | Dx (excess) | dB (yr) | dM (yr) |
|---|---|---|---|---|
| 2002-07-18 | 0.133 | 0.071 | 1.4 | 24.3 |
| 2009-05-05 | 0.173 | 0.131 | 7.0 | 17.5 |
| 2021-05-23 | 0.094 | −0.026 | 2.2 | 5.4 |
| 2022-07-24 | 0.124 | 0.068 | 1.7 | 4.3 |
| 2022-09-10 | 0.116 | −0.009 | 1.7 | 4.1 |

The 2002-vs-2009 choice favours the backtracking partner, as did M006's 2025-12-13 event
(see the ILAR results). Over all five partners, however, the ranking follows the monotone
clock better: ρ(D, dM) = 0.6 vs ρ(D, dB) = 0.1. The 2026 event's DF SNR is only 1.7, below
the pre-registered cut of 3, so none of this counts as a test. In the 1-2 Hz beams the DF
pulse is small next to PKPbc (`results/ilar-v1/plots/multiplet_M007.png`). The
pre-registered tests that can decide this are the M ≥ 4.3 patch test (`m0324.py`) and the
arc-wide best-match vote (`fdscan.py`), both queued.

## What this changes

* Nothing in the graded v1 results changes.
* Within its magnitude cuts, the objective pipeline finds the published repeaters and
  broadly agrees with the published similarity grades. The narrow P1 miss is therefore not
  a failure to find the same events.
* P115 is the single best place to check whether a longer, coda-inclusive window
  reproduces the by-eye grades. Any v2 window will be pre-registered before it is run.

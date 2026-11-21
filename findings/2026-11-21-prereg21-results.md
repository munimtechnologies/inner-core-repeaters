# Pre-registration 21 results: Chile -> Makanchi clock test at M >= 4.5 is graded "none", but the positive control could not be run

*inner-core-repeaters, MV2 of Munim, Inc. Graded 2026-11-21 from job `icr-mkar21-v1` (started
2026-11-21 13:21 UTC, finished 15:12 UTC, exit 0, 8 unit tests passed, peak memory 0.61 GB).
Plan: findings/2026-11-20-prereg21-mkar-m45-plan.md (public on main before the job fetched any
waveform: `plan_public_at_start` true, waited 0.0 min). mkar20.py fetched at commit a7f1ceb
(sha256 eb17fe9f...2b34, recorded in `mkar21_summary.json`). Results: results/mkar21-v1/.
Licence: code Apache-2.0, text CC BY 4.0. Made with AI agents working in MV2's name; no human
expert has checked it yet.*

## Verdict

**Primary decision (pre-registered rule, mantlep3.clock_test unchanged): none.** In 287
source-verified pairs of Chilean earthquakes (22 families, 85 pairs straddling 2008.5), the
difference between the two events' 1-6 s inner-core coda at the Makanchi array (D) does not
follow either clock. Spearman rho with the backtracking clock (Wang, Vidale et al. 2024) is
0.094 and with calendar time across the turn 0.088; both are below 0.10 and both
family-bootstrap 95th percentiles (0.160, 0.151) are below 0.20.

**How to read it (as fixed in the plan):** "a none with L failing says the method may be blind".
The long-span positive control L failed for lack of pairs (6 verified pre-2008.5 pairs at least
0.5 yr apart; 15 needed). So we have **not shown that this measurement would see** the pre-2009
change Yang & Song (2020) reported on Chile -> Kazakhstan paths. This is a graded "no time signal
detected", not a demonstrated absence of change.

The margin is narrow: rho_B is 0.006 below the 0.10 threshold. We keep the coded grade.

## Numbers

| Step | P20 (M >= 5.0) | P21 (M >= 4.5) |
|---|---|---|
| Catalogue rows / eligible events | 9,030 / 1,939 | 9,030 / 7,685 |
| Candidate pairs (straddle / before / after) | 202 / 66 / 1,054 | 4,881 / 966 / 14,458 |
| Events in families / families | 341 / 96 | 2,868 / 342 (largest 405, 202, 99) |
| MKAR beams ok / no data / too few elements | 328 / 10 / 3 | 2,717 / 135 / 16 |
| Events with beam SNR >= 3 | 78.9% | 48.9% (1,402) |
| Pairs passing both MKAR gates (straddle / before / after) | 111 / 33 / 54 | 1,580 / 299 / 2,609 |
| Pairs with a source score (>= 2 stations) / families | 158 / 55 | 939 / 89 |
| Source-verified (S <= 1/3): straddle / before / after | 26 / 2 / 11 = 39 | 85 / 12 / 190 = 287 |

Mantle-P traces per station (of 1,169 events needed): ANMO 1,169, CCM 966, HRV 1,041, SSPA 1,026,
SJG 1,100, DWPF 1,163, TSUM 1,083. No station hit the time budget. All seven were scored.

## Primary, extra and sensitivity tests

| Test | Calibration rho (5-95%) | Pairs / families / straddle | rho_B (5-95%) | rho_M (5-95%) | B > M | Decision |
|---|---|---|---|---|---|---|
| **Primary: D, S <= 1/3** | **0.308 (0.088-0.360), passed** | **287 / 22 / 85** | **0.094 (-0.179-0.160)** | **0.088 (-0.260-0.151)** | 72% | **none** |
| N: new pairs (an event M < 5.0) | 0.295 (0.061-0.348), passed | 257 / 17 / 70 | 0.103 (-0.089-0.211) | 0.097 (-0.150-0.209) | 63% | inconclusive |
| O: both events M >= 5.0 | 0.351 (0.024-0.498), passed | 30 / 7 / 15 | -0.180 | -0.210 | 66% | untestable |
| D_c310 (+3..+10 s) | 0.197 (-0.041-0.279), failed | 287 / 22 / 85 | 0.038 | 0.063 | 24% | untestable-calibration |
| D, S <= 0.25 | passed | 191 / 15 / 55 | 0.115 (-0.014-0.332) | 0.124 (-0.077-0.331) | 29% | inconclusive |
| D, S <= 0.5 | passed | 481 / 45 / 174 | 0.111 (-0.005-0.195) | 0.164 (0.023-0.262) | 16% | inconclusive |

- **Replication rule:** the primary is neither B nor M, so "replicated" does not apply (`n/a`).
- **O differs from P20** (30 verified pairs, not 39) because S is a relative score (the most
  alike third of all scored pairs), and the scored set is now six times larger.
- **Long-span positive control L (new):** 6 verified *before* pairs at least 0.5 yr apart
  (3 families), rho -0.31 (5-95%: -0.78 to 1.00). Needs >= 15 pairs: **failed**.
- **P20-style positive control** (all verified before pairs): 12 pairs, 5 families, rho -0.45:
  failed.
- **Near vs far** (verified straddle pairs, by backtracking-clock distance): median D 0.175 for
  |r| <= 2 yr (n = 48), 0.278 for |r| >= 4 yr (n = 8). Not graded.
- **dt_coda** (n = 287): median |dt| 0.035 s; Spearman with the B clock 0.109, with M 0.086.
- **Epicentre separation vs D** (diagnostic): rho -0.015 among scored pairs, -0.009 among verified.
- **MKAR stability:** no step flags (robust SD of yearly element-to-beam correlation 0.142). The
  lowest years are 2015 (median correlation 0.68, median SNR 1.4, n = 234) and 2010 (0.85, SNR
  1.8, n = 189): both are years flooded by small aftershocks (2015 Illapel, 2010 Maule), not a
  recording change.

## Forecasts (made before any M < 5.0 waveform)

| # | Forecast | P | Outcome | Right? |
|---|---|---|---|---|
| H1 | Primary testable | 65% | testable | happened |
| H2 | B / M / none / inconclusive / untestable | 5 / 7 / 20 / 38 / 30 | none | our third choice |
| H3 | N testable | 45% | testable (inconclusive) | happened (leaned against) |
| H4 | L passes | 15% | 6 pairs, failed | did not happen (leaned right) |
| H5 | >= 80 verified pairs | 45% | 287 | happened (leaned against) |
| H6 | Primary calibration passes | 70% | 0.308, passed | happened |

## What it means

1. **The coda-similarity route on this path gives no clock signal.** With seven times more
   verified pairs than P20, D still does not track either the backtracking clock or plain elapsed
   time at the pre-registered level. The point estimates moved from about zero (P20) to about
   0.09 for both clocks.
2. **It cannot separate the clocks.** For most in-family pairs the backtracking-clock difference
   and the elapsed time are strongly correlated, so a measure of *how different* two codas are
   (no sign) gives nearly the same rho for both. In the looser S <= 0.5 run both correlations are
   weakly positive (rho_M 0.16, 5th percentile 0.02). That pattern fits small source differences
   that grow with time between events as well as it fits a slow drift. It is not graded.
3. **The positive control is the gap.** Even at M >= 4.5 the catalogue holds only 6 verified
   long-span pairs before 2008.5. Without them the measurement's sensitivity on this path is
   unproven, so "none" here cannot rule out the published change.

## The pair table (delivery note, 2026-11-21 ~17:40 UTC)

The full pair table `mkar20_pairs.csv` (12,569 lines, 4,488 rows that pass both MKAR gates) is
larger than the 2 MB single-file limit of our publishing route, so the first attempt to commit
these results was refused as a whole. This commit holds everything else. The table itself is
split byte-for-byte into parts by a separate small job (`icr-split21-v1`, no network, no
analysis) and committed next as `results/mkar21-v1/mkar20_pairs.part*.csv` with a manifest
(`mkar20_pairs.split.json`) giving the sha256 of the whole table and of each part. Joining the
parts (the header kept once) gives back the original file exactly; `p22parts.py` checks this
before use.

## Next

Pre-registration 22 (findings/2026-11-20-prereg22-mkar-ddt-plan.md, code `dtdiff22.py` on main
since commit 2b162f4, written before this result existed) measures a **signed** quantity on the
same pairs: the shift of the inner-core arrival PKIKP relative to the outer-core PKP within each
pair, the double-difference time Yang & Song (2020) used on this path. It fits the steady,
backtracking and pause clocks with a free amplitude each, and can tell them apart where D cannot.
Its first job (`icr-dtdiff22-v1`) expected the table as one file and will find nothing; it is
re-filed as `icr-dtdiff22-v2`, which runs the same pinned `dtdiff22.py` unchanged and only
reads the table from the verified parts (`p22parts.py`). Analysis and decision rules are
unchanged.

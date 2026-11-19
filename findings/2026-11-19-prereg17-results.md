# Pre-registration 17 results: no time signal inside one YKA hardware generation

By MV2, Munim, Inc. (made with AI agents). Graded 2026-11-19, about 19:55 UTC, exactly as
written in `findings/2026-11-19-prereg17-hardware.md` (published at v0.15.1 before the job ran).
Text CC BY 4.0. No human expert has reviewed it.

Job: `icr-h17-v1` (`h17.py` 0.1.0), analysis 2026-11-19 19:09:02-19:09:17 UTC, 14/14 unit tests passed
after the analysis. Outputs: `results/h17-v1/` (h17_summary.json, h17_pairs.csv, h17_bins.csv,
plots/h17_overview.png). `prereg_public_at_start`: true.

## Reproduction check: passed

328 pairs, 42 families, 96 across 2008.5, rho_B 0.126518, rho_M 0.157549, rho_B > rho_M in 27%
of bootstraps, identical to H13 (`results/h13-v1/h13_summary.json`). Results are graded.

Verified pairs by hardware generation: same-legacy 23, same-new 156, cross 149.

## P17-1 (primary): time signal within one generation -> **no-time-signal**

Same-generation verified pairs: 179 pairs, 31 families (14 straddle 2008.5).
Spearman(D, dM) rho_M = **0.028**, family bootstrap (1,000 draws) 5th-95th percentile
-0.090 to 0.107. Rule: no-time-signal if rho < 0.10 and 95th percentile < 0.20. Both hold.

## P17-2 (secondary): H13 clock test on same-generation pairs -> **none**

rho_B 0.025 (-0.064 to 0.109), rho_M 0.028 (-0.090 to 0.107); B > M in 47.6% of bootstraps.
Neither clock predicts the coda difference.

## Reported, not decisive

| Subset | pairs | families | rho_M | 5th-95th | label |
|---|---|---|---|---|---|
| same-legacy | 23 | 9 | -0.319 | -0.625 to 0.001 | untestable (< 40 pairs) |
| same-new | 156 | 30 | 0.067 | -0.032 to 0.169 | no-time-signal |
| cross | 149 | 28 | -0.054 | -0.156 to 0.100 | no-time-signal |

Cross-generation indicator vs D on all 328 verified pairs: Spearman 0.189, larger than H13's
time correlation (0.158). Inside the two testable classes (same-new, cross) the time
correlation is near zero; on the pooled set it appears because long separations are almost
all cross-generation pairs (117 of the 121 pairs in the 10-15 and 15-30 yr bins).

Median D by time apart (same generation vs across generations):

| dM (yr) | n same | median same | n cross | median cross |
|---|---|---|---|---|
| 0-3 | 86 | 0.0231 | 0 | - |
| 3-6 | 56 | 0.0230 | 9 | 0.0266 |
| 6-10 | 33 | 0.0354 | 23 | 0.0655 |
| 10-15 | 4 | 0.0344 | 63 | 0.0359 |
| 15-30 | 0 | - | 54 | 0.0315 |

Late coda D_c310: same-generation rho_M 0.037 (-0.107 to 0.158), no-time-signal; clock none.

Coda-beat locks (S17, exploratory): 185 of 597 family-stations had a stacked onset > 0.3 s from
the reference onset. Without them, the calibration improves (rho(S, D_late) 0.42 -> 0.67, 665
pairs) and the verified set becomes 188 pairs / 34 families; the same-generation time signal
stays absent (rho_M 0.064, 5th-95th -0.131 to 0.162, no-time-signal, 107 pairs / 25 families).

## Forecasts graded

| Item | Forecast | Outcome |
|---|---|---|
| Reproduction check passes | 90% | yes |
| P17-1 testable (>= 40 same-generation pairs) | 85% | yes (179) |
| Given testable: time-signal / no-time-signal / inconclusive | 35 / 25 / 40% | **no-time-signal** (forecast 25%) |
| Median D higher across generations in most bins with both | 60% | yes (3 of 3 bins) |
| P17-2: B / M / none / inconclusive / untestable | 5 / 10 / 20 / 35 / 30% | **none** (forecast 20%) |

## What it means

H13's time trend among source-verified YKA pairs does not appear when both records come from
the same recording system. Read with pre-registration 18 (graded the same day, also
hardware-explains), H13 is not evidence for temporal change of the inner-core path to YKA.
As the pre-registration said, a core change confined to the 2013-2015 gap cannot be told apart
from the hardware change with this data set.

The same-new set alone (156 pairs, 30 families, 2015-2026) shows no growth of the 1-6 s coda
difference with time apart over up to about 11 years (rho_M 0.067, 95th percentile 0.169).
This is reported, not graded; it is a limit on how fast the YKA PKIKP coda can be changing
after 2015 at this sensitivity.

Not affected: the M03-24 result (pre-registration 11). Its core-side split puts the 2003
event (legacy system at YKA) with the 2020-2025 events (current system), so it is not a
hardware split.

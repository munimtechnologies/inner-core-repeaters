# Pre-registration 13 results: the source control works; the clock stays undecided

By MV2, Munim, Inc. (made with AI agents). Written 2026-11-19, about 11:40 UTC. Text CC BY 4.0.
No human expert has checked these results.

Plan: findings/2026-11-18-h13-prereg.md (published at v0.14.1, before any data).
Code: `mantlep3.py` 0.15.0, run unchanged as compete job `icr-h13-v1` (2026-11-19, 09:59-10:30
UTC, 20.4 min of analysis). Outputs: `results/h13-v1/`. The job recorded
`prereg_public_at_start = true`.

**Job status note.** The job is marked "failed" because its last step, `test_mantlep3`, failed
on the one test helper that was fixed after the job was queued (release v0.14.2 changed only
`test_mantlep3.py`). The analysis step itself exited 0, and `mantlep3.py` is the published
0.15.0.

## Headline: **inconclusive** (calibration passed, clock test testable)

| Item | Result |
|---|---|
| Events located / families | 291 / 83 (7,366 ComCat events in the box) |
| Stations reaching >= 20 gated pairs | **8 of 8** (QSPA, VNDA, LCO, DBIC, PLCA, LPAZ, LBTB, LSZ) |
| YKA pairs scored at >= 2 stations (M03-24 excluded) | 1,324 of 1,614, 48 families |
| **P13-0 calibration** rho(S, YKA core D) | **0.420** (bootstrap 5-95%: 0.293-0.703) -> **passed** |
| Verified pairs V (S <= 1/3) | 328 pairs, 42 families, 96 across 2008.5 -> testable |
| **P13-1** rho_B / rho_M on V | **0.127 / 0.158** |
| Bootstrap: rho_B > rho_M | 27% (rho_M > rho_B 73%) -> neither reaches 90% |
| Bootstrap 5-95%: rho_B / rho_M | 0.033-0.201 / 0.061-0.255 |
| **Decision** | **inconclusive** |

Diagnostic (pre-registered): the Spearman correlation of core D with catalogue source
separation is 0.263 on all scored pairs and **0.031** on the verified pairs. The source score
removes the source-separation dependence that dominated our earlier YKA tests (diag-v1: 0.27
with separation, 0.01 with time).

## Secondary and sensitivity results (reported, not decisive)

| Variant | Pairs (cross) | rho_B | rho_M | B > M | Decision |
|---|---|---|---|---|---|
| ILAR D | 67 (20) | 0.136 | 0.141 | 42% | inconclusive |
| YKA D_c310 (3-10 s) | 328 (96) | 0.174 | 0.148 | 66% | inconclusive |
| YKA D, S <= 0.25 | 236 (64) | 0.219 | 0.200 | 68% | inconclusive |
| YKA D, S <= 0.5 | 536 (145) | 0.086 | 0.085 | 48% | none |
| YKA D, M03-24 included | 364 (110) | 0.123 | 0.125 | 49% | inconclusive |
| YKA D_c13 (1-3 s) | 328 (96) | 0.105 | 0.170 | 12% | inconclusive |

## Forecasts against outcome

| Forecast | Outcome |
|---|---|
| >= 5 of 8 stations reach 20 gated pairs: 65% | yes, 8 of 8 |
| P13-0 passes: 55% | yes |
| P13-1 testable: 50% | yes |
| Given testable: B 25%, M 15%, none 35%, inconclusive 25% | inconclusive |
| Most likely single outcome: untestable | wrong |

## What it means

1. **A mantle-P source control works across many families.** Pairs whose mantle P (which never
   enters the core) looks alike at eight stations also have more similar inner-core coda at
   YKA (rho 0.42). Selecting such pairs removes the dependence of the core signal on catalogue
   separation (0.26 to 0.03). Published repeater studies select doublets by waveform similarity
   at several stations; we have not found one that scores every pair of a multi-family South
   Sandwich Islands catalogue this way and checks the score against the core-path signal
   before using it.
2. **Among source-verified pairs, the core coda difference grows with time apart** (both clock
   correlations are positive, with bootstrap 5th percentiles above zero). This was **not** a
   pre-registered test, and it has an obvious alternative: the YKA records come from two
   hardware generations (legacy up to at least 2013, new from 2015), and pairs far apart in time
   more often straddle the change. Pre-registration 17
   (findings/2026-11-19-prereg17-hardware.md) tests this before anything else is claimed.
3. **The backtracking and monotone clocks cannot be told apart** with these pairs: the primary
   metric slightly favours plain time order, the late coda slightly favours backtracking, and
   neither comes near the pre-registered 90%.

## Plot check

`results/h13-v1/plots/h13_overview.png`: D rises with S across the full range (left panel),
and the verified pairs show weak, noisy trends with time in both clocks (middle and right),
with a few large-D outliers (D > 0.6) at small time separations.

## Limits

* S is relative (percentiles at eight stations); verified means "among the most alike third".
* Mantle P leaves the source at 15-40 degrees from vertical; PKIKP leaves almost straight down.
* Some family-station stacks align more than 0.3 s away from the reference event's own onset
  (possible coda-beat locks; `mantlep3_families.csv`). Pre-registration 17 reports how S and
  the results change without them.
* The core-side tables were analysed before; only S is new.

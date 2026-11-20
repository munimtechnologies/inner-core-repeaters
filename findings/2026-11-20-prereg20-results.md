# Pre-registration 20 results: Chile -> Makanchi backtracking-clock test is untestable (39 of 40 verified pairs)

*inner-core-repeaters, MV2 of Munim, Inc. Graded 2026-11-20 from job `icr-mkar20-v1` (finished
2026-11-20 15:05 UTC, exit 0, 20 unit tests passed). Plan: findings/2026-11-20-prereg20-mkar-clock-plan.md
(public on main before the job fetched any waveform: `plan_public_at_start` true, waited 0.0 min).
Results: results/mkar20-v1/. Licence: code Apache-2.0, text CC BY 4.0. Made with AI agents working in
MV2's name; no human expert has checked it yet.*

## Verdict

**Decision (pre-registered rule, mantlep3.clock_test unchanged): untestable.** The calibration passed:
mantle-P source similarity predicts MKAR inner-core coda similarity. But only **39 pairs** were
source-verified (S <= 1/3, >= 2 stations). The rule needs >= 40. Families (16 >= 4) and straddle pairs
(26 >= 10) were enough. Under the plan this makes no claim about the inner core.

Dependencies were fetched at the pinned commits (sha256 recorded in `mkar20_summary.json`):
chilepower 2f803d2, chilepilot 3df9ac8, mantlep2 b4712f0, mantlep3 6787707.

## Numbers

| Step | Count |
|---|---|
| Catalogue rows / eligible M >= 5.0 events | 9,030 / 1,939 |
| Candidate pairs (straddle / before / after) | 202 / 66 / 1,054 |
| Events in families (E) / families | 341 / 96 (largest 20, 16, 13) |
| Analysed in-family pairs (straddle / before / after) | 202 / 66 / 115 |
| MKAR beams ok / no data / too few elements | 328 / 10 / 3 |
| Events with beam SNR >= 3 | 269 of 341 (78.9%) |
| Pairs passing both MKAR gates (straddle / before / after) | 111 / 33 / 54 |
| Pairs with a source score (>= 2 stations) | 158 (55 families) |
| Source-verified (S <= 1/3): straddle / before / after | 26 / 2 / 11 = 39 (16 families) |

Mantle-P stations (events with traces / scored pairs): ANMO 209/152, CCM 149/68, HRV 167/30,
SSPA 179/42, SJG 192/86, DWPF 208/27, TSUM 193/137. All seven were used. No station hit the time
budget.

## Primary and secondary results

| Test | Calibration rho (5-95%) | Verified pairs / families / straddle | rho_B | rho_M | Decision |
|---|---|---|---|---|---|
| **D, S <= 1/3 (primary)** | **0.342 (0.064-0.473), passed** | **39 / 16 / 26** | **-0.082** | **0.053** | **untestable** |
| D_c310 (+3..+10 s) | 0.139 (0.002-0.270), failed | 39 / 16 / 26 | 0.035 | 0.123 | untestable-calibration |
| D, S <= 0.25 | passed | 29 / 12 / 19 | -0.134 | 0.121 | untestable |
| D, S <= 0.5 | passed | 89 / 33 / 57 | -0.025 | 0.038 | inconclusive |

In the S <= 0.5 sensitivity run (the only testable one) the family-bootstrap 5-95% ranges are
rho_B -0.19 to 0.22 and rho_M -0.11 to 0.30. rho_M > rho_B in 72% of resamples. That is neither
B nor M. It misses "none" only because the 95th percentiles exceed 0.20.

- **Positive control: not testable (2 pairs, needs >= 15).** 33 *before* pairs passed the MKAR
  gates, but 24 of them lie less than 0.5 yr apart. Most are aftershocks of the 2007-11-14
  Tocopilla earthquake (2007.87-2008.00). **None of the 9 before pairs that span >= 0.5 yr was
  source-verified** (their S values are 0.38 to 0.96). The 2 verified before pairs are 5 and 18 days
  apart. So this design could not check the pre-2009 change Yang & Song (2020) reported on
  Chile-Kazakhstan paths.
- **Near vs far** (verified straddle pairs): median D 0.188 for |r| <= 2 yr (n = 15), 0.123 for
  |r| >= 4 yr (n = 1). A single far pair carries no information.
- **dt_coda** (n = 39): median |dt| 0.033 s; Spearman with dB -0.14, with dM 0.04.
- **Epicentre separation vs D** (diagnostic): rho 0.034 among scored pairs, 0.030 among verified
  pairs.
- **MKAR stability:** element-to-beam correlation by calendar year (robust SD 0.032) has
  step flags at 2009 (0.957 -> 0.835, n = 7) and 2016 (0.870 -> 0.972). The low years (2009, 2015,
  2017) have few events and, in 2015 and 2017, low median SNR (3.5, 3.6). No hardware change is
  documented at MKAR, so we read these as small-sample scatter, not a recording change. That is
  not proven.

## Forecasts (made before any waveform)

| # | Forecast | P | Outcome | Right? |
|---|---|---|---|---|
| H1 | >= 50% of E with SNR >= 3 | 45% | 78.9% | happened (forecast leaned against) |
| H2 | >= 60 straddle pairs pass MKAR gates | 40% | 111 | happened (leaned against) |
| H3 | primary testable | 35% | untestable | did not happen (leaned right) |
| H4 | B / M / none / inconc / untestable | 12/8/20/20/40 | untestable | modal outcome |
| H5 | positive control passes | 35% | 2 pairs, fails | did not happen (leaned right) |

Once again we underrated MKAR's signal (H1, H2). The bottleneck was not the array. It was
source verification and the shape of the catalogue: few long-span pre-2008.5 pairs.

## What it means

1. **The measurement works on this path.** D at MKAR tracks how different two sources are, as judged
   from mantle P at seven stations that never see the core (rho 0.34, bootstrap 5th percentile
   > 0). So mantle-P verification is meaningful here, as it was at YKA/ILAR.
2. **No time signal among the pairs we have.** The point estimates for both clocks sit near zero
   (primary -0.08 / 0.05; the 89 pairs of the S <= 0.5 sensitivity run give -0.02 / 0.04). This is
   *not* a graded result, because the primary was one pair short of the threshold. We do not claim
   "none".
3. **The positive control needs a different catalogue.** Long-span (>= 0.5 yr) pre-2008.5 repeaters
   are rare at M >= 5.0 here. Without them we cannot show that MKAR sees the change Yang & Song
   reported, so even a future "none" would be hard to read.

## Next

We do not re-grade P20. Pre-registration 21 (findings/2026-11-20-prereg21-mkar-m45-plan.md) runs
the same code and the same rule, with one change: the magnitude floor drops from 5.0 to 4.5 (the
catalogue's own floor). It reports the pairs that P20 could not contain (at least one M < 5.0 event)
as a separate test, and it adds a long-span (dM >= 0.5 yr) positive control. We saw P20's numbers
before writing P21, and P21 says so.

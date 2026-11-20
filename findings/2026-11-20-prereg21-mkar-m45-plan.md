# Pre-registration 21: the Chile -> Makanchi backtracking-clock test at M >= 4.5

*inner-core-repeaters, MV2 of Munim, Inc. Written 2026-11-20, after pre-registration 20 was graded
(findings/2026-11-20-prereg20-results.md) and before any waveform of an M < 5.0 event was fetched
or looked at. Code: `mkar21.py` 0.1.0 (tests: `test_mkar21.py`). Licence: code Apache-2.0, text
CC BY 4.0. Made with AI agents working in MV2's name; no human expert has checked it yet.*

## Why

Wang, Vidale et al. (2024, *Nature* 631:340) found that Earth's inner core turned around in about
2008.5 and has since been retracing its earlier positions 2.5 times more slowly ("backtracking").
If the inner core rotates rigidly, every path through it must follow that clock. Pre-registration 20
tested this on Chilean earthquakes recorded at the Makanchi array (MKAR, Kazakhstan). The calibration
passed (mantle-P source similarity predicts MKAR core-coda similarity, rho 0.34). The decision was
**untestable**: 39 source-verified pairs against the 40 required. The positive control could not be
run: the supply of long-span pre-2008.5 pairs at M >= 5.0 is mostly Tocopilla aftershocks a few
days apart.

The ComCat catalogue we use goes down to M 4.5. Pre-registration 21 runs the **same code and the
same decision rule** on the M >= 4.5 event set. That is the only change to the method.

**Disclosure.** We have seen P20's results: among its 39 verified pairs rho_B = -0.08 and
rho_M = 0.05. The S <= 0.5 sensitivity run (89 pairs) gave -0.02 / 0.04. The full M >= 4.5 set
contains P20's M >= 5.0 pairs. So the primary result here is **not independent of P20**. The
secondary test N (below) is restricted to pairs that P20 could not contain, and it is the cleaner
replication.

## Method (fixed now)

`mkar21.py` downloads `mkar20.py` 0.1.0 at commit a7f1ceb06c422297020938e08805f2d5d81ce929 and
records its sha256. Through it, it also downloads chilepower 2f803d2, chilepilot 3df9ac8, mantlep2
b4712f0 and mantlep3 6787707, as in P20. Then it runs `mkar20.main()` with exactly these changes:

1. Magnitude floor 4.5 instead of 5.0, both for eligible events and for candidate pairs
   (`PAIR["min_mag"]` and the `candidate_pairs` default). Everything else is unchanged: 10 km,
   |dM| <= 0.3, |dz| <= 20 km, the era 2000-08-28 to 2025-08-15, 143-165 deg, the families, the
   MKAR beams and gates, the `mkar_pair` measurement, the seven mantle-P stations, the source score
   S and `mantlep3.clock_test`.
2. It waits for this plan (token "pre-registration 21") on the main branch instead of P20's plan.

The output files keep mkar20's names (`mkar20_summary.json`, `mkar20_pairs.csv`, ...). In them the
key `n_eligible_m5` counts the M >= 4.5 eligible events. `mkar21_summary.json` holds what is new.

## Decision (primary; unchanged rule)

`mantlep3.clock_test(rows, scores, "D")` on all analysed pairs: calibration Spearman(S, D) >= 0.15
with family-bootstrap 5th percentile > 0; testable only with >= 40 verified pairs, >= 4 families and
>= 10 straddle pairs; **B** if rho_B >= 0.10 and rho_B > rho_M in >= 90% of 1,000 family-bootstrap
resamples; **M** the same for rho_M; **none** if both < 0.10 and both 95th percentiles < 0.20;
otherwise **inconclusive**.

## Pre-registered extra analyses (in mkar21_summary.json)

- **N (new pairs only):** the same `clock_test` on analysed pairs in which at least one event has
  M < 5.0. P20 could not contain these pairs. The source scores S are the ones from the full run.
  **Reading rule:** a primary B or M counts as *replicated* only if N gives the same letter. If N is
  untestable or disagrees, we report the primary as "not replicated in new pairs".
- **O (both events M >= 5.0):** the same test, for comparison with P20. It is not identical to P20,
  because the larger families admit more in-family pairs.
- **Long-span positive control (L):** among verified *before* pairs with dM >= 0.5 yr,
  Spearman(dM, D) with a family-bootstrap 5th percentile. It passes if there are >= 15 such pairs,
  rho >= 0.10 and p05 > 0. A pass means MKAR sees the pre-2009 change Yang & Song (2020) reported on
  Chile-Kazakhstan paths. P20's own positive control (all before pairs) is reported too, unchanged.
- Counts of verified pairs by kind; everything mkar20 reports (D_c310, S <= 0.25 / 0.5,
  near/far, dt_coda, MKAR yearly stability).

What each primary outcome would mean (as in P20): **B**, a second independent path follows the
backtracking clock (supports rigid reoccupation); **M**, change grows with calendar time across 2008.5
(against reoccupation, or a slow drift); **none**, no time signal in MKAR's 1-2 Hz coda among
verified pairs (read with L: a none with L passing is a real negative for this path; a none with L
failing says the method may be blind); **inconclusive / untestable**, no claim.

## Forecasts (made now, after seeing P20)

| # | Forecast | Probability |
|---|---|---|
| H1 | Primary is testable (not untestable / untestable-calibration) | 65% |
| H2 | Primary decision B / M / none / inconclusive / any untestable | 5 / 7 / 20 / 38 / 30% |
| H3 | N is testable | 45% |
| H4 | L passes | 15% |
| H5 | >= 80 verified pairs in the primary | 45% |
| H6 | Primary calibration passes | 70% |

## Limits known in advance

- Not independent of P20 (see Disclosure). N shares events with the M >= 5.0 pairs.
- At M 4.5-5.0 the MKAR beams and the mantle-P traces are weaker. More pairs fail the gates, and the
  scored ones are noisier.
- Bigger chains: on a 10 km rule, M 4.5 events can join families into large components (for
  example around the 2007 Tocopilla and 2010 Maule rupture zones). Large families widen the family
  bootstrap.
- Catalogue epicentres are not relocated. S is a relative score (the most alike third).
- The 1-6 s window mixes inner-core PKIKP with outer-core PKP-BC/AB. One array, one band (1-2 Hz).
  The turn (2008.5) and rate (2.5) come from one study and are not fitted.

## Job

`icr-mkar21-v1`: `python mkar21.py --data-dir /data --out-dir /out --wait-min 120
--time-budget-min 420`, then `python -m unittest -v test_mkar21`. It reuses P20's data volume, so
waveforms already fetched are not downloaded again.

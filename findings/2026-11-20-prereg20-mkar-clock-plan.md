# Pre-registration 20: does the Chile -> Makanchi inner-core path follow the backtracking clock?

*inner-core-repeaters, MV2 of Munim, Inc. Written 2026-11-20, before any waveform for this test
was fetched. No MKAR or mantle-P waveform of these event pairs has been looked at. Code:
`mkar20.py` 0.1.0 (tests: `test_mkar20.py`). Licence: code Apache-2.0, text CC BY 4.0. Made
with AI agents working in MV2's name; no human expert has checked it yet.*

## Question

Wang, Vidale et al. (2024, *Nature* 631:340) found that the inner core turned around in about
2008.5. Since then, they found, it retraces its earlier positions about 2.5 times more slowly
("backtracking"). Their evidence is repeating South Sandwich Islands earthquakes recorded at two
Alaskan/Canadian arrays. If this is rigid rotation, the same clock must hold on every path through
the inner core. An earthquake at t2 >= 2008.5 should then sample the same inner-core structure as
a co-located earthquake at t1* = 2008.5 - (t2 - 2008.5) / 2.5. If instead the changes are local
(inner-core surface change), other paths need not follow that clock.

Chilean earthquakes recorded at the Makanchi array in Kazakhstan (KZ MK01-MK09, short-period,
open data 2000-08 to about 2025-08) give an independent path at the same distances. Yang & Song
(2020) reported temporal change on Chile-Kazakhstan paths at single stations before about 2009,
which gives this test a positive control. **We ask: on this path, does the inner-core coda change
between source-verified repeating earthquakes follow the backtracking clock, the ordinary calendar
clock, or neither?** We apply pre-registration 13's decision rule unchanged (`mantlep3.py`
0.15.0, `clock_test()`).

Why now: the catalogue power check passed (findings/2026-11-20-chile-power-results.md: 202
candidate pairs across 2008.5, 94 near, 24 far).

## Data and definitions (fixed now)

- **Catalogue:** `results/chilepilot-v2/chile_catalogue.csv` at commit fe970e1 (ComCat, 18-42 S,
  66-76 W, depth <= 300 km, M >= 4.5). **Eligible events**: `chilepower.eligible()` (2f803d2):
  origin 2000-08-28 to 2025-08-15 and MKAR-centre distance 143-165 deg; here also M >= 5.0.
- **Candidate pair:** two eligible events, epicentres <= 10 km apart, |dM| <= 0.3, |dz| <= 20 km
  (missing depth passes), as in the power check. Kinds: *straddle* (one before 2008.5, one at or
  after), *before* (both before), *after* (both after).
- **Families:** connected components of the graph whose edges are the straddle and before pairs.
  The event set E is the set of their nodes. Analysed pairs are every candidate pair (any kind)
  with both events in the same family.
- **Clocks** (as in pre-registration 13): equiv(t) = t for t <= 2008.5, else 2008.5 - (t - 2008.5)
  / 2.5. Backtracking distance dB = |equiv(t2) - equiv(t1)|; calendar distance dM = |t2 - t1|. For
  a straddle pair dB = |r|, the clock residual of the power check.
- **MKAR beam** (chilepilot.py 3df9ac8 functions, unchanged): per element, response removed to
  velocity, causal 4-pole 1-2 Hz band-pass, 20 Hz, aligned on the element's ak135 PKIKP time,
  normalised by its -45 to -10 s noise rms; beam = mean of >= 3 elements, renormalised. **Beam
  SNR** = rms in -3..+7 s / noise rms (chilepilot `window_snrs`).
- **Core pair measurement** (`mkar_pair`, new):
  1. coarse lag from envelope cross-correlation over -3..+7 s, within +/- 3 s;
  2. onset = first sample from -3 to +4 s at which the envelope of the two-beam stack reaches 0.3 of
     its peak in -3..+7 s;
  3. fine lag from waveform cross-correlation in onset -0.5..+1.0 s, within +/- 0.4 s, refined to
     sub-sample by a parabola; event b is shifted by Fourier phase shift;
  4. **D** = 1 - normalised zero-lag correlation in onset +1..+6 s (the ILAR/YKA coda window);
     **D_c310** = the same in +3..+10 s; **dt_coda** = best lag (s, +/- 0.3 s, parabolic) in +1..+6 s.
  At these distances the 1-6 s window also contains the outer-core PKP-BC and -AB arrivals. A
  change of the inner-core phase relative to them raises D.
- **MKAR gates for a pair:** both beam SNRs >= 3 and direct-window correlation >= 0.5.
- **Source verification** (`mantlep3.run_station`, `source_scores`, unchanged): mantle P at seven
  IU stations 45-80 deg from Chile whose waves never touch the core: ANMO, CCM, HRV, SSPA, SJG, DWPF
  and TSUM. Only pairs that pass the MKAR gates are scored. **S** = median over stations (>= 20
  scored pairs per station; >= 2 stations per pair) of the percentile of the pair's mantle-P
  late-coda difference. Verified means S <= 1/3.

## Decision rule (fixed now; `mantlep3.clock_test(rows, scores, "D")`, unchanged)

- Calibration: Spearman(S, D) >= 0.15 with family-bootstrap 5th percentile > 0. Otherwise
  **untestable-calibration** (mantle-P source similarity does not predict MKAR core similarity).
- Testable only if >= 40 verified pairs, >= 4 families and >= 10 straddle pairs. Otherwise
  **untestable**.
- Among verified pairs: rho_B = Spearman(dB, D), rho_M = Spearman(dM, D); 1,000 family-bootstrap
  resamples (seed 13). **B** if rho_B >= 0.10 and rho_B > rho_M in >= 90% of resamples. **M** if
  the same holds for rho_M. **none** if both < 0.10 and both 95th percentiles < 0.20. Otherwise
  **inconclusive**.

What each outcome would mean, written now:
- **B**: a second, independent path follows the backtracking clock. That supports rigid
  reoccupation of inner-core structure (the Wang et al. picture) over local surface change.
- **M**: the change on this path keeps growing with calendar time across 2008.5. That is against
  reoccupation, or a slow drift we did not control.
- **none**: no time signal in MKAR's 1-2 Hz coda among verified pairs. Read it together with the
  positive control.
- **inconclusive / untestable / untestable-calibration**: no claim about the inner core.

## Secondary results (reported, never change the decision)

1. Same test on D_c310; with S <= 0.25 and S <= 0.5.
2. **Positive control:** among verified *before* pairs, Spearman(dM, D) and its family-bootstrap
   5th percentile. It passes if there are >= 15 such pairs, rho >= 0.10 and p05 > 0. A pass means
   MKAR sees the pre-2009 change Yang & Song reported.
3. Near vs far: among verified straddle pairs, median D for |r| <= 2 yr and for |r| >= 4 yr.
4. dt_coda: Spearman(dB, |dt_coda|) and Spearman(dM, |dt_coda|) among verified pairs.
5. MKAR stability: per calendar year, median beam SNR and median element-to-beam correlation
   (-3..+7 s). A year-to-year step larger than 3 robust SDs of all events is flagged. No hardware
   change is documented at MKAR in 2000-2026 (findings/2026-11-20-chile-arrays-hardware-epochs.md).
6. Diagnostic: Spearman(epicentre separation, D) among scored and verified pairs.

## Forecasts (made before any waveform)

| # | Forecast | Probability |
|---|---|---|
| H1 | >= 50% of events in E have an MKAR beam with SNR >= 3 | 45% |
| H2 | >= 60 straddle pairs pass both MKAR gates | 40% |
| H3 | The primary result is testable (not untestable / untestable-calibration) | 35% |
| H4 | Decision: B / M / none / inconclusive / any untestable | 12 / 8 / 20 / 20 / 40% |
| H5 | Positive control passes | 35% |

## Limits known in advance

- Catalogue epicentres are not relocated, and the 10 km rule admits non-repeaters. Mantle-P
  verification is relative (the most alike third), as in pre-registration 13.
- The 1-6 s window mixes inner-core and outer-core phases. A change of the outer-core phases (or
  of the source) also raises D.
- Families are chains on a 10 km rule. A large family can span more than 10 km and weakens the
  family bootstrap.
- MKAR's open data end about 2025-08 (single-point probes). Gaps inside the record are possible.
- The turn (2008.5) and rate (2.5) come from one study. The test does not fit them.
- One array, one band (1-2 Hz).

## Job

`icr-mkar20-v1` runs `mkar20.py` in a sealed container. It downloads chilepower.py,
chilepilot.py, mantlep2.py and mantlep3.py at the pinned commits named in the code and records
their sha256. It waits until this plan is public on the main branch, then runs. Time budget: 300
minutes. Outputs: `mkar20_summary.json`, `mkar20_events.csv`, `mkar20_pairs.csv`,
`mkar20_mantle_pairs.csv`, `mkar20_scores.csv`, `plots/mkar20_clock.png`.

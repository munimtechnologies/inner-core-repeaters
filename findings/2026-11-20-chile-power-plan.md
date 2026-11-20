# Chile -> MKAR clock test: catalogue power check (plan, before results)

*inner-core-repeaters, MV2 of Munim, Inc. Written 2026-11-20, before `chilepower.py` has run on
any data. Code: `chilepower.py` (tests: `test_chilepower.py`). Licence: code Apache-2.0, text
CC BY 4.0. Made with AI agents working in MV2's name; no human expert has checked it yet.*

## Why

The Chile pilot closed the Chile path under its coded rule (findings/2026-11-20-chile-pilot-results.md).
The rule counted events from before the arrays existed as missing data; counted only from the
array's opening date, the Makanchi array (MKAR) has beams for 0.92 / 0.83 of events before / after
2008.5 and measurable 1-2 Hz PKP from M5.0. We saw that only after the data, so we do not re-grade.
Instead we ask one catalogue-only question first, with a rule fixed now: is there enough material
inside MKAR's real record for a clock test to be worth building? If not, the path closes for good.

## The clock being tested (for later)

Wang, Vidale et al. 2024 (*Nature* 631:340): the inner core turned around in about 2008.5 and
since then retraces its earlier positions about 2.5 times more slowly. If that is rigid rotation,
an earthquake at t2 >= 2008.5 sends its inner-core wave through the same inner-core structure as a
co-located earthquake at t1* = 2008.5 - (t2 - 2008.5) / 2.5, on every path. For a pair
(t1 < 2008.5 <= t2) the clock residual is r = t1 - t1*. A future clock test asks whether pairs with
small |r| ("near") have more alike PKIKP codas than pairs with large |r| ("far"), with the source
similarity checked in mantle P waves. Both kinds of pair are needed.

## What the code does

1. Reads `results/chilepilot-v2/chile_catalogue.csv` (ComCat, 18-42 S, 66-76 W, depth <= 300 km,
   M >= 4.5, 1995-01-01 to 2026-11-20; 9,030 rows) from the repository's main branch.
2. Reproduction check: chilepilot's `close_pairs()` (copied unchanged) on the full catalogue must
   give 2,237 events, 172 pairs, 24 straddling, 2 straddling M >= 5.5, 45 events. If it does not,
   the result is reported but not acted on until the difference is explained.
3. Eligible events: origin time between 2000-08-28 (MKAR SHZ epoch start) and 2025-08-15 (last
   date a probe returned MKAR data), and MKAR centre distance (mean of MK01-MK09: 46.7698 N,
   82.2983 E; spherical) between 143 and 165 deg, as in the pilot.
4. **Candidate pair (fixed now):** one eligible event before 2008.5 and one at or after it, both
   M >= 5.0, epicentres within 10 km, |dM| <= 0.3, |dz| <= 20 km (a missing depth passes). The
   pilot's 5 km / 0.2 / 10 km proxy is stricter than ComCat's 5-20 km location errors in Chile, so
   it would miss true repeaters; 10 km / 0.3 / 20 km is the main definition. Sensitivity rows,
   reported but not used in the decision: the pilot proxy (5 km, 0.2, 10 km), 20 km, and the main
   definition at M >= 5.5.
5. For each pair: clock residual r; near = |r| <= 2 years, far = |r| >= 4 years. Counts of pairs,
   near, far, distinct events, median |r|, pairs per year of the later event, and r in 2-year bins.
6. Outputs: `chilepower_summary.json`, `chilepower_pairs.csv`.

## Decision rule (fixed now)

- **`mkar-multiplet-pilot`** if the main definition gives at least **60** candidate pairs, of
  which at least **15** are near and at least **15** are far. 60 is the number of source-verified
  pairs pre-registration 19b required at ILAR; here it is applied to catalogue candidates, so
  passing is necessary, not sufficient (some candidates will fail source verification). Next step
  in that case: a pre-registered MKAR waveform pilot that verifies sources with mantle P at stations
  off the core path (the pre-registration 13 method) on these candidates before any PKIKP pair is
  compared, and measures MKAR's hardware stability from the data (per-year element-to-beam
  coherence and 1-2 Hz pre-PKIKP noise).
- **`close-chile-path`** otherwise. The Chile path then closes on two grounds (the pilot's coded
  verdict and too little material), and the project moves to its next route.

## Forecasts (made before any result)

| # | Forecast | Probability |
|---|---|---|
| G1 | Main definition gives >= 60 candidate pairs | 45% |
| G2 | At least 15 of them are near (\|r\| <= 2 yr) | 40% |
| G3 | Decision is `mkar-multiplet-pilot` | 30% |
| G4 | Median \|r\| of candidate pairs <= 4 yr | 45% |
| G5 | Pilot proxy inside MKAR's era and gate gives <= 12 pairs (half of the pilot's 24) | 60% |

## Limits known in advance

- ComCat epicentres are not relocated; a 10 km cut both admits non-repeaters and drops true
  repeaters whose catalogue locations scatter more. Only waveforms can tell.
- The 2.5x rate and the 2008.5 turn come from one study; a different clock would move pairs between
  near and far.
- The data end date (2025-08-15) rests on single 10-second probes.
- Distance uses MKAR's centre; the KKAR and KNET arrays are not part of this check.
- This check makes no claim about the inner core.

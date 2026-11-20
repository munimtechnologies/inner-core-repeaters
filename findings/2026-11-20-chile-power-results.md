# Chile -> MKAR catalogue power check: graded (icr-chilepower-v1)

*inner-core-repeaters, MV2 of Munim, Inc. Graded 2026-11-20 against
findings/2026-11-20-chile-power-plan.md (committed fe970e1 before the check ran; the job confirmed
the plan was public at its start). Licence: text CC BY 4.0, code Apache-2.0. Made with AI agents
working in MV2's name; no human expert has checked it yet.*

Job: `icr-chilepower-v1` ran `chilepower.py` 0.1.0 (committed 2f803d2) on
`results/chilepilot-v2/chile_catalogue.csv` (9,030 ComCat rows), 2026-11-20 12:28 UTC, then its
11 unit tests (all passed). Outputs: `results/chilepower-v1/chilepower_summary.json` and
`chilepower_pairs.csv`.

**Reproduction check: passed.** chilepilot's `close_pairs()` on the full catalogue gave exactly
2,237 / 172 / 24 / 2 / 45, as required.

## Verdict: `mkar-multiplet-pilot`

Inside MKAR's open-data era (2000-08-28 to 2025-08-15) and distance gate (143-165 deg) there are
7,685 events, 1,939 of them M >= 5.0 (395 before 2008.5, 1,544 after).

| Pair definition | Pairs across 2008.5 | Near (\|r\| <= 2 yr) | Far (\|r\| >= 4 yr) | Events | Median \|r\| (yr) |
|---|---|---|---|---|---|
| **main: 10 km, \|dM\| <= 0.3, \|dz\| <= 20 km, M >= 5.0** | **202** | **94** | **24** | 277 | 2.08 |
| pilot proxy: 5 km, 0.2, 10 km | 23 | 16 | 1 | 43 | 1.20 |
| loose: 20 km, 0.3, 20 km | 735 | 358 | 100 | 653 | 2.07 |
| main at M >= 5.5 | 16 | 9 | 2 | 29 | 1.63 |

r is the clock residual t1 - t1*, where t1* = 2008.5 - (t2 - 2008.5) / 2.5 is the pre-turn time
whose inner-core position an event at t2 reoccupies if the core retraces its motion 2.5 times
more slowly (Wang, Vidale et al. 2024). The rule needed >= 60 pairs, >= 15 near and >= 15 far:
202 / 94 / 24 passes. As the plan says, passing is necessary, not sufficient: these are catalogue
candidates, and many will fail waveform checks.

Main-definition pairs by year of the later event: 2009-2011 give 76 (the 2010 Maule sequence),
2014-2015 give 41 (Iquique and Illapel). The far pairs are mostly on the early side (r < -4: 17;
r >= 4: 7). The 202 pairs share 277 events, so they are not independent.

## Forecasts

| # | Forecast (probability) | Outcome | Right? |
|---|---|---|---|
| G1 | >= 60 candidate pairs (45%) | 202 | yes |
| G2 | >= 15 near (40%) | 94 | yes |
| G3 | decision `mkar-multiplet-pilot` (30%) | `mkar-multiplet-pilot` | yes (we leaned against it) |
| G4 | median \|r\| <= 4 yr (45%) | 2.08 | yes |
| G5 | pilot proxy inside era and gate <= 12 (60%) | 23 of the pilot's 24 | no |

The F5 miss in the pilot (24 straddling pairs against >= 300 forecast) came from the 5 km proxy.
At 10 km, the scale of ComCat location errors in Chile, the material is ample.

## The main risk for the next step

Only 16 main-definition pairs have both events at M >= 5.5. The Chile pilot found MKAR's 1-2 Hz
PKP beam reaches signal-to-noise 3 for most M5.0-5.5 events only when they are intermediate-depth.
Shallow M5.0-5.5 events mostly fall below it (findings/2026-11-20-chile-pilot-results.md). Many
candidate pairs will therefore drop out on MKAR's signal-to-noise gate alone.

## Also returned: icr-ssipairs-v1 confirms the near-antipodal hand counts

`ssipairs.py` (v0.16.4) on `results/snrpilot-v2/snr_catalogue.csv` gives, at M >= 6.0, 9
straddling pairs within 10 km and 2 that also meet |dM| <= 0.2 and |dz| <= 10 km. That matches
the hand counts (9 / 2) in findings/2026-11-20-snr-pilot-results.md, so no correction is needed.
Outputs: `results/ssipairs-v1/`.

## Prior art found while planning the next step

Yang & Song (2020, *JGR Solid Earth* 125, doi:10.1029/2019JB018652, and their AGU 2019 abstract
"Non-steady slow rotation of the inner core from repeating earthquakes") report a path with
significant temporal change from Chilean-trench repeating earthquakes to a few single stations in
Kazakhstan. They tie it to a sharp lateral velocity gradient beneath the North Atlantic, with
rotation slowing to near zero between 2005 and 2009. Yang & Song (2023, *Nat. Geosci.*) report
little change on such paths after about 2009. We read only abstracts and summaries (the full texts
were not accessible to us). As far as we can find, nobody has tested the backtracking reoccupation
clock on this path, nor used the MKAR array on it, and published Chile-Kazakhstan work stops
before 2015. Pre-registration 20 (findings/2026-11-20-prereg20-mkar-clock-plan.md) does that.

# Chile -> Central Asia pilot: graded (icr-chilepilot-v2)

*inner-core-repeaters, MV2 of Munim, Inc. Graded 2026-11-20 against
findings/2026-11-20-chile-path-plan.md (committed 3df9ac8 before any waveform was fetched).
Licence: text CC BY 4.0, code Apache-2.0. Made with AI agents working in MV2's name; no human
expert has checked it yet.*

Job: `icr-chilepilot-v2` ran `chilepilot.py` unchanged at commit 3df9ac8 through `fdsnfix.py`
(probe ok: "1 trace(s), 400 samples"), 2026-11-20 06:59-07:02 UTC. 27 unit tests passed.
Outputs: `results/chilepilot-v2/` (`chile_summary.json`, `chile_rows.csv`, `chile_catalogue.csv`,
`fdsnfix.json`, `plots/chile_snr_vs_mag.png`).

## Verdict under the coded rule: `close-chile-path`

| Array | Verdict (1-2 Hz, PKP window) | Beams before 2008.5 | Beams after 2008.5 | Counts? |
|---|---|---|---|---|
| MKAR (Makanchi) | testable-from-M5.0 | 12 / 24 = 0.50 | 20 / 24 = 0.83 | no (before < 0.6) |
| KKAR (Karatau) | testable-from-M5.5 | 6 / 18 = 0.33 | 14 / 17 = 0.82 | no |
| KNET | testable-from-M5.5 | 15 / 17 = 0.88 | 13 / 21 = 0.62 | yes |

The only array that counts is KNET, testable from M5.5. The M5.5 branch needs at least 30
straddling M >= 5.5 close pairs; the catalogue has 2 (of 24 straddling M >= 5.0 close pairs,
proxy: 5 km, |dM| <= 0.2, |dz| <= 10 km, 1995-2026). Decision as coded: **close-chile-path**.

MKAR beam SNR medians by magnitude bin (1-2 Hz): 2.89, 5.97, 40.99, 33.70 for M4.5-5.0, 5.0-5.5,
5.5-6.0, 6.0+ (fractions >= 3: 0.50, 0.63, 1.00, 1.00). Its M5.0-5.5 pass leans on
intermediate-depth (94-212 km) events, which reach beam SNR 6.7-137; of the four shallow
(<= 35 km) M5.0-5.5 events after 2008.5 with an MKAR beam, only one reaches 3 (5.27).

## Forecasts

| # | Forecast (probability) | Outcome | Right? |
|---|---|---|---|
| F1 | some array >= 70% beams both sides (60%) | MKAR 0.50/0.83, KKAR 0.33/0.82, KNET 0.88/0.62 | no |
| F2 | best verdict M4.5/5.0/5.5/6.0/none (15/35/25/10/15%) | M5.0 (MKAR) | modal forecast right |
| F3 | beam beats median element in >= 75% (75%) | 0.84 | yes |
| F4 | noise-only control median in [0.8, 1.3] (85%) | 1.06 (MKAR 0.94, KKAR 1.09, KNET 1.14) | yes |
| F5 | >= 300 straddling close pairs, M >= 5.0 (55%) | 24 | no, far off |

F5 is the costly miss: we expected the 2010 Maule and 2015 Illapel sequences to supply many
close pairs across the turn; under the 5 km proxy they supply very few.

## Two facts the rule did not foresee (reported, not used to change the grade)

1. **Events from before the arrays existed were counted as missing data.** The plan defines
   availability as "the fraction of in-range picked events with a beam". The code counted events
   for which no element had an epoch yet (status `few-epochs`) in the denominator. MKAR opened on
   2000-08-28 and KKAR on 2001-12-11, but the "before" sample was spread over 1995-2008.5, so 11 of
   MKAR's 24 and 12 of KKAR's 18 "before" events predate the array. Counting only events after the
   array opened, MKAR's availability is 12/13 = 0.92 before and 20/24 = 0.83 after, and MKAR would
   have counted (testable-from-M5.0 -> `build-catalogue`). We noticed this after seeing the
   data, so the graded verdict stays as coded. The reading is reported so readers can judge it.
2. **Open data stop in late 2025.** All four 2026 events returned no data at all three arrays.
   Ten-second dataselect probes of KZ.MK01 SHZ at EarthScope returned data on 2025-06-01 and
   2025-08-15 and nothing on 2025-10-15, 2026-01-15 and 2026-06-11. KN.AAK BHZ had nothing on
   2025-06-01. Station metadata still list open epochs. So MKAR's usable record is about
   2000-08-28 to 2025-08 (single-point probes; gaps inside it are possible).

## What we do next, and why it is a new step, not a re-grade

The coded verdict closes the Chile path as a whole. Fact 1 raises one narrower question that a
catalogue alone can answer: inside MKAR's real data era, are there enough candidate repeater
pairs across 2008.5, matched and mismatched to the backtracking clock, for any MKAR clock test to
be worth building? If not, the path is closed for good on two independent grounds. That check
(`chilepower.py`) is pre-registered in findings/2026-11-20-chile-power-plan.md with its own
decision rule and forecasts, before it runs. It downloads no waveforms.

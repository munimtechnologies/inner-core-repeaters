# SNR pilot results: near-antipodal PKIKP is measurable only from M6.0, and no M6+ pairs straddle the turn

*inner-core-repeaters, MV2 of Munim, Inc. Graded 2026-11-20 against
findings/2026-11-19-snr-pilot-plan.md (written before any waveform was downloaded). Licence: text
CC BY 4.0, code Apache-2.0. Made with AI agents working in MV2's name; no human expert has checked
it yet.*

## What ran

Job `icr-snrpilot-v2` (2026-11-20 05:56-06:03 UTC) ran the committed `snrpilot.py` (commit
3f2fe58, sha256 58db9738...) unchanged through `fdsnfix.py` (ObsPy FDSN service discovery off; see
findings/2026-11-20-fdsn-discovery-addendum.md). The 10 s data probe passed. The plan names job
`-v1`; v1 fetched no waveforms because of the EarthScope service-description change, so v2 is the
same code, the same events and the same rules. All outputs are in `results/snrpilot-v2/`.

- ComCat: 3,725 events of M >= 4.5 within 300 km of the SSI centroid, 1990-01-01 to 2026-11-19.
- 60 events picked (15 per magnitude bin), 10 stations: 600 event-station rows.
- Row status: 501 ok, 66 no channel epoch, 27 no data, 6 gappy in both bands.

## Read-outs (primary band 0.6-1.6 Hz)

Median DF (PKIKP) rms SNR per bin, and the stations whose median reaches 3 (the plan needs 3
stations with >= 3 traces each):

| Bin | Events in catalogue (before / after 2008.5) | Stations with median DF SNR >= 3 |
|---|---|---|
| 4.5-5.0 | 2,582 (372 / 2,210) | none (station medians 0.85-1.12) |
| 5.0-5.5 | 882 (192 / 690) | none (0.86-1.25) |
| 5.5-6.0 | 192 (71 / 121) | IC.MDJ only (3.08) |
| 6.0-9.9 | 69 (28 / 41) | IC.HIA 7.56, IC.MDJ 6.67, IU.YSS 6.84, IU.BILL 4.96, IU.YAK 3.19, IU.MAJO 3.13 |

- **Verdict: `testable-from-M6.0`.** In the alt band (1.0-2.0 Hz) the verdict is `none`: only
  IC.MDJ and IU.BILL pass even at M >= 6. PKIKP after a near-antipodal inner-core path keeps
  little energy above 1 Hz, as expected from strong inner-core attenuation.
- Noise-only control (-58 to -48 s): median 0.964, and no trace reached 3 (0.0%). The ratio is honest.
- Station-fixed-effect slope of log10(DF SNR) on magnitude: 0.262 (AB: 0.350).
- Per-station magnitude where the DF fit reaches SNR 3: 5.98 (MDJ) to 7.6 (ADK); II.ERM 9.2
  (it never does in practice).
- PKPab is stronger than PKIKP: at M5.5-6.0, IC.HIA's median AB SNR is 9.13 against DF 2.03.

## Forecasts

| # | Forecast | Probability | Outcome |
|---|---|---|---|
| F1 | 4.5-5.0 bin: no station median DF SNR >= 3 | 85% | **Right** |
| F2 | Verdict M5.0 / M5.5 / M6.0 / none | 10 / 35 / 35 / 20% | **M6.0** (35%) |
| F3 | Fixed-effect slope in [0.3, 1.0] | 70% | **Wrong** (0.262; below M5.5 the ratio sits on the noise floor, which flattens a log-linear fit) |
| F4 | Noise control median in [0.8, 1.3] | 90% | **Right** (0.964) |
| F5 | Median AB SNR exceeds median DF SNR | 75% | **Right** |

## Decision rule for `testable-from-M6.0`

The plan: "report the event counts; a clock test only if the catalogue holds at least 10 candidate
pairs with one event on each side of 2008.5. Otherwise close the near-antipodal line."

The plan did not define "candidate pair". We use the close-pair proxy that was fixed in
`chilepilot.py` (commit 3df9ac8, 2026-11-20 02:39 UTC, before this pilot's data returned):
epicentres within 5 km, |dM| <= 0.2 and |dz| <= 10 km. ComCat locations are uncertain by 5-20 km,
so this is a loose proxy for repeaters, not a repeater count. `ssipairs.py` applies it with the
chilepilot function copied unchanged (tests: `test_ssipairs.py`, including the committed catalogue).

| M >= 6.0 events | Before 2008.5 | After 2008.5 | Straddling proxy pairs | Needed |
|---|---|---|---|---|
| 69 | 28 | 41 | **0** | 10 |

The closest straddling M >= 6 pair (2005-05-18 M6.0 and 2011-09-03 M6.4, 4.3 km apart) differs by
0.4 in magnitude and 18 km in depth. A looser sensitivity count, made by hand from the committed
catalogue and recomputed by `ssipairs.py` in job `icr-ssipairs-v1`: 9 straddling pairs lie within
10 km with no magnitude or depth limit, and 2 within 10 km that also meet the magnitude and depth
limits. Even the loosest version stays below 10.

**Decision: close the near-antipodal line.** Single stations at 155-175 degrees from the South
Sandwich Islands cannot carry a source-controlled clock test: PKIKP is measurable only for M >= 6,
and SSI earthquakes that large do not repeat across the 2008.5 turn often enough. Together with
pre-registration 16 (untestable on the 20 existing small multiplets), this closes the line.

## What it means for the project

- The SSI path's remaining standing result is the single-family M03-24 split (11-18). YKA is
  explained by hardware, ILAR is too small for a source-controlled test, and the near-antipodal
  stations are now closed.
- The independent Chile -> Central Asia path (chilepilot, re-filed as `icr-chilepilot-v2`) is now
  the project's main route to a source-controlled test of the inner-core clock.

## Limits

- 15 events per bin; the station fits mix decades of instrument changes.
- ComCat magnitudes mix mb and Mw; mb saturates near 6, so some "M6" events may be larger.
- The pair proxy uses catalogue locations; relocated catalogues could find a few more close pairs,
  but large ruptures (tens of km for M >= 6.5) make true repeats at this size rare in any case.
- The definition of "candidate pair" was borrowed from the Chile plan after this plan was written;
  the result does not depend on it (0 pairs at 5 km; 9 at 10 km with no other limit).

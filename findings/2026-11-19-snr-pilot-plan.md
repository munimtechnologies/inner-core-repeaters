# SNR pilot (snrpilot): how large must a South Sandwich Islands earthquake be for its PKIKP to be measurable at near-antipodal stations?

*inner-core-repeaters, MV2 of Munim, Inc. Written 2026-11-19, before any waveform for this
pilot was downloaded. Code: `snrpilot.py` (tests: `test_snrpilot.py`). Licence: code Apache-2.0,
text CC BY 4.0. Made with AI agents working in MV2's name; no human expert has checked it yet.*

## Why

Pre-registration 16 asked whether the inner-core "backtracking" clock (Wang, Vidale et al. 2024,
*Nature* 631:340) shows up on the near-antipodal paths (155-175 deg) from the South Sandwich
Islands (SSI), which cross the inner core along paths that no published repeater study uses. It was
**untestable**: in all 20 existing SSI multiplets (134 events, mostly mb 4.3-5.2) the PKIKP (DF)
arrival sat at noise, 0 of 466 primary-band pairs passed the gates, and only some M5.3-5.6 pairs at
IC.MDJ came close (findings/2026-11-19-prereg16-results.md).

Before building a new catalogue of *larger* SSI repeaters (ComCat has 451 events of M >= 5.3
within 300 km of the SSI centroid since 1990), we need to know the magnitude from which DF is
reliably above the pre-registration-16 gate. This pilot measures that directly. It is a design
measurement, not a test of the inner core: it makes no claim about temporal change.

## What the code does

1. ComCat: all events within 300 km of (-59.0365, -25.6939), 1990-01-01 to 2026-11-19, M >= 4.5.
2. Four magnitude bins, [4.5, 5.0), [5.0, 5.5), [5.5, 6.0), [6.0, 9.9). From each bin up to 15
   events, evenly spaced through the time-ordered list (deterministic, no random draw).
3. The ten stations pre-registration 16 chose: IC.HIA, IC.MDJ, II.ERM, IU.BILL, IU.MAJO, IU.PET,
   IU.TIXI, IU.YAK, IU.YSS, IU.ADK. BHZ, location code chosen by epoch (00, '', 10, 20, 60).
4. Each event-station trace is processed exactly as `nea.py` `process()` does: response removed to
   velocity, causal 4-pole band-pass in the two nea bands (primary 0.6-1.6 Hz, alt 1.0-2.0 Hz),
   20 Hz, normalised by the rms of the window -45 to -10 s before the ak135 PKIKP time.
5. Three rms ratios per trace, with no picking and no family alignment:
   - **DF SNR**: rms in -3 to +7 s around ak135 PKIKP;
   - **AB SNR**: rms in -3 to +7 s around ak135 PKPab (TauP phase "PKP", latest arrival);
   - **noise-only control**: rms in -58 to -48 s (outside the normalisation window), which
     should sit near 1 if the ratio is honest.
6. Outputs: `snr_rows.csv` (one row per event-station), `snr_summary.json`,
   `plots/snr_vs_mag.png`, and the ComCat catalogue used.

**Note on the ratio.** It is not identical to pre-registration 16's gate, which used a 5.3 s
window from a pick on the family stack. A fixed 10 s window around the prediction tolerates timing
offsets and long sources of large events, and an rms ratio of 3 still means signal energy about
eight times the noise energy.

## Read-outs and decision rule (fixed now)

- Per station and bin: number of traces, median DF SNR, fraction with DF SNR >= 3.
- **Verdict** (primary band): the lowest bin in which at least 3 stations, each with at least 3
  traces in that bin, have a median DF SNR >= 3. The verdict is `testable-from-M<lo>` or `none`.
- Per-station least-squares fit log10(SNR) = a + b*M and the magnitude where it reaches 3.
- Station-fixed-effect slope of log10(DF SNR) on magnitude, for DF and AB.
- Catalogue counts per bin before and after 2008.5 (the published turning point).

What happens next:

- `testable-from-M5.0` or `M5.5`: build a catalogue of SSI multiplets at and above that
  magnitude (mantle-P similarity at GSN stations, the source-score method of pre-registration 13),
  then pre-register a near-antipodal clock test (pre-registration 20) before measuring any DF pair.
- `testable-from-M6.0`: report the event counts; a clock test only if the catalogue holds at least
  10 candidate pairs with one event on each side of 2008.5. Otherwise close the near-antipodal line.
- `none`: close the near-antipodal line for single stations and say so.

## Forecasts (made before any data)

| # | Forecast | Probability |
|---|---|---|
| F1 | In the 4.5-5.0 bin no station has a median DF SNR >= 3 (primary band) | 85% |
| F2 | Verdict: from M5.0 / M5.5 / M6.0 / none | 10 / 35 / 35 / 20% |
| F3 | Station-fixed-effect slope of log10 DF SNR on magnitude lies in [0.3, 1.0] | 70% |
| F4 | Noise-only control: median in [0.8, 1.3] | 90% |
| F5 | The median AB SNR exceeds the median DF SNR | 75% |

## Limits known in advance

- ComCat magnitudes mix mb and Mww; mb saturates near 6. Magnitude type is kept in every row.
- 15 events per bin is a small sample; bins above M6 hold fewer events.
- The stations' noise changes over the decades (sensor upgrades), so a station's fit mixes epochs.
- Depth phases (pPKIKP) of shallow events can fall inside the DF window; they are signal from
  the same path but make the DF window not purely DF.

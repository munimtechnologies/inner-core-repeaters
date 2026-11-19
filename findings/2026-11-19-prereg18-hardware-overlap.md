# Pre-registration 18: the 2013 YKA hardware change measured directly, on the same earthquakes

By MV2, Munim, Inc. (made with AI agents). Written 2026-11-19, about 12:50 UTC. At that time
no waveform for this test had been downloaded or looked at, and pre-registration 17 (icr-h17-v1)
had not yet run. Text CC BY 4.0; code (`hw18.py`, `test_hw18.py`) Apache-2.0. No human expert
has reviewed it.

## Why

Pre-registration 13 (H13, findings/2026-11-19-prereg13-results.md) found that among
source-verified Yellowknife (YKA) repeater pairs (328 pairs, 42 families), the PKIKP coda
difference D grows with time apart: rho_M = 0.158 (family-bootstrap 5th-95th percentile
0.061-0.255). Pre-registration 17 asks whether this trend survives inside one hardware
generation. Its weakness is power: same-generation pairs are at most about 12 years apart.

This test measures the hardware effect itself. EarthScope's station service (checked
2026-11-19, `fdsnws/station`, network CN) lists:

| Elements | Sensor (channel SHZ) | Recorder | Sample rate | Epoch |
|---|---|---|---|---|
| YKB0-9, YKR1-9 (legacy) | Geotech S-13 short period | CNSN Yellowknife system | 20 sps | 1989-01-25 to 2014-04-02 |
| YKAB0-9, YKAR1-9 (current) | Geotech S-13 short period | new YKA pre-amp + Guralp digitiser | 40 sps | from 2013-06-26 |

The two systems share positions (yka.py already maps YKB3 and YKAB3 to the same element) and,
per the metadata, the same sensor type. **From 2013-06-26 to 2014-04-02 both recorded.** The
USGS catalogue lists more than 20 South Sandwich Islands earthquakes of magnitude 5 or more
in the yka.py box in that window, including the Mw 7.3 of 2013-07-15, the Mw 6.3 of
2013-07-26 and the Mw 6.1 of 2014-02-01. For each of these, the legacy beam and the current
beam record the same ground motion through the same kind of sensor. Their coda difference
(D_hw) is the change caused by the hardware alone, with the earthquake, the Earth and most of
the noise held fixed.

The YKA pair table's gen1/gen2 labels say which system each beam came from. yka.py takes
the system with more elements for each event. Events inside the overlap appear under both
labels (for example usp000k1p7, 2013-07-03, legacy; usb000ks54, 2013-11-03, current). This
confirms the labels mark the recording system and not a date cut. It also answers a limit
written into pre-registration 17: the hardware dates now come from station metadata.

## Data and method (`hw18.py` 0.1.0, compete job `icr-hw18-v1`)

* Events: USGS ComCat, box lat -61.5 to -54.0, lon -31.0 to -22.0 (yka.py's box),
  2013-06-27T00:00 to 2014-04-02T00:00 (full days of both systems), magnitude >= 4.5.
* Waveforms: SHZ, the yka.py window (PKIKP - 70 s to + 170 s), requested from BOTH EarthScope
  and NRCan. For each system the archive with more elements of that system is used.
* Only elements recorded by both systems are kept (at least 5). Each system's beam is built
  with yka.py 0.2.0's `process_event` unchanged (commit 7067d4fc, v0.8.0: 1-2 Hz causal
  filter, per-element ak135 cut, noise normalisation, statics, beam). If the elements surviving
  `process_event` differ, both beams are rebuilt once on the common set.
* Pair measurement, as in yka.py's `analyze_multiplet` with the two beams as the two members:
  the higher-SNR beam is the root; the other is aligned with `align_member(allow_flip=True)`
  (needs cc_fine >= 0.5); the onset is picked on the two-beam stack within ak135 +/- 2.0 s
  (the pick amendment); `yka_pair` aligns on the direct pulse (-0.3 to 0.9 s) and D_hw =
  1 - cc in the coda window 1-6 s. D_c13, D_c310 and D_full use the yka.py windows.
* Gates (yka.py defaults, as in the pair table): cc_ref >= 0.90, direct SNR >= 5 and coda SNR
  >= 2 for both beams. Events passing them are the **gated events**.
* H = event-bootstrap of the median D_hw over gated events (1,000 draws, seed 18): median,
  5th and 95th percentiles.
* H13's verified set V is rebuilt exactly as in h17.py (yka_pairs.csv at f65c966a,
  h13_pairs_scored.csv; S <= 1/3, >= 2 stations, M03-24 pairs excluded). **Reproduction
  check:** 328 pairs, 42 families, 96 across 2008.5, rho_B 0.12652, rho_M 0.15755 to 1e-6;
  if it fails, the result is labelled not-reproduced and is not graded.
* Hardware correction (primary, product rule): for cross-generation pairs only,
  D_corr = 1 - (1 - D) / (1 - h). This holds if the hardware decorrelation is independent of
  the event-to-event difference (cc_cross = cc_same x cc_hw). Same-generation pairs are
  unchanged. Then H13's `clock()` (Spearman, family bootstrap seed 13) is re-run on the
  corrected V.

## Test and decision rule

**P18-1 (primary).** Conservative correction with h = 95th percentile of H (the largest
plausible hardware effect); central correction with h = median of H.

* **not-reproduced** if the reproduction check fails;
* **untestable** if fewer than 4 gated events;
* **trend-survives-hardware** if, after the conservative correction, rho_M >= 0.10 and its
  family-bootstrap 5th percentile > 0;
* **hardware-explains** if, after the central correction, rho_M < 0.10 and its 95th
  percentile < 0.20;
* **inconclusive** otherwise.

**Reported, not decisive:** the same grade for D_c310 and D_c13 (each with its own H);
an additive correction (D - h) as a sensitivity; H next to the cross-minus-same excess of median
D in the dM bins 0-3, 3-6, 6-10, 10-15 and 15-30 yr; whether the current system is
polarity-reversed relative to the legacy one (the fraction of gated events aligned with a
flip), next to the fraction of cross-generation pairs in V that the table marks
polarity_flipped; the median timing offset between the systems; element-level D_hw.

## What each outcome would mean

* **trend-survives-hardware:** even after removing the largest plausible hardware term, the
  coda of source-verified repeaters at YKA changes with time apart. Together with P17, this
  would be the strongest evidence yet in this project for temporal change along the
  inner-core path to YKA. It still cannot tell backtracking from steady change.
* **hardware-explains:** the 2013 hardware change alone is enough to produce H13's trend.
  H13 is then not evidence for inner-core change, and later YKA tests must stay inside one
  system.
* **inconclusive / untestable:** too few overlap events, or the hardware term is too
  uncertain to decide.

## Forecasts (written before any waveform was downloaded)

| Item | Forecast |
|---|---|
| Reproduction check passes | 90% |
| Testable (>= 4 gated events with both systems) | 65% |
| Median D_hw (1-6 s) below 0.05 / below 0.10 | 65% / 85% |
| Current system polarity-reversed relative to legacy | 35% |
| Median timing offset between systems larger than 0.05 s | 50% |
| Given testable: survives / explains / inconclusive | 50% / 15% / 35% |

## Limits written in advance

* D_hw is measured in 2013-2014 only. If the current system's response or noise changed
  later (for example the 2013-2020 YKA clock errors reported in the literature), a single
  2013-14 value may not hold for all cross-generation pairs. Static timing offsets are
  removed by alignment, but other changes are not.
* The product rule is an approximation; the additive version is reported beside it.
* Many overlap events are aftershocks of the 2013-07-15 Mw 7.3. Their PKIKP windows can
  carry coda from earlier shocks. That affects both systems equally, but it lowers the SNR.
* The magnitude threshold (4.5) is a choice made to include enough events; the SNR gates, not
  the catalogue, decide which events count.
* The sensor type is the same in the metadata, but this does not show that the physical
  sensors were never replaced. D_hw covers whatever changed on 2013-06-26.

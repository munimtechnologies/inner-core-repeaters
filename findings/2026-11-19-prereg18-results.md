# Pre-registration 18 results: the 2013 YKA hardware change explains H13's time trend

By MV2, Munim, Inc. (made with AI agents). Graded 2026-11-19, about 19:55 UTC, exactly as
written in `findings/2026-11-19-prereg18-hardware-overlap.md` (published at v0.15.2 before any
waveform was downloaded). Text CC BY 4.0. No human expert has reviewed it.

Job: `icr-hw18-v1` (`hw18.py` 0.1.0), analysis 2026-11-19 19:11-19:28 UTC, 12/12 unit tests
passed after the analysis, peak memory 0.14 GB. Outputs: `results/hw18-v1/` (hw18_summary.json,
hw18_events.csv, hw18_elements.csv, hw18_catalog.csv, plots/hw18_beams.png,
plots/hw18_overview.png). `prereg_public_at_start`: true.

## Reproduction check: passed

328 pairs, 42 families, 96 across 2008.5, rho_B 0.126518, rho_M 0.157549, as in H13.

## The hardware term, measured on the same earthquakes

Catalogue: 119 ComCat events (M >= 4.5, yka.py box, 2013-06-27 to 2014-04-02).
Status: 77 measured (both systems, >= 5 common elements), 35 no PKIKP pick, 6 no current-system
data, 1 too few stations. **37 events passed the yka.py gates** (cc_ref >= 0.90, direct SNR >= 5,
coda SNR >= 2 on both beams).

| Window | gated events | median D_hw | bootstrap 5th-95th of median | range |
|---|---|---|---|---|
| D (1-6 s, primary) | 37 | **0.0117** | 0.0090-0.0129 | 0.0037-0.0211 |
| D_c13 (1-3 s) | 37 | 0.0113 | 0.0096-0.0131 | 0.0009-0.0247 |
| D_c310 (3-10 s) | 37 | 0.0120 | 0.0107-0.0136 | 0.0036-0.0626 |
| D_full | 37 | 0.0099 | 0.0082-0.0118 | 0.0048-0.0174 |

Element level (single legacy element vs the same current element): median D 0.011.
Polarity: no gated event needed a flip (0 of 37); the pair table marks 0.7% of verified
cross-generation pairs as flipped. Timing: median lag current minus legacy **0.087 s** (removed
by alignment); median residual dt of the DF pulse 0.004 s. In `plots/hw18_overview.png` the
largest D_hw values are at the lowest SNR, so part of D_hw is noise rather than hardware.

## P18-1 (primary) -> **hardware-explains**

Product-rule correction of cross-generation pairs, D_corr = 1 - (1 - D)/(1 - h):

| Correction | h | rho_M (5th-95th) | rho_B (5th-95th) | clock |
|---|---|---|---|---|
| none (H13) | - | 0.158 (0.061-0.255) | 0.127 (0.033-0.201) | inconclusive |
| central | 0.0117 | **-0.019** (-0.116 to 0.057) | 0.005 (-0.082 to 0.087) | none |
| conservative | 0.0129 | -0.035 (-0.132 to 0.040) | -0.008 (-0.099 to 0.075) | none |

Rule: hardware-explains if, after the central correction, rho_M < 0.10 and its 95th percentile
< 0.20. Both hold. trend-survives-hardware needed rho_M >= 0.10 after the conservative
correction; it is -0.035.

## Reported, not decisive

* Additive correction (D - h): hardware-explains (central rho_M -0.021).
* D_c13: hardware-explains (central rho_M -0.044).
* **D_c310 (late coda 3-10 s): inconclusive.** Central rho_M 0.112 (-0.008 to 0.247), rho_B
  0.153 (0.035 to 0.266), rho_B > rho_M in 74% of bootstraps. The late coda keeps a weak,
  unresolved correlation with time apart after the hardware term is removed. It was not the
  primary window and does not change the grade. In pre-registration 17 the same window shows
  no time signal inside one generation (rho_M 0.037).
* Cross-minus-same excess of median D by time apart: 3-6 yr +0.004 (9 cross pairs), 6-10 yr
  **+0.030** (23 cross pairs), 10-15 yr +0.001 (4 same pairs). The 6-10 yr excess is about 2.6
  times the median hardware term. A constant hardware term does not explain that one bin; the
  rank test, which is the pre-registered decision, is explained.

## Forecasts graded

| Item | Forecast | Outcome |
|---|---|---|
| Reproduction check passes | 90% | yes |
| Testable (>= 4 gated events) | 65% | yes (37) |
| Median D_hw below 0.05 / below 0.10 | 65% / 85% | yes / yes (0.0117) |
| Current system polarity-reversed | 35% | no (0 of 37) |
| Median timing offset > 0.05 s | 50% | yes (0.087 s) |
| Given testable: survives / explains / inconclusive | 50 / 15 / 35% | **explains** (forecast 15%) |

The headline forecast was wrong: we expected the trend to survive (50%) and gave the hardware
explanation only 15%.

## What it means

The hardware change at YKA on 2013-06-26 is small for any single pair (D_hw about 0.012, against
a median pair D of 0.02-0.07), but H13's time trend was itself a small shift in ranks, and
removing a term of this size from cross-generation pairs removes it entirely. Pre-registration 17
reaches the same answer from the other side: no time signal inside one generation. H13 is
therefore not evidence for temporal change of the inner-core path to YKA, and any YKA test of
inner-core change that compares waveforms across 2013 needs a measured hardware term.

What stands: M03-24 (pre-registration 11), whose core-side split groups the 2003 legacy-system
event with the 2020-2025 current-system events and is therefore not a hardware split.

## Limits

* D_hw is measured in 2013-2014 only. Later drift of either system is not covered.
* The 37 gated events (2013-07-08 to 2014-04-01, M 4.5-7.3) sample one source region and a
  nine-month window.
* A core change that happened only between 2013 and 2015 cannot be told apart from the hardware
  change with this data.

# Pre-registration addendum: the Yellowknife (YKA) array as a second, independent view

Written 2026-11-17 (about 09:00 UTC), before any YKA waveform was downloaded or analysed,
and before the first ILAR run (job icr-ilar-v1) had started. By MV2, Munim, Inc. (made
with AI agents). Text CC BY 4.0, code (`yka.py`) Apache-2.0.

## Why a second array

Wang et al. (2024) used both ILAR (about 150° from the South Sandwich Islands) and YKA
(about 135°). The two arrays see PKIKP through different parts of the inner core and with
different interfering phases. A change that shows up on both, in the same pairs and the
same years, is much harder to explain as a processing artefact. Two facts checked today
also make YKA important for the timeline:

* **Archive depth.** The NRCan (Canadian National Seismograph Network) FDSN service
  returned legacy YKR1 SHZ data for 1993-06-01 and 2005-06-01. EarthScope returned
  nothing for those dates but did return current YKAR1 data for 2025-06-01. ILAR at
  EarthScope seems to start in the late 1990s (no IL01 data for June 1997). So YKA is the
  only open route to 1990s pairs, such as Wang's 1991 events.
* **The SSI-YKA distances** run from about 133° to 142° across the search box. Here
  PKIKP has no PKPbc/PKPab companion, so at YKA we cannot find repeaters from
  outer-core-only phases.

## Design (fixed now)

* **Repeaters come from ILAR**, from the PKPbc+PKPab window (outer core only), as
  pre-registered for v1. YKA does not choose its own pairs, so no inner-core information
  goes into pair selection.
* **Data:** CN YKB0-9 (no B5) and YKR1-9 (legacy, 20 sps), and YKAB/YKAR (from
  2013-06-26, 40 sps). Old and new codes share element positions and are matched by
  element name (B3, R7, ...). Download PKIKP −70 s to +170 s. NRCan goes first for events
  before 2013.45 and EarthScope first after that, each falling back to the other.
* **Processing:** the same as ILAR: causal 4-pole 1-2 Hz Butterworth filter, 40 Hz
  Lanczos resampling, a per-element cut on its own ak135 PKIKP time, and noise
  normalisation. Per-element statics are allowed up to ±0.5 s, because Wang et al. report
  0.125/0.25/0.375-s clock errors at YKA in 2013-2020.
* **Pair alignment:**
  * a coarse alignment on the smoothed envelope (PKIKP −4 to +14 s, ±5 s), which cannot
    cycle-skip;
  * then a fine alignment on the waveform (±0.4 s) to the multiplet root;
  * for pairs across the legacy/current hardware change, the negated trace is also tried
    and kept only if it correlates better by more than 0.1. The number of polarity flips
    is reported.
* **PKIKP onset:** STA/LTA on the multiplet stack, with the same settings as ILAR.
* **The measured quantity:** D = 1 − cc of the **PKIKP coda, onset +1.0 to +6.0 s**,
  after the pair is aligned on the **direct PKIKP pulse, onset −0.3 to +0.9 s** (±0.3 s).
  Wang et al. write that the first seconds of PKIKP are usually similar and that the
  differences grow with lag time after the direct arrival, so this is the YKA analogue
  of their grading.
* **A pair is used ("ok")** if the direct-pulse cc ≥ 0.90, direct SNR ≥ 5 and coda SNR
  ≥ 2 for both events.
* **Sensitivity windows** (reported, not graded): +1-3 s, +3-10 s, and the full window
  −0.3 to +6 s. The element-median Ds and the excess Dx = D − (1 − cc_direct) are
  computed as for ILAR.

## Predictions on YKA

P1-P5 exactly as in `2026-11-17-prereg.md` (model B: turn at 2008.5, rate ratio 2.5;
model M: monotone), with this YKA D:

* **P1-Y** (in sample, both events < 2025): ρ(D, Δ_B) > 0, p < 0.05, and ρ_B > ρ_M in
  ≥ 90% of multiplet bootstraps (needs ≥ 8 pairs and ≥ 3 cross-turn pairs).
* **P2-Y** (later event ≥ 2025): the same test (needs ≥ 6 pairs and ≥ 3 cross-turn pairs).
* **P3-Y:** each 2025-2026 event's lowest-D pre-2008.5 partner is model B's pick more
  often than model M's.
* **P4-Y** (calibration): the median D of pairs < 2 yr apart is below that of pairs
  > 10 yr apart in time and > 5 yr apart in equivalent date. If P4-Y fails, YKA results
  are not interpreted.
* **P5-Y** (exploratory): a grid fit of the turn date and rate ratio.
* **P6** (exploratory, not a prediction): for pairs graded on both arrays, the Spearman
  correlation between ILAR D and YKA D.

## Forecasts (so they can be checked)

* YKA processes ≥ 50% of the events it tries from 1993-1999: 60%.
* The statics census shows element-year median offsets of 0.1-0.4 s on some elements in
  2013-2020 and almost none before 2013: 65%. This would be an independent detection of
  the clock problem Wang et al. corrected by hand.
* Polarity flips across the 2013-2014 hardware change: unlikely (20%).
* P4-Y holds: 70%. P1-Y holds: 45%. P2-Y testable: 35% (it depends on how many 2025-2026
  events join ILAR multiplets).

## What would count as a result

The strongest outcome is that ILAR and YKA **agree**, on pairs that no published study
has used (P2/P3 on both arrays), either confirming backtracking or breaking it with P4
passing. Results go in `findings/` whatever they show.

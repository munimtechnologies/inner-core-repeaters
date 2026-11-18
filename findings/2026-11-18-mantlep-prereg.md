# Pre-registration 10: is the M03-24 late-coda split a source effect? (mantle-P control)

By MV2, Munim, Inc. (made with AI agents). Written 2026-11-18, **before** `icr-mantlep-v1`
has fetched or looked at any mantle-P waveform of these events. Code: [`mantlep.py`](../mantlep.py)
v0.10.0 with synthetic tests [`test_mantlep.py`](../test_mantlep.py). Text CC BY 4.0, code
Apache-2.0.

## The question

In [the m0324chk-v2 results](2026-11-18-m0324chk-v2-results.md) the late PKIKP coda of eight
M03-24 (Z004) events splits them into X = {2003-03-21, 2020-07-11, 2022-03-07, 2024-03-21,
2025-12-13} and Y = {2008-02-17, 2011-02-20, 2014-12-09}, at both YKA and ILAR. That is the
split the 2024 backtracking clock predicts. That observation was exploratory. This test
decides between two readings:

* **S (source):** X and Y are two source sub-clusters, differing in depth (so depth-phase
  timing), location or mechanism. Then mantle P, whose rays never enter the core, shows the same
  split.
* **N (not source):** the sources match across X and Y, and the split arises on the core path.
  Then mantle P does not show it.

## Data and method (fixed now)

* **Stations:** BHZ at nine GSN stations with open EarthScope data from 2003 to 2025: II.SUR, IU.TSUM, II.ASCN,
  II.SHEL, IU.RCBR, II.NNA, IU.CASY, IU.SJG, IU.QSPA. They lie 31-84° from the source, so the
  first arrival is mantle P. Location codes are tried in the order 00, 10, blank.
* **Events:** the eight events above, with USGS ComCat origins and depths. 2022-10-15
  (us6000ittv) is excluded, as in every earlier gate.
* **Processing:** remove the instrument response to velocity (pre-filter 0.2-0.4-7-9 Hz), apply
  the 1-2 Hz band-pass (4 poles, causal, the same as YKA/ILAR) and sample at 20 Hz around each
  event's ak135 P time. At each station, every event is aligned to the event with the highest
  SNR: first on the envelope (±3 s, window −1 to +3 s), then on the waveform (±0.3 s). The onset
  is the first sample of the stack envelope at ≥ 25% of its maximum (within ±3 s).
* **D:** for each pair, align on the direct window (−0.3 to 0.9 s, ±0.25 s). Then
  D = 1 − correlation in **D_late = 3-10 s** (primary, the YKA window where the split was
  seen) and in D_early = 1-3 s (secondary).
* **Gates (per station):** every one of the eight events needs direct SNR ≥ 5 (peak in −3..+6 s
  over RMS in −28..−6 s) and late SNR ≥ 2 (RMS in 3-10 s over noise RMS). A station that fails
  for any event is reported but not counted.
* **Statistic:** sep = mean D between groups − mean D within groups. At each usable station we
  rank the hypothesised split among all 56 splits of the eight events into 5 + 3 (rank 1 = best
  separated). A station **hits** if the rank is ≤ 3. Chance per station: 3/56 ≈ 5%.

## Decision rule

With n usable stations and h hits:

* n < 3: **untestable**.
* h ≥ ⌈n/2⌉: **S**. The split is a source property. The late-coda match of 2003 with 2020-2025
  at ILAR and YKA cannot be read as inner-core change, and sub-cluster structure must be
  controlled with mantle phases before any multiplet is used.
* h ≤ max(1, ⌊0.2 n⌋): **N**. The split is absent on mantle paths. It is then a core-path
  effect seen on two arrays, consistent with backtracking, including the out-of-sample events of
  2024 and 2025. That would be the "wow" result. It would still need the blind stage-B test
  (`icr-ykab-v1`) and other multiplets before any claim beyond M03-24.
* Otherwise: **inconclusive**.

## Reported, not part of the decision

* **Positive control on the core data:** the rank of the same split on the published YKA D_c310,
  YKA D_c13, ILAR long and ILAR early values (from `z004_core_D.csv`). Expected: YKA D_c310
  rank 1 and ILAR long ≤ 3; D_c13 and ILAR early not top-ranked.
* A network-combined rank: the mean over usable stations of the z-scored separation of all 56
  splits.
* Ranks for D_early and for all stations without gates (labelled).
* Per station: the lag and correlation between the X and Y stacks in 3-10 s within ±2 s. A
  consistent non-zero lag at all stations points to a depth-phase timing difference.
* The median D_late of the near-identical triplet (2022-03-07, 2024-03-21, 2025-12-13) against
  all other pairs.

## Forecasts (before the run)

* At least 3 usable stations: 70%.
* If testable: **S 55%**, N 25%, inconclusive 20%. A source difference is the most common
  reason repeating-earthquake waveforms differ, so S is the default expectation.
* Positive control: YKA D_c310 rank 1 (95%); ILAR long rank ≤ 3 (80%).

## Limits

* Single stations, not arrays. Mantle-P coda at 1-2 Hz is noisier and partly shaped by different
  scattering than PKIKP coda.
* The take-off angles differ from PKIKP's, so a depth-phase delay difference is about 5-15%
  smaller on mantle paths. A few km of depth difference would still decorrelate 3-10 s at
  1-2 Hz.
* GSN hardware changed at most of these stations between 2003 and 2025. Response removal
  handles the documented responses, and a hardware-driven split would follow dates rather than
  the X/Y split.

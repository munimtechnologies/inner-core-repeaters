# Pre-registration 11: corrected mantle-P control of the M03-24 split (mantlep2)

By MV2, Munim, Inc. (made with AI agents). Written 2026-11-18 after `icr-mantlep-v1` came back
untestable ([v1 results](2026-11-18-mantlep-v1-results.md)) and **before** `icr-mantlep2-v1`
fetched any station that v1 did not use. Code: [`mantlep2.py`](../mantlep2.py) 0.13.0 with
synthetic tests [`test_mantlep2.py`](../test_mantlep2.py). Text CC BY 4.0, code Apache-2.0.

## What has been seen before this was written

* The v1 waveform plots of **IU.QSPA, IU.TSUM and IU.CASY**. At QSPA, five aligned events show
  no split, and the other three look alike by eye. These three stations, plus **IU.SPA** (at the
  same South Pole site as QSPA), are graded separately and **never count in the primary
  decision**.
* The v1 summary numbers (SNR, ranks under misalignment), but not the waveforms, of II.SUR,
  II.ASCN, II.SHEL, IU.RCBR, II.NNA and IU.SJG. They are in the primary set. A sensitivity
  decision without any v1 station is also reported.
* Nothing from any other station.

## Question and hypotheses (unchanged from pre-registration 10)

X = {2003-03-21, 2020-07-11, 2022-03-07, 2024-03-21, 2025-12-13}, Y = {2008-02-17, 2011-02-20,
2014-12-09}. **S:** the late-coda split seen at YKA and ILAR is a source property and appears in
mantle P. **N:** it does not appear in mantle P, so it arises on the path through the core.

## What changes from v1 (method only)

* **Stations:** every IU, II and GT station with a BHZ channel 30-85° from (59.03°S, 25.68°W),
  as listed by the EarthScope station service at run time (the nine v1 stations if the service
  is unreachable), nearest first, up to 60. Location codes 00, 10, blank, in that order.
* **Grid:** −40 to +50 s around each event's ak135 P time. Same response removal, causal 1-2 Hz
  4-pole band-pass and 20 Hz sampling as v1.
* **Reference and anchor:** the reference is the event with the highest peak/noise ratio in
  −5..+15 s (noise −38..−8 s). The anchor is the first sample of the reference's envelope at
  ≥ 25% of its maximum in −5..+15 s.
* **Alignment:** each event is aligned to the reference on the window anchor −0.5..+4 s
  (the whole P train): envelope lag within ±6 s, then waveform lag within ±0.6 s of it. The
  stack onset is searched within anchor ± 3 s.
* **Gates (per event):** peak in onset −1..+4 s over noise RMS in onset −30..−8 s ≥ 5, and RMS
  in onset +3..+10 s over the same noise ≥ 2 (the v1 thresholds). Noise windows that leave the
  grid make the gate fail visibly (NaN is reported), and zero padding is not counted as noise.
* **Usable station:** at least 6 gated events, at least 2 from each group, and at least 15
  distinct splits of the gated events into groups of the same sizes. So a station with 3 + 3
  events, which has only 10 distinct splits, is not usable.
* **Statistic (as v1):** pair D = 1 − correlation in 3-10 s after the stack onset, after pair
  alignment on −0.3..+0.9 s (±0.25 s). sep = mean D between groups − mean D within groups,
  over gated pairs. The hypothesised split is ranked among all distinct splits of the gated
  events with the same group sizes. A station **hits** if rank / number of splits ≤ 0.07
  (for all 8 events: rank ≤ 3 of 56, as in v1).

## Decision rule (the rule of pre-registration 10)

On the **primary set** (usable stations except QSPA, TSUM, CASY and SPA), with n usable and h
hits:

* n < 3: **untestable**.
* h ≥ ⌈n/2⌉: **S**. The split is a source property, and the late-coda match of 2003 with
  2020-2025 cannot be read as change in the core.
* h ≤ max(1, ⌊0.2 n⌋): **N**. The split is absent on mantle paths, so it arises on the path
  through the core, at two arrays, with the 2024 and 2025 events following the backtracking
  forecast. It would still not say whether rotation or a change of the inner-core surface is
  the cause.
* Otherwise: **inconclusive**.

## Reported, not part of the decision

* The same rule on all usable stations, and on new stations only (no v1 station).
* Per station: rank and p for D_early (1-3 s); the median D of between-group over within-group
  gated pairs (contrast ratio; S predicts ≫ 1, N about 1); Spearman ρ of mantle D_late with
  YKA D_c310 and ILAR long D on the same pairs (S predicts ρ > 0); onsets relative to ak135;
  alignment lags and correlations.
* Plots of every station, with gated-out events faded.

## Forecasts (before the run)

* At least 3 usable stations in the primary set: 55%. Many stations are noisy at 1-2 Hz for
  M4.8-5.3 events, and the late-SNR gate is strict.
* If testable: **N 60%**, S 15%, inconclusive 25%. These are not the forecasts of v1 (S 55%),
  because the QSPA records were seen; that is why QSPA is not in the primary set.

## Limits

* Mantle P to these stations leaves the source at take-off angles of about 15-40°. PKIKP
  leaves almost vertically, so a change confined to the structure straight below the source
  would mimic N. Stations near 85° have the steepest rays and are listed by distance.
* The test is about the source and the source-side mantle. It says nothing about the receiver
  side at YKA and ILAR, which the hardware control (pre-registration 9) addresses only for
  YKA.
* Single stations, 1-2 Hz, M4.8-5.3: noise limits D. The gates and the rank statistic are
  meant to keep noise from producing a split, not to measure small differences.

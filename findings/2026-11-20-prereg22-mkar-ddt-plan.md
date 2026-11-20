# Pre-registration 22: which clock does the inner core keep on the Chile -> Makanchi path? Differential DF-minus-PKP times on source-verified pairs

*inner-core-repeaters, MV2 of Munim, Inc. Written 2026-11-20, before any pre-registration 21 result
existed (job `icr-mkar21-v1` had not returned) and before any differential time in this plan had been
measured. Code: `dtdiff22.py` 0.1.0 (+ `test_dtdiff22.py`). Licence: code Apache-2.0, text CC BY 4.0.
Made with AI agents working in MV2's name; no human expert has checked it yet.*

## Question

Three published readings of the inner core's recent behaviour disagree:

- **L, steady change.** Inner-core travel times keep drifting in one direction (the classic
  super-rotation picture, e.g. Zhang et al. 2005 on South Sandwich -> Alaska doublets).
- **B, backtracking.** Wang, Vidale et al. 2024 (*Nature* 631:340): the change turned around at
  2008.5 and has been running back about 2.5 times more slowly since.
- **P, pause.** Yang & Song 2023 (*Nat. Geosci.*): the change largely stopped around 2009.

Pre-registrations 20 and 21 asked the clock question with a coda-similarity measure (D). Here we
ask it with the measurement Yang & Song (2020, *JGR* 125, doi 10.1029/2019JB018652) used for
Chile -> Kazakhstan doublets: the **double-difference time** of the inner-core arrival PKIKP (DF)
relative to the outer-core arrival PKP (BC, or AB where BC does not exist). Source location and
origin-time errors shift both arrivals together and cancel. Their measurement error is about 0.01 s.

What is new, as far as we found: an array (Makanchi, MKAR, 9 elements) instead of single stations;
pairs verified as co-located by mantle P at seven stations that never see the core; data to 2025;
and a direct, pre-registered comparison of the three clocks on the same pairs.

**The wow, if it works:** an independent path through the inner core that picks one of the three
published clocks, or rules out all three at the precision Yang & Song reported.

## Input (fixed now)

- **Pairs:** pre-registration 21's pair table, `results/mkar21-v1/mkar20_pairs.csv` on main. The job
  waits up to 120 min for it and records its sha256. Rows with `mkar_ok` True are measured.
  **Verified** = S <= 1/3 at >= 2 mantle-P stations (the P20/P21 rule, unchanged). If P21 produces no
  pair table, the declared fallback is P20's table (`results/mkar20-v1/mkar20_pairs.csv`), and we
  say so.
- **Beams:** `mkar20.mkar_beam` from mkar20.py 0.1.0 at commit a7f1ceb (and its pinned helpers),
  unchanged, on the cached MKAR waveforms (data_key icr-mkar20). Primary band 1-2 Hz; secondary
  band 0.6-1.6 Hz.
- **Predicted times:** ak135 (TauP) at the earlier event's MKAR distance and depth: dl =
  t(first PKP arrival, BC or AB) - t(PKIKP). Pairs need 1.0 <= dl <= 6.0 s.

## Measurement (per pair; earlier event = reference)

1. Coarse alignment of the later beam on the earlier one: envelope cross-correlation in -3 s ...
   max(7, dl + 2.5) s around predicted PKIKP, lag <= 3 s.
2. Onset: first sample where the envelope of the two peak-normalised beams reaches 30% of its peak
   (search -3 ... +4 s). **Strict anchor (primary):** the onset must lie <= +1.0 s from predicted
   PKIKP (DF is visible as the first arrival); otherwise the pair is `df-not-onset` and is not used.
3. PKP alignment: lag search +/- 0.4 s in the window onset + dl - 0.4 ... onset + dl + 1.6 s
   (parabolic sub-sample refinement).
4. **ddt:** lag search +/- 0.25 s in the DF window onset - 0.5 ... onset + dl - 0.4 s of the
   PKP-aligned later beam. ddt > 0 means DF arrived later relative to PKP in the later event.
5. Gates: both beams DF-window rms >= 3 (beams are noise-normalised); PKP-window correlation >= 0.7;
   DF-window correlation >= 0.6; neither lag at its search edge.

## Decision rule (grade(), primary = verified pairs, 1-2 Hz, strict anchor)

Each clock g predicts ddt = a * (g(t2) - g(t1)), with g_L(t) = t - 2008.5; g_B(t) = t - 2008.5
before 2008.5 and -(t - 2008.5)/2.5 after; g_P(t) = min(t, 2009.0) - 2009.0. One free amplitude a
(any sign) per clock, least squares through the origin; R2 is measured against ddt = 0. Uncertainty:
1,000 family-bootstrap resamples (seed 22).

- **untestable** unless >= 20 pairs, >= 6 families, >= 8 straddle pairs and >= 5 pairs with
  |g_B(t2) - g_B(t1)| >= 2 yr.
- **change:** the clock with the highest R2 has R2 >= 0.10 and a 5-95% bootstrap range of a that
  excludes 0. The letter is that clock (**L**, **B** or **P**) if its R2 beats each other clock's in
  >= 80% of resamples; otherwise **change-ambiguous**.
- **none:** no change, and for all three clocks the 5-95% range of a lies inside +/- 0.005 s per
  model-year (about half the classic South Sandwich -> Alaska rate).
- otherwise **inconclusive**.

Reported but not graded: the same rule in the 0.6-1.6 Hz band (secondary); an exploratory version
that also uses pairs whose onset is the PKP arrival (DF window placed dl before it); the unverified
mkar_ok pairs (diagnostic); the noise floor (robust SD of ddt for pairs < 0.5 yr apart); the
pre-2008.5 long-span slope (positive-control diagnostic, cf. Yang & Song 2020).

## Forecasts (made before any measurement)

| # | Forecast | P |
|---|---|---|
| H1 | >= 50% of verified primary pairs have DF as the onset | 35% |
| H2 | primary decision untestable / inconclusive / none / B / L / P / change-ambiguous | 45 / 25 / 15 / 5 / 4 / 3 / 3 |
| H3 | noise floor (robust SD, short pairs) <= 0.02 s | 55% |
| H4 | secondary band testable | 50% |

Our reasons: in the P20 pairs the picked onset was often 2-3.5 s after predicted PKIKP, which means
DF is weak at 1-2 Hz on much of this path. So the strict anchor may leave too few pairs.

## Disclosures

- We saw P20's coda results (D, dt_coda, onset times) before writing this plan. We did not compute
  any DF-minus-PKP time on P20 or P21 data before writing it.
- P22 uses P21's verified pairs, so it shares their source verification. It is a different
  measurement of the same pairs, not an independent sample.
- A rigid-rotation signal on this path needs lateral structure in the sampled inner core. "none"
  means no change at this precision on this path. It does not by itself refute rotation.

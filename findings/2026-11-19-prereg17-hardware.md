# Pre-registration 17: is the time trend among source-verified YKA pairs a hardware artefact?

By MV2, Munim, Inc. (made with AI agents). Written 2026-11-19, about 11:45 UTC, after the
result of pre-registration 13 (H13) and before any of the splits below were computed. Text
CC BY 4.0; code (`h17.py`, `test_h17.py`) Apache-2.0. No human expert has reviewed it.

## Why

H13 (findings/2026-11-19-prereg13-results.md) applied a mantle-P source control to every
published Yellowknife (YKA) repeater pair. Among the 328 pairs whose sources are verified alike
(S <= 1/3, 42 families), the PKIKP coda difference D grows with time apart: Spearman
rho_M = 0.158 (family-bootstrap 5th-95th percentile 0.061-0.255) and, on the backtracking clock,
rho_B = 0.127 (0.033-0.201). The pre-registered clock test was inconclusive (rho_B > rho_M in
only 27% of bootstraps). Before source control, the same table showed no time dependence at
all (diag-v1: rho with time 0.01 over 1,650 pairs).

A positive rho with time apart is what a changing inner core predicts. It is also what an
instrument change predicts. The YKA pair table marks each record as **legacy** (older
waveform source; events up to at least mid-2013 in the table) or **new** (events from
September 2015 on). Pairs far apart in time are more likely
to straddle that change, and if the two systems record the same wave slightly differently,
cross-generation pairs would have larger D for a reason that has nothing to do with the core.
Nobody has looked at this split in our data: `h17.py` is the first code that reads gen1/gen2
next to the source score.

## Data (`h17.py` 0.1.0, compete job `icr-h17-v1`; no waveforms fetched)

* `results/ykab-v2/yka_pairs.csv` at commit f65c966a (the table H13 used), rows with ok = True,
  with gen1, gen2. Clocks recomputed from t1, t2 exactly as in mantlep3.py 0.15.0.
* `results/h13-v1/h13_pairs_scored.csv` (S and the number of stations per pair, as published).
* Statistics copied unchanged from mantlep3.py 0.15.0 (Spearman with average ranks, family
  bootstrap with 1,000 draws). **Reproduction check:** the script must reproduce H13's primary
  clock numbers (328 pairs, 42 families, 96 pairs across 2008.5, rho_B 0.12652, rho_M 0.15755,
  to 1e-6). If it does not, every result below is labelled not-reproduced and not graded.

Verified set V: as in H13 (YKA, M03-24 pairs excluded, scored at >= 2 stations, S <= 1/3).
Generation class: same-legacy, same-new, cross (one record of each), unknown (excluded).

## Tests and decision rules

**P17-1 (primary): a time signal within one hardware generation.** On V restricted to
same-generation pairs (same-legacy and same-new together): rho = Spearman(D, dM), family
bootstrap (seed 17), 5th and 95th percentiles.

* **untestable** if fewer than 40 pairs or fewer than 4 families;
* **time-signal** if rho >= 0.10 and the 5th percentile > 0;
* **no-time-signal** if rho < 0.10 and the 95th percentile < 0.20;
* **inconclusive** otherwise.

**P17-2 (secondary, the H13 clock test on same-generation pairs).** mantlep3's `clock()` and
pre-registered decision rule on the same set (needs >= 40 pairs, >= 4 families, >= 10 pairs
across 2008.5; only legacy-legacy pairs can straddle 2008.5): B / M / none / inconclusive /
untestable.

**Reported, not decisive.**
* The time signal separately for same-legacy, same-new and cross pairs.
* The cross-generation indicator's Spearman correlation with D on all of V, and median D of
  same- vs cross-generation pairs in bins of time apart (0-3, 3-6, 6-10, 10-15, 15-30 yr), so
  a hardware step can be seen at matched time separation.
* Everything above for the late coda D_c310.
* **Coda-beat locks (S17).** In H13, a family's stacked mantle-P onset at a station sometimes
  sits more than 0.3 s from the reference event's own onset (`mantlep3_families.csv`), a sign
  that the alignment locked onto a later wiggle. S is recomputed from `mantlep3_pairs.csv`
  (D_late as published, 5 decimals) with and without those family-stations; calibration,
  clock and P17-1 are reported for both. Exploratory.

## What each outcome would mean

* **time-signal:** the core coda of source-verified repeaters changes with time even when both
  records come from the same hardware. H13's time trend is not an instrument artefact. That
  is the first multi-family, source-controlled evidence for temporal change of the inner-core
  path at YKA from this data set. It would still not tell backtracking from steady change.
* **no-time-signal:** within one hardware generation the trend vanishes, so H13's trend most
  likely comes from the 2013-2015 hardware change, or from a core change that happened in
  that same gap; this data set cannot tell those apart. Either way, H13 would not be evidence
  for a time trend.
* **inconclusive / untestable:** too few same-generation verified pairs.

## Forecasts (written before the split was computed)

| Item | Forecast |
|---|---|
| Reproduction check passes | 90% |
| P17-1 testable (>= 40 same-generation verified pairs) | 85% |
| Given testable: time-signal / no-time-signal / inconclusive | 35% / 25% / 40% |
| Median D higher across generations than within, in most dM bins with both | 60% |
| P17-2: B / M / none / inconclusive / untestable | 5% / 10% / 20% / 35% / 30% |

## Limits written in advance

* The YKA hardware dates come from the table's gen labels, not from station metadata.
  A core change between 2012 and 2015 and the hardware change cannot be separated.
* Same-generation pairs are shorter in time (at most about 12 years), which lowers the power to
  see a slow trend.
* S is relative (percentiles at eight stations); "verified" means among the most alike third.
* The core-side table and H13's overall rho are already known; only the split is new.

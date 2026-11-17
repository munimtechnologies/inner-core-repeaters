# ILAR results (icr-ilar-v1): graded against the pre-registration

By MV2, Munim, Inc. (made with AI agents). Written 2026-11-17 after the first full run.
Text and tables CC BY 4.0; code Apache-2.0. Pre-registration:
[`2026-11-17-prereg.md`](2026-11-17-prereg.md). Raw outputs: [`../results/ilar-v1/`](../results/ilar-v1/).

## Short version

* The open, objective pipeline works end to end on 28 years of public ILAR data. It found
  **75 South Sandwich Islands (SSI) multiplets** and measured **874 pairs**, of which 136 passed
  the quality cuts.
* **P4 (calibration) holds.** Pairs less than 2 years apart are more alike in PKIKP than
  pairs that differ under both rotation models (median D 0.080 vs 0.105).
* **P1 (in-sample replication of backtracking) does not hold, narrowly.** The D ordering
  prefers Wang et al.'s backtracking clock over a monotone clock (ρ_B = 0.19, p = 0.032;
  ρ_M = 0.07, p = 0.47). But only **88%** of multiplet bootstraps favour it, and the
  pre-registered bar was 90%.
* **P2 and P3 (the out-of-sample tests on 2025-2026 earthquakes) are untestable in v1.**
  Of the 13 usable new pairs, only one crosses the 2008.5 turn, and no new event has two or
  more partners from before 2008.5. That was the most likely outcome we forecast. The
  decisive tests are the queued M03-24 patch test and the arc-wide scan, which go down to
  M4.3 where the 2001-2002 partners exist.
* **P5 (exploratory):** the turn date that best orders the in-sample pairs is **2011.75**
  (rate ratio 2.0), not 2008.5. It sits outside our pre-stated 2008 ± 2 window, so it does
  not count as an independent recovery of Wang et al.'s turn. No uncertainty is attached
  yet (see "Next").

This is a partial, honest result: the metric is calibrated and leans toward the
backtracking ordering in-sample, but not at the confidence we set in advance. The 2025-2026
test is still open.

## Run

* 2026-11-17 16:14-16:56 UTC (42 min), 2 CPU, 0.52 GB peak. All 19 unit tests passed.
* Events: USGS ComCat, region lat −61.5 to −54, lon −31 to −22. M ≥ 5.0 to mid-2024 and
  M ≥ 4.7 after. 1,559 events were made into beams.
* **ILAR coverage at EarthScope starts in 2000.** Beams made: 0 of 231 events in 1990-1997,
  4 of 31 in 1998, 0 of 16 in 1999, then 11/11 (2000), 20/20 (2001), 26/27 (2002) and
  close to 100% every year after.

## Grading (pre-registered rules, unchanged)

| Test | Result | Numbers |
|---|---|---|
| P4 calibration | **holds** | median D: 0.080 (n = 26, < 2 yr apart) vs 0.105 (n = 14, > 10 yr and > 5 yr equivalent) |
| P1 in-sample | **does not hold** (bootstrap below bar) | n = 123 pairs, 38 cross-turn, 28 multiplets; ρ_B = 0.194 (p = 0.032), ρ_M = 0.065 (p = 0.47); ρ_B > ρ_M in 88.0% of bootstraps (needed 90%) |
| P2 out of sample | untestable | 13 pairs, 1 cross-turn (needed 3) |
| P3 best match | untestable | 0 discriminating new events (needed 2) |
| P5 exploratory | best T = 2011.75, ratio 2.0 | ρ = 0.268 in-sample (ρ_M = 0.065); the same on all pairs |

Robustness versions of P1 (also pre-registered as checks, not as the decision):

| Variant | ρ_B (p) | ρ_M | bootstrap ρ_B > ρ_M |
|---|---|---|---|
| D (decision) | 0.194 (0.032) | 0.065 | 88.0% |
| Dx = D − (1 − cc_ref) | 0.240 (0.007) | 0.090 | 93.8% |
| Ds, element median | 0.199 (0.027) | 0.074 | 89.5% |
| pair threshold cc_ref ≥ 0.90 (175 pairs) | 0.165 (0.029) | 0.025 | 95.3% |
| pair threshold cc_ref ≥ 0.95 (75 pairs) | 0.194 (0.095) | 0.118 | 66.0% |

Two of four variants clear the 90% bar and two do not. We read this as a consistent
lean toward the backtracking ordering that the v1 sample cannot establish.

### Our forecasts, scored

Stated before any waveform was downloaded: P4 holds 75% (it held); P1 holds 45% (it did
not, narrowly); P2 testable 35% (it was not); P3 testable 35% (it was not). The forecast
"most likely v1 outcome" (a working catalogue and calibration, too few cross-turn
2025-2026 pairs) is what happened.

## The 13 usable pairs with a 2025-2026 event (descriptive only)

Δ_M = |t2 − t1| (years); Δ_B = |e(t2) − e(t1)| under Wang et al.'s backtracking clock.

| Multiplet | Earlier event | New event | D | Δ_M | Δ_B |
|---|---|---|---|---|---|
| M043 | 2003-03-21 | 2026-01-31 | 0.244 | 22.9 | 1.75 (cross-turn) |
| M047 | 2010-02-03 | 2026-11-12 | 0.135 | 16.8 | 6.71 |
| M028 | 2009-02-05 | 2025-02-20 | 0.165 | 16.0 | 6.42 |
| M006 | 2011-02-20 | 2025-12-13 | 0.020 | 14.8 | 5.92 |
| M008 | 2015-09-09 | 2025-12-12 | 0.072 | 10.3 | 4.10 |
| M008 | 2020-03-02 | 2025-12-12 | 0.117 | 5.8 | 2.31 |
| M006 | 2020-07-11 | 2025-12-13 | 0.073 | 5.4 | 2.17 |
| M037 | 2021-08-31 | 2026-05-25 | 0.043 | 4.7 | 1.89 |
| M015 | 2021-12-19 | 2026-06-20 | 0.175 | 4.5 | 1.80 |
| M037 | 2022-03-30 | 2026-05-25 | 0.049 | 4.2 | 1.66 |
| M005 | 2022-10-27 | 2026-07-11 | 0.027 | 3.7 | 1.48 |
| M071 | 2022-08-12 | 2025-03-26 | 0.150 | 2.6 | 1.05 |
| M006 | 2024-03-21 | 2025-12-13 | 0.043 | 1.7 | 0.69 |

Two cases stand out, but each is a single low-SNR pair (DF SNR 3.2-3.6) and proves nothing:

* **M043:** the only cross-turn new pair. Backtracking says a 2026 event should look like
  a 2001-2002 one, so its 2003 partner should be close (Δ_B = 1.75 yr). It is in fact
  quite dissimilar (D = 0.24, against a median of 0.08 for pairs < 2 yr apart).
* **M006** contains four of the six events of Zhang & Wen's (2026) multiplet M03-24
  (2011-02-20, 2020-07-11, 2022-03-07, 2024-03-21), re-found here by the objective
  pipeline at M ≥ 5.0. Its 2025-12-13 event matches the **2011** event best (D = 0.020),
  better than the 2024 (0.043) or 2020 (0.073) events. Both rotation clocks predict the
  2011 partner should match *worst*. This fits episodic, partly reversible change better
  than either smooth clock, but it is three pairs from one patch.

## Sanity checks

**Phase timing on this path.** The DF onset picked on each multiplet stack fell **before**
the ak135 PKIKP time in 72 of 75 multiplets: median −4.35 s, range −6.05 to −0.73 s. Three
picks were late (M031, M057, M064; probably not DF). Five sat at the −6 s edge of the
search window (M006, M020, M039, M054, M056), so their true onset may be earlier. The large
arrival in the plots (`multiplet_M001.png`, `multiplet_M006.png`) is PKPbc, about on its
ak135 time, roughly 11 s after the picked DF onset where ak135 predicts about 7 s. DF itself
is weak at 1-2 Hz, yet it repeats across events: usable pairs reach DF correlation above
0.95. Noise from different earthquakes cannot do that. This matches the SSI → Alaska path
being the most anomalous PKP path known: DF arrives 1-4 s early, the scatter in BC−DF is
about 6 s, and Frost & Romanowicz (2019, PEPI) attribute about 2 s to the slab beneath
Alaska. Consequences:

* The 2.8 s DF window does contain PKIKP (and, at edge picks, early PKIKP coda).
* The reference window starts at the *predicted* BC time. On this path that is about 4 s
  before the real BC, so it also contains about 4 s of weak late-DF coda. BC and AB
  dominate the energy, so alignment is driven by the outer-core phases. Still, Dx (which
  subtracts 1 − cc_ref) removes a little inner-core signal too.
* **The planned BC−DF job (`bcdf.py` v0.4.0) would have failed on this path.** It searched
  BC at the DF onset plus the ak135 BC−DF time (window −0.6 to +1.5 s), about 4 s too early.
  It was fixed in v0.5.0 before running, and so was `cola.py` (College, Alaska, 40 km from
  ILAR on the same path, which calls the same measurement); see
  [`2026-11-17-bcdf-amendment.md`](2026-11-17-bcdf-amendment.md).
* For a v2 there is room for a much longer DF-plus-coda window (DF to BC is about 11 s
  here, not 2.8 s). That is closer to what Wang et al. graded by eye. It needs its own
  pre-registration.

**Instrument epochs.** The ILAR SHZ response metadata change on 2018-09-25 (new digitiser
chain; the sensor is unchanged and the gain changes by 1.9%; the FIR filters are
linear-phase with corrected delay). D is a normalised correlation, so the gain does not
matter. Whether cross-epoch pairs carry extra decorrelation is checked by `ilar_diag.py`
(below). Note that P4 "near" pairs are often both post-2018 (the 2021 sequence), and "far"
pairs often span the break.

## Next

* `ilar_diag.py` (exploratory, written after these results were read): instrument-epoch
  split, SNR partial correlations, leave-one-multiplet-out, the P1/P4 tests without the
  2021 sequence, a within-multiplet sign test, and a multiplet bootstrap of the P5 turn
  date. It shows whether the 2011.75 optimum is sharp or just one point on a flat ridge.
* `bcdf.py` v0.5.0 (single-event BC−DF on the same cached beams, 1998-2026).
* The decisive out-of-sample tests are already queued: the M03-24 patch test (M ≥ 4.3,
  with partners on 2001-01-29 and 2002-06-29) and the arc-wide forecast-date scan.

## References

* Wang, W., Vidale, J. E. et al. 2024, *Nature* 631, 340. doi:10.1038/s41586-024-07536-4
* Vidale, J. E. et al. 2025, *Nature Geoscience* 18, 267. doi:10.1038/s41561-025-01642-2
* Zhang, Y. & Wen, L. 2026, *JGR Solid Earth* 131, e2026JB034595. doi:10.1029/2026JB034595
* Wu, Song & Yang 2026, *GRL*, doi:10.1029/2025GL119295
* Tkalčić, H. et al. 2010, *GRL* (large travel-time variations on SSI paths); Leykam et al.
  2010, *GJI* 180, 1329 (SSI-Alaska PKP scatter compared with Antarctic paths)
* Frost, D. A. & Romanowicz, B. 2019, *PEPI* (Alaska slab contribution to SSI-Alaska PKP anomalies)

# Results (exploratory): why YKA's calibration failed, and what the ILAR BC−DF "branches" are

MV2, Munim, Inc. (made with AI agents). Written 2026-11-18. Text CC BY 4.0, code Apache-2.0.
These are the results of `icr-diag-v1` (`icrdiag.py` 0.1.0), read against the
[diagnostic plan](2026-11-18-diag2-plan.md). That plan was published before the job ran. Tables and plots are in
[`results/diag-v1/`](../results/diag-v1/).

**Status: exploratory.** Nothing here changes a pre-registered verdict. YKA stage B is
still uninterpretable (P4-Y failed), and the pre-registered ILAR BC−DF result (+0.0014 s/yr,
not decisive) stands.

## Short version

* **YKA: what limits the coda measurement is the source, not time.** Across 1,650 YKA pairs,
  coda dissimilarity D (+1 to +6 s after PKIKP) has **no detectable dependence on the time
  between the two events**. The rank-regression coefficient is −0.0004 (5-95%: −0.026 to
  +0.027; Spearman ρ = 0.013). D rises with catalogue epicentral separation (coefficient
  +0.26, 5-95%: +0.11 to +0.42) and with magnitude difference (+0.16, +0.10 to +0.20). It
  falls with signal-to-noise (−0.16, −0.30 to −0.08). Of the four explanations in the
  plan, only **location** meets its pre-stated criterion.
* So families defined by ILAR's outer-core phases are not tight enough for YKA's inner-core
  coda. On the PKIKP coda path to YKA, the source differences between "repeaters" that are
  ~15 km apart (median) are larger than any change over time. Under the plan, this leads to a
  **prospective** YKA test with a selection rule built from location
  ([pre-registration 12](2026-11-18-yka-prospective-prereg.md)). That test can only use
  earthquakes after this note.
* **ILAR BC−DF: the "branches" are failed onset picks, not geography or cycle skips.** The
  two-component fit found a narrow core (σ 0.57 s, 70% of events) plus a wide tail (σ 2.1 s,
  30%). It did not find two branches 0.6 s apart. The tail is mostly DF onsets picked in the
  wrong place: 20 of the 22 events with residual below −2.5 s had DF picked late (after
  −2.1 s, where the core sits at −3 to −5.5 s). 11 of the 20 events above +2.5 s are
  intermediate-depth (100-131 km) events from one northern-arc cluster, picked 5-8 s early.
  The tail's share is flat in time (27%, 31%, 30% and 32% across the four periods), so it adds
  noise but no trend.
* **The ILAR single-event BC−DF lane is closed.** Even the core has σ ≈ 0.57 s. Its
  pre-2025 slope is known only to ±0.004 s/yr. Telling the clocks apart on the 2025-26 events
  would need a precision of a few hundredths of a second, more than 10 times better than this
  picking gives. The plan's bar for a new version (σ below about 0.2 s) is not within reach.

## Part Y: YKA stage-B pairs (`results/ykab-v2/yka_pairs.csv` + USGS ComCat)

2,971 pairs were located and 1,650 are ok. Near = less than 2 years apart (321 pairs, 21
families). Far = more than 10 years apart and more than 5 years apart on the backtracking
clock (92 pairs, 15 families).

| Explanation | Pre-stated criterion | Result | Met? |
|---|---|---|---|
| Simpson's paradox | within-family sign test p < 0.05 **and** family-demeaned P4 holds | D: near < far in 6 of 8 families, p = 0.14; demeaned P4 holds (bootstrap 75%) | **no** (half) |
| The 2021-22 sequence | P4 holds without it; most near pairs come from it | 58% of near pairs are from it, but without it P4 gets *worse* (0.117 vs 0.086) | **no** |
| Noise | min-SNR coefficient dominates; P4 holds at high SNR | coefficient −0.16 (excludes 0) but smaller than separation's; high-SNR half is a tie (0.0950 vs 0.0951) | **no** (partial) |
| Location | separation coefficient > 0 with 5-95% excluding 0 | +0.26 (+0.11 to +0.42) | **yes** |

Also reported:

* **The 1-3 s window (D_c13).** Near pairs are more alike than far pairs in **all 8
  families** (sign test p = 0.004), and family-demeaned P4 holds in 95% of bootstraps. So the
  calibration that fails in the pre-registered +1 to +6 s window holds in the window nearest
  the direct pulse. D_c13 was not pre-registered for P4 and is not promoted.
* **The close half alone does not fix P4 for D** (0.072 vs 0.066). Only the high-SNR, close
  quarter has near < far (0.065 vs 0.091, bootstrap 69%). This stratum was picked after
  looking, which is why the next test is prospective.
* **Same-generation pairs** (both events recorded by the same YKA hardware): near < far
  (0.097 vs 0.146), but only 10 far pairs in 3 families.

Plot: [`ykadiag_D.png`](../results/diag-v1/plots/ykadiag_D.png) (D against time separation
and against SNR, with the 2021-22 sequence marked).

## Part B: ILAR BC−DF residuals (`results/bcdf-v2/bcdf_events.csv`, 477 quality-passed events)

| Check | Result |
|---|---|
| Mixture | μ −0.15 / +0.16 s, σ 0.57 / 2.07 s, weights 0.59 / 0.41 (posterior labels: 334 core, 143 tail); ΔBIC −109 against one Gaussian |
| Neighbour agreement | 1 neighbour: p = 0.10; 5 neighbours: p = 0.002 |
| Tail vs core | more southern (median lat −57.7 vs −56.6, p = 0.003); lower correlation (0.90 vs 0.93, p = 10⁻⁷); lower DF SNR (4.4 vs 5.0, p = 0.01); DF onset median −4.7 vs −4.2 s |
| Share of the tail in time | pre-2008.5 27%, 2008.5-2021.6 31%, 2021.6-2025 30%, 2025+ 32% |
| Core slope | pre-2025 +0.0092 ± 0.0039 s/yr (n 304); all +0.0049 ± 0.0037 s/yr (n 334) |

**Reading.** By the plan's table, the neighbour agreement and the latitude difference point
to "geography". The per-event table shows what lies under that: **DF pick failures**, which
cluster in space because some source regions give weaker or more complex DF arrivals.
Residuals below −2.5 s go with DF onsets picked at −2.1 s or later, up to +0.7 s, where the
true early DF should be near −4 s. Residuals above +2.5 s go with onsets at −5 to −8 s. More
than half of them are the 100-131 km deep events near 56.0-56.3°S, 27.5-28.3°W. A real BC−DF
anomaly of +3 s that depends on source depth is not plausible, because the two phases share
the source. The cycle-skip question (two modes 0.6 s apart *inside* the core) was not tested
by this fit, which modelled the tail instead. With a core σ of 0.57 s it would not change
the conclusion below.

**Not interpreted:** the core's median for 2025-26 events sits 0.30 s below the pre-2025
median. Here the spatial terms were fitted on all events, outliers included. The
pre-registered v2 analysis used a ±1 s screen and fitted space and time jointly, and there the
constant model predicted the 2025-26 events as well as either clock. The core slope
(+0.009 s/yr) and the pre-registered slope (+0.0014 s/yr) differ by about 2σ, so the trend
depends on how outliers are handled. Neither is evidence for a change.

Plot: [`bcdf_branches.png`](../results/diag-v1/plots/bcdf_branches.png).

## What this means for the project

1. YKA's coda cannot carry a retrospective test at M4.3-5.3: catalogue-level "repeaters"
   differ at the source by more than the inner core changes. This does **not** apply to the
   M03-24 family, where both arrays show the same split. There the queued mantle-P control
   (`icr-mantlep2-v1`, [pre-registration 11](2026-11-18-mantlep2-prereg.md)) separates source
   from path directly.
2. Single-event BC−DF at ILAR has reached its precision limit with this picking. The lane is
   closed.
3. The YKA lane now continues only prospectively ([pre-registration 12](2026-11-18-yka-prospective-prereg.md)).

## Files (`results/diag-v1/`)

`icrdiag_summary.json` (every number above), `ykadiag_pairs.csv` (pairs with ComCat
separation, depth difference, min SNR, Δmag, 2021-sequence flag), `bcdf_branches.csv`
(events with spatially corrected residual, onset, correlation, SNR, branch), `plots/`. The job
ran on 2026-11-18 (233 s, all steps exit 0, 13 unit tests pass).

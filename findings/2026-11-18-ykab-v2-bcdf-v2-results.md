# Results: YKA stage B v2 (pre-registration 8) and ILAR BC−DF times (pre-registration 3)

MV2, Munim, Inc. (made with AI agents). Written 2026-11-18. Text CC BY 4.0, code Apache-2.0.
Graded against [pre-registration 8](2026-11-17-yka-stageb-prereg.md), the
[pick amendment](2026-11-18-yka-pick-amendment.md), [pre-registration 9](2026-11-17-yka-hardware-prereg.md)
(H4), [pre-registration 3](2026-11-17-bcdf-prereg.md) and its [amendment](2026-11-17-bcdf-amendment.md).
All of these were written before the jobs ran.

## Short version

* **YKA stage B cannot be interpreted.** Its own calibration check (P4-Y) failed. Pairs recorded
  less than 2 years apart are *less* alike at YKA than pairs recorded far apart (median D 0.098 vs
  0.084). Under the pre-registered rule, nothing at YKA is interpreted.
* For the first time, the out-of-sample vote on 2025-2026 earthquakes (F1-Y) **could be run**: 6
  discriminating cases in 4 families. It is **inconclusive**: 2 cases match the backtracking pick,
  1 matches the monotone pick and 3 match neither (binomial p = 0.50). Even if P4 had held, this
  would not confirm or break the 2024 backtracking model.
* **The YKA hardware control is complete.** H4 holds: same-generation near pairs (median D 0.097)
  sit 8.6 times above the hardware-only floor (0.011). The 2013-14 change of electronics is not
  what limits YKA.
* **ILAR BC−DF is not decisive, as forecast.** Over 280 events (2000-2024) there is no measurable
  time trend in the BC−DF residual. The slope is +0.0014 s/yr, and the monotone, backtracking and
  constant models fit equally well (residual sums of squares 60.28, 60.31 and 60.31). The 24 new
  2025-2026 events could not separate the models: the predicted gap is 0.012 s, against the 0.19 s
  the rule required.
* **What the null says.** Using open data and code, out of sample, on 2025-2026 earthquakes, our
  objective pipeline finds **no detectable signal either way** at YKA or in ILAR BC−DF times. That
  is consistent with the 2024 model and does not support it. The YKA failure has a likely cause
  that can be tested (below).

## 1. YKA stage B v2 (`icr-ykab-v2`, `yka_b.py` 0.10.1 `--repick`)

78 merged families, 410 members (35 families with a 2025+ member). 2,971 YKA pairs, 1,650 ok,
187 out of sample. Run on 2026-11-18, 10.5 min, peak 0.2 GB. All steps exited 0.

### Pick validity (Rule A, for v2 itself)

The constrained picker fixed the v1 failure. 62 of 75 picked families are within ±2 s of ak135,
none is before −4 s, and 13 are between −2.0 and −2.55 s. That is 17% invalid, below the 25%
limit, so Rule A does not trigger again. Most picks are 0.7-2.0 s before ak135. PKIKP is
early on this quasi-polar path, as expected. The 13 slightly-over picks come from the picker's
backward refinement, not from an edge pile-up. On the plots
([Y001](../results/ykab-v2/plots/yka_multiplet_Y001.png), [Y000](../results/ykab-v2/plots/yka_multiplet_Y000.png)),
the alignment window now sits on the first swing of PKIKP. The largest PKIKP peak is about 1.2 s
later, inside the measured +1 to +6 s window.

### Graded tests (in the pre-registered order)

| Test | Rule | Result | Forecast | Verdict |
|---|---|---|---|---|
| **P4-Y** calibration | median D near (< 2 yr) < far (> 10 yr and > 5 yr on the backtracking clock) | 0.0983 (321 pairs) vs 0.0838 (92) | holds 70% | **fails → nothing at YKA is interpreted** |
| F1-Y new-event vote | ≥ 6 discriminating cases in ≥ 3 families; supported needs n_B > n_M, p < 0.05, bootstrap ≥ 95% | 6 cases, 4 families: B 2, M 1, neither 3; p_B 0.50; bootstrap B ahead 63% | testable 35% | **testable, inconclusive** |
| F1x-Y (excess D) | must not point the other way | B 2, M 2, neither 2 | | inconclusive (F1-Y not robust) |
| P2-Y out of sample | rank test, ≥ 6 pairs, ≥ 3 cross-turn | 187 pairs, 29 cross-turn; ρ_B 0.060 (p 0.41), ρ_M 0.029; B better in 85% of bootstraps | testable 55% | testable, **does not hold** |
| P1-Y in sample | | ρ_B −0.016 vs ρ_M 0.017 (1,463 pairs) | holds 45% | does not hold |
| P6 ILAR vs YKA D | ρ > 0, p < 0.05 | ρ 0.446, p 9×10⁻⁸, n 132 | 50% | **holds** |
| H4 (pre-reg 9) | D_hw below same-generation near-pair D | 0.011 vs 0.097 (318 pairs, 21 families), ratio 8.6 | 75% | **holds** |

The F1-Y cases (best partner by D in the +1 to +6 s window):

| New event | Partners | B pick | M pick | Best match | Matches |
|---|---|---|---|---|---|
| us6000rii9 (2025.81) | 2 | 2002.49 | 2006.96 | 2006.96 | M |
| us6000rk5n (2025.83) | 5 | 2001.21 | 2007.22 | 2001.21 | B |
| us6000runq (2025.95) | 2 | 2003.22 | 2008.13 | 2003.22 | B |
| us6000t7ng (2026.48) | 4 | 2001.08 | 2006.26 | 1998.27 | neither |
| us6000thwz (2026.59) | 4 | 2001.08 | 2006.26 | 2003.40 | neither |
| us6000tmzw (2026.64) | 5 | 2001.21 | 2007.22 | 2005.44 | neither |

Sensitivity windows (reported, not graded): F1-Y in D_c13 is B 3 / M 1; in D_c310 B 2 / M 2; in
D_full B 3 / M 0 (p 0.125). None is significant. The same counts hold when only families with a
valid pick are used (`ykab_grade.json`, `v2_valid_only`). In D_c13, P4 passes narrowly
(0.044 < 0.047). That was not the pre-registered window and is not promoted.

The v1 tables were not on the public repo when the v2 grader ran (that publish failed a file
check and is re-filed with this release). So the v1 variants of H4 were not computed. They are
superseded by v2 and will not be chased.

### The 2025-12-13 event at YKA (reported, not interpreted)

At ILAR, the 2025-12-13 event (us6000runq) looked like the 2011 event rather than the 2003 one.
At YKA it is the other way round: D to 2003 = 0.004, to 2008 = 0.014, to 2011 = 0.013. Its YKA
D is lowest to the 2003, 2020, 2022-03 and 2024 events (0.0001-0.0055) and higher to 2008, 2011
and 2014 (0.006-0.014). That is the same X/Y split of the M03-24 family that the late coda
showed on both arrays ([m0324chk-v2](2026-11-18-m0324chk-v2-results.md)). Whether the split is
inner-core change or two source sub-clusters is what the queued mantle-P control
(`icr-mantlep-v1`, [pre-registration 10](2026-11-18-mantlep-prereg.md)) decides.

### Why P4 might fail, and what we do about it

The 2021 M8.1 South Sandwich sequence supplies most near-in-time pairs: 110 of the YKA events
are from 2021-22. Its aftershocks spread along a long rupture, so pooled near pairs may mix
true repeaters with merely nearby events. Pooled medians can then reverse even when every
family behaves (Simpson's paradox). P6 adds to the suspicion: ILAR and YKA dissimilarity are
strongly correlated pair by pair (ρ 0.45). Their inner-core paths are different, so what they
share is mostly the source side. An exploratory diagnostic was written down before it ran
([plan](2026-11-18-diag2-plan.md), `icrdiag.py`, job `icr-diag-v1`). Any new confirmatory YKA
test will be prospective, using only earthquakes after 2026-11-18, because every 2025-2026
pair has now been seen.

## 2. ILAR single-event BC−DF (`icr-bcdf-v2`, `bcdf.py` 0.5.0)

ILAR beams were rebuilt from the pinned v0.8.0 code (1,559 cached events), and BC−DF was
measured on 1,269. 477 passed the quality cuts and 304 the ±1 s outlier screen: 280 before 2025
and 24 from 2025-2026. DF arrives a median 4.2 s before ak135. The BC−DF residual is +3.5 s
(10-90%: 3.0-4.4 s), the same before and after the 2018-09-25 ILAR electronics change.

| Question | Rule | Result | Forecast | Verdict |
|---|---|---|---|---|
| Q1 in sample | lower RSS of L vs Bk, sign of b | RSS L 60.283, Bk 60.309, constant 60.311; b_L = +0.0014 s/yr | testable 70% | testable; **L marginally lower, no meaningful difference from a constant** |
| Q2 out of sample | decisive if ≥ 8 new events and \|Bk − L\| ≥ 2σ/√n | 24 new events; gap 0.012 s vs required 0.191 s | decisive 25% | **not decisive** |
| Q3 exploratory turn scan | best T, ratio | T 2017.0, ratio 1.0; RSS 59.98 (0.5% lower) | | no echo of 2008.5 |

Out of sample, the constant model predicts the 24 new events as well as either clock (MSE 0.186
vs 0.188 for L and 0.186 for Bk). The element-median measurement gives the same verdicts.

**Phase check and a measurement caveat.** On the example beams
([`bcdf_examples.png`](../results/bcdf-v2/bcdf_examples.png)), the DF pick falls on a small
emergent arrival about 4 s before ak135, and BC falls on the strong arrival near its ak135 time.
After the spatial terms, the residuals form **two branches about 0.6 s apart**
([`bcdf_vs_time.png`](../results/bcdf-v2/bcdf_vs_time.png)). That is about one period at 1-2 Hz,
so some BC picks may be one cycle off. The ±1 s screen cannot remove that, so it inflates the
scatter (σ = 0.47 s). `icrdiag.py` tests whether the branches are cycle skips or geography.

**What it constrains.** A slope of +0.0014 s/yr over 2000-2024, with the scatter above, means
that, by our rough estimate, ILAR BC−DF times rule out a steady change much larger than
about 0.01 s/yr since 2000. The diagnostic computes the exact standard error. That is the scale
reported for this path in the 1990s, so ILAR alone cannot say whether the change slowed. It
is consistent with Wu, Song & Yang 2026's slow post-2000 change at COLA. With the present
picking, the series cannot resolve backtracking: that would need a gap of a few hundredths of a
second.

## Forecast scorecard

Pre-registration 8: P4 holds (70%) ✗; P2 testable (55%) ✓; F1 testable (35%) ✓; inconclusive
if testable (50%) ✓; P1 holds (45%) ✗; P6 positive (50%) ✓. Pre-registration 9: H4 (75%) ✓.
Pre-registration 3: Q1 testable (70%) ✓; Q2 not decisive (75%) ✓.

## Files

* [`results/ykab-v2/`](../results/ykab-v2/): `yka_b_summary.json`, `ykab_grade.json`,
  `yka_pick_diagnostics.csv`, `yka_pairs.csv` (every pair), `yka_b_cases.csv`,
  `yka_multiplets.csv`, `merged_multiplets.csv`, `yka_event_status.csv`, `yka_statics.csv`,
  `yka_predictions.json`, `yka_report.md`, plots.
* [`results/bcdf-v2/`](../results/bcdf-v2/): `bcdf_summary.json`, `bcdf_events.csv`, plots.

# Mantle-P source control of the M03-24 split: results (pre-registration 11)

By MV2, Munim, Inc. (made with AI agents). Written 2026-11-18 after `icr-mantlep2-v1` finished
(17:57-18:01 UTC, 4.6 min, all steps exit 0, 17 unit tests passed). Text CC BY 4.0.
Tables and plots: [`results/mantlep2-v1/`](../results/mantlep2-v1/). Pre-registration:
[`2026-11-18-mantlep2-prereg.md`](2026-11-18-mantlep2-prereg.md). Code: [`mantlep2.py`](../mantlep2.py) 0.13.0.

## The question

Eight earthquakes of one South Sandwich Islands repeating source, M03-24 (59.03°S, 25.68°W),
split into two groups by the late (3-10 s) coda of their waves through the inner core. The
split is the same at two arrays, ILAR in Alaska and YKA in Canada:

* X = {2003-03-21, 2020-07-11, 2022-03-07, 2024-03-21, 2025-12-13}
* Y = {2008-02-17, 2011-02-20, 2014-12-09}

So the inner-core waves changed after 2003 and changed back by 2020. This is the pattern the
2024 "backtracking" clock predicts (Wang, Vidale et al., *Nature* 631:340). Such a pattern could
also come from the earthquakes themselves, if the 2008-2014 events sat slightly apart from the
others (a source sub-cluster). Mantle P never enters the core. **S:** the split shows in mantle
P too, so it is a source effect. **N:** it does not, so it arises on the path through the core.

## Pre-registered result: **N**

| Item | Forecast (before the run) | Result |
|---|---|---|
| Stations found by the station service | – | 42 (IU, II, GT, 30-85°); station_source = station-service |
| At least 3 usable stations in the primary set | 55% | **Yes: 7** |
| Decision | N 60%, S 15%, inconclusive 25% | **N: 0 hits among 7** (rule: h ≤ max(1, ⌊0.2n⌋) = 1) |
| All usable stations (secondary; includes QSPA) | – | N: 0 of 8 |
| New stations only (sensitivity) | – | N: 0 of 7 |

Primary stations, D_late (3-10 s after the P onset). A hit needed p ≤ 0.07.

| Station | Distance | Gated events (X+Y) | Rank of the X/Y split | p | Contrast ratio* |
|---|---|---|---|---|---|
| GT.PLCA (Argentina) | 33.3° | 8 (5+3) | 34 of 56 | 0.61 | 1.20 |
| IU.LCO (Chile) | 42.8° | 7 (5+2) | 6 of 21 | 0.29 | 1.89 |
| GT.VNDA (Antarctica) | 43.4° | 8 (5+3) | 38 of 56 | 0.68 | 0.88 |
| GT.LBTB (Botswana) | 49.1° | 7 (5+2) | 3 of 21 | 0.14 | 1.15 |
| GT.LPAZ (Bolivia) | 52.8° | 7 (5+2) | 7 of 21 | 0.33 | 1.19 |
| IU.LSZ (Zambia) | 58.8° | 7 (4+3) | 27 of 35 | 0.77 | 1.11 |
| GT.DBIC (Côte d'Ivoire) | 67.8° | 8 (5+3) | 11 of 56 | 0.20 | 1.29 |
| *IU.QSPA (South Pole, seen in v1; secondary)* | 31.0° | 8 (5+3) | 23 of 56 | 0.41 | 1.12 |

\*Median D of X-Y pairs over median D of within-group pairs. S predicts a ratio well above 1;
N predicts about 1.

The early window (1-3 s) gives no hit either (best p = 0.095 at LBTB).

## The same pairs through the core

Computed on all 28 pairs of the eight events from the published core table (`z004_core_D.csv`):

| Path | Median D, X-Y pairs (15) | Median D, within groups (13) | Ratio | Rank of the X/Y split |
|---|---|---|---|---|
| YKA PKIKP coda, 3-10 s | 0.132 | 0.0165 | **8.0** | 1 of 56 |
| ILAR PKIKP, long window | 0.211 | 0.063 | **3.4** | 1 of 56 |
| Mantle P, 7 primary stations | – | – | **median 1.19** (0.88-1.89) | never better than 3 of 21 |

The Spearman correlation between mantle D_late and YKA D_c310 on the same pairs has a median of
0.15 over the primary stations (S predicts clearly positive values).

## Checks made before reading the decision

* **Alignment.** Every usable station's plot was inspected. The events are aligned on the
  whole P train, and the stack onset is not at the search edge: at all usable stations the
  stack onset lies within 0.05 s of the reference onset. At VNDA, QSPA, DBIC and LCO the eight
  records are nearly identical wiggle for wiggle out to 10 s.
* **Stations that failed.** Of the 42 stations, 8 were usable, 10 kept some events but too
  few, 19 lost every event at the gates and 5 had 0-1 records (SPA, BDFB, PTGA, MSKU, KOWA).
  At most of the 19, the reference onset landed at or near the first sample of the search
  window (-5 s), and the plots (for example `II.SHEL`) show noise with no clear P at 1-2 Hz.
  These are M4.8-5.3 events. The gates are per event and blind to the groups; they reject
  noisy stations as intended and do not decide which way the test goes.
* **Power.** The late windows of both paths contain the depth phases (pP, sP). A depth or
  location difference between the sub-groups large enough to produce YKA's 8-fold contrast would
  change these phases at mantle stations too. At the quiet stations (QSPA, VNDA) within-group
  D is only 0.011-0.017, so a source difference of even a tenth of the core-side contrast would
  stand out. The split ranks 1 of 56 on both core paths and never better than 3 of 21 on any
  mantle path.

## What this means

By the pre-registered rule, the M03-24 split is **not a property of the source**. The same
repeating earthquake sent nearly identical mantle P waves to seven independent stations on
three continents and to the South Pole in 2003, 2008, 2011, 2014, 2020, 2022, 2024 and 2025.
Over the same years, the waves it sent through the inner core to Alaska and Canada changed
after 2003 and changed back by 2020, and the 2024 and 2025 events stay on the "back" side.
Two arrays, whose paths share only the source region and the inner core, see the same split.

This removes the main alternative explanation, that the 2008-2014 events were a separate
source. It also covers the 2025-12-13 event, which no published study of this multiplet has
used.

**What it does not say:**

* **Rotation or a changing surface.** Both an inner core that turned back (Wang et al. 2024;
  Vidale et al. 2025) and an inner-core surface that changed and recovered (Zhang & Wen 2026,
  who use this same multiplet through 2024) predict a core-path change. This test cannot
  separate them.
* **Structure straight under the source.** Mantle P leaves the source at about 15-40° from
  vertical; PKIKP leaves almost straight down. A time-varying structure confined to a narrow
  cone beneath the source would mimic N. It would need to change and change back over 20
  years in the uppermost mantle, which seems less likely than a change at depth, but it is not
  excluded here.
* **The receiver side.** ILAR and YKA would both need to change in the same way at the same
  times. YKA's 2013-14 hardware change was measured (D ≈ 0.01, pre-registration 9) and the
  split cuts across it.
* **One multiplet.** This is one source and eight events. Whether other South Sandwich
  families show the same thing, with split dates on the same clock, is the next test.

## Where the project stands (synthesis)

| Lane | Result |
|---|---|
| ILAR, 75 multiplets, objective D (v0.7) | Calibration holds; weak lean towards backtracking (ρ_B 0.19 vs ρ_M 0.07; 88% of bootstraps, 90% needed). Depends on the 2021 sequence. Turn ~2012 (2008.5-2018.5). |
| YKA, merged families (v0.12-0.13.1) | Calibration fails: catalogue "repeaters" differ at the source more than the core changes. Prospective test only (pre-registration 12, looks in 2027 and 2028). |
| BC−DF times, ILAR and COLA | Pick scatter (σ ≈ 0.5 s) is ten times too large; closed. |
| **M03-24, two arrays + mantle-P control (v0.14)** | **Core-path change confirmed against the source alternative, with 2024-2025 events following the "back" side (pre-registered, N).** |

Across many families the signal is weak and limited by how alike the sources really are. In
the one family where the source is shown to be a near-perfect repeater, the inner-core change
and change back is clear at two arrays.

## Next

* Apply the same two-step design (late-coda split at both arrays, then a mantle-P source
  control on unseen stations) to other South Sandwich families with at least five members from
  2003 to 2025, and ask whether their split dates follow one clock. That is what separates
  rotation, which moves all paths together on one time-line, from local surface change.
* The prospective YKA test (pre-registration 12) runs unchanged.

# Pre-registration 34: results of the M03-24 mirror test (untestable: no new members exist)

MV2, Munim, Inc. - graded 2026-11-27. Text CC BY 4.0; code Apache-2.0. Made with AI agents working in
MV2's name; no human expert has checked it.

Graded exactly as written in `findings/2026-11-24-prereg34-m0324-mirror-plan.md` and amendments 1
(`...-amendment1-clock-families.md`) and 1b (`...-amendment1b-wsy-rates.md`). All three were public
before either job produced a result (`prereg_public_at_start: true` in the job output).

- Job `m34-mirror-v1` (m34.py 0.1.0; unit tests passed as its first step): Dell, 2026-11-26 17:29-17:48Z.
- Job `clocks34-v1` (clocks34.py 0.1.0; unit tests passed): Dell, 2026-11-26 17:22-17:28Z.
- Files: `results/m34-v1/`, `results/clocks34-v1/`, hand-copied grade `results/p34-grade/p34_grade.json`.

## Headline

**None of the 31 candidate earthquakes is a member of the M03-24 repeating multiplet.** USGS ComCat lists
31 events of mb >= 4.6 within 25 km of the multiplet between March 2003 and February 2008 and between
December 2014 and July 2020, the two windows where a member would test whether the multiplet's states
retrace in time order, as a single-turn rotation requires. Every candidate that could be measured is
clearly *not* a repeat: mantle P says "different" for all 22 that reached three stations, and YKA's late
PKIKP coda puts all 24 measurable candidates far outside both known states. So **P34 is untestable**. The
contrast between one-turn rotation (Wang, Vidale et al. 2024) and episodic surface change (Zhang & Wen 2026)
cannot be decided from this multiplet. The core's state in 2003-2008 and 2015-2020 on this path rests on
interpolation between members, in our work and in every published study. The gap is in nature (no repeats
happened, among catalogued events), not in what earlier studies chose to use.

## Results by item

| Item | Pre-registered rule | Result |
|---|---|---|
| R34 reproduction | old-nine YKA pairs match published `m0324chk-v2` to < 0.001 | **reproduced** (36 pairs, max abs diff 0.0 in D and D_c310) |
| Membership (mantle P) | same if median q <= 0.75, different if >= 0.999, >= 3 stations | 0 same, 0 inconclusive, **22 different**, 9 untestable (0-2 stations gated) |
| YKA state | X / Y / intermediate / new-state / ambiguous | **24 new-state**, 7 untestable (fewer than 3 X / 2 Y ok partners), 0 X, 0 Y, 0 intermediate |
| R1 (Wang clock) | needs >= 1 member used | **untestable** |
| R2 (2,346 one-turn clocks) | needs >= 1 member used | **untestable** (old eight alone: 1,742 compatible) |
| Brackets | reported | unchanged: outbound 2003.2177-2008.1306, return 2014.9392-2020.5267 (both clean) |
| Amendment 1: R3, R4, R5 | on old eight plus usable new members | **untestable** (no usable new member) |
| Amendment 1b: R6a, R6b | same | **untestable** |

Every decided candidate has median q = 1.0: at each counted station its late mantle-P coda is farther from
the old members than every old member is from the others. At YKA the smallest distance from any candidate
to the nearer state is D_c310 = 0.58 (2020-01-08, us70006w5y), against a new-state threshold of 0.098 and
an X-Y separation near 0.1; the largest is 1.34 (2006-09-29). The two lines of evidence are independent
(mantle-only path vs core path) and agree. A candidate that repeated M03-24 while the core changed would
sit near the X-Y scale (~0.1) at YKA and stay "same" in mantle P. None does.

`clocks34-v2` (planned to add new members) was **not filed**: with no usable member its input would equal
the old eight and its output would equal `clocks34-v1`.

## Forecasts

| Forecast (before data) | p | Happened |
|---|---|---|
| R34 reproduces | 85% | yes |
| >= 1 candidate confirmed a member | 40% | no |
| >= 1 member with a usable state | 30% | no |
| R1 testable | 30% | no |
| Amendment 1 / 1b forecasts | - | void (all were conditional on a usable new member) |

## clocks34-v1 on the old eight (exploratory; the labels were known when the readings were written)

Every value derived by hand in amendments 1 and 1b was reproduced by the code:

- **Drift-free oscillation:** 8,269 of 9,801 periods (2.00-100.00 yr) can order the eight states. Passing
  intervals: 2.07-2.11, 3.27-3.51, 3.71-3.90 and **17.82-100 yr** (hand: 17.819). The 5.5-6.5 yr band
  (Wang & Vidale 2022) passes at 0% of periods; the 20-30 yr band at 100%. By chance 14.3% of phase
  arrangements pass, so the short isolated windows mean little.
- **Equal speeds both ways:** turn window 2009.079-2014.328 (hand 2009.078-2014.329; the code's 0.001-yr
  grid sits inside). A symmetric turn at 2009.0 fails, 2009.5 passes. Without 2014-12-09 the lower edge
  drops to 2007.179 (hand 2007.18); without 2020-07-11 the upper edge rises to 2015.155 (hand 2015.16).
- **Wide one-turn grid:** 11,361 of 39,349 clocks pass. The Wang et al. (2024) clock (turn 2008.5, ratio
  2.5) passes.
- **Amendment 1b check:** the earliest turn with a compatible ratio <= 2.0 is **2007.20** (ratio 1.9613;
  2007.15 needs 2.0594), exactly the hand value. A 2010 turn allows ratios 0.7392-5.4641 on the grid
  (hand range 0.73-5.63; the grid points nearest inside). So amendment 1b's exploratory reading stands as
  written: Wu, Song & Yang's (2026) rates fit these states only with a reversal starting around 2007-2010
  or later, with forward drift of at least about 7-22% of the 1980s rate continuing after 2003.

## What it means

- The mirror (backtracking-order) test needs repeats inside the leg windows. For M03-24, with catalogued
  events of mb >= 4.6, they do not exist. That is an honest untestable result, not a pass or fail for
  either reading.
- Being near a multiplet in a catalogue is a poor proxy for being a member: 0 of 22 decided candidates
  within 25 km repeat. Future designs should start from published repeater lists (Wang et al. 2024's 121
  events in 42 multiplets) or from waveform-detected members, not from catalogue distance.
- The one clear statement this multiplet still makes is the exploratory one above: under rotation alone
  with one boundary, the eight known states need a turn no earlier than about 2007 for a return no faster
  than twice the outbound speed, and a drift-free oscillation period of 17.8 yr or more.

## Limits

One multiplet; one array (YKA) for the state. Membership rests on mantle P at up to eight stations and
smaller events were often gated out (9 untestable). Events below mb 4.6 or missing from ComCat were not
searched. The X/Y labels of the old members came from exploratory work later confirmed as a core-path
effect (pre-registration 11). The 2014-12-09 Y label is the only one after the YKA 2013-14 hardware
change. The work was carried out by AI agents; no human expert has checked it.

## Credit

Wang, Vidale et al. (2024) for the multiplet and the one-turn clock; Zhang & Wen (2026) for the episodic
reading and open data; Wu, Song & Yang (2026) and Yang & Song (2023) for the multidecadal picture;
Wang & Vidale (2022) for the 6-year oscillation; USGS ComCat, EarthScope (IRIS DMC) and NRCan for the
catalogue and waveforms.

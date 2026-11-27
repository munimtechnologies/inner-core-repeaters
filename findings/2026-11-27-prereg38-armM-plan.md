# Pre-registration 38, arm M: does one rotation clock order every multiplet, at both arrays?

MV2, Munim, Inc. (mv2@munimtech.com), 2026-11-27. Text CC BY 4.0, code Apache-2.0.
Made by AI agents working in MV2's name; no human expert has checked it.

**Status: FROZEN (decision rules).** Decision code: `m2core/m2core.py` v0.1.0 with
`m2core/test_m2core.py` (expected values derived by hand, listed below). No M2 waveform has been
fetched or measured. The measurement front end (`m2measure.py`, which writes m2core's input file)
is written and frozen separately before its job is filed; it may not change anything in this
note. The job waits until zw33-v1, zw35g2-v1, zw36-v1 and dprk37probe-v2 have run.
Arm W (the published hand grades, `findings/2026-11-27-prereg38-gradeclock-plan.md`) is separate
and unaffected. This plan freezes `findings/2026-11-27-m2-design-draft.md` and lists every change
from it at the end.

## Question

Wang, Vidale et al. 2024 (Nature 631:340, https://doi.org/10.1038/s41586-024-07536-4) read
South Sandwich Islands (SSI) repeating earthquakes, recorded at the YKA and ILAR arrays, as
the inner core rotating forward until about 2008.5 and then back about 2.5 times more slowly.
Rival readings: a multidecadal oscillation (Wu, Song & Yang 2026 GRL,
https://doi.org/10.1029/2025GL119295); episodic local change of the inner-core surface (Zhang & Wen 2026
JGR SE, https://doi.org/10.1029/2026JB034595); and artefacts (Yang et al. 2026 CEE,
https://doi.org/10.1038/s43247-026-03961-7; Ackerley 2026 SRL, https://doi.org/10.1785/0220260111).

Rigid rotation makes the inner core's angle one function of time. Two consequences can be
tested together on measured waveforms, and as far as our searches on 2026-11-27 found, nobody has
published either test:

- **M2a, one clock.** One clock (a turn date T and a back-speed ratio r) must order the late-coda
  states of **every** multiplet **at both arrays at once**. Patch-by-patch surface change does
  not require this.
- **M2b, cross-array closure.** ILAR decorrelates more than YKA (Wang, Pang & Vidale 2026 EPSL
  692:120274, https://doi.org/10.1016/j.epsl.2026.120274). So a pair that re-matches at ILAR must
  also re-match at YKA, outside YKA's timing-error eras.

## Input (fixed here; produced by m2measure)

One row per event pair, multiplet and array: `array, mult, ev1, ev2, t1, t2, D, ok`.
- Multiplets: Wang et al. 2024 Supplementary Data Table 3, except I ("not good repeaters") and M
  ("too noisy"). No catalogue-proximity additions.
- YKA: D = 1 - cc in the c310 late-PKIKP-coda window. It comes from the same calls m0324_check
  v0.8.0's `run_yka` makes: yka.py v0.2.0 `fetch_event` / `process_event` / `analyze_multiplet`,
  with the yka_repick 0.9.2 constrained picker, as in m34.py. The functions are unchanged; only
  the event list differs, one multiplet at a time. `ok` = that pipeline's pair flag.
- ILAR: the events are beamed by the ilar-v1 pipeline with the parameters in
  results/ilar-v1/run_summary.json, then D comes from m0324_check v0.8.0's `ilar_frame` + `all_D`,
  unchanged and called once per multiplet. (m0324_check itself reads cached beams and has a
  fixed event list, so m2measure imports its functions rather than running its `main`.) `ok` =
  both events beamed and both pass ilar-v1's SNR gates (min_snr_group 3.0, min_snr_ref 5.0,
  min_snr_df 3.0).
- m2measure must reproduce our published multiplet-M03-24 rows (results/m0324chk-v2/yka_z004_pairs.csv
  and ilar_z004_matrix.csv, pairs common to both) within 0.001 in D. If it does not, nothing is
  graded until a recount note explains why.

## Rules (m2core v0.1.0)

- **YKA timing-error eras** (Ackerley 2026): events from 2014-12-01 up to 2020-04-08, and from
  2021-12-01 up to 2022-04-01. In M2a every YKA pair touching these eras is dropped. In M2b such
  pairs are reported but never decide the grade.
- **Same/changed cut.** YKA: same if D <= 0.098 (P34). ILAR: same if D <= the 75th percentile
  (nearest rank) of usable within-multiplet ILAR pairs less than 2 years apart, if there are
  at least 8 such pairs; otherwise 0.098. The cut is set once from the full data.
- **States.** For each multiplet at each array, members are events with at least one usable pair.
  Complete-linkage clustering: the linkage of two clusters is the largest D over their
  *measured* cross pairs (infinite if none was measured). Clusters merge while the smallest
  linkage is <= cut. Ties go to the first pair in sorted order.
- **Testable unit** (multiplet x array): >= 3 members, >= 2 states, and members on both sides of
  2009.0.
- **Clocks.** One-turn equiv(t, T, r) as in clocks34.py, plus the no-turn clock (position = date).
  - Test grid: T from 2000.00 to 2020.00 in 0.25-yr steps x 28 ratios (every 4th of WIDE_RATIOS,
    0.1-19.4), plus no-turn: 2,269 clocks.
  - Fine grid: T in 0.05-yr steps x all 109 WIDE_RATIOS, plus no-turn: 43,710 clocks. It contains
    the test grid.
  - A clock *orders* a unit if, with members sorted by clock position (ties by date), every state
    forms one contiguous run.
- **I** = the clocks that order every testable unit. I_test is I on the test grid; I_fine is I on
  the fine grid.
- **Null.** 1,000 shuffles (seed 3813) of event dates within each multiplet. States and
  testability stay fixed, and I_test is recomputed on each shuffle.
  p_consist = (1 + shuffles with I_test non-empty) / 1,001.
- **Drop-one.** Remove one event at a time (each event in a testable unit), then rebuild states and
  testability. A drop that leaves fewer than 3 testable multiplets counts as not robust.

**M2a grade** (needs >= 3 multiplets with a testable unit, else *untestable*):
- *broken*: I_fine is empty, and stays empty after every single drop.
- *monotone*: I_test contains the no-turn clock (no turn is needed).
- *holds*: I_test is non-empty, it does not contain no-turn, and p_consist <= 0.05.
- *mixed*: anything else.

**M2b grade.** Over within-multiplet pairs usable at both arrays and outside the YKA eras:
n = pairs with ILAR same, and v = those among them with YKA changed.
- *broken*: v >= 2 and v >= 0.25 n.
- *holds*: n >= 5 and v = 0.
- *untestable*: n < 5 (and not broken).
- *mixed*: anything else.

**Report only:**
- For each unit: its states, its share of the test grid and of the fine grid, and whether
  no-turn and each named clock order it.
- Named clocks in I: Wang (2008.5, 2.5) and Yang & Song (2009.0 and 2009.5, r 1).
- Variants: YKA only, ILAR only, YKA era pairs kept, multiplet A left out, and both cuts x0.8 and
  x1.2.
- The drop-one breakers, and the count of pairs that are YKA same but ILAR changed (rotation
  allows these).

## How we will read it

- **M2a broken:** no single one-turn clock orders all multiplets at both arrays, even with any one
  event removed. This counts against rigid one-turn rotation and for patch-dependent change or
  artefacts. It does not rule out rotation with several turns.
- **M2a holds:** one clock with a turn orders all of them, and shuffled dates rarely allow that. This
  is the first objective, multi-multiplet, two-array support for the backtracking reading.
- **M2a monotone:** the measured states need no turn.
- **M2b broken:** ILAR re-matches that YKA does not share, outside the hardware eras. Rotation does
  not predict these.

## Not blind, said plainly

We have already seen the published hand grades (Wang Table 2; 4 "YKA different, ILAR similar"
pairs: P30, P36, P38, P100). We have also seen our own YKA and ILAR D for multiplet M03-24 (~
Wang's A; P34, m0324chk-v2). The "without A" variant is reported for this reason. No other
multiplet has been measured by us.

## Hand checks (test_m2core.py)

Toy events a 2003, b 2006, c 2011, d 2013, with D = 0.05 within a group and 0.2 across.
- Pattern "mirror" {a, d} | {b, c}: two states. The no-turn clock does not order it. (2009.0, 1)
  does: positions a 2003, d 2005, b 2006, c 2007. Wang's (2008.5, 2.5) does not: a 2003,
  b 2006, d 2006.7, c 2007.5.
- Pattern "date" {a, b} | {c, d}: no-turn orders it. (2009.0, 1) does not.
- Case analysis over T (< 2003, 2003-2006, 2006-2011, 2011-2013, >= 2013), with ties broken by
  date, shows that no clock orders both patterns. So I_fine is empty for {mirror A, date B,
  mirror C}.
- Dropping c from B frees it: (2009.0, r 0.5003) orders A, C and B-without-c, so drop-one is not
  robust. The grade is *mixed*.
- Three "date" multiplets give *monotone*. Two multiplets give *untestable*.
- Complete-linkage chain: e1-e2 0.05, e2-e3 0.05, e1-e3 0.15 gives {e1, e2} | {e3}.
- ILAR cut: eight short pairs 0.05-0.12 give 0.10. Seven give the default 0.098.
- Closure toy: n = 3, v = 2 gives *broken*. One era pair is reported, and one reverse pair is
  counted.

## Forecasts (made before any M2 measurement)

- Testable multiplets >= 3: 60%.
- M2a: holds 20%, monotone 10%, broken 25%, mixed 30%, untestable 15%.
- M2b: broken 35%, holds 20%, mixed 15%, untestable 30%.

## Changes from the design draft

1. Units are multiplet x array, and the primary test is joint across both arrays. The draft
   used YKA states only for M2a.
2. Linkage uses measured pairs only, so a missing pair no longer forces a split.
3. The ILAR cut rule is now fixed (the draft left it open).
4. YKA eras now include Ackerley's Dec 2021-Mar 2022 window. Era pairs are dropped from M2a.
5. The null and the "holds" condition use the 0.25-yr test grid; "broken" needs the 0.05-yr fine
   grid. The draft used 0.01-yr steps. "Holds" is now judged by how often shuffles allow any
   common clock (p_consist). The draft compared grid shares, which fails when most shuffles
   allow none.
6. 1,000 shuffles instead of 2,000 (CPU budget).
7. A new outcome, *monotone*, is split out of *mixed*.

## Limits fixed in advance

- States come from a fixed cut on a noisy D.
- Multiplets have 3-7 events.
- The clock family has one turn only.
- The same events appear at both arrays, so the units are not independent.
- Missing pairs make states coarser.
- The YKA era rule removes post-2014 information at YKA.
- A *broken* M2a refutes one-turn rigid rotation of the sampled region only, under these
  measurements.

# M2 design draft: one rotation clock for every South Sandwich multiplet? (not yet frozen)

MV2, Munim, Inc. (mv2@munimtech.com), 2026-11-27. Text CC BY 4.0, code Apache-2.0.
Made by AI agents working in MV2's name; no human expert has checked it.

**Status: DESIGN DRAFT.** This is not a pre-registration yet. No waveform has been fetched or
measured for M2. Thresholds below are proposals. The frozen plan will be committed as
`findings/<date>-prereg38-m2-plan.md` before any M2 job runs, and it will say what changed
from this draft.

## Why

The three rival readings of the South Sandwich Islands (SSI) inner-core repeaters make
different predictions once **all** multiplets are looked at together:

- **Rigid rotation with one turn** (Wang, Vidale et al. 2024, Nature 631:340,
  https://doi.org/10.1038/s41586-024-07536-4): the inner core's angle relative to the mantle is one
  function of time. If it returns to the same angle at dates a and b, every path, multiplet
  and array sees the same structure again. So one clock (one turn date T, one speed ratio r)
  must order the coda states of every multiplet at once, and a re-match at one array must
  also be a re-match at the other.
- **Local episodic surface change** (Zhang & Wen 2026, JGR SE, https://doi.org/10.1029/2026JB034595):
  each patch of the inner-core surface changes on its own schedule. Different multiplets
  (and the two arrays, which sample different patches) may need different clocks.
- **Artifacts** (Yang et al. 2026 CEE, https://doi.org/10.1038/s43247-026-03961-7; Ackerley 2026 SRL,
  https://doi.org/10.1785/0220260111): source differences and station timing/hardware problems
  (YKA element timing errors 2014-2020) can make "change" where there is none.

Our P34 test of a single multiplet (A, the 2003-03-21 / M03-24 family) ended **untestable**:
catalogue proximity found no new members (findings/2026-11-27-prereg34-results.md). M2 therefore
starts from the published membership itself: the 16 multiplets of Wang et al. 2024 Supplementary
Data Table 3, which we fetched from the publisher in job icr-wang-v1 (we do not redistribute the
tables).

## Exploratory look at the published hand grades (not a test)

Wang et al. 2024 Supplementary Table 2 grades 142 pairs by eye at YKA and ILAR
(0 similar, 1 somewhat similar, 2 different, -1 noisy, -2 no data). We tallied it by hand:

- 70 pairs have a grade of 0-2 at **both** arrays. 31 agree exactly.
- 11 are fully discordant (0 at one array, 2 at the other):
  - YKA similar, ILAR different (7): P86 2000-09-11/2020-06-06; P25 2002-07-18/2009-05-05 (J);
    P60 2003-05-26/2005-04-01 (E); P27 2003-05-26/2017-01-18 (E); P150 2004-03-23/2023-01-07 (B);
    P137 2005-06-10/2022-05-19 (F); P44 2010-12-11/2014-12-11.
  - YKA different, ILAR similar (4): P30 2004-03-23/2009-12-02 (B); P36 2005-06-10/2012-12-03 (F);
    P38 2006-05-29/2015-05-24; P100 2013-08-23/2018-09-04.
- ILAR decorrelates more than YKA (Wang, Pang & Vidale 2026 EPSL 692:120274,
  https://doi.org/10.1016/j.epsl.2026.120274), so "YKA similar, ILAR different" fits rotation with a
  small residual offset. **"YKA different, ILAR similar" is the direction rotation does not explain.**
  Two of these four (P38, P100) have a second event inside the YKA timing-error era (2014-2020;
  both events are in Wang's own Table 4 of YKA clock corrections), so hardware is a live
  explanation for them. P30 and P36 are earlier and are not covered by that explanation.
- In multiplets B and F, the same first event (2004-03-23 in B, 2005-06-10 in F) shows both
  directions of discordance with two different partners.

These are hand grades by the original authors, read off their table. They motivate M2; they are
not evidence by themselves.

## Proposed test

**Members.** Wang Table 3 multiplets except I ("not good repeaters") and M ("too noisy"): 14
multiplets, 3-7 events each. No catalogue-proximity additions.

**Measurement (pinned, already public in the repo).** YKA: m0324_check v0.8.0 + yka_repick 0.9.2,
D = 1 - cc in the c310 coda window, the P34 thresholds (same if D <= 0.098; changed otherwise).
YKA 2014-12..2020-04 events get the hardware flag; a pair is graded at YKA only if it passes the
P34 SNR gates. ILAR: the ilar-v1 pipeline (D on PKPbc+ab, 2.8 s), plus the longer v2 window as
report-only.

**States per multiplet.** Complete-linkage clustering of members on YKA D with cut 0.098. A
multiplet is **testable** if it has >= 3 graded members, >= 2 states, and members on both sides of
2009.0.

**Clock family.** One-turn clock as in clocks34.py: effective position equiv(t, T, r), turn T from
2000.00 to 2020.00 in 0.01-yr steps, speed ratio r on clocks34's WIDE_RATIOS grid; plus the
no-turn (monotone) clock. A clock **orders** a multiplet if, sorted by position, its states do not
interleave (each state occupies one contiguous run).

**M2a, one clock (primary).** Let S_m be the set of (T, r) grid points that order multiplet m.
- *One-clock broken:* >= 3 testable multiplets, every S_m non-empty, the intersection of all S_m
  empty, and still empty after dropping any single member of any multiplet (drop-one).
- *One-clock holds:* the intersection is non-empty, the monotone clock fails in >= 1 testable
  multiplet (a turn is needed), and the intersection's grid share is at or below the 5th percentile
  of a null that shuffles member dates within each multiplet (2,000 shuffles, fixed seed).
- *Mixed:* anything else with >= 3 testable multiplets. *Untestable:* fewer.
- Report-only: whether Wang's own clock (T 2008.5, r 2.5) and Yang & Song's symmetric turn
  (T 2009.0 / 2009.5, r 1) lie in the intersection.

**M2b, cross-array closure.** For every within-multiplet pair graded "same" at one array
and measurable at the other: rotation predicts "same" at ILAR implies "same" at YKA.
- *Closure broken:* >= 2 pairs with ILAR same and YKA changed outside the YKA hardware era, and
  they make up >= 25% of the ILAR-same pairs measurable at YKA.
- *Closure holds:* >= 5 such ILAR-same pairs and none changed at YKA outside the hardware era.
- Hardware-era pairs are reported separately and never decide the grade.

**Forecasts (to be frozen with the plan).** Testable multiplets >= 3: 60%. M2a holds 30%,
broken 25%, mixed 30%, untestable 15%. M2b broken 35%.

## Compute (when the queue drains)

One job, network on, 2 GB, 1 CPU, data_key m2 (YKA + ILAR windows for about 60 events). It
waits until zw33-v1, zw35g2-v1, zw36-v1 and dprk37probe-v2 have run; it adds nothing to the
queue before then.

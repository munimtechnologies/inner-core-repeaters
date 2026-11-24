# Pre-registration 34: a mirror test of inner-core backtracking with unused members of M03-24

MV2, Munim, Inc. - 2026-11-24. Text CC BY 4.0; code Apache-2.0. Made with AI agents working in
MV2's name; no human expert has checked it.

Written and published **before any waveform of the 31 candidate earthquakes below has been
fetched, plotted or measured** by this project. Code: `m34.py` 0.1.0 (committed separately, with
`test_m34.py`); job `m34-mirror-v1`.

## The question

Two explanations now compete for the same repeating-earthquake records from the South Sandwich
Islands:

- **Rotation with one turn** (Wang, Vidale et al. 2024, *Nature* 631:340): the inner core
  turned back around 2008-2010 and has since been retracing its path, more slowly.
- **Local, episodic change of the inner-core surface** (Zhang & Wen 2026, *JGR Solid Earth*,
  doi 10.1029/2026JB034595): patches change, rest and recover on their own timetables.

For the multiplet M03-24 we already know (pre-registrations 11 and 14, results/mantlep2-v1,
results/m0324chk-v2) that YKA's late PKIKP coda (3-10 s, `D_c310`) splits the eight good members
into **X** = {2003-03-21, 2020-07-11, 2022-03-07, 2024-03-21, 2025-12-13} and **Y** =
{2008-02-17, 2011-02-20, 2014-12-09}, and that mantle P shows no split, so the split is on the
core path. Both explanations allow this.

**What separates them.** On a one-turn rotation clock every moment maps to an "effective date"
(the inner core's position): forward in time until the turn, then backward. If the coda state
depends only on the inner core's position, the members must form **contiguous blocks in
effective date**: all X, then any intermediate states, then all Y. Episodic surface change has
no such constraint; it allows, for example, X-Y-X flips during the recovery. The eight known
members cannot test this: the boundary sits somewhere between effective dates 2003.69 and 2005.92
on the Wang clock (turn 2008.5, return 2.5 times slower). New members inside 2003-2008 and
2015-2020 can. No published study of this multiplet uses any event from those two windows
(Zhang & Wen list 1999, 2003, 2008, 2011, 2020, 2022-03, 2022-10 and 2024).

## Candidates (fixed now)

USGS ComCat, queried 2026-11-24: within 25 km of 59.03 S, 25.68 W, mb >= 4.6, every event in
2003-03-22..2008-02-16 (outbound leg, 11) and 2014-12-10..2020-07-10 (return leg, 20). The Wang
clock with the old eight members **requires** the state shown (X or Y) or leaves it free.

| Leg | Date (UTC) | ComCat id | mb | Effective date (Wang) | Wang requires |
|---|---|---|---|---|---|
| out | 2003-05-16 | usp000bxap | 4.8 | 2003.37 | X |
| out | 2003-05-26 | usp000by6r | 4.8 | 2003.40 | X |
| out | 2003-08-03 | usp000c3yx | 5.0 | 2003.59 | X |
| out | 2003-12-13 | usp000cf4p | 4.6 | 2003.95 | free |
| out | 2005-08-21 | usp000dxs7 | 5.1 | 2005.64 | free |
| out | 2005-09-17 | usp000dzku | 4.7 | 2005.71 | free |
| out | 2006-01-02 | usp000e7e5 | 4.6 | 2006.00 | Y |
| out | 2006-09-29 | usp000eu1v | 4.8 | 2006.74 | Y |
| out | 2006-12-17 | usp000f0ed | 4.6 | 2006.96 | Y |
| out | 2007-10-10 | usp000fq8f | 4.7 | 2007.77 | Y |
| out | 2008-01-12 | usp000fwhh | 4.6 | 2008.03 | Y |
| return | 2015-08-14 .. 2020-07-06 | 20 events (list in `m34.py`) | 4.6-5.2 | 2003.69-2005.65 | free (but must be monotone: later = closer to X) |

## Method (all code already published and pinned)

For each candidate separately (old members plus that one candidate, so one candidate cannot
change another's result):

1. **Membership (mantle P).** `mantlep2.py` v0.13.0 (commit b4712f0) at the eight stations of
   pre-registration 13 (QSPA, VNDA, LCO, DBIC, PLCA, LPAZ, LBTB, LSZ), D_late 3-10 s after P,
   1-2 Hz. Per station q = share of the gated old members whose median D_late to the others is
   below the candidate's (pre-registration 14's rule; at least 5 old members gated). Decision
   over counted stations (at least 3): **same** if median q <= 0.75, **different** if >= 0.999,
   otherwise inconclusive. Only "same" candidates are members.
2. **State (YKA).** `m0324_check.run_yka` v0.8.0 (commit 7067d4f) with `yka_repick` 0.9.2
   (+/-2 s pick), on the nine Z004 events plus the candidate, exactly as icr-m0324chk-v2.
   dX, dY = median `D_c310` from the candidate to the X and Y members over `ok` pairs (at least
   3 X and 2 Y partners). Labels:
   - **X**: dX <= 0.07 and dY >= 0.07 and dY >= 3 dX; **Y**: the same with X and Y swapped;
   - **intermediate**: otherwise, if max(dX, dY) < 0.098 (the smallest published X-Y value);
   - **new-state**: min(dX, dY) >= 0.098; **ambiguous**: anything else; **untestable**: too few
     partners.
   The continuous coordinate s = dX / (dX + dY) is reported, not graded.
3. Order of work: YKA for all candidates (45% of the time budget), then mantle P, candidates
   with a testable YKA state first.

## Tests (on members with state X, Y or intermediate; old eight keep their published labels)

- **R1, Wang clock** (turn 2008.5, ratio 2.5). Violation = a pair whose state order (X <
  intermediate < Y) disagrees with its effective-date order. **consistent** (0 violations),
  **violated** (violations remain after dropping any one new member), **violated-fragile**
  (otherwise), **untestable** (no new member used).
- **R2, any one-turn clock**: turn 2007.0-2012.0 (step 0.1) x return ratio 0.5-5.0 (step 0.1),
  2,346 clocks. **compatible** (at least one clock with 0 violations), **contradicted** (none,
  even after dropping any one new member), **contradicted-fragile** (otherwise).
- **Brackets** (reported): outbound = last X before first Y among members before 2011.0;
  return = last Y before first X after 2011.0; "clean" or "overlap", with the number out of order.
  Today they are 2003.22-2008.13 and 2014.94-2020.53.
- **R34 reproduction**: the old-nine YKA run must reproduce the published
  `results/m0324chk-v2/yka_z004_pairs.csv` (max |difference| in D and D_c310 < 0.001).

## What each outcome means

- R1 consistent with at least one member whose state the Wang clock **required** (the "X"/"Y"
  rows above): the clock passed a risky prediction on unseen data.
- R1 violated but R2 compatible: the published clock is wrong in detail, but a one-turn
  rotation still fits; the compatible range of turn and ratio is reported.
- R2 contradicted: no one-turn rotation explains this multiplet's states; local episodic change
  (or a more complex rotation) is needed. That would be the first direct test against the
  one-turn picture that uses members no published model was fitted to.
- Few or no members: an honest untestable result; the 2003-2008 and 2015-2020 states then rest
  on interpolation in every published study.

## Forecasts (before any data)

- R34 reproduces: 85%.
- At least one candidate confirmed a member (source "same"): 40%.
- At least one member with a usable state (X, Y or intermediate): 30%.
- R1 testable: 30%. If testable: consistent 60%, violated or violated-fragile 40%.
- R2 contradicted (either kind), if testable: 15%.
- At least one member labelled intermediate, if any are usable: 30%.

## Limits fixed now

One multiplet, one array for the state (ILAR is noisier for this family and is not used).
Membership rests on mantle P at up to eight stations; M4.6-4.8 events may be gated out at most of
them. Catalogue depths are partly fixed (10, 33, 35 km), so ComCat depth is not used. The X/Y
labels of the old members were found exploratorily (icr-m0324chk-v2) and then confirmed as a
core-path effect (pre-registration 11); this test takes them as given.

## Credit

Wang, Vidale et al. (2024) for the multiplet and the clock; Zhang & Wen (2026) for the episodic
surface-change reading and their open data (Zenodo 10.5281/zenodo.19625711); USGS ComCat and
EarthScope (IRIS DMC) for the catalogue and waveforms.

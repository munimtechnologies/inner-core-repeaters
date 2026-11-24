# Pre-registration 34, amendment 1: which rotation clocks can order the M03-24 states?

MV2, Munim, Inc. - 2026-11-24. Text CC BY 4.0; code Apache-2.0. Made with AI agents working in
MV2's name; no human expert has checked it.

Written and published **before job `m34-mirror-v1` produced any result** (it is still queued; no
candidate waveform has been fetched or looked at). Code: `clocks34.py` 0.1.0 with
`test_clocks34.py` (committed separately; every expected test value was worked out by hand first).
Job: `clocks34-v1` (old eight members only); after m34-mirror-v1 lands, `clocks34-v2` reruns it
with `--candidates m34_candidates.csv`. Nothing in pre-registration 34 changes. This amendment only
adds readings.

## Why amend

Pre-registration 34 tests one family of rotation clocks: one turn, constant speed each way, with the
turn in 2007.0-2012.0 (R2). Since it was written we found a published clock outside that family.
**Wu, Song & Yang (2026, GRL 53:e2025GL119295, doi 10.1029/2025GL119295)** use individual events of
long earthquake sequences on the same South Sandwich-Alaska path (1982-2024) and on Kuril-Argentina
(1994-2024). They find the inner core rotated faster than the mantle by about 0.10 deg/yr from the
early 1980s, "decelerated to a near-zero rate around 2000", and was "perhaps slower than the mantle
after about 2010". They support the multidecadal oscillation model (Yang & Song 2023, Nature Geoscience
16:182: a halt around 2009, turning points about 35 years apart) and "do not support shorter-term
oscillations or bursts". Another published picture is a ~6-year oscillation (Wang & Vidale 2022,
Science Advances, doi 10.1126/sciadv.abm9916, from 1969-1974 nuclear tests). R2 covers none of these.

## The scale-free condition (exploratory, on the eight labels already known)

The YKA late-coda states of M03-24 (pre-registrations 11 and 14; results/m0324chk-v2) are
X = {2003-03-21, 2020-07-11, 2022-03-07, 2024-03-21, 2025-12-13} and
Y = {2008-02-17, 2011-02-20, 2014-12-09}. Suppose the state depends only on the inner core's angle
phi(t) relative to the mantle, with one boundary between the X and Y regions. Then for **any**
history that goes forward once and then back (any speeds, not just constant ones), the eight
states fit if and only if

- **C1:** phi(2014-12-09) is still on the Y side of phi(2003-03-21). By December 2014 the inner core
  had not retraced back to where it was in March 2003.
- **C2:** phi(2020-07-11) is beyond phi(2008-02-17) on the X side. By July 2020 it had retraced
  past where it was in February 2008.

(The other members are implied: on a rise-then-fall history a time between two members is never
below both of them.) These two conditions use no wave speed or structure size. They hold on the
Wang et al. (2024) clock (turn 2008.5, return 2.5 times slower). There the old members put the X-Y
boundary at effective dates 2003.69-2005.92, as recorded in pre-registration 34.

Consequences, worked by hand and checked by the unit tests:

1. **Equal speeds both ways.** This covers the turn of any long-period oscillation: the turn must lie
   in **2009.078-2014.329**. That is (2003.2177 + 2014.9392)/2 and (2008.1306 + 2020.5268)/2. A
   symmetric turn at exactly 2009.0 fails by 0.08 yr, and one at 2009.5 passes. The lower edge rests
   on one member, 2014-12-09. Without it the bound drops to (2003.2177 + 2011.1390)/2 = 2007.18. The
   upper edge rests on 2020-07-11. Without it the bound rises to 2015.16.
2. **Unequal speeds.** With the turn at 2005.0, the return must be more than 5.58 times slower than
   the outbound leg. With the turn at 2012.0, the return-to-outbound time ratio must lie in
   0.335-2.204. A turn before 2003.22 is impossible.
3. **A near-standstill around 2000-2010.** This is the reading of Wu, Song & Yang (2026). C1 then
   requires that the core's forward motion between March 2003 and the turn exceeded its backward
   motion between the turn and December 2014. If the rotation in 2003-2010 was "near-zero" and then
   reversed, that inequality is the whole question. We have asked for their numbers (below) and
   will evaluate it only after this amendment is public.
4. **Pure oscillation without drift** (phase free). Periods 6.0, 10, 12, 17.0, 17.5 and 17.8 yr
   cannot order the eight states. The three Y phases are always split by X phases. **Every period
   from 17.82 yr upward can.** The edge is where the phase of the 2025-12-13 member
   (25.950 - P) drops below that of 2008-02-17 (8.131), at P = 17.819 yr. At P = 6.0 yr the exclusion survives dropping any single
   member. By chance alone, 8/56 = 14.3% of random phase arrangements of 5 X and 3 Y would pass, so
   isolated passing windows at short periods mean little. The full scan (2-100 yr, step 0.01 yr)
   comes from clocks34-v1.

These are statements about the eight labels we already had. They are exploratory, not tests.

## Added readings for pre-registration 34 (fixed now, before any m34 result)

Run on the old eight plus the confirmed new members that m34-mirror-v1 gives a usable state
(`used_in_mirror` = True, state X, Y or intermediate), with the exact grids in `clocks34.py`:

- **R3, periodic clocks** (position = A sin(2 pi t/P + c), free phase, no drift; P = 2.00-100.00 yr
  in steps of 0.01 yr; intermediate states must sit between X and Y). Labels: **compatible** (at
  least one period passes), **excluded** (none passes, even after dropping any one new member),
  **excluded-fragile** (otherwise). Reported, not graded: `continuum_from`, the shortest period
  from which every longer one passes (17.82 yr today); band fractions; the 5.5-6.5 and 20-30 yr
  bands.
- **R4, wide one-turn grid**: turn 2003.00-2021.00 in steps of 0.05, times return ratio 0.1 x 1.05^k
  (k = 0..108, 0.1-19.4), 39,349 clocks. Labels: **compatible**, **contradicted** (none passes,
  even after dropping any one new member), **contradicted-fragile**.
- **R5, symmetric window**: the turn interval for equal speeds both ways, on a 0.001-yr grid.
  Labels: **unchanged**, **narrowed** or **closed**, against 2009.078-2014.329.

## Forecasts (before any m34 result; conditional on at least one usable new member, 30% per P34)

- R3 compatible: 90%. continuum_from rises above 17.82 yr: 30%.
- R4 contradicted (either kind): 12%.
- R5 narrowed: 45%. Closed: 10%.

## What would matter

If R4 is contradicted, no one-turn rotation history of any speed profile on the grid orders this
multiplet's states. Rotation would then need more than one turn, or the change is not rotation.
If R5 narrows, the members constrain when an oscillation turned, independently of any published
fit. Wu, Song & Yang's near-standstill and Yang & Song's 2009 halt can then be checked against
data they did not use.

## Limits fixed now

All readings assume that the coda state depends only on the inner core's angle and that one boundary
separates the X and Y regions along the sampled track. Episodic surface change (Zhang & Wen 2026)
breaks the first assumption by design. Small-scale structure with several boundaries breaks the
second. The labels come from one array (YKA). The 2014-12-09 member is the only one after the YKA
2013-14 hardware change with a Y label. That change has a D floor near 0.01, against an X-Y
separation near 0.1. "Turn" means the extremum of position, i.e. when the differential rotation
rate passes through zero.

## Credit

Wang, Vidale et al. (2024) for the multiplet and the one-turn clock; Yang & Song (2023) and Wu,
Song & Yang (2026) for the multidecadal picture; Wang & Vidale (2022) for the 6-year oscillation;
Zhang & Wen (2026) for the episodic-change reading.

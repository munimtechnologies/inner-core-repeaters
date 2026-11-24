# Pre-registration 34, amendment 1b: do Wu, Song & Yang's rates allow the M03-24 states?

MV2, Munim, Inc. - 2026-11-24. Text CC BY 4.0; code Apache-2.0. Made with AI agents working in
MV2's name; no human expert has checked it.

Written and published **before jobs `m34-mirror-v1` and `clocks34-v1` produced any result**. Both
are still queued, and no candidate waveform has been fetched. Code: `wsy1b.py` 0.1.0 and
`test_wsy1b.py`. Tables: `results/wsy1b/`. Every value was worked out by hand first and checked
by the unit tests. Nothing in pre-registration 34 or amendment 1 changes. This note does two
things. It evaluates point 3 of amendment 1, which is exploratory because the eight labels are
already known. It also adds one reading, R6, for the new members.

## Their numbers

Source: Wu, Song & Yang (2026), *GRL* 53:e2025GL119295,
[doi 10.1029/2025GL119295](https://doi.org/10.1029/2025GL119295). The numbers were copied from
the article by a read-only browser session on 2026-11-24.

- Section 4: "about 3 deg in about 30 years, with an average rotation rate of 0.10 deg/yr" from
  the early 1980s. Also: a "westward reversal in the recent 20 years, average rate of
  0.01 deg/yr", "and a maximum of 0.03 deg/yr".
- Abstract: "decelerated to a near-zero rate around 2000". Discussion: the rate was much reduced
  after about 2005.
- Reversal: "2010 from the best-fit spline". No uncertainty is given. It is "well supported
  within the 68% CI but less certain under the 95% CI".
- No sinusoid or other parametric model was fitted.

The text supports two readings of when the turn happened. One is **2010**, the spline zero
crossing. The other is **about 2004-2005**, where a "reversal in the recent 20 years" (to 2024)
would begin.

## The test (from amendment 1)

Assume the YKA late-coda state of M03-24 depends only on the inner core's angle, with one X/Y
boundary. The states are X: 2003-03-21, 2020-07-11, 2022-03-07, 2024-03-21, 2025-12-13 and
Y: 2008-02-17, 2011-02-20, 2014-12-09. Assume the core moves forward at a constant mean rate f
until the turn T, then back at a constant mean rate b. Then the eight states fit if and only if:

- **C1:** f (T - 2003.2177) > b (2014.9392 - T)
- **C2:** b (2020.5268 - T) > f (T - 2008.1306)

Only f/b matters. This is clocks34's return `ratio`.

## Result (exploratory; the labels were known)

Minimum mean forward rate in 2003.22-T that C1 needs, and the maximum that C2 allows, in deg/yr
(`results/wsy1b/wsy1b_turn_table.csv`):

| turn | b = 0.01: f must exceed | f may not exceed | b = 0.03: f must exceed | f may not exceed |
|---|---|---|---|---|
| 2004 | 0.140 | - | 0.419 | - |
| 2005 | 0.0558 | - | 0.167 | - |
| 2006 | 0.0321 | - | 0.0964 | - |
| 2007 | 0.0210 | - | 0.0630 | - |
| 2008 | 0.0145 | - | 0.0435 | - |
| 2009 | 0.0103 | 0.133 | 0.0308 | 0.398 |
| 2010 | 0.0073 | 0.0563 | 0.0218 | 0.169 |
| 2011 | 0.0051 | 0.0332 | 0.0152 | 0.0996 |
| 2012 | 0.0033 | 0.0220 | 0.0100 | 0.0661 |

1. **Their 2010 spline turn is compatible.** With the stated average westward rate of 0.01 deg/yr,
   the mean forward rate over March 2003 to 2010 must lie between 0.0073 and 0.0563 deg/yr. That
   is 7-56% of their 1980s rate. A modest residual forward drift is enough. At their stated
   maximum westward rate (0.03), it must be at least 0.0218 deg/yr (22%). A smooth turn, with
   equal speeds just either side of it, also passes at 2010, inside the window 2009.078-2014.329.
2. **A turn around 2004-2005 is not compatible with their own rates.** At b = 0.01, a turn in
   2005 needs more than 0.0558 deg/yr of forward motion in 2003-2005, and a turn in 2004 needs
   more than 0.140. That is 56% and 140% of their peak early-1980s rate. Both are at odds with
   their "decelerated to a near-zero rate around 2000". So if M03-24 changes by rotation alone,
   the "reversal in the recent 20 years" cannot have started before about 2006 at their rates.
3. **How late the turn must be depends on what "near-zero" means**
   (`results/wsy1b/wsy1b_earliest_turn.csv`):

| backward b | forward f (deg/yr) | earliest turn |
|---|---|---|
| 0.01 | 0.01 | 2009.08 |
| 0.01 | 0.02 | 2007.12 |
| 0.01 | 0.03 | 2006.15 |
| 0.03 | 0.01 | 2012.01 |
| 0.03 | 0.02 | 2010.25 |
| 0.03 | 0.03 | 2009.08 |

   Suppose the forward rate after 2003 was at most 0.02 deg/yr and the westward rate reached
   their stated maximum of 0.03. Then the turn must come after 2010.25, later than their spline's
   2010.
4. **Scale.** On the 2010 turn at b = 0.01, the core moves 0.036-0.28 deg between the March 2003 (X)
   and February 2008 (Y) events. That is at most about 0.8-5.9 km at the inner-core surface
   (21.3 km per degree there; less at depth). This is the whole X-to-Y change, and it is
   comparable to or smaller than one P wavelength in the inner core at 1-2 Hz (about 5-11 km).
   Under these rates, the rotation reading needs structure that changes the coda over a
   sub-wavelength move. We note this and do not grade it.

**In short:** the M03-24 states (Wang et al. 2024 multiplet; our YKA labels) are compatible with
Wu, Song & Yang's picture only if their reversal began around 2007-2010 or later, not in
2004-2005, and only if the core kept a small forward drift (at least about 7-22% of its 1980s
rate) after 2003.

## Added reading R6 (pre-registered now, for m34/clocks34-v2)

The source is `clocks34_oneturn.csv` from clocks34-v2 (the wide one-turn grid, the per-turn
ratio_min, sets "old8" and "all"). It is read by hand from that file, with no new code.

- **R6a:** the earliest grid turn with a compatible ratio of 2.0 or less. This is the WSY
  "near-zero" reading: forward at most 0.02 and backward 0.01 deg/yr. Our hand value for old8 is
  **2007.20**: at 2007.15 the smallest admissible grid ratio is 2.0594, and at 2007.20 it is
  1.9613. Labels: **later** if the "all" value is after old8's value, **unchanged** if equal.
- **R6b:** whether turn 2010.00 still has a compatible ratio of 2.0 or less in "all". Labels:
  **wsy-2010-compatible** or **wsy-2010-excluded**.
- Forecasts, if at least one new member is usable: R6a later 35%. R6b compatible 80%.
- Check: clocks34-v1's old8 rows must give 2007.20 for R6a. If they do not, the discrepancy is
  reported before any grade.

## Limits

These are mean rates, not rate histories, and only the ratio of the forward and backward means
over the member spans matters. WSY's rates are averages over many paths and years. They come from
travel times of different events (not M03-24), and their stated 0.03 is a maximum, not a bound
with uncertainty. Everything rests on rotation alone and one X/Y boundary. Episodic surface change
(Zhang & Wen 2026) does not need any of this. The 2014-12-09 Y label is the only one after the YKA
2013-14 hardware change. The hands copy of their text reached us only in part: what is quoted
above is all we used.

## Credit

Wu, Song & Yang (2026) for the rates. Wang, Vidale et al. (2024) for the multiplet. Yang & Song
(2023) for the multidecadal picture.

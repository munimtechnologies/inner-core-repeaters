# Pre-registration 38, arm W: results

MV2, Munim, Inc. (mv2@munimtech.com), 2026-11-27. Text CC BY 4.0, code Apache-2.0.
Made by AI agents working in MV2's name; no human expert has checked it.

Plan (frozen before the run): [findings/2026-11-27-prereg38-gradeclock-plan.md](2026-11-27-prereg38-gradeclock-plan.md).
Code: `gradeclock38/gradeclock38.py` v0.1.0 (unchanged since the freeze). Job `gradeclock38-v1` ran once on
the compete lane (Blade, 2026-11-27 18:17-18:34Z; unit tests passed, then the real run).
Outputs committed unchanged: `results/gradeclock38-v1/gradeclock38_summary.json`,
`gradeclock38_turn_profile.csv`, `gradeclock38_null.csv`. Wang et al.'s table itself is read at run time
and is not redistributed.

## Hand checks

`checks_ok: true`. 143 pair rows; 193 YKA comparisons and 141 ILAR comparisons; multiplet A at YKA on
Wang's clock 4 discordant of 72. All equal the values derived by hand before the run, so the run is graded.

## Grades (exactly as pre-registered)

| | YKA | ILAR |
|---|---|---|
| graded pairs in multiplets / comparisons | 95 / 193 | 62 / 141 |
| s, no-turn clock (date only) | 0.4715 | 0.4965 |
| s*, best test-grid clock | 0.1244 | 0.1206 |
| null s* (500 date shuffles): 5th / 50th / 95th pct | 0.259 / 0.358 / 0.440 | 0.199 / 0.355 / 0.440 |
| p1 (structure) | 0.002 (no shuffle reached s*) | 0.002 (no shuffle reached s*) |
| **W1** | **structure** | **structure** |
| improvement I = s(no-turn) - s* | 0.347 | 0.376 |
| p2 (turn) | 0.006 | 0.006 |
| **W2** | **turn-needed** | **turn-needed** |
| s on Wang's clock (turn 2008.5, back 2.5x slower) | 0.1399 | 0.1915 |
| share of 43,710 fine-grid clocks scoring strictly better | 0.34% | 2.23% |
| **W3** | **top5** | **top5** |

**W4: consistent.** The best YKA clock sits at the 2.23rd percentile of the ILAR scores, and the best
ILAR clocks at the 2.58th percentile of the YKA scores (both <= 10%).

p = 0.002 is the floor for 500 shuffles (1/501).

## Report only (pre-registered as report only)

- Best fine-grid clocks: YKA one clock, turn 2009.20, back-ratio 2.27 (s 0.1192). ILAR 33 tied clocks,
  turns 2007.25-2008.10, ratios 3.52-4.96 (s 0.1206).
- Symmetric turns (back at the same speed, ratio 1), which we use as a crude stand-in for a turn in an
  oscillation: turn 2009.0 gives s 0.4404 at YKA (56% of fine clocks better) and 0.3901 at ILAR (23%);
  turn 2009.5 gives 0.4041 (46%) and 0.3688 (20%). The grades prefer a slower return than the forward
  motion.
- The no-turn clock is no better than shuffled dates (p 0.33 at YKA, 0.49 at ILAR): on the calendar,
  pairs further apart in time are not graded less similar than chance would give.
- Turn profile (best s over ratios per turn date; `gradeclock38_turn_profile.csv`): at ILAR, s stays
  at 0.46-0.48 for every turn before 2003.5 and drops to 0.156 by 2004.25; at YKA, s falls from 0.34 (2003)
  to 0.24 (2004) and stays at 0.119-0.150 for turns 2008.0-2012.9. Turn date and ratio trade off along a ridge.
- Per multiplet (discordant comparisons, no-turn / Wang / best): YKA A 40/4/0 of 72, B 20/1/0 of 47,
  D 2/7/7 of 9, E 6/4/5 of 16, F 6/4/3 of 7, G 7/4/4 of 24, H 6/3/4 of 8. ILAR A 30/6/0 of 74,
  B 7/4/5 of 21, E 20/15/9 of 24, H 7/1/2 of 8, J 4/1/1 of 7. Multiplets A and B carry most of the YKA
  gain (60 of the 68 comparisons that improve), and A carries 30 of 53 at ILAR. D at YKA gets worse on the
  best clock (2 -> 7 of 9), and E at ILAR stays poor (9 of 24).
- Pooled with single pairs (report only): YKA s no-turn 0.442, best 0.177, Wang 0.186 (4,252 comparisons);
  ILAR 0.467, 0.166, 0.180 (1,704).

## Post-hoc reading (not pre-registered; computed by hand from the turn profile)

Along the low-s ridge, the clock position of an event in 2020.0 is T - (2020 - T)/r. YKA: 2004.4 (best clock),
2004.1 (Wang), 2004.7 (turn 2010.0, r 1.87), 2004.9 (turn 2011.0, r 1.46). ILAR best set: 2004.6-2004.7.
So the grades pin, at both arrays, roughly *which earlier year a 2020 event lines up with* (about 2004),
much better than they pin the turn date itself.

## Forecasts made before the run

W1 structure YKA 65%, ILAR 50%; W2 turn-needed 40% / 35%; W3 top5 50% / 30%; W4 consistent 40%. Every
outcome came out on the side we gave less than or about even odds for W2-W4; our forecasts were too
sceptical of the grades' internal consistency.

## How to read it (as fixed in the plan)

W1 structure + W2 turn-needed at both arrays + W4 consistent: **Wang et al.'s published grades carry a
one-clock signal with a turn, beyond what a flexible clock family finds in date-shuffled grades, and the
two arrays agree on the clock.** This supports the backtracking reading as an *internal consistency of
the published grades*. It is not independent evidence: the grades were assigned by eye by the authors
who proposed the reading, knowing the dates, and the pairs they chose to grade were not random.

Zhang & Wen (2026, https://doi.org/10.1029/2026JB034595) report, with their own quantitative similarity
definitions, that the time spans of similar pairs cannot be fit by any backtracking model. Arm W does not
test their grades, so the two statements are about different grade sets. Next: apply the same frozen
statistic to an independent grade set (Zhang & Wen's classification, or our own measured D in arm M).

## Limits

Categorical by-eye grades; 193 and 141 comparisons that share events (the shuffle null, not a formula,
carries the p-values); one-turn clocks only; multiplet quality varies; YKA has documented timing errors in
2014-2020 (Ackerley 2026, https://doi.org/10.1785/0220260111), which arm W does not correct.

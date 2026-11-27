# Pre-registration 38, arm W: do the published hand grades single out one rotation clock?

MV2, Munim, Inc. (mv2@munimtech.com), 2026-11-27. Text CC BY 4.0, code Apache-2.0.
Made by AI agents working in MV2's name; no human expert has checked it.

**Status: FROZEN before any run on the real table.** Code: `gradeclock38/gradeclock38.py` v0.1.0 with
`gradeclock38/test_gradeclock38.py` (expected values derived by hand). Job: `gradeclock38-v1`
(compete lane, 0.5 GB, 1 CPU). This is arm W of pre-registration 38. The measured-waveform arms
(M2a one clock on our own YKA/ILAR D, M2b cross-array closure; findings/2026-11-27-m2-design-draft.md)
stay in design and will be frozen separately before their job is filed.

## Question

Wang, Vidale et al. 2024 (Nature 631:340, https://doi.org/10.1038/s41586-024-07536-4) graded 143
South Sandwich Islands repeating-earthquake pairs by eye at the YKA and ILAR arrays (Supplementary
Data Table 2: 0 similar, 1 somewhat similar, 2 different; -1 noisy, -2 no data) and read them as the
inner core rotating forward until about 2008.5, then back about 2.5 times more slowly. Rival readings:
a multidecadal oscillation (Yang & Song 2023; Wu, Song & Yang 2026 GRL, https://doi.org/10.1029/2025GL119295),
episodic local change of the inner-core surface (Zhang & Wen 2026 JGR SE, https://doi.org/10.1029/2026JB034595),
and artefacts (Yang et al. 2026 CEE, https://doi.org/10.1038/s43247-026-03961-7; Ackerley 2026 SRL,
https://doi.org/10.1785/0220260111).

Arm W asks a narrow question that, as far as our searches on 2026-11-27 found, nobody has asked in
public: **taken alone, do the authors' own grades prefer one rotation clock, and a turn, more than
the same grades do after their dates are shuffled within each multiplet?** And do the YKA grades
and the ILAR grades prefer the same clock?

## Statistic

If the late PKIKP coda depends only on the inner core's angle, then on the right clock a pair's
grade can only rise with the distance between its two events' clock positions. Within each
multiplet, every comparison of two pairs with different grades is *discordant* if the less similar
pair is closer on the clock (a tie counts one half).

s(clock) = discordant comparisons / all comparisons, summed over multiplets (0 perfect, about 0.5 no relation).

- Pairs: grade 0, 1 or 2 at the array; multiplet letter from Table 2; multiplets I ("not good
  repeaters") and M ("too noisy") excluded as in Wang Table 3; pairs outside multiplets ("-") are
  used only in the report-only pooled score.
- Event dates: the table's dates as decimal years (day resolution).
- Clocks: one-turn `equiv(t, T, r)` (forward until T, then back r times slower; clocks34.py), plus
  the no-turn clock (position = date). Test grid: T 2000.00-2020.00 in 0.25-yr steps x 28 ratios
  (every 4th of clocks34's WIDE_RATIOS, 0.1-19.4), plus no-turn = 2,269 clocks. Fine grid (report
  and W3/W4): T in 0.05-yr steps x all 109 WIDE_RATIOS, plus no-turn = 43,710 clocks.
- Null: 500 shuffles (seed 38) of event dates within each multiplet; the full test-grid search is
  repeated on each shuffle. One shuffle is used for both arrays.

## Hand checks (the job flags `checks_ok: false` if any differs; then nothing is graded until a recount note explains it)

- 143 pair rows parsed.
- Comparisons: YKA 193 (A 72, B 47, D 9, E 16, F 7, G 24, H 8, K/L/N/O/P 2 each; C and J have one
  grade only); ILAR 141 (A 74, B 21, E 24, F 3, G 2, H 8, J 7, O 2).
- Multiplet A at YKA on Wang's clock (T 2008.5, r 2.5): 4 discordant of 72. (All four are grade-1
  vs grade-2 comparisons; the grade-0 pairs are all closer than every other pair.)

## Pre-registered outcomes (each array separately; an array with fewer than 30 comparisons is untestable)

- **W1, clock structure.** s* = min over the test grid. p1 = share of shuffles with s* at or below
  the real s* (with the +1 correction). *structure* if p1 <= 0.05; *no-structure* if p1 > 0.20;
  otherwise *weak*.
- **W2, turn needed.** Improvement I = s(no-turn) - s*. p2 = share of shuffles whose I is at or
  above the real I. *turn-needed* if p2 <= 0.05 and I >= 0.02; *not-needed* if I < 0.02 or
  p2 > 0.20; otherwise *weak*.
- **W3, Wang's clock.** Share of fine-grid clocks scoring strictly better than (T 2008.5, r 2.5).
  *top5* if <= 0.05, else *not-top5*.
- **W4, the two arrays agree on the clock.** For each array, the best fine-grid percentile, at the
  other array, of any clock in its own best set (all clocks tied at the minimum). *consistent* if
  both <= 0.10; *inconsistent* if both > 0.25; otherwise *mixed*.
- Report only: s at Yang & Song's symmetric turns (2009.0 and 2009.5, r 1); per-multiplet
  discordance at no-turn, at Wang's clock and at the best clock; the turn profile (best s per turn
  date); the pooled score including single pairs; p of the no-turn score against its own null.

## How we will read it

- W1 *structure* + W2 *turn-needed* at both arrays + W4 *consistent*: the grades carry a one-clock
  signal with a turn beyond chance. That supports the backtracking reading **as an internal
  consistency of the published grades**; it is not independent evidence, because the same people who
  proposed the reading graded the waveforms by eye, knowing the dates.
- W1 *no-structure* or W2 *not-needed* at an array: the grades at that array do not single out a turn
  any better than a flexible clock does for shuffled dates. Then the backtracking inference from that
  array's grades depends on more than the grades (for example on the waveform alignments shown in
  the paper's figures).
- W4 *inconsistent*: YKA and ILAR prefer different clocks. One rigid rotation predicts the same clock
  at both arrays, so this would favour array- or patch-dependent change (or array-dependent artefacts,
  such as the YKA 2014-2020 timing errors).

## Forecasts (made before the run)

W1 structure: YKA 65%, ILAR 50%. W2 turn-needed: YKA 40%, ILAR 35%. W3 top5: YKA 50%, ILAR 30%.
W4 consistent 40%, inconsistent 20%.

## Limits fixed in advance

Grades are categorical and by eye; the data are 193 and 141 comparisons that share events, so the
shuffle null (not a formula) carries the p-values; multiplets differ in quality; the clock family is
one turn only (no sinusoid in this arm). The table is the publisher's supplementary file; it is read
at run time and not redistributed here.

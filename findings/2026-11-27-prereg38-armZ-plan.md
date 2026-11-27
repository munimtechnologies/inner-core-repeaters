# Pre-registration 38, arm Z: do the surface-change camp's own calls order on one rotation clock?

MV2, Munim, Inc. (mv2@munimtech.com), 2026-11-27 (written 23:45-00:30Z). Text CC BY 4.0, code Apache-2.0.
Made with AI agents working in MV2's name; no human expert has checked it.

**Frozen before the grades are transcribed or read.** Code: `gradez/gradez.py` 0.1.0, tests
`gradez/test_gradez.py` (expected values derived by hand). Job (to be filed after the transcription
is committed): `gradez-v1`, compete lane, 0.5 GB, 1 CPU, network (reads the arm W summary from this
repo and Wang et al.'s workbook from the publisher; data_key icr-wang).

## Why

Arm W ([plan](2026-11-27-prereg38-gradeclock-plan.md), graded 27 Nov) found that Wang, Vidale et al.
(2024, Nature 631:340, https://doi.org/10.1038/s41586-024-07536-4)'s own by-eye grades order on one
one-turn rotation clock far better than date-shuffled grades at YKA and ILAR (p = 0.002 each; a turn is
needed, p = 0.006), and that the two arrays prefer the same clock. We said then that this is internal
consistency of one camp's grades, not independent evidence.

Zhang & Wen (2026, JGR Solid Earth 131, e2026JB034595, https://doi.org/10.1029/2026JB034595, "ZW26")
classify waveform pairs of the same kind of South Sandwich repeaters (multiplet M03-24, 2003-2024, and
multiplet F, 2005/2012/2022) at ILAR and YKAR with their own quantitative definitions, as change,
quiescence and recovery of the inner-core surface, and state that the time spans of similar pairs
cannot be fit by any backtracking model. Nobody has run a common, calibrated statistic on both grade
sets. Arm Z runs arm W's frozen statistic, unchanged, on ZW26's calls.

The question: **do ZW26's own calls order on a one-turn clock better than date-shuffled calls, and is
it the same clock that Wang's grades pick?** Either answer is informative. If ZW26's calls pick Wang's
clock, the disagreement between the camps is about interpretation, not about the waveforms. If they
pick no clock, or another one, the two camps' grade sets genuinely conflict and the dispute moves to how
pairs are graded, which arm M (our own measured differences) then decides.

## Data and transcription rules (fixed now, before reading)

1. Source: ZW26 main text, figures, figure captions and Supporting Information, read in the browser
   (open access, CC BY 4.0). A read-only hands copy job transcribes; it is told these rules.
2. A row is one explicit statement or figure mark that two **named events** (dates) at one array, in one
   region (sub-array, sub-region or "whole array"), are similar / unchanged (label **S**) or different /
   changed (label **C**). Quiescence between A and B = S. Change between A and B = C. A recovery between
   A and B = C for (A, B); a pair (earlier event, B) is S only if the paper names that earlier event as
   the one B matches again.
3. Nothing is inferred: no pair is graded from chains of episodes, from a model, or from our own data.
   Unstated pairs are absent.
4. Columns: multiplet,array,date1,date2,region,label,source (source = figure/table/section).
   The transcription is committed as `data/zw26_pair_calls.csv` (CC BY 4.0, credit ZW26) before the job.
5. Merging (code): all region labels of one pair at one array: all S = 0, all C = 2, mixed = 1.
6. Hand tallies of n_pairs and n_comp per array are made from the committed CSV and passed to the job
   with `--expect`; any mismatch sets checks_ok false and nothing is graded until a recount note.

## Statistic, clocks, null (identical to arm W)

Within each multiplet, every two pairs with different grades form a comparison; discordant if the less
similar pair is closer on the clock (ties 1/2); s = discordant / comparisons. One-turn clocks
equiv(t, T, r) plus the no-turn clock; test grid 2,269 clocks (T 2000-2020 step 0.25 x 28 ratios + no
turn), fine grid 43,710 (T step 0.05 x 109 ratios + no turn). Null: 500 shuffles of event dates within
each multiplet, full test-grid search on each, **seed 40**. An array with < 30 comparisons is untestable.

## Outcomes (each array separately)

- **Z1 structure:** as W1 (p <= 0.05 structure; > 0.20 no-structure; else weak).
- **Z2 turn needed:** as W2 (turn-needed if p <= 0.05 and improvement >= 0.02; not-needed if
  improvement < 0.02 or p > 0.20; else weak).
- **Z3 Wang's clock (2008.5, 2.5):** top5 if <= 5% of fine-grid clocks score strictly better.
- **Z4 same clock as arm W:** lowest fine-grid percentile, on ZW26's calls at that array, of any clock in
  arm W's best set for the same array (read from results/gradeclock38-v1/gradeclock38_summary.json:
  YKA (2009.2, 2.2705); ILAR 33 clocks, T 2007.25-2008.1, r 3.52-4.96). agrees <= 0.10; disagrees > 0.25;
  else mixed.
- **Headline per array:** same-clock (Z1 structure, Z2 turn-needed, Z4 agrees); other-clock (structure,
  turn-needed, Z4 disagrees); turn-unclear (structure, turn-needed, Z4 mixed); calendar-only (structure,
  not-needed); no-clock (Z1 no-structure); weak otherwise.
- **Z5 shared pairs:** pairs with the same two event dates in ZW26's calls and Wang's Table 2 at the same
  array, both graded 0 or 2: agree if >= 80% match, disagree if <= 50%, else mixed; < 5 shared pairs
  untestable.

Report-only: the variant without multiplet F's 2012-12-03 event (prereg 35 G2 found it is not a
mantle-P repeat); ZW26's own cross-array agreement (as W4); turn profile; s at Wang's clock.

## Reading

- same-clock at both arrays: the surface-change camp's own calls order like backtracking; the camps
  disagree on interpretation, not on which waveforms match. Still internal to by-eye or rule-based calls.
- other-clock or no-clock at an array, with Z5 disagree: the two published grade sets conflict on the
  same waveforms; the question is decided only by independent measurement (arm M).
- calendar-only: ZW26's calls behave like progressive change with no return, contrary to both papers.
- untestable (too few explicit calls): reported as such; arm Z then closes.

## Forecasts (before reading)

Testable at >= 1 array: 55%. Given testable: Z1 structure 60%; Z2 turn-needed 35%; headline same-clock
20%, other-clock 15%, calendar-only 20%, no-clock 25%, weak/turn-unclear 20%. Z5 testable 40%; given
testable, agree 35%, mixed 35%, disagree 30%.

## Limits

ZW26 grades mostly one multiplet (M03-24), so the shuffle null has few events to permute and the clock
search has much freedom relative to the data; the null pays for that, but power is low. Calls stated in
text are a selected subset. Regions differ between the papers. Only one-turn clocks are tested.

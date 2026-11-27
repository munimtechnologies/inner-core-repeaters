# Pre-registration 38, arm W: robustness checks (rules fixed before they run)

MV2, Munim, Inc. (mv2@munimtech.com), 2026-11-27. Text CC BY 4.0, code Apache-2.0.
Made by AI agents working in MV2's name; no human expert has checked it.

**Status: written after the arm W results were known** ([findings/2026-11-27-prereg38-W-results.md](2026-11-27-prereg38-W-results.md)),
and committed before the robustness job runs. These are labelled robustness checks of a graded result,
not a new confirmatory test. Code: `gradeclock38r/gradeclock38r.py` v0.1.0 with hand-derived tests
`gradeclock38r/test_gradeclock38r.py`; it imports the frozen `gradeclock38.py` v0.1.0 unchanged
(same grids, 500 shuffles, seed 38, same W1-W3 rules). Job: `gradeclock38r-v1` (0.5 GB, 1 CPU).

## Why

Arm W found that Wang et al. (2024)'s own YKA and ILAR grades favour one turn clock far beyond chance at
both arrays. Two obvious objections: (1) multiplet A, the largest family, supplies most of the gain
(60 of 68 improving comparisons at YKA with B; 30 of 53 at ILAR), so the result might be one family;
(2) YKA had timing errors in 2014-2020 and Dec 2021-Mar 2022 (Ackerley 2026, SRL,
https://doi.org/10.1785/0220260111), which could make pairs look different for non-core reasons.

## Checks and rules

- **Gate:** the job first reproduces arm W's counts and best test-grid scores (143 rows; 193 and 141
  comparisons; s* 0.124352 and 0.120567). If `checks_ok` is false, nothing is graded.
- **R1, drop multiplet A.** W1 and W2 exactly as in arm W, per array (under 30 comparisons: untestable).
- **R2, YKA without its timing-error eras.** YKA pairs with either event in 2014-12-01..2020-04-08 or
  2021-12-01..2022-04-01 count as ungraded. W1, W2 and W3 exactly as in arm W.
- **R3, does A's clock transfer?** Take every fine-grid clock tied for best on multiplet A alone; score the
  other multiplets on each. Graded quantity: share of fine-grid clocks scoring strictly better on the
  others than the upper median of those scores. *transfers* if <= 0.10; *does-not-transfer* if > 0.25;
  otherwise *mixed*. Per array. Report only: the most favourable member's percentile and Wang's clock on
  the others.

## How we will read it

- R1 *structure* + *turn-needed* at an array: the arm-W signal there does not rest on multiplet A alone.
  R1 *no-structure* or *not-needed* at both: the arm-W result is essentially one family, and we will say
  so wherever arm W is cited.
- R2 *turn-needed*: the YKA signal survives removing the bad-clock years. *not-needed*: it does not, and
  the YKA half of arm W is not evidence.
- R3 *transfers*: a clock learned on one family orders the others, which a single rigid rotation predicts.

## Forecasts (before the run)

R1 structure: YKA 55%, ILAR 35%. R1 turn-needed: YKA 40%, ILAR 25%. R2 structure 85%, turn-needed 70%,
top5 55%. R3 transfers: YKA 40%, ILAR 30%.

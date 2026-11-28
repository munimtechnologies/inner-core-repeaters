# At YKA, the two inner-core camps read the same waveform ordering and differ only in where they draw the line

MV2, Munim, Inc. (mv2@munimtech.com), 2026-11-28. Text and data CC BY 4.0. Made with AI agents working in MV2's
name; no human expert has checked it.

**Status: exploratory (post hoc).** This note compares two published sets of calls with numbers we published on
2026-11-21 ([`results/m0324chk-v2/yka_z004_pairs.csv`](../results/m0324chk-v2/yka_z004_pairs.csv)), before
either set of calls had been copied. The measure is D_c310: 1 minus the zero-lag correlation of the YKA array
beams in the 3-10 s late PKIKP coda, at 1-2 Hz. It has been our graded YKA measure since pre-registration 11.
Nothing here was tuned. But the comparison was not pre-registered, so its confirmatory test is frozen separately
([arm T](2026-11-28-prereg38-armT-plan.md)).

## The calls

- **ZW26:** Zhang & Wen (2026, JGR Solid Earth, https://doi.org/10.1029/2026JB034595) call each M03-24 pair
  at YKAR similar (S) or changed (C) using their sliding S(t) rule. We copied the nine calls word for word and
  confirmed them ([check](2026-11-28-zw26-disputed-calls-check.md)).
- **Wang:** Wang, Vidale et al. (2024, Nature 631:340, https://doi.org/10.1038/s41586-024-07536-4), Table 2,
  graded by eye: 0 similar, 1 intermediate, 2 different.

## YKA: one ordering, two cuts

| pair | our D_c310 | ZW26 | Wang |
|---|---|---|---|
| 2022-03-07 / 2024-03-21 | 0.0018 | S | (not graded) |
| 2003-03-21 / 2024-03-21 | 0.0130 | S | (not graded) |
| 2003-03-21 / 2022-03-07 | 0.0147 | S | 0 |
| 2008-02-17 / 2011-02-20 | 0.0179 | **C** | **0** |
| 2020-07-11 / 2022-03-07 | 0.0209 | **C** | **0** |
| 2003-03-21 / 2020-07-11 | 0.0228 | **C** | **0** |
| 2003-03-21 / 2008-02-17 | 0.1273 | C | 2 |
| 2003-03-21 / 2011-02-20 | 0.1365 | C | 2 |
| 2011-02-20 / 2020-07-11 | 0.1556 | C | 1 |

- **ZW26's nine YKAR calls are perfectly ordered by D_c310.** Every S pair has a lower D than every C pair: all
  18 S-versus-C comparisons are concordant. Their implied cut lies between 0.0147 and 0.0179, a gap of only a
  factor 1.22. If the nine values were in random order, a split this clean would happen with probability
  1 / C(9,3) = 1/84.
- **Wang et al.'s binary YKA grades on the same pairs are also perfectly ordered**, with all 8 comparisons
  concordant. Their implied cut lies between 0.0228 and 0.1273, a factor 5.6 gap. Their intermediate grade (1,
  for 2011/2020 at 0.156) does not fit the order.
- **The three pairs the camps dispute at YKA are exactly the three pairs with D between the two cuts**
  (0.018-0.023). On every pair below 0.015 and above 0.12 the camps agree.

So at YKA the dispute over these pairs is not about what the waveforms look like. Both camps' calls match one
measured ordering. The dispute is about whether a late-coda decorrelation of about 2% (cc of 0.977-0.982, against
0.985-0.998 for the pairs both camps or ZW26 call similar) counts as a change of the inner core.

## The wider M03-24 picture at YKA (all 36 published pairs)

The nine events fall into two tight groups and one noisy event. Group X is 2003, 2020-07, 2022-03, 2024 and
2025; within it D ranges from 0.0012 to 0.0228. Group Y is 2008, 2011 and 2014; within it D ranges from 0.018
to 0.033. Every X-Y pair lies between 0.098 and 0.156. The 2022-10-15 event has a direct SNR of 7.8 and D of 0.069
or more against every other event. The pre-registered arm M cut (0.098) sits at the edge of the X-Y gap, so arm
M calls all three disputed pairs "same". That is Wang et al.'s side, and it follows from the cut alone. This is
why arm M cannot settle these pairs neutrally and arm T is needed.

## ILAR (weaker; report-only)

At ILAR, our published D (results/m0324chk-v2/ilar_z004_matrix.csv, column D_reproduced) does **not** order
ZW26's ILAR calls: 12 of 20 S-versus-C comparisons are concordant. For example, 2003/2020 is called C with D
0.038, while 2003/2022 is called S with 0.075. Wang et al.'s binary ILAR grades are ordered (6 of 6; implied cut
between 0.075 and 0.160). Most of these ILAR rows fail our own ILAR quality gate (ok = False), so this is weak
evidence.

## Limits

- Post hoc: we saw both sets of calls before comparing them, though our D values were published first.
- One measure, one band and one window. In the earlier 1-3 s window (D_c13) and over the full window (D_full),
  ZW26's calls are not ordered. For example, 2008/2011 has D_c13 = 0.0005. The separation belongs to the late
  coda, which is also where ZW26's dissimilar windows lie (about 9.5-21 s on their S(t) axis).
- Wang et al.'s grades and ZW26's S(t) are different measures from ours. Agreement in ordering does not mean
  the same numbers.
- 2022-03-07 falls inside a documented YKA timing-error interval (Ackerley 2026, SRL,
  https://doi.org/10.1785/0220260111: B4/R4, Dec 2021-Mar 2022). Its pair with 2024-03-21 is nonetheless the
  third tightest of the 36 (D 0.0018).

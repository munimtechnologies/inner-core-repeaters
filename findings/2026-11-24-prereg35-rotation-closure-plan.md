# Pre-registration 35: the rotation-closure test (one angle or many patches?)

MV2, Munim, Inc. - 2026-11-24. Text CC BY 4.0; code Apache-2.0. Made with AI agents working in
MV2's name; no human expert has checked it.

Written **before any waveform of Zhang & Wen's Zenodo bundles has been read or measured**
(the header-only inventory job zw26-inv-v1 has not run yet). Code: `closure35.py` 0.1.0, tests
`test_closure35.py` (standard library only).

## The question

Two readings of the same repeating-earthquake records compete:

- **Rigid differential rotation** of the inner core (Wang, Vidale et al. 2024, *Nature* 631:340;
  Wu, Song & Yang 2026, *GRL* 53, e2025GL119295): between repeats, only the inner core's
  orientation changes.
- **Local, episodic surface change** (Zhang & Wen 2026, *JGR Solid Earth*,
  [10.1029/2026JB034595](https://doi.org/10.1029/2026JB034595)): patches change, rest and
  recover on their own timetables.

## The closure rule

If source and mantle stay fixed and only the orientation changes, every array's waveform
depends on **one** angle phi(t). Call (a, b, c), with dates a < b < c, a *recovery triple* at
array A when A sees a and b as different and a and c as the same. The plain rotation reading of
a recovery is phi(c) = phi(a), and then **every** array must see a and c as the same. Another
array B that sees a and c as different breaks the closure: rotation then survives only if A's
waveform happens to repeat at a different angle (a coincidence). Local change predicts no such
constraint.

This rule needs no clock, no rate and no turn year, so it does not depend on any of the clock
families in pre-registration 34.

## Why multiplet F is the sharp case

Zhang & Wen's own labels for their multiplet F (2005-06-10, 2012-12-03, 2022-05-19), as copied
from the paper on 2026-11-24 (`labels/closure35/zw26_multiF_paper_labels.csv`):

- YKAR, region R1: changed 2005 -> 2012, then "recovered to the state on 2005/06/10" by 2022.
  That is a recovery triple.
- ILAR, region R2: "no changes" 2005 -> 2012, then "exhibited new changes" 2012 -> 2022.

By transitivity, the ILAR labels imply 2005 differs from 2022 at ILAR. With the YKAR recovery,
the labels taken together break closure. `closure35.py` grades this `violated-implied`: the
paper's text does not state the ILAR 2005-vs-2022 comparison itself. **The decisive
measurement is that one pair, ILAR 2005-06-10 vs 2022-05-19, together with YKAR 2005 vs 2022,
measured the same way on the authors' own traces.**

## Measurement (fixed now)

- Data: Data_multiF.zip from Zenodo 10.5281/zenodo.19625711, the authors' processed ILAR and
  YKAR traces. If the bundle lacks either array for any of the three events, the test is
  **untestable** and is reported as such.
- Dissimilarity D and its floor D_floor: D is the one pre-registration 33 (E1) fixes for the
  M03-24 traces, applied unchanged to multiplet F. D_floor is that pair's noise-predicted
  dissimilarity under the same rule. Pre-registration 33 may choose window and band from the
  inventory (copying the paper's band), but it may not change anything below.
- Labels: r = D / max(D_floor, 0.005); **same if r <= 2.0, changed if r >= 3.0**, otherwise nd.
- Gates: (G1) YKAR must show changed(2005, 2012) and same(2005, 2022) under this rule, or there
  is no recovery triple and the grade is **untestable (no recovery)**. (G2) mantle-P source
  check (rule of pre-registration 11) if Data_ForRelocation or EarthScope has GSN mantle P for
  the three events. A pair graded S (different source) makes the grade **untestable (source)**.

## Grades (per triple; `closure35.grade`)

- **violated**: ILAR measured changed(2005, 2022). Rotation needs a coincidence; this supports
  local change.
- **holds**: ILAR measured same(2005, 2022). Then the paper's ILAR quiescence 2005 -> 2012 and
  change 2012 -> 2022 cannot both hold under our D, and we report which one fails.
- **untestable**: ILAR 2005 vs 2022 is nd, or a gate fails.

The same grader also runs on every M03-24 recovery triple once pre-registration 33 has measured
it. As a check that it reads labels correctly, our own pre-registration 11 split of M03-24 (X =
2003, 2020, 2022, 2024, 2025; Y = 2008, 2011, 2014, seen at both ILAR and YKA) gives 24 triples
that all hold (unit test). At that coarse level, our M03-24 labels are closure-compatible.

## Forecasts (stated now)

- Data_multiF has ILAR and YKAR traces for all three events: 65%.
- G1 passes (YKAR recovery reproduced under our D): 55%.
- If testable: **violated** 55%, **holds** 30%, **untestable (nd)** 15%.

## What a result would mean

A measured violation would be the first rule-based, clock-free check showing that one
repeating-earthquake data set cannot come from rigid rotation alone. A pass would tell against
a central label of the local-change reading. Either way it is graded on the authors' own data,
with rules published first.

## Prior work

Zhang & Wen (2026) for the labels and data. They conclude that rotation is "an unnecessary
add-on" to surface change. Arguments that one rotation cannot fit all the data are not new:
Yao, Tian, Sun & Wen (2019) argued rotation would need rates differing by a factor of 27 or
more across data subsets, and Vidale et al. (2025, *Nature Geoscience*) found one path steady
while another changed. What this pre-registration adds is narrower: a clock-free closure rule
with fixed thresholds, graded independently on the authors' own multiplet-F traces at the one
pair (ILAR 2005 vs 2022) their text does not state. Searches on 2026-11-24 (Crossref, Semantic
Scholar, web) found no comment, reply or replication of the 2026 paper. Wang,
Vidale et al. (2024) and Wu, Song & Yang (2026) for the rotation reading. Yao, Tian, Sun & Wen
(2019, *JGR*, 10.1029/2019JB017532) for the earlier surface-change argument.

# Correction to pre-registration 14: the 2026-11-13 event had already been processed at ILAR

By MV2, Munim, Inc. (made with AI agents). Written 2026-11-19, about 04:00 UTC, **before any
result of icr-m26-v1 exists**. Text CC BY 4.0. No human expert has reviewed it.

## What was wrong

Pre-registration 14 (findings/2026-11-19-prereg14-m0324-2026-event.md) says: "nobody working
on this project had fetched, plotted or measured any waveform of the events named below."
That is **not true for the ILAR array**. On 2026-11-17 (18:46-18:59 UTC) the automated patch
run `icr-m0324-v1` (m0324.py 0.5.0) queried USGS for every M >= 4.3 event within 45 km of the
M03-24 centroid up to 2026-11-17. That included us7000ty3q (2026-11-13) and us6000thwz
(2026-08-04). It fetched their ILAR records, built beams, linked events into multiplets and
wrote the results to `results/m0324-v1/`, which is public in this repo. The rebuild inside
`icr-m0324chk-v2` (2026-11-18) did the same from the cached beams.

The agents that wrote pre-registration 14 did not use or mention these rows. But they were
public, so the blinding statement was too strong. The statement holds for what
pre-registration 14 actually measures: **no YKA record and no mantle-P record of us7000ty3q
or us6000thwz had been fetched** by any job of this project when it was written. Its rules,
thresholds and forecasts stay exactly as published, and icr-m26-v1 runs unchanged.

## What the published ILAR rows already say about us7000ty3q

From `results/m0324-v1/m0324_multiplets.csv` and `m0324_pairs.csv` (found while preparing the
next test; exploratory):

* At ILAR, us7000ty3q does **not** join multiplet Z004 (M03-24). Its PKPbc+PKPab beam
  correlates at cc >= 0.90 with events of a different multiplet, **Z006**: 2000-10-08,
  2007-11-21, 2012-08-23, 2019-05-01, 2021-09-04 and 2023-10-15. Its highest beam SNR makes it
  Z006's reference event.
* us6000thwz (2026-08-04, the exploratory extra event of pre-registration 14) links at ILAR
  to ten events from 2003-2026, none of which is an M03-24 member.
* Other 2026 events in the same patch are also in ILAR multiplets: us6000ty0j (2026-09-27,
  Z007) and us7000txqr (2026-11-12, Z009).

**What this means for pre-registration 14.** Its YKA quality gate (reference-phase cc >= 0.90
with the Z004 members) and its mantle-P source score are the membership tests. Since ILAR
already links the event to a different family, **untestable** or **not-a-repeat** is more
likely than the 20% + 10% we forecast. This is a note made after seeing ILAR data. The
pre-registered forecasts are not changed, and the result will be graded against them as
written.

## ILAR values for us7000ty3q inside Z006 (exploratory, all pairs fail the DF-SNR gate)

`D` = 1 - zero-lag PKIKP (DF) correlation after reference alignment, ILAR 1-2 Hz beams. Every
Z006 pair is flagged `ok = False` because the DF-window SNR of these mb 4.3-5.2 events is below
the pipeline's gate of 3 (it is 1.2-2.8). Effective distance uses the backtracking clock of
Wang, Vidale et al. (2024): turn at 2008.5, motion back 2.5 times slower.

| Partner | D at ILAR | Years apart | Effective years apart (backtracking) |
|---|---|---|---|
| 2000-10-08 (usp000a1fh) | 0.118 | 26.10 | 0.38 |
| 2007-11-21 (usp000ft0f) | 0.159 | 18.98 | 6.74 |
| 2023-10-15 (us6000lfpq) | 0.163 | 3.08 | 1.23 |
| 2021-09-04 (us7000f7ug) | 0.288 | 5.19 | 2.08 |
| 2019-05-01 (us70003jaz) | 0.370 | 7.54 | 3.01 |
| 2012-08-23 (usp000jqn5) | 0.534 | 14.22 | 5.69 |

The partner most like the 2026 event at ILAR is the 2000 event, the one the backtracking
clock picks. The 2023 event, which the plain time order picks, is third. This is suggestive
but weak: the DF SNR is below the gate, the 2007 event (which backtracking puts far away) is
second, and inside Z006 the 2007 and 2023 events are the most similar old pair (D = 0.066),
which backtracking does not predict. The 2000-versus-2023 contrast at ILAR is only 1.38.

Because these ILAR values are already known, they cannot be a pre-registered test. The same
question at YKA, which nobody has looked at, is pre-registered in
findings/2026-11-19-prereg15-z006-2026-event.md.

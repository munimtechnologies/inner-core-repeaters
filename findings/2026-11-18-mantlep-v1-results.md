# Mantle-P control of the M03-24 late-coda split: v1 results (pre-registration 10)

By MV2, Munim, Inc. (made with AI agents). Written 2026-11-18 after `icr-mantlep-v1` finished
(10:43-11:15 UTC, 29.5 min). Text CC BY 4.0. Tables and plots: [`results/mantlep-v1/`](../results/mantlep-v1/).
Pre-registration: [`2026-11-18-mantlep-prereg.md`](2026-11-18-mantlep-prereg.md).

## The question

At YKA and ILAR the late PKIKP coda of eight M03-24 (Z004) earthquakes splits them into
X = {2003-03-21, 2020-07-11, 2022-03-07, 2024-03-21, 2025-12-13} and Y = {2008-02-17, 2011-02-20,
2014-12-09}, which is the split the 2024 backtracking clock predicts. Mantle P never enters the
core. If the split is a **source** effect (S), mantle P shows it too. If it arises on the **core
path** (N), mantle P does not.

## Pre-registered result: **untestable**

| Item | Forecast | Result |
|---|---|---|
| At least 3 usable stations | 70% | **No: 1 usable (IU.QSPA)** |
| Decision | S 55%, N 25%, inconclusive 20% | **untestable** (n < 3) |
| Positive control: YKA D_c310 rank 1 of 56 | 95% | **rank 1** ✓ |
| Positive control: ILAR long rank ≤ 3 | 80% | **rank 1** ✓ |
| YKA D_c13 and ILAR early not top-ranked | expected | ranks 8 and 10 ✓ |
| QSPA (the only usable station), D_late | – | rank 19 of 56 (not a hit) |

So by the pre-registered rule nothing is decided. The core-side positive controls behave as
forecast: the split is the best of all 56 possible 5+3 splits in both arrays' late coda, and
it is not top-ranked in the early windows.

## Why eight of nine stations failed: three faults in v1, found after the run

1. **Windows anchored on ak135 P.** The alignment window (−1 to +3 s) and the stack-onset search
   (±3 s) were centred on each event's ak135 P time. At QSPA the observed P comes about 2.5 s
   late and at CASY 3-7 s late, so the alignment window held little or none of the P wave.
   Events were then misaligned by whole cycles: at QSPA 2003 and 2020 sit about 1 s, and 2024
   about 2 s, away from the other five (see `plots/mantlep_IU_QSPA.png`). At CASY the P waves are
   scattered 6-10 s after the "onset" (`plots/mantlep_IU_CASY.png`).
2. **A silent gate failure.** The late-SNR noise window began 28 s before the stack onset, but
   the grid began only 30 s before ak135 P. Where the onset landed at the −3 s search edge
   (SUR −3.0 s, SHEL −2.85, RCBR −2.8, CASY −3.0) the noise slice ran off the grid, late SNR
   became NaN, and every event failed the gate. The onset at the search edge is itself the
   symptom of fault 1.
3. **All eight events required.** TSUM lacks 2003 and 2022, ASCN three events, NNA two, CASY and
   SHEL one each. Most long-running GSN stations miss one or two of these dates.

Some stations are also simply too noisy at 1-2 Hz for M4.8-5.3 events: direct SNR was 1.6-5 at
SUR, ASCN, SHEL, RCBR and SJG.

## Exploratory reading (not the pre-registered decision)

Read these as hypotheses, made after looking at the data.

**At the South Pole (QSPA) the five events that were aligned are indistinguishable in mantle P.**
2008-02-17, 2011-02-20, 2014-12-09, 2022-03-07 and 2025-12-13 have direct-P correlation
0.987-0.9999 and late (3-10 s) coda dissimilarity D = 0.002-0.013 in all ten pairs. The other
three events (2003, 2020, 2024) look the same by eye but are shifted by whole cycles (fault 1).

| Same pairs (five aligned events) | between X and Y (6 pairs) | within groups (4 pairs) | ratio |
|---|---|---|---|
| QSPA mantle P, D 3-10 s | mean 0.0083 | mean 0.0080 | **1.05** |
| YKA PKIKP coda, D_c310 | mean 0.126 | mean 0.020 | **6.4** |
| ILAR PKIKP, long window | mean 0.208 | mean 0.103 | **2.0** |

On the path that avoids the core, 2008-2014 and 2022-2025 look like the same source to within
D ≈ 0.01. On the path through the inner core, the same pairs differ about six times more than
pairs within a group at YKA, and twice as much at ILAR.

**At TSUM (Namibia)** all six recorded events align well (direct cc 0.95-0.995). Among the five
that pass the gates (2008, 2011, 2020, 2024, 2025) the mean between-group D (0.119) and the
within-group D (0.110) are nearly the same: no split.

**A side note on catalogue depths.** ComCat lists depths of 10 to 72.6 km for these events. A
depth difference of tens of km would change the depth-phase timing and decorrelate 1-2 Hz
waveforms, but the QSPA records of 2008 (66.5 km), 2011 (10 km) and 2025 (72.6 km) correlate
at 0.99-0.9997. The catalogue depths cannot all be right. This is common for small South
Sandwich events, and a reason not to use catalogue depths to split repeaters.

**What this does and does not mean.** If it holds up, the late-coda split of M03-24 is not a
source sub-cluster effect. The PKIKP coda of one repeating source would then have changed
between 2003 and 2008-2014 and changed back by 2020-2025, at two arrays. That is the pattern
the backtracking clock predicts, and also what an episodic change of the inner-core surface
(Zhang & Wen 2026) could produce. But it rests on one station, five events and a reading made
after the fact. The QSPA ray leaves the source towards the south at a take-off angle of about
30-40°. PKIKP leaves almost straight down, so a change in the structure right beneath the source
along that steep ray cannot be excluded by QSPA. Receiver-side changes at YKA and ILAR are not
tested here; the split cuts across YKA's 2013-14 hardware change, which the hardware control
showed is small (D ≈ 0.01).

## Next

A corrected test was pre-registered before any new station was fetched:
[pre-registration 11](2026-11-18-mantlep2-prereg.md) (`mantlep2.py`, job `icr-mantlep2-v1`). It
anchors on the observed P onset and allows a station with 6 or more good events. It adds every
IU, II and GT station 30-85° from the source. Its primary decision excludes the stations whose
v1 waveforms were looked at (QSPA, TSUM, CASY, and SPA at the QSPA site).

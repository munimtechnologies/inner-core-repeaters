# Results: ILAR diagnostics (icr-ilardiag-v1) and the M03-24 check (icr-m0324chk-v1)

By MV2, Munim, Inc. (made with AI agents). Written 2026-11-18. Text CC BY 4.0. Raw outputs:
[`../results/ilardiag-v1/`](../results/ilardiag-v1/) and
[`../results/m0324chk-v1/`](../results/m0324chk-v1/).

## Short version

* **ILAR's lean toward the backtracking ordering is consistent but weak, and it rests on the
  2021 sequence.** Leaving out any single multiplet never flips it, and controlling for SNR
  keeps it. But without the 2021 aftershock sequence it falls to p = 0.15, and the
  calibration P4 then **fails**: pairs less than 2 years apart are *less* alike than pairs
  that differ under both clocks (median D 0.116 vs 0.093). Inside the pre-2018 instrument
  epoch alone, the monotone clock orders the pairs better (ρ_M 0.21 vs ρ_B 0.10, neither
  significant). The best turn date is 2012.0, with a bootstrap 90% range of 2008.5-2018.5.
  That range includes Wang et al.'s 2008.5 but does not pin it.
* **M03-24 check:** the ILAR half (R0, C1, C2) **did not run**. Its cached beams are on
  another machine, because compete dataset volumes are per host. The YKA half ran, but
  **Y1 and Y0 are untestable**, because the onset picker triggered 4.8 s before PKIKP on
  weak early energy (see the [amendment](2026-11-18-yka-pick-amendment.md)). Ungated, the
  values point the same way as ILAR (D 0.044 to 2003 vs 0.017 to 2011), but they are not a
  valid measurement.
* `icr-bcdf-v1` also found no cached ILAR data on its host and measured nothing. All three
  are re-filed, with each job rebuilding its own cache first.

## 1. ILAR diagnostics (exploratory; written after the ILAR results were read)

123 usable pairs in 28 multiplets (the P1 set of [the ILAR results](2026-11-17-ilar-results.md)).
ρ_B and ρ_M are Spearman correlations of D with the time separation under Wang et al.'s
backtracking clock and under a monotone clock.

| Check | Result |
|---|---|
| Instrument epochs (ILAR digitiser change 2018-09-25) | median D: both pre-2018 0.187 (n 33), cross-epoch 0.080 (n 61), both post-2018 0.048 (n 29). Crossing the change does **not** raise D. |
| P1 inside the pre-2018 epoch only | n 33, 11 multiplets: ρ_B 0.10 (p 0.57), ρ_M 0.21 (p 0.24); B better in 37% of bootstraps |
| SNR control | ρ(D, SNR) = −0.24; partial ρ_B 0.21, partial ρ_M 0.09 |
| Leave one multiplet out | ρ_B − ρ_M stays between 0.05 and 0.17 (full 0.13); no flips |
| Without the 2021 sequence | n 104: ρ_B 0.14 (p 0.15), ρ_M 0.02; B better in 87% of bootstraps. **P4 fails**: near 0.116 (n 17) vs far 0.093 (n 9) |
| Within-multiplet sign test | B better in 4, M better in 1, tied in 4; p = 0.19 |
| Turn-date scan, multiplet bootstrap (300) | best T 2012.0 (ratio 2.0, ρ 0.263); bootstrap median 2011.75, 5-95% 2008.45-2018.5; 32% of draws in 2006-2011 |

The scan never prefers "no turn" (0% of draws). That is **not** evidence for a turn: the
maximum over about 45 turn dates is compared with a single monotone value, so it is biased
toward finding one.

**Reading.** The in-sample ILAR signal is an honest lean, not a detection. The out-of-sample
2025-2026 tests remain the real test. The P4 failure without 2021 means that, outside one
energetic sequence, D at ILAR does not clearly separate short from long time gaps, which
limits what ILAR alone can say. This is why YKA, where PKIKP is strong, carries the main test.

## 2. M03-24 check (pre-registered: [prereg](2026-11-17-m0324-check-prereg.md))

| Test | Forecast | Result |
|---|---|---|
| R0 reproduction (ILAR) | holds 90% | **not run**: no cached ILAR beams on the host |
| C1 element bootstrap (ILAR) | robust 55% | **not run** |
| C2 windows (ILAR) | robust 50% | **not run** |
| Y1 YKA replication | testable 25% | **untestable**: 0 of 36 pairs pass the gate |
| Y0 YKA in-sample pattern | testable 30% | **untestable** |

YKA processed all 9 events (12-18 elements each; 2003-2011 legacy system, 2014-2025 current), with
whole-window beam SNR 49-115, except 2022-10-15 (us6000ittv) at SNR 6. That event
anti-correlates with all the others (D > 1). Either it is not a repeater at YKA, or its record
is bad. It is excluded by the gates anyway.

**Why Y1/Y0 failed.** The 9-event stack pick was 4.8 s before ak135. In the stack plot
([`yka_z004.png`](../results/m0324chk-v1/plots/yka_z004.png)) the strong impulsive arrival sits
4.9-6.5 s after the pick, i.e. at the ak135 PKIKP time. Direct-pulse SNR was 0.8-4.5, below
the gate of 5, for every event. The [amendment](2026-11-18-yka-pick-amendment.md) fixes a
pick-validity rule for the queued stage B run and re-runs this check with the onset searched
within ak135 ± 2 s.

Ungated values, which are **not a valid measurement**: D(2003, 2025-12-13) = 0.044,
D(2011, 2025-12-13) = 0.017, D(2003, 2011) = 0.009, D(2003, 2020) = 0.065,
D(2003, 2024) = 0.004.

## 3. Lesson recorded

Compete dataset volumes belong to the machine that wrote them. A job that reads another
job's cache can land on a machine where the volume is empty. From now on, each job first runs
the published fetch/processing script, which is a no-op when the cache exists
(`fetch_pinned.py` + `m0324.py` / `run.py`).

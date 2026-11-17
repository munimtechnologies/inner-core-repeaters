# Results: the arc-wide forecast-date scan (ILAR) and YKA stage A

MV2, Munim, Inc. (made with AI agents), 2026-11-17. Text CC BY 4.0. Both runs were graded
against their pre-registrations:
[`2026-11-17-fdscan-prereg.md`](2026-11-17-fdscan-prereg.md) and
[`2026-11-17-yka-prereg.md`](2026-11-17-yka-prereg.md). Raw outputs are in
[`results/fdscan-v1/`](../results/fdscan-v1/) and [`results/yka-v1/`](../results/yka-v1/).

## In one paragraph

The arc-wide scan of 2025-2026 South Sandwich Islands earthquakes at Alaska's ILAR array
is **untestable**, and the reason is now clear. ILAR sees the outer-core phases strongly,
but its PKIKP (inner-core) signal from M4.5-5.3 events barely rises above noise (SNR about
1-5). So no 2025+ event kept two usable pre-2008 partners. Stage A at Canada's
Yellowknife array (YKA) shows the opposite:

* its archive is almost complete from 1990 to 2026 (1,861 of 1,881 events processed);
* its PKIKP is strong (131 events from 2025-2026 have a beam SNR of 5 or more);
* it independently found the 2014-2019 timing offsets that Wang et al. (2024) corrected by
  hand.

The next step follows directly: ILAR's outer-core phases choose the repeating pairs, and
YKA measures the inner core. It is pre-registered
([`2026-11-17-yka-stageb-prereg.md`](2026-11-17-yka-stageb-prereg.md)) and queued as
icr-ykab-v1.

## 1. Forecast-date scan (icr-fdscan-v1, fdscan.py 0.6.0, ILAR)

| Item | Value |
|---|---|
| Catalogue | 1,381 events: 378 events from 2025+ at M ≥ 4.5, plus the 1999-2008.5 pool at M ≥ 4.3 |
| Candidate cases in the catalogue | **296** (forecast: at least 50 at 70%, so this held) |
| Events processed / beamed | 666 / 640 with reference SNR ≥ 3 |
| Repeater edges (PKPbc+PKPab cc ≥ 0.90) | 155 |
| Families with a 2025+ member and a pre-turn member | 14 (F000-F013) |
| Pairs measured / passing all cuts | 736 / **16** |
| Pairs from 2025+ events passing all cuts | 3: F002 (2026-01-31 vs 2003-03-21), F009 (2025-12-13 vs 2003-03-21), F011 (2026-01-04 vs 2005-06-10) |

**Grades.**

* **F1, F1x and F1e (the best-match vote): untestable.** There were 0 cases, because no
  2025+ event had two ok pre-turn partners. The forecast had given "testable" 35%.
* **F3 (calibration): untestable.** There were 7 near pairs and 0 far pairs.
* **The P3 cross-check: untestable.**
* **Exploratory F5:** the turn-date fit on the 16 ok pairs gives T = 2009.25 and ratio 5
  (ρ 0.52, against ρ_M 0.10). 13 of the 16 pairs come from one 2002-2006 family (F013),
  so this says almost nothing about 2025-2026. It is not interpreted.

**Why it failed.** Most pairs fail one gate, PKIKP SNR ≥ 3 at ILAR. Two examples:

* In family F009 (2003-03-21 M4.8, 2008-02-17 M4.9 and 2025-12-13 M5.0), the reference-phase
  correlation is 0.99-1.00 and the reference SNR is 24-30. But the PKIKP SNR is only
  2.7-3.6.
* In the plot `fd_F009.png`, the PKIKP window is a small wiggle in front of large
  PKPbc/PKPab arrivals.

At 146-152° PKIKP is weak relative to PKP, so ILAR alone can grade only the larger events,
which is what Wang et al. did by eye. No below-gate vote was computed. Those D values are
dominated by noise, and reading them would not be a test.

**Lesson.** Because ILAR picks the pairs so well (the reference correlation is high), it
remains the right detector. It is the wrong array to measure PKIKP on for M < 5.5.

## 2. YKA stage A (icr-yka-v1, yka.py 0.2.0)

| Item | Value |
|---|---|
| Events tried / processed | 1,881 / **1,861** (1990-2026, M ≥ 5.0 before mid-2024, M ≥ 4.7 after) |
| 1990s coverage | 273 of 278 tried events processed (forecast: at least 50% at 60%, so this held) |
| Events with beam SNR ≥ 5 | 412 from 1990-2008; 793 of 1,150 current-hardware (2013+) events; **131 from 2025-2026** |
| Example | 2025-12-15 M5.3 (`yka_quicklook_2025-12-15_us6000ruvy.png`): beam SNR 88.7. The anomaly event 2025-12-13 M5.0 (us6000runq) has 56.3 at YKA, against 3.3 at ILAR |
| Source | NRCan for legacy elements, EarthScope for current ones, as pre-registered |

**Statics census.** The forecast gave 65% that element-year median offsets of 0.1-0.4 s
would appear on some elements in 2013-2020 and almost none before 2013. **This held.**
13 element-years have a median static of at least 0.1 s:

* all of them are in 2014-2019, between 0.108 and 0.203 s;
* they are on elements B3, B4, B6, B7, B9, R3 and R4;
* none is before 2014, and none is from 2020 to 2026.

`yka_statics.png` shows the same thing: flat rows for 1990-2013 and 2021-2026, and
step-like excursions of an eighth to a quarter of a second in 2014-2020. This is an
independent, automatic detection of the YKA clock errors that Wang et al. (2024) report
and corrected by hand (0.125, 0.25 and 0.375 s in 2013-2020). Stage B aligns every pair on
the direct PKIKP pulse, so these offsets cannot change its coda similarity.

**Polarity flips** across the 2013 hardware change can only be counted in stage B.

## What this means for the question

The 2024 *Nature* backtracking model makes a sharp, dated forecast for 2025-2026 South
Sandwich earthquakes. It has not been tested yet, and the reason is now understood. The
combined design (ILAR families, YKA measurement) uses each array where it is strong. It is
pre-registered (pre-registration 8) and queued.

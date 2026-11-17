# Pre-registration 4: the College, Alaska BC−DF series, 1982-2026, tested out of sample

Written 2026-11-17 (about 10:30 UTC), before any COL or COLA waveform was downloaded and
before any of this project's earlier jobs (`icr-ilar-v1`, `icr-yka-v1`, `icr-wang-v1`) had
run. By MV2, Munim, Inc. (made with AI agents). Text CC BY 4.0, code (`cola.py`,
`test_cola.py`, `bcdf.py` 0.4.0) Apache-2.0.

## Why this station

The first evidence for inner-core rotation came from this path: South Sandwich Islands (SSI)
earthquakes recorded at College, Alaska. Song & Richards (1996, *Nature* 382:221) found that
the BC−DF differential time grew by about 0.3 s over three decades. The newest study of it,
**Wu, Song & Yang 2026** (*GRL* 53, e2025GL119295, doi:10.1029/2025GL119295), measured 151
single-event SSI→COLA residuals from 1982 to **2024**. It found a fast rise before about
2000, a slower change after, and agreement with the multidecadal-oscillation model of
Yang & Song (2023, *Nat. Geosci.*): a pause around 2009, then a gradual turning back. They
found no support for short bursts or 6-year oscillations.

As far as we found (searches on 2026-11-17), no published BC−DF or waveform study of this
path uses events from 2025 on.

The EarthScope archive holds one record of the College site with no break except a few
weeks in 1996 (station service, checked 2026-11-17):

| years | channel | rate |
|---|---|---|
| 1982-01-08 to 1991-03-04 | DW.COL SHZ | 20 Hz |
| 1991-03-04 to 1996-05-11 | IU.COL BHZ | 20 Hz |
| 1996-06-14 to 1999-02-05 | IU.COLA BHZ | 20 Hz |
| 1999-02-05 on | IU.COLA 00 BHZ (primary) | 20, then 40 Hz |
| 1999-02-05 on | IU.COLA 10 BHZ (second sensor, cross-check only) | 40 Hz |

BC−DF is a time difference on a single trace, so station clock errors cancel. The COL to
COLA move is about 4 km, against ~16,800 km of path.

## Data and measurement (fixed)

* **Events:** USGS ComCat, box lat −61.5..−54, lon −31..−22, M ≥ 5.0 from 1982-01-08 to
  2024-07-01 and M ≥ 4.7 after (the same cuts as the ILAR pipeline). 2024-07-01 onward is
  fetched first; then the rest in order of decreasing magnitude, within a 190-minute budget.
* **Trace:** 140 s around the ak135 PKIKP time. Demean, detrend, 5% taper, causal 4-pole
  Butterworth **1-2 Hz** (primary), Lanczos resampling to 40 Hz, normalised by the rms noise
  35-5 s before PKIKP. Gappy, flat or short traces are dropped.
* **BC−DF:** exactly `bcdf.measure_bcdf` from pre-registration 3 (STA/LTA DF onset within
  ±4 s; DF template from onset − 0.2 s to onset + clip(BC−DF_pred − 0.4, 0.8, 1.6) s, correlated
  at lags BC−DF_pred − 0.6 to + 1.5 s, capped at PKPab_pred − 0.8 s). Events with ak135 BC−DF
  < 1.0 s are not measured.
* **Residual:** res = observed − ak135 BC−DF at the event's own distance and depth.
* **Usable:** cc ≥ 0.75, DF rms ≥ 3× noise, |DF onset − ak135| ≤ 3 s, and |res − median| ≤ 1.0 s.
* **Robustness, not decisive:**
  * the same measurement in a 0.6-1.6 Hz band, reusing the 1-2 Hz DF onset;
  * IU.COLA 10 measured for every usable event from 1999 on, with the 00 − 10 difference
    reported.

## Models (fixed)

res = c + a1·Δlat + a2·Δlon·cos(lat0) + a3·Δdepth/100 + b·g(t), with:

* **L (monotone):** g = t − 2008.5.
* **Bk (backtracking, Wang et al. 2024):** g = e(t) − 2008.5, where e(t) = t up to 2008.5 and
  2008.5 − (t − 2008.5)/2.5 after.
* **P (pause):** g = min(t, 2008.5) − 2008.5, so all change stops at the turn. A simple
  stand-in for "quiescence"; the gradual turning-back of Yang & Song lies between P and Bk.

All models are fitted by least squares on usable events before 2025.0.

## Questions and decision rules

* **C1 (calibration, must hold for the others to count as more than a measurement):**
  on usable events before 2008.5, the L slope b is **positive** in at least 90% of 1,000
  event bootstrap resamples. This is the classic Song & Richards signal; if our objective
  measurement cannot see it, the out-of-sample verdicts are reported but flagged as unreliable.
* **Q1 (in sample, descriptive):** which of L, Bk and P has the lowest residual sum of squares
  before 2025. This re-measures Wu et al.'s series, so it is a replication, not a discovery.
* **Q2 (the test): L vs Bk on 2025-2026 events.** Q2 is decisive only if there are at least
  **8** usable new events and the mean predicted gap |Bk − L| is at least **2×** the in-sample
  scatter / √n_new. If it is decisive, the model with the lower mean squared error in at
  least 90% of 2,000 bootstrap resamples of the new events wins; otherwise **neither** wins.
  This is the same rule as pre-registration 3.
* **Q2b: P vs Bk** and **Q2c: L vs P**, with the same rule. These are secondary.
* **Q3 (exploratory):** the turn date (1998-2022) and rate ratio (1-5) of the Bk-type model
  that best fit events before 2025.

## Our forecasts (stated now)

* **Coverage:** DW.COL 1982-1991 data are served and usable: **60%**. Usable events overall
  (all years): **150-500**. Usable 2025-2026 events: **10-40**, because single-station
  PKIKP at 150° is weak for M < 5.3.
* **C1 holds: 65%.** The early rise is the best-established part of the record, but our cuts
  are not Wu et al.'s and we fit a spatial gradient at the same time.
* **Q1 best:** P 40%, Bk 35%, L 25%.
* **Q2 decisive: 45%.** If decisive: **Bk 40%, L 25%, neither 35%**.
* **Q2b decisive: 20%.** Q2c decisive: 40%.
* **Cross-check:** the median IU.COLA 00 − 10 difference is within ±0.02 s, and at least 70%
  of events agree within 0.05 s.

## What would count

* **Wow:** C1 holds and Q2 is decisive. Then the station that founded the field would say,
  on earthquakes nobody has analysed, whether the inner core is retracing its path, and this
  could be checked against the ILAR waveform and BC−DF tests (pre-registrations 1-3).
* **Honest null:** C1 fails, or Q2 is not decisive. We would publish the series, the
  measurement code, and how many more events (or how much better picking) would make Q2
  decisive.

Either way, `cola_events.csv` (one row per event: ComCat ID, channel, BC−DF, residual,
quality) will be published under CC BY 4.0. It is an open, event-by-event series that can
be compared with Wu et al. 2026 and rebuilt from the code alone.

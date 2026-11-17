# Pre-registration 3: single-event BC−DF times at ILAR, out of sample in 2025-2026

Written 2026-11-17 (09:30 UTC), before any BC−DF time was measured, and before the first ILAR
run (`icr-ilar-v1`, queued) had downloaded a waveform. By MV2, Munim, Inc. (made with AI
agents). Text CC BY 4.0, code (`bcdf.py`, `test_bcdf.py`) Apache-2.0.

## Why a second observable

Pre-registrations 1 and 2 test inner-core "backtracking" (Wang, Vidale et al. 2024, *Nature*
631:340) with waveform similarity of repeating earthquakes. That needs repeaters, and
2025-2026 repeaters of old events may be few. The classic observable needs none: the
differential travel time of PKPbc (BC, outer core only) minus PKIKP (DF, through the
inner core) for South Sandwich Islands (SSI) earthquakes recorded in Alaska. Its change with
time was the first evidence for inner-core rotation (Song & Richards 1996, *Nature*), and
it underlies the multidecadal-oscillation model (Yang & Song 2023, *Nat. Geosci.*).

The newest study of this path, **Wu, Song & Yang 2026** (*GRL* 53, doi:10.1029/2025GL119295),
measured BC−DF residuals for 151 single SSI events at COLA (40 km from ILAR) from 1982 to
**2024**. It reports a fast rise before about 2000 and a slower one after, and argues for a
multidecadal oscillation with no burst or reversal at a few-year scale. As far as we
found (searches on 2026-11-17), no published BC−DF or waveform study of this path uses
events from 2025 on.

## Data and measurement (fixed)

* Events: every SSI event processed by `run.py` v0.1.0 (USGS M ≥ 5.0 to mid-2024 and
  M ≥ 4.7 from 2024-07-01; box lat −61.5..−54, lon −31..−22) with an ak135 BC−DF of at least
  **1.0 s** at the array centre (about 147° and beyond).
* Waveforms: ILAR (IM IL01-IL19 SHZ) 1-2 Hz causal beams, as cached by `run.py`
  (`/data/proc_v1/*.npz`, same `data_key` icr-ilar, no new downloads).
* DF onset: STA/LTA on the beam (0.5 s / 8 s, threshold 5) within ±4 s of the ak135 DF time.
* BC−DF: correlate the DF template (onset − 0.2 s to onset + clip(BC−DF_pred − 0.4, 0.8, 1.6) s)
  with the beam at lags from BC−DF_pred − 0.6 s to BC−DF_pred + 1.5 s. The upper bound is
  capped at PKPab_pred − 0.8 s. Parabolic sub-sample peak.
* Residual: res = observed BC−DF − ak135 BC−DF (at the event's distance and depth).
* Usable: correlation ≥ 0.75, DF window rms ≥ 3× noise, |onset − DF_pred| ≤ 3 s, and
  |res − median| ≤ 1.0 s (cycle-skip screen).
* Robustness: the same measurement element by element (median of elements with cc ≥ 0.6,
  at least 5), run through the same tests.

## Models (fixed)

res = c + a1·Δlat + a2·Δlon·cos(lat0) + a3·Δdepth/100 + g(t), with

* **L (monotone rotation):** g = b·(t − 2008.5)
* **Bk (backtracking):** g = b·(e(t) − 2008.5), where e(t) = t before 2008.5 and
  2008.5 − (t − 2008.5)/2.5 after (Wang et al. 2024's turn date and rate ratio).

Both models are fitted by least squares on events before 2025.0.

## Questions and decision rules

* **Q1 (in sample, descriptive):** which of L and Bk has the lower residual sum of squares on
  events before 2025, and the sign of b. This partly re-measures what Wu et al. measured at
  COLA, so it is a replication, not a discovery.
* **Q2 (out of sample, the test):** both fits predict each 2025-2026 event. Q2 is **decisive**
  only if there are at least **8** new usable events and the mean predicted gap
  |Bk − L| is at least **2×** the in-sample scatter divided by √(n_new). If it is decisive:
  * **Bk wins** if, in at least 90% of 2,000 bootstrap resamples of the new events, Bk's mean
    squared error is lower;
  * **L wins** if that happens in at most 10%;
  * otherwise **neither** wins.
* **Q3 (exploratory):** the turn date T (1998-2022) and rate ratio (1-5) that best fit events
  before 2025. Recovering T ≈ 2008 from travel times alone would be an independent echo of
  Wang et al.

## Our forecasts (stated now)

* ILAR SHZ coverage starts around 1999-2000 (an IL01 request for June 1997 returned no
  data), so the pre-turn baseline is short, about 8 years. **Q1 testable: 70%.**
* Wu et al. found only a slow change after 2000. The expected BC−DF signal from 2008.5 to
  2026 is therefore a few hundredths of a second, against a scatter of roughly 0.1 s from
  location and picking. **Q2 decisive: 25%.** If decisive: Bk 35%, L 35%, neither 30%.
* The most likely outcome is an honest "not decisive" result, together with the first
  public BC−DF series for this path that includes 2025-2026 events. A run on the nearby
  long-running station COLA (IU), the direct continuation of Wu et al.'s series, would be
  the next step.

## What would count

* **Wow:** a decisive Q2 (either way), agreeing with the waveform test of pre-registrations
  1-2. Then two independent observables, travel time and waveform shape, would agree out of
  sample on whether the inner core is retracing its path.
* **Honest null:** not decisive. We would still publish the series, the measurement code and
  the needed precision (how many 2025-2027 events, or which picking accuracy, would make
  Q2 decisive).

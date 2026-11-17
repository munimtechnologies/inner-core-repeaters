# inner-core-repeaters

**An open, objective and out-of-sample test of whether Earth's inner core is still
"backtracking".**

By MV2, Munim, Inc. (mv2@munimtech.com). The code and text were made with AI agents
working in MV2's name. Code: Apache-2.0. Findings text and derived tables: CC BY 4.0.

## Why

In 2024, Wang, Vidale and colleagues reported in *Nature* (doi:10.1038/s41586-024-07536-4)
that waveforms of repeating South Sandwich Islands earthquakes, after passing through
the inner core to arrays in Alaska (ILAR) and Canada (YKA), changed and then changed
back. They read this as the inner core reversing its rotation relative to the mantle
around 2008.5 and then moving back along the same path about 2-3 times more slowly.

Later work:

* Vidale et al. 2025 (*Nature Geoscience*) extended the data to 2024.
* Wu, Song & Yang 2026 (*GRL*, doi:10.1029/2025GL119295) used single events at COLA
  through 2024 and found a multidecadal oscillation with no sudden reversal.
* Zhang & Wen 2026 (*JGR Solid Earth*) argue that episodic changes of the inner-core
  surface explain the same records without any rotation.

Two things are still missing:

1. **An objective, open similarity measure.** The published similarity grades were
   given by eye, and the code is "available on request".
2. **A test on earthquakes that no published analysis has seen.** The backtracking
   model predicts that 2025-2026 events should resemble, through the inner core, their
   repeating partners from about 2001-2002, and not those from 2004-2008. The USGS
   catalogue lists 226 South Sandwich Islands earthquakes of magnitude 4.7 or more in
   2025-2026 (to 13 November 2026).

This repository supplies both. Each pre-registration was written before the data it
uses were analysed:

* [`findings/2026-11-17-prereg.md`](findings/2026-11-17-prereg.md): ILAR waveform test.
* [`findings/2026-11-17-yka-prereg.md`](findings/2026-11-17-yka-prereg.md): YKA second
  array.
* [`findings/2026-11-17-bcdf-prereg.md`](findings/2026-11-17-bcdf-prereg.md): single-event
  BC−DF travel times at ILAR, which need no repeaters.
* [`findings/2026-11-17-cola-prereg.md`](findings/2026-11-17-cola-prereg.md): the
  College, Alaska BC−DF series from 1982 to 2026 (DW.COL, IU.COL, IU.COLA), the station
  where inner-core rotation was first seen.

## What the ILAR pipeline does (`run.py`)

1. Gets SSI earthquakes from the USGS ComCat FDSN event service.
2. Downloads ILAR (network IM, IL01-IL19 SHZ) records around PKIKP from EarthScope FDSN.
3. Applies a causal 1-2 Hz filter, cuts each element on its own ak135 PKIKP time,
   normalises by noise, applies per-event static corrections and makes a beam.
4. Finds repeating earthquakes by cross-correlating the PKPbc+PKPab window. These
   phases stay in the outer core, so finding repeaters does not depend on the inner core.
5. For every pair in a multiplet: aligns on PKPbc+PKPab only, then measures the PKIKP
   window. It records D = 1 − cc, the residual PKIKP delay, and element-by-element
   medians.
6. Grades the pre-registered predictions P1-P4 (`model.py`) and an exploratory fit of
   the turn date and rate ratio (P5). It writes CSV tables, `predictions.json`,
   `report.md` and plots.

## Single-event BC−DF times (`bcdf.py`, `cola.py`)

For each event, the PKIKP onset is picked and the PKIKP pulse is correlated with PKPbc to
get the BC−DF time and its ak135 residual.

* `bcdf.py` uses the ILAR beams cached by `run.py`, with no new downloads.
* `cola.py` uses the single broadband station at College, Alaska. Its record runs
  from 1982 (DW.COL short-period, then IU.COL and IU.COLA broadband), the same span as
  Wu et al. 2026. It adds 2025-2026, a second band (0.6-1.6 Hz), and a cross-check on
  the second COLA sensor (location 10).

A spatial gradient (latitude, longitude, depth) is fitted together with one time term:

* a monotone rotation, L;
* Wang et al.'s backtracking, Bk;
* in `cola.py` also a pause, P, where change stops at 2008.5.

The terms are fitted on events before 2025. Each then predicts the 2025-2026 events, and
the pre-registered rule decides whether one prediction wins or the data are not decisive.
`cola.py` first checks that it recovers the classic rise of BC−DF before 2008.5
(calibration C1).

## The published repeater list (`wang_table.py`)

This script:

1. Downloads Supplementary Table 1 of Wang et al. 2024 (xlsx) and writes every sheet to
   CSV.
2. Finds each of their 121 events in USGS ComCat.
3. Flags which events fall inside this pipeline's magnitude cuts.
4. Lists every ComCat event since February 2023 near each of their 42 locations. These
   are candidates for new, out-of-sample members of the published multiplets.

## The second array: Yellowknife (`yka.py`)

YKA sits about 133-142° from the SSI, where PKIKP has no outer-core companion phase, so
repeaters are taken from ILAR's outer-core detection. Then:

* **Stage A** downloads YKA records for every catalogue event and makes beams.
  * Legacy YKB/YKR elements, 1989-2014, come from the Canadian National Seismograph
    Network (NRCan) FDSN service.
  * Current YKAB/YKAR elements, from 2013, come from EarthScope.
  * It writes archive coverage by year, a census of per-element static (clock) offsets,
    and quick-look plots with ak135 phases marked.
* **Stage B** (`--ilar-multiplets multiplets.csv`) aligns each pair on the direct PKIKP
  pulse and measures the similarity of the PKIKP coda (+1 to +6 s). It then grades the
  same predictions on YKA alone, and compares ILAR and YKA pair by pair.

```
pip install -r requirements.txt openpyxl
python run.py --data-dir ./data --out-dir ./out --time-budget-min 200
python bcdf.py --data-dir ./data --out-dir ./out_bcdf                     # after run.py
python cola.py --data-dir ./data_cola --out-dir ./out_cola
python wang_table.py --data-dir ./data --out-dir ./out_wang
python yka.py --data-dir ./data --out-dir ./out_yka                       # stage A
python yka.py --data-dir ./data --out-dir ./out_yka \
    --ilar-multiplets ./out/multiplets.csv --ilar-pairs ./out/pairs.csv   # stage B
python -m unittest -v test_xcorr test_model test_cluster_phases test_yka test_bcdf test_wang_table test_cola
```

## Status

* v0.1.0, 2026-11-17: ILAR code and pre-registration published; the first full ILAR run
  is queued.
* v0.2.0, 2026-11-17: the YKA second array (`yka.py`, `test_yka.py`) and its
  pre-registration were added before any YKA data were downloaded. The YKA stage A run
  is queued.
* v0.3.0, 2026-11-17: added single-event BC−DF travel times (`bcdf.py`) with their
  pre-registration, and the published-repeater cross-match (`wang_table.py`). Both were
  written before any ILAR waveform had been downloaded.
* v0.4.0, 2026-11-17: added the College, Alaska BC−DF series (`cola.py`,
  `test_cola.py`) and pre-registration 4, before any COL/COLA waveform was downloaded.
  `bcdf.py` gained the pause model P. The COLA run is queued.

Results go in `findings/` whatever they show, including null or negative results.

## Data credits

* Seismic data, through the EarthScope Consortium FDSN web services:
  * the Eielson array ILAR (FDSN network code IM);
  * College, Alaska: DW.COL (DWWSSN), and IU.COL and IU.COLA (Global Seismograph
    Network, IRIS/USGS, doi:10.7914/SN/IU);
  * the Yellowknife array YKA (network CN), from the Canadian National Seismograph
    Network (Natural Resources Canada) and EarthScope.
* Earthquake catalogue: the USGS ComCat.
* Travel times: ak135, via ObsPy TauP (ObsPy is LGPL-3.0, used as a library).
* Event list: Supplementary Table 1 of Wang et al. 2024 (*Nature*), fetched from the
  publisher at run time and not redistributed here.

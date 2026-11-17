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
around 2008.5 and then moving back along the same path about 2-3 times more slowly. A
2025 follow-up (Vidale et al., *Nature Geoscience*) extended the data to 2024. A 2026
paper (Zhang & Wen, *JGR Solid Earth*) argues that episodic changes of the inner-core
surface explain the same records without any rotation.

Two things are still missing:

1. **An objective, open similarity measure.** The published similarity grades were
   given by eye, and the code is "available on request".
2. **A test on earthquakes that no published analysis has seen.** The backtracking
   model predicts that 2025-2026 events should resemble, through the inner core, their
   repeating partners from about 2001-2002, and not those from 2004-2008. The USGS
   catalogue lists 226 South Sandwich Islands earthquakes of magnitude 4.7 or more in
   2025-2026 (to 13 November 2026).

This repository supplies both. The pre-registrations, written before any data were
analysed, are in [`findings/2026-11-17-prereg.md`](findings/2026-11-17-prereg.md) (ILAR)
and [`findings/2026-11-17-yka-prereg.md`](findings/2026-11-17-yka-prereg.md) (YKA).

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
pip install -r requirements.txt
python run.py --data-dir ./data --out-dir ./out --time-budget-min 200
python yka.py --data-dir ./data --out-dir ./out_yka                       # stage A
python yka.py --data-dir ./data --out-dir ./out_yka \
    --ilar-multiplets ./out/multiplets.csv --ilar-pairs ./out/pairs.csv   # stage B
python -m unittest -v test_xcorr test_model test_cluster_phases test_yka
```

## Status

* v0.1.0, 2026-11-17: ILAR code and pre-registration published; the first full ILAR run
  is queued.
* v0.2.0, 2026-11-17: the YKA second array (`yka.py`, `test_yka.py`) and its
  pre-registration were added before any YKA data were downloaded. The YKA stage A run
  is queued.

Results go in `findings/` whatever they show, including null or negative results.

## Data credits

* Seismic data:
  * the Eielson array ILAR (FDSN network code IM), through the EarthScope Consortium
    FDSN web services;
  * the Yellowknife array YKA (network CN), from the Canadian National Seismograph
    Network (Natural Resources Canada) and EarthScope.
* Earthquake catalogue: the USGS ComCat.
* Travel times: ak135, via ObsPy TauP (ObsPy is LGPL-3.0, used as a library).

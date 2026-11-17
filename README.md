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
   repeating partners from about 2001-2002, and not those from 2004-2008.

This repository supplies both. The pre-registration, written before any data were
analysed, is in [`findings/2026-11-17-prereg.md`](findings/2026-11-17-prereg.md).

## What the pipeline does (`run.py`)

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

```
pip install -r requirements.txt
python run.py --data-dir ./data --out-dir ./out --time-budget-min 200
python -m unittest -v test_xcorr test_model test_cluster_phases
```

## Status

v0.1.0, 2026-11-17: the code and the pre-registration are published. The first run on
the full catalogue is queued. Results will go in `findings/` whatever they show,
including null or negative results.

## Data credits

* Seismic data: the Eielson array ILAR (FDSN network code IM), served through the
  EarthScope Consortium FDSN web services.
* Earthquake catalogue: the USGS ComCat.
* Travel times: ak135, via ObsPy TauP (ObsPy is LGPL-3.0, used as a library).

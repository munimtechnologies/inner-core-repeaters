# inner-core-repeaters

**An open, objective and out-of-sample test of whether Earth's inner core is still
"backtracking".**

By MV2, Munim, Inc. (mv2@munimtech.com). The code and text were made with AI agents
working in MV2's name. Code: Apache-2.0. Findings text and derived tables: CC BY 4.0.

## Results so far (2026-11-17)

**Latest: the arc-wide scan and the second array**
([`findings/2026-11-17-fdscan-yka-a-results.md`](findings/2026-11-17-fdscan-yka-a-results.md),
tables in [`results/fdscan-v1/`](results/fdscan-v1/) and [`results/yka-v1/`](results/yka-v1/)).

* The arc-wide forecast-date scan at ILAR found 296 candidate cases and 14 families that
  join a 2025-2026 event to pre-2008 partners. **Its tests were untestable.** Only 16 of 736
  pairs passed the cuts, because ILAR's PKIKP signal from M4.5-5.3 events is barely above
  noise (SNR about 1-5). ILAR's outer-core phases, which pick the repeaters, are strong.
* **The Yellowknife array (YKA) is the opposite.** It processed 1,861 of 1,881 events from
  1990-2026, and 131 events from 2025-2026 have a PKIKP beam SNR of 5 or more. Its automatic
  statics census independently finds the 2014-2019 timing offsets (0.11-0.20 s) that Wang
  et al. corrected by hand, with none before 2014 or after 2019.
* **Next (queued, pre-registered before any measurement):**
  * ILAR families measured at YKA, which is the first test that can put 2025-2026 events
    against their pre-2008 partners
    ([`findings/2026-11-17-yka-stageb-prereg.md`](findings/2026-11-17-yka-stageb-prereg.md),
    `yka_b.py`);
  * a **hardware control**. Every one of those pairs crosses YKA's 2013-2014 change of
    electronics (same Geotech S-13 sensors, new pre-amp and digitiser). The two systems
    recorded side by side for nine months, so the same earthquakes can be compared on both
    to measure what the change alone does to the similarity metric
    ([`findings/2026-11-17-yka-hardware-prereg.md`](findings/2026-11-17-yka-hardware-prereg.md),
    `ykagen.py`).

**The M03-24 patch test**
([`findings/2026-11-17-m0324-results.md`](findings/2026-11-17-m0324-results.md), tables in
[`results/m0324-v1/`](results/m0324-v1/)).

* The objective pipeline re-found Zhang & Wen's multiplet M03-24 exactly: all 6 published
  events land in one multiplet, plus 2014-12-09, 2022-10-15 and 2025-12-13. Calibration
  holds.
* Before 2025 the patch looks like backtracking. The 2003 event's PKIKP resembles
  2020-2024 and not 2008-2014. The objective D also agrees with Wang et al.'s by-eye ILAR
  grades for this family (Spearman 0.67, 15 pairs).
* The pre-registered 2025-2026 tests were **untestable**. Of 46 new events near the patch,
  only one has a PKIKP signal above the pre-set gate.
* That one event, 2025-12-13, **goes against the backtracking forecast** (exploratory,
  one event). Its PKIKP matches the 2011 event (D 0.033), not the 2003 event (D 0.167).
  A pre-registered check (element bootstrap, other windows, and the YKA array in Canada)
  was filed before being run:
  [`findings/2026-11-17-m0324-check-prereg.md`](findings/2026-11-17-m0324-check-prereg.md).
* College, Alaska single-station BC−DF (`cola.py`) kept only 27 of 2,166 events (one from
  2025). Not decisive; that lane is paused (same findings file, section 2).

**The first full ILAR run**
([`findings/2026-11-17-ilar-results.md`](findings/2026-11-17-ilar-results.md), tables in
[`results/ilar-v1/`](results/ilar-v1/)).

* 1,559 South Sandwich Islands earthquakes were beamed at ILAR. The pipeline found 75
  repeating multiplets and measured 874 pairs, 136 of which passed the quality cuts.
* **Calibration (P4) holds.** **In-sample replication (P1) narrowly fails** the
  pre-registered bar: ρ_B = 0.19 (p = 0.03) against ρ_M = 0.07, but 88% of multiplet
  bootstraps favour backtracking where 90% was required.
* **The 2025-2026 tests (P2, P3) were untestable at M ≥ 5.0.**
* On this path PKIKP arrives about 4 s before ak135, the known SSI-Alaska anomaly. The
  BC−DF code was corrected for it before it ever ran
  ([amendment](findings/2026-11-17-bcdf-amendment.md)).
* **Cross-check against the published grades** (exploratory,
  [`findings/2026-11-17-wang-crossmatch.md`](findings/2026-11-17-wang-crossmatch.md)).
  120 of Wang et al.'s 121 events were matched to ComCat, but only 59 are inside the v1
  magnitude cuts. The pipeline re-found 17 of the 19 ILAR-graded published pairs it
  could see. On those 17 pairs, the objective D agrees with the by-eye ILAR grades
  (Spearman 0.62; 16 of 17 separated at D ≈ 0.07). The one clear disagreement is
  2009-05-05 vs 2022-09-10, graded "different" by eye but nearly identical in D.

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
* Zhang & Wen 2026 (*JGR Solid Earth*, doi:10.1029/2026JB034595) argue that episodic
  changes of the inner-core surface explain the same records without any rotation. They
  base this on one multiplet, M03-24, observed 2003-2024.

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
  BC−DF travel times at ILAR, which need no repeaters, and its
  [amendment](findings/2026-11-17-bcdf-amendment.md) (search windows only, written before
  any BC−DF time was measured).
* [`findings/2026-11-17-cola-prereg.md`](findings/2026-11-17-cola-prereg.md): the
  College, Alaska BC−DF series from 1982 to 2026 (DW.COL, IU.COL, IU.COLA), the station
  where inner-core rotation was first seen.
* [`findings/2026-11-17-m0324-prereg.md`](findings/2026-11-17-m0324-prereg.md): the
  M03-24 patch test, which goes down to M4.3 at Zhang & Wen's multiplet.
* [`findings/2026-11-17-fdscan-prereg.md`](findings/2026-11-17-fdscan-prereg.md): the
  forecast-date scan, which runs the same best-match test on every South Sandwich patch
  where it can tell the two models apart.
* [`findings/2026-11-17-m0324-check-prereg.md`](findings/2026-11-17-m0324-check-prereg.md):
  robustness and second-array check of the 2025-12-13 M03-24 anomaly (written after the
  anomaly was seen, before the checks ran).
* [`findings/2026-11-17-yka-stageb-prereg.md`](findings/2026-11-17-yka-stageb-prereg.md):
  YKA stage B on the merged ILAR families, with the new-event vote F1-Y (written after
  YKA stage A and before any stage-B YKA measurement).
* [`findings/2026-11-17-yka-hardware-prereg.md`](findings/2026-11-17-yka-hardware-prereg.md):
  the YKA hardware control on the 2013-2014 overlap of the legacy and current systems
  (written before any overlap waveform or stage-B measurement).

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

`ilar_diag.py` runs exploratory robustness checks on the published pair table: the
instrument epoch, SNR, leave-one-multiplet-out, results without the 2021 sequence, a
within-multiplet sign test, and a bootstrap of the turn date. It was written after the
ILAR results were read and is labelled as such.

## The M03-24 patch test (`m0324.py`) and its follow-up check (`m0324_check.py`)

Zhang & Wen's multiplet M03-24 has six events: 2003-03-21, 2008-02-17, 2011-02-20,
2020-07-11, 2022-03-07 and 2024-03-21, centred near 59.04°S, 25.69°W. The same patch has
two features no published study has used together:

* **46 USGS events of M ≥ 4.3 within 40 km in 2025-2026.**
* **Events on 2001-01-29 and 2002-06-29.** Under Wang et al.'s model these are the partners
  a 2025-2026 event should match best, because e(2025-2026) = 2001.9-2001.5.

`m0324.py` runs the unchanged `run.py` pipeline on every M ≥ 4.3 event within 35 km,
1991-2026. It then grades four things:

* Z0: does the objective pipeline re-find the published multiplet?
* Z1, the main test: does each new event's best PKIKP match pick the backtracking partner
  or the latest pre-turn partner?
* Z2: an out-of-sample rank test.
* Z3: calibration.

`m0324_check.py` takes the nine-event family and does the following:

* recomputes its D matrix from the cached ILAR beams;
* resamples ILAR's elements (bootstrap and leave-one-out);
* tries three other PKIKP windows;
* measures the same nine events at YKA with the published `yka.py` (fetched at a pinned
  commit).

Inputs `z004_events.csv`, `z004_edges.csv` and `z004_published_D.csv` are copied from the
icr-m0324-v1 output.

## The forecast-date scan (`fdscan.py`)

The scan runs the M03-24 best-match test on the whole South Sandwich arc.

* **Which events.** It takes every 2025+ event of M ≥ 4.5 that has catalogue neighbours
  (within 25 km, M ≥ 4.3) in two time windows:
  * the backtracking forecast window, |t − e(t)| ≤ 1 yr, about 2000-2003;
  * later, up to 2008.5, the monotone pick.
* **Waveforms.** It processes those events through the unchanged `run.py` pipeline, at
  most 700 events, largest new events first.
* **The vote.** For each new event it records whether the best PKIKP match is the
  backtracking pick or the monotone pick.
* **Decision rule.** The vote is graded with a one-sided binomial test and a bootstrap
  over whole multiplets. Robustness: excess dissimilarity Dx.
* **Independence.** A variant drops the M03-24 patch.
* **Controls.** It also reports the chance level and an SNR control.

The USGS catalogue, counted on 2026-11-17, holds 378 SSI events of M ≥ 4.5 since 2025
and 1,001 of M ≥ 4.3 from 1999 to mid-2008.

## Single-event BC−DF times (`bcdf.py`, `cola.py`)

For each event, the PKIKP onset is picked and the PKIKP pulse is correlated with PKPbc to
get the BC−DF time and its ak135 residual.

* `bcdf.py` uses the ILAR beams cached by `run.py`, with no new downloads. Since v0.5.0
  it searches DF from 9 s before to 1 s after ak135 and BC around its own ak135 time,
  because on the SSI-Alaska path DF arrives about 4 s early.
* `cola.py` uses the single broadband station at College, Alaska. Its record runs
  from 1982 (DW.COL short-period, then IU.COL and IU.COLA broadband), the same span as
  Wu et al. 2026. It adds 2025-2026, a second band (0.6-1.6 Hz), and a cross-check on
  the second COLA sensor (location 10). Its first run kept only 27 usable events, and one
  of its own synthetic unit tests fails (see the M03-24 results, section 2).

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
4. Lists every ComCat event since February 2023 near each published source location
   (12 after merging locations within 20 km). These are candidates for new,
   out-of-sample members of the published multiplets.

Its first run (2026-11-17) and the comparison with the ILAR v1 pairs are in
[`findings/2026-11-17-wang-crossmatch.md`](findings/2026-11-17-wang-crossmatch.md) and
[`results/wang-v1/`](results/wang-v1/). The pair-by-pair tables quote only the published
ILAR grades they compare against. The publisher's tables are not redistributed.

## The second array: Yellowknife (`yka.py`, `yka_b.py`, `ykagen.py`)

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
* **`yka_b.py`** (v0.9.0) runs stage B on the ILAR families of ilar-v1, m0324-v1 and
  fdscan-v1, merged where they share an event. It lowers the YKA catalogue to M4.3 so
  ILAR's small pre-2008 partners are included, and grades the new-event vote F1-Y of
  pre-registration 8. It fetches `yka.py` and its helpers at a pinned commit.
* **`ykagen.py`** (v0.9.1) is the hardware control.
  * **Which events:** every M ≥ 4.7 event 128.5-145.5° from YKA between 2013-06-27 and
    2014-04-01, when the legacy (CNSN digitiser, 20 sps) and current (YKA pre-amp and
    Guralp digitiser, 40 sps) systems both recorded from Geotech S-13 sensors on the same
    piers.
  * **How:** each system is processed separately with `yka.py`'s unchanged steps, and the
    two beams of the same earthquake are compared exactly as a stage-B pair.
  * **What it reports:** the hardware-only floor of D, the survival of the cc ≥ 0.90
    gate, polarity, and the timing offset between the systems (pre-registration 9).

```
pip install -r requirements.txt openpyxl
python run.py --data-dir ./data --out-dir ./out --time-budget-min 200
python ilar_diag.py --pairs ./out/pairs.csv --out-dir ./out_diag            # after run.py
python m0324.py --data-dir ./data_m0324 --out-dir ./out_m0324
python m0324_check.py --data-dir ./data_m0324 --out-dir ./out_m0324_check  # after m0324.py
python fdscan.py --data-dir ./data_fd --out-dir ./out_fd
python bcdf.py --data-dir ./data --out-dir ./out_bcdf                     # after run.py
python cola.py --data-dir ./data_cola --out-dir ./out_cola
python wang_table.py --data-dir ./data --out-dir ./out_wang
python yka.py --data-dir ./data --out-dir ./out_yka                       # stage A
python yka.py --data-dir ./data --out-dir ./out_yka \
    --ilar-multiplets ./out/multiplets.csv --ilar-pairs ./out/pairs.csv   # stage B
python yka_b.py --data-dir ./data --out-dir ./out_ykab --tables-dir .     # stage B, merged families
python ykagen.py --data-dir ./data_gen --out-dir ./out_ykagen             # hardware control
python -m unittest -v test_xcorr test_model test_cluster_phases test_yka test_bcdf test_wang_table test_cola test_m0324 test_fdscan test_ilar_diag test_m0324_check test_yka_b test_ykagen
```

## Data access notes (2026-11-17)

* The EarthScope FDSN dataselect service returned data for ILAR in 2002 and 2026, for
  DW.COL in 1985 and for IU.COLA.
* ILAR waveforms at EarthScope effectively start in 2000: the first full run beamed 0 of
  231 SSI events in 1990-1997, 4 of 31 in 1998 and 0 of 16 in 1999.
* YKA, by contrast, is nearly complete from 1990 (NRCan legacy elements) to 2026
  (EarthScope): 1,861 of 1,881 events.
* EarthScope has retired fdsnws-event and irisws-timeseries, and fdsnws-availability has
  been offline since July 2026. This code uses USGS ComCat for events and asks dataselect
  directly.
* ILAR SHZ response metadata change on 2018-09-25 (new digitiser chain, same sensor).
* YKA SHZ station metadata (EarthScope, 2026-11-17): the legacy YKB*/YKR* channels
  (1989-01-25 to 2014-04-02, 20 sps) and the current YKAB*/YKAR* channels (from
  2013-06-26, 40 sps) both use Geotech S-13 sensors at the same coordinates. The
  electronics differ (CNSN digitiser vs the YKA pre-amp and a Guralp digitiser).
* ObsPy cannot re-write DW.COL records stored in the DWWSSN encoding, so `cola.py` measures
  them in memory without caching them.

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
* v0.5.0, 2026-11-17: added the M03-24 patch test (`m0324.py`, `test_m0324.py`) and
  pre-registration 5, before any waveform from the patch was processed, plus the data
  access notes.
* v0.6.0, 2026-11-17: added the forecast-date scan (`fdscan.py`, `test_fdscan.py`) and
  pre-registration 6, before any waveform for it was processed.
* v0.7.0, 2026-11-17: **first results.** The ILAR run was graded against its
  pre-registration and published with its raw tables. Also added: `bcdf.py` and `cola.py`
  v0.5.0 (search windows fixed for the SSI-Alaska anomaly, before any BC−DF or COLA
  measurement) and `ilar_diag.py` (exploratory robustness checks).
* v0.7.1, 2026-11-17: cross-match with Wang et al. 2024 and an objective check of their
  ILAR grades (exploratory), plus a new 2026-11-03 member of their multiplet J (an
  anecdote below the SNR cut).
* v0.8.0, 2026-11-17: **M03-24 and COLA results** graded and published with raw tables.
  Added `m0324_check.py` and pre-registration 7 for the 2025-12-13 anomaly, before it ran.
* v0.9.0, 2026-11-17: **forecast-date scan and YKA stage A results** (scan untestable,
  limited by ILAR's PKIKP SNR; YKA strong and complete, and its statics census finds the
  2014-2019 clock offsets on its own). Added `yka_b.py` and pre-registration 8 (ILAR
  families measured at YKA) before any stage-B measurement.
* v0.9.1, 2026-11-17: added the YKA hardware control (`ykagen.py`, `test_ykagen.py`) and
  pre-registration 9, before any overlap-window waveform or stage-B measurement.

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
* Event list and the quoted ILAR grades: Supplementary Tables 1-2 of Wang et al. 2024
  (*Nature*), fetched from the publisher at run time and not redistributed here.
* M03-24 membership: Zhang & Wen 2026, from the file names of their CC BY 4.0 Zenodo
  dataset (doi:10.5281/zenodo.19625711). Event locations: USGS ComCat.

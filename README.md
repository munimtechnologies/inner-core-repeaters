# inner-core-repeaters

**An open, objective and out-of-sample test of whether Earth's inner core is still
"backtracking".**

By MV2, Munim, Inc. (mv2@munimtech.com). The code and text were made with AI agents
working in MV2's name. Code: Apache-2.0. Findings text and derived tables: CC BY 4.0.

## Results so far (2026-11-18)

**Latest: the inner-core change at M03-24 is on the path through the core, not at the
source** (pre-registered, decision **N**;
[`findings/2026-11-18-mantlep2-results.md`](findings/2026-11-18-mantlep2-results.md),
tables and plots in [`results/mantlep2-v1/`](results/mantlep2-v1/)).

* Eight earthquakes of one repeating South Sandwich Islands source (2003-2025) split into
  {2003, 2020, 2022, 2024, 2025} and {2008, 2011, 2014} by their inner-core waves at both
  Alaska (ILAR) and Canada (YKA). The inner-core waves changed after 2003 and changed back by
  2020, as the backtracking clock predicts.
* The mantle P waves of the same earthquakes, which never enter the core, show **no split
  at any of 7 usable global stations** that had not been looked at before (0 hits; the
  pre-registered rule needed at most 1 for N). Through the core, pairs across the groups differ
  8.0 times (YKA) and 3.4 times (ILAR) more than pairs within a group; in mantle P the median
  ratio is 1.19.
* So the earthquakes did not change; the path through the core did. The test cannot tell an
  inner core that turned back from an inner-core surface that changed and recovered, and a
  change confined to the structure straight below the source is not excluded. It is one
  multiplet; other families are next.

**Earlier on 2026-11-18: at Yellowknife the source sets the limit, not time; the ILAR BC−DF
lane is closed** (exploratory, [`findings/2026-11-18-diag2-results.md`](findings/2026-11-18-diag2-results.md),
tables in [`results/diag-v1/`](results/diag-v1/)).

* At YKA, coda dissimilarity between "repeaters" shows **no dependence on the time between
  them**. It grows with catalogue source separation (median ~15 km) and with magnitude
  difference, and it falls with SNR. Of the four explanations written down before the
  diagnostic ran, only location is supported. Families defined by ILAR's outer-core phases
  are not tight enough for YKA's inner-core coda at M4.3-5.3.
* In the 1-3 s window, calibration holds in all 8 families (sign p = 0.004). That window was
  not pre-registered for this purpose.
* The "two branches" in ILAR BC−DF times are failed DF onset picks (late picks, and early
  picks on an intermediate-depth cluster), not geography or cycle skips. Their share is flat
  in time. The core scatter (σ ≈ 0.57 s) is about 10 times too large to separate the clocks,
  so **single-event BC−DF at ILAR is closed**.
* The YKA lane continues only **prospectively**:
  [pre-registration 12](findings/2026-11-18-yka-prospective-prereg.md) fixes a
  location-and-SNR selection rule now, for earthquakes after 2026-11-19 only. It has two
  looks, in 2027 and 2028.

**Earlier on 2026-11-18: the first source control for M03-24 (mantle P) was untestable as
run** ([`findings/2026-11-18-mantlep-v1-results.md`](findings/2026-11-18-mantlep-v1-results.md),
tables and plots in [`results/mantlep-v1/`](results/mantlep-v1/)).

* Only 1 of 9 stations passed the gates (forecast: at least 3), so the pre-registered decision
  was **untestable**. The positive controls held: the hypothesised split is the best of 56 in
  both YKA's and ILAR's late PKIKP coda.
* v1 had three faults, found after the run. Windows were anchored on the ak135 P time, but P
  arrives 2.5-7 s late at some stations. A noise window ran off the grid, which made a gate
  fail silently. And every station needed all eight events.
* Exploratory: at the South Pole (QSPA), the five events that were aligned looked
  indistinguishable in mantle P. The corrected test,
  [pre-registration 11](findings/2026-11-18-mantlep2-prereg.md) (`mantlep2.py`), was published
  before it fetched any new station and left the stations already seen out of its primary
  decision. Its result is the "Latest" entry above.

**Earlier on 2026-11-18: the out-of-sample tests on 2025-2026 earthquakes return a null**
([`findings/2026-11-18-ykab-v2-bcdf-v2-results.md`](findings/2026-11-18-ykab-v2-bcdf-v2-results.md),
tables in [`results/ykab-v2/`](results/ykab-v2/) and [`results/bcdf-v2/`](results/bcdf-v2/)).

* **Yellowknife (YKA) stage B cannot be interpreted.** Its pre-registered calibration
  (P4-Y) failed: pairs less than 2 years apart were not more alike than pairs far apart
  (median D 0.098 vs 0.084). The new-event vote on 2025-2026 earthquakes could be run for the
  first time (6 cases in 4 families) and was **inconclusive**: 2 backtracking, 1 monotone,
  3 neither.
* **The YKA hardware control is complete and passes.** YKA's 2013-14 change of electronics
  alters the similarity measure by about 0.01, an order of magnitude below real pair
  differences (H1-H4 hold;
  [`findings/2026-11-18-ykab-v1-ykagen-results.md`](findings/2026-11-18-ykab-v1-ykagen-results.md)).
* **ILAR BC−DF travel times show no time trend** over 2000-2024 (+0.0014 s/yr). They cannot
  separate a monotone clock from backtracking on the 24 new 2025-2026 events (not decisive,
  as forecast).
* ILAR diagnostics: the lean towards backtracking at ILAR is consistent but weak. It
  depends on the 2021 sequence, and the best turn date is 2012 (bootstrap 2008.5-2018.5)
  ([`findings/2026-11-18-diag-m0324chk-results.md`](findings/2026-11-18-diag-m0324chk-results.md)).
* The 2025-12-13 anomaly at M03-24 is not robust at ILAR. Two arrays' late coda split
  M03-24 into {2003, 2020, 2022-03, 2024, 2025} and {2008, 2011, 2014}
  ([`findings/2026-11-18-m0324chk-v2-results.md`](findings/2026-11-18-m0324chk-v2-results.md)).

**The arc-wide scan and the second array (2026-11-17)**
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
* That one event, 2025-12-13, went against the backtracking forecast at ILAR (exploratory,
  one event). The follow-up check found this not robust, and at YKA the same event matches
  2003 best.
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
  YKA stage A and before any stage-B YKA measurement), and its
  [pick amendment](findings/2026-11-18-yka-pick-amendment.md).
* [`findings/2026-11-17-yka-hardware-prereg.md`](findings/2026-11-17-yka-hardware-prereg.md):
  the YKA hardware control on the 2013-2014 overlap of the legacy and current systems
  (written before any overlap waveform or stage-B measurement).
* [`findings/2026-11-18-mantlep-prereg.md`](findings/2026-11-18-mantlep-prereg.md):
  mantle-P control of the M03-24 late-coda split (source sub-clusters or not).
* [`findings/2026-11-18-mantlep2-prereg.md`](findings/2026-11-18-mantlep2-prereg.md):
  the corrected mantle-P control (written after v1 was untestable, before any new station
  was fetched; stations already seen are left out of the primary decision).
* [`findings/2026-11-18-yka-prospective-prereg.md`](findings/2026-11-18-yka-prospective-prereg.md):
  the prospective YKA test, which uses only earthquakes from 2026-11-19 on, with a
  location-and-SNR selection rule fixed in advance.

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
ILAR results were read and is labelled as such. `icrdiag.py` does the same kind of
exploratory work for the YKA stage-B pairs and the ILAR BC−DF series.

## The M03-24 patch test (`m0324.py`), its check (`m0324_check.py`) and the source control (`mantlep.py`, `mantlep2.py`)

Zhang & Wen's multiplet M03-24 has six events: 2003-03-21, 2008-02-17, 2011-02-20,
2020-07-11, 2022-03-07 and 2024-03-21, centred near 59.04°S, 25.69°W. The same patch has
two features no published study has used together:

* **46 USGS events of M ≥ 4.3 within 40 km in 2025-2026.**
* **Events on 2001-01-29 and 2002-06-29.** Under Wang et al.'s model these are the partners
  a 2025-2026 event should match best, because e(2025-2026) = 2001.9-2001.5.

`m0324.py` runs the unchanged `run.py` pipeline on every M ≥ 4.3 event within 35 km,
1991-2026, and grades Z0 (re-finding the published multiplet), Z1 (the new-event vote),
Z2 (an out-of-sample rank test) and Z3 (calibration). `m0324_check.py` re-measures the
nine-event family with an element bootstrap, three other windows and the YKA array.
`mantlep.py` (v1) and `mantlep2.py` (corrected) measure mantle P of the same events at
global stations 30-85° away. Mantle P does not touch the core, so it tells source changes
from path changes. `mantlep2.py` found 42 IU, II and GT stations, of which 8 were usable;
the decision is N (no source split).

## The forecast-date scan (`fdscan.py`)

The scan runs the M03-24 best-match test on the whole South Sandwich arc: every 2025+
event of M ≥ 4.5 with catalogue neighbours (within 25 km, M ≥ 4.3) both in the
backtracking forecast window (about 2000-2003) and later, up to 2008.5. The vote is graded
with a one-sided binomial test and a bootstrap over whole multiplets.

## Single-event BC−DF times (`bcdf.py`, `cola.py`), closed

For each event, the PKIKP onset is picked and the PKIKP pulse is correlated with PKPbc to
get the BC−DF time and its ak135 residual. `bcdf.py` uses the ILAR beams cached by
`run.py`; since v0.5.0 it searches DF from 9 s before to 1 s after ak135 and BC around its
own ak135 time. `cola.py` uses the single station at College, Alaska (1982-2026). A
spatial gradient is fitted together with one time term (monotone L, backtracking Bk, and in
`cola.py` a pause P) on events before 2025, and each predicts the 2025-2026 events. Both
lanes are closed: the onset-pick scatter (σ ≈ 0.5 s) is far above the few hundredths of a
second needed ([diag2 results](findings/2026-11-18-diag2-results.md)).

## The published repeater list (`wang_table.py`)

Downloads Supplementary Table 1 of Wang et al. 2024 at run time, finds each of their 121
events in USGS ComCat and lists new catalogue events near each published source. The
publisher's tables are not redistributed.

## The second array: Yellowknife (`yka.py`, `yka_b.py`, `ykagen.py`)

YKA sits about 133-142° from the SSI, where PKIKP has no outer-core companion phase, so
repeaters are taken from ILAR's outer-core detection.

* **Stage A** downloads YKA records (legacy YKB/YKR 1989-2014 from NRCan, current
  YKAB/YKAR from 2013 via EarthScope), makes beams, and writes coverage and a statics census.
* **Stage B** aligns each pair on the direct PKIKP pulse and measures the similarity of the
  PKIKP coda (+1 to +6 s).
* **`yka_b.py`** (0.10.1) runs stage B on the merged ILAR families and grades the new-event
  vote F1-Y. `--repick` searches the stack onset within ak135 ± 2 s (`yka_repick.py`).
  `ykab_grade.py` grades the pick-validity rule and H4.
* **`ykagen.py`** is the hardware control on the 2013-14 overlap of the two systems.

```
pip install -r requirements.txt openpyxl
python run.py --data-dir ./data --out-dir ./out --time-budget-min 200
python ilar_diag.py --pairs ./out/pairs.csv --out-dir ./out_diag            # after run.py
python m0324.py --data-dir ./data_m0324 --out-dir ./out_m0324
python m0324_check.py --data-dir ./data_m0324 --out-dir ./out_m0324_check  # after m0324.py
python mantlep2.py --data-dir ./data_mp --out-dir ./out_mantlep2           # source control
python fdscan.py --data-dir ./data_fd --out-dir ./out_fd
python bcdf.py --data-dir ./data --out-dir ./out_bcdf                     # after run.py
python cola.py --data-dir ./data_cola --out-dir ./out_cola
python wang_table.py --data-dir ./data --out-dir ./out_wang
python yka.py --data-dir ./data --out-dir ./out_yka                       # stage A
python yka_b.py --data-dir ./data --out-dir ./out_ykab --tables-dir . --repick   # stage B
python ykagen.py --data-dir ./data_gen --out-dir ./out_ykagen             # hardware control
python icrdiag.py --out-dir ./out_diag2                                    # exploratory
python -m unittest -v test_xcorr test_model test_cluster_phases test_yka test_bcdf test_wang_table test_cola test_m0324 test_fdscan test_ilar_diag test_m0324_check test_yka_b test_ykagen test_ykab_grade test_icrdiag test_mantlep test_mantlep2
```

## Data access notes

* ILAR waveforms at EarthScope effectively start in 2000. YKA is nearly complete from 1990
  (NRCan legacy elements) to 2026 (EarthScope): 1,861 of 1,881 events.
* EarthScope has retired fdsnws-event and irisws-timeseries, and fdsnws-availability has
  been offline since July 2026. This code uses USGS ComCat for events and asks dataselect
  directly.
* ILAR SHZ response metadata change on 2018-09-25 (new digitiser chain, same sensor).
* YKA's legacy (20 sps) and current (40 sps) channels use Geotech S-13 sensors at the same
  coordinates; the electronics differ.
* ObsPy cannot re-write DW.COL records stored in the DWWSSN encoding, so `cola.py` measures
  them in memory without caching them.
* Mantle P from the South Sandwich Islands arrives 2.5-7 s after ak135 at some global
  stations (QSPA, CASY), so windows must be anchored on the observed onset (`mantlep2.py`).

## Status

* v0.1.0-v0.6.0, 2026-11-17: code and pre-registrations 1-6 (ILAR, YKA, BC−DF, COLA,
  M03-24, forecast-date scan), each published before its data were analysed.
* v0.7.0-v0.7.1, 2026-11-17: first ILAR results and the cross-match with Wang et al.
* v0.8.0, 2026-11-17: M03-24 and COLA results; pre-registration 7.
* v0.9.0-v0.9.1, 2026-11-17: forecast-date scan and YKA stage A results; pre-registrations
  8 (YKA stage B) and 9 (hardware control).
* v0.9.2, 2026-11-18: ILAR diagnostics and the YKA pick amendment.
* v0.10.0, 2026-11-18: M03-24 check v2 results; pre-registration 10 (mantle P).
* v0.11.0, 2026-11-18: YKA stage B v1 (pick rule triggered) and hardware-control results.
* v0.12.0, 2026-11-18: YKA stage B v2 and ILAR BC−DF results (calibration fails at YKA,
  vote inconclusive, BC−DF not decisive); `icrdiag.py` and its plan.
* v0.13.0, 2026-11-18: mantle-P control v1 results (untestable; exploratory QSPA reading
  points to a path effect) and pre-registration 11 (`mantlep2.py`).
* v0.13.1, 2026-11-18: diagnostic results (YKA limited by source separation, not time;
  ILAR BC−DF closed) and pre-registration 12 (prospective YKA test).
* v0.14.0, 2026-11-18: **mantle-P control results (pre-registration 11): N.** The M03-24
  split is absent in mantle P at 7 unseen stations, so it arises on the path through the
  core; synthesis of all lanes.

Results go in `findings/` whatever they show, including null or negative results.

## Data credits

* Seismic data, through the EarthScope Consortium FDSN web services:
  * the Eielson array ILAR (FDSN network code IM);
  * College, Alaska: DW.COL (DWWSSN), and IU.COL and IU.COLA (Global Seismograph
    Network, IRIS/USGS, doi:10.7914/SN/IU);
  * the Yellowknife array YKA (network CN), from the Canadian National Seismograph
    Network (Natural Resources Canada) and EarthScope;
  * Global Seismograph Network stations of networks IU (doi:10.7914/SN/IU) and II
    (IRIS/IDA, doi:10.7914/SN/II), and the GT network, for the mantle-P control.
* Earthquake catalogue: the USGS ComCat.
* Travel times: ak135, via ObsPy TauP (ObsPy is LGPL-3.0, used as a library).
* Event list and the quoted ILAR grades: Supplementary Tables 1-2 of Wang et al. 2024
  (*Nature*), fetched from the publisher at run time and not redistributed here.
* M03-24 membership: Zhang & Wen 2026, from the file names of their CC BY 4.0 Zenodo
  dataset (doi:10.5281/zenodo.19625711). Event locations: USGS ComCat.

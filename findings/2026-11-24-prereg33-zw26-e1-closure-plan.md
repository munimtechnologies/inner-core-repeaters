# Pre-registration 33: do Zhang & Wen's own records give their labels, and does rotation close? (E1 + the P35 measurement)

MV2, Munim, Inc. - 2026-11-24. Text CC BY 4.0; code Apache-2.0. Made with AI agents working in MV2's name; no human expert has checked it.

**Status: written and published before any waveform in Zhang & Wen's Zenodo bundles was read.** The only things looked at are the file trees that Zenodo's preview shows (folder and file names and sizes; listed below). Code: `zw33.py` 0.1.0 with `test_zw33.py` (stdlib-only tests). Job: `zw33-v1`. It also grades pre-registration 35 with `closure35.py` 0.1.0, unchanged.

## The question

Zhang & Wen (2026, *JGR Solid Earth*, [doi 10.1029/2026JB034595](https://doi.org/10.1029/2026JB034595)) read the repeating South Sandwich multiplet M03-24 (2003-2024) at the ILAR and YKA arrays as "irregular episodes of change, quiescence and recovery" of the inner-core surface. They also read multiplet F (2005-06-10, 2012-12-03, 2022-05-19): YKA changed and then recovered, while ILAR stayed quiet and then changed. They call rotation "an unnecessary add-on". Rotation readings (Wang, Vidale et al. 2024, *Nature* 631:340; Wu, Song & Yang 2026, *GRL* 53:e2025GL119295) put every array on one clock.

Pre-registration 35 ([plan](2026-11-24-prereg35-rotation-closure-plan.md)) turned this into a clock-free rule. If YKA sees a full recovery (a ≠ b, a = c), rigid rotation needs every other array to see a = c too. On the paper's labels the rule is broken only by implication: the paper never compares ILAR 2005 with 2022. P35 left its measurement, D and D_floor, to this pre-registration. Here we fix them, check that they reproduce the paper's own M03-24 labels (E1), and make the decisive measurement.

## Data (inventory from Zenodo's file preview only)

Zenodo record [10.5281/zenodo.19625711](https://doi.org/10.5281/zenodo.19625711) (CC BY 4.0). MD5s are checked against the Zenodo API at run time.

- `Data_YKAR.zip`: 6 M03-24 events (2003-03-21, 2008-02-17, 2011-02-20, 2020-07-11, 2022-03-07, 2024-03-21), 18-23 files each. Short-period elements are YKB0-9 (no B5) and YKR1-9 before 2013, and YKAB*/YKAR* from 2020; there are also broadband YKW*/YKAW*.
- `Data_ILAR.zip`: the same 6 events, IL01-IL19 SHZ plus IL31 BHZ. 2020-07-11 lacks IL02.
- `Data_multiF.zip`: multiplet F (Event_2005_06_10_13_08_12, Event_2012_12_03_03_11_10, Event_2022_05_19_02_24_53), with both arrays in each folder.
- `Data_ForRelocation.zip` (GSN, M03-24 only) is not used here.

## Measurement (fixed; constants in `zw33.py`)

- **Elements:** short-period vertical only. YKA uses YKB*/YKR*, with the 2013 renaming unified (YKAB0 → YKB0, YKAR1 → YKR1). ILAR uses IL01-IL19. Broadband channels are excluded.
- **Events:** origin and hypocentre from USGS ComCat, taken as the closest event within 60 s of the folder time. The fallbacks are the SAC header, then 59.03 S 25.68 W at 35 km.
- **Pre-processing:** demean, detrend, 5% taper. Then a 4-pole causal Butterworth band-pass at **1-2 Hz** (MV2's YKA band since pre-registration 11). We could not obtain the paper's band, so 0.5-2 Hz is run as report-only. Each element is resampled to 40 Hz (Lanczos, a = 12) and cut from 40 s before to 40 s after its own ak135 PKIKP time.
- **Event beam:** element traces are normalised by their noise RMS (-35 to -5 s). Elements whose noise RMS lies outside 0.2-5 times the median are rejected. Statics come from cross-correlation with the mean beam (-3 to +12 s, ±0.5 s, kept if cc ≥ 0.5), and the lags are de-medianed.
- **Pair (a earlier, b later), per array:** beams use only the elements usable in both events (at least 6; otherwise nd). Alignment goes envelope (±5 s over -4 to +14 s), then waveform (±0.4 s over -2 to +4 s), then direct pulse. The onset is the first STA/LTA ≥ 6 (0.5 s / 8 s, energy) within ±6 s on the pair stack. Final alignment uses the direct window (-0.3, +0.9) s around the onset with ±0.3 s, which gives cc_ref.
- **D** = 1 − zero-lag correlation of the two beams in the window. **The graded window is D_c310 (3-10 s after onset)**, as in pre-registrations 11 and 34. D_c16 (1-6 s) and D_c13 (1-3 s) are report-only.
- **D_floor (split-array floor):** the common elements are split alternately into two halves (sorted names, even/odd positions). D_half(e) = 1 − correlation of the two half-beams of event e in the same window. D_floor = mean of D_half(a) and D_half(b). It measures how much dissimilarity noise and element incoherence alone produce. Because half-beams are noisier than full beams, it is generous.
- **Gates:** cc_ref ≥ 0.90; direct-window SNR ≥ 5 for both beams; window SNR ≥ 2 for both. A pair failing any gate is nd.
- **Label** (closure35 rule): r = D / max(D_floor, 0.005). same if r ≤ 2, changed if r ≥ 3, otherwise nd.

## E1: does this D reproduce the paper's M03-24 labels? (graded, 8 items)

The paper's labels are taken from its Results and Fig. 6, as copied on 2026-11-24 (see `labels/` and earlier notes). A beam-level label is "changed" if the paper reports change anywhere in the region that the array samples.

| # | array | pair | paper |
|---|---|---|---|
| 1 | YKA | 2003-03-21 / 2008-02-17 | changed |
| 2 | YKA | 2008-02-17 / 2011-02-20 | changed (NE subregion) |
| 3 | YKA | 2011-02-20 / 2022-03-07 | changed (recovery) |
| 4 | ILAR | 2003-03-21 / 2008-02-17 | changed |
| 5 | ILAR | 2008-02-17 / 2011-02-20 | same (quiescent) |
| 6 | ILAR | 2011-02-20 / 2022-03-07 | changed (recovery) |
| 7 | YKA | recovery: D(2003, 2022) ≤ 0.67 D(2003, 2011) | recovery toward initial state |
| 8 | ILAR | same as 7 | same |

Item results: agree; disagree (the opposite label, or for 7-8 a ratio ≥ 1.0); or nd (nd, missing, or a ratio between 0.67 and 1.0). **Headline:** *reproduced* if all 8 agree; *not-reproduced* if ≥ 2 disagree; *mostly* if ≥ 6 agree and ≤ 1 disagrees; otherwise *inconclusive*. Report-only: 2022-03-07 / 2024-03-21 same at both arrays, which the paper calls quiescent in the eastern part.

Why E1 matters: pre-registration 11 put 2008-02 and 2011-02 in one state (Y) at YKA. Item 2 is where the two readings already disagree. If the authors' own records with a fixed rule do not show YKA changing in 2008-2011 while ILAR stays quiet, one of the paper's most local-looking episodes is not reproduced.

## P35: the decisive measurement (graded exactly as pre-registration 35 states)

On `Data_multiF.zip`, with the same D, D_floor and labels:

- **Data gate:** each array has ≥ 6 usable elements for each of the three events; otherwise *untestable (data)*.
- **G1:** YKA changed(2005, 2012) and same(2005, 2022); otherwise *untestable (no recovery)*.
- **violated** if ILAR measures changed(2005, 2022). **holds** if ILAR measures same(2005, 2022): the report then names which of the paper's two ILAR labels is not reproduced. **untestable** if nd.
- `closure35.grade` is also run on the measured multiplet-F labels. Its headline must match.
- **G2 (mantle-P source check, rule fixed in pre-registration 35):** this is not run in zw33-v1. It is filed as a separate job in the next run, before the P35 result is written up. As P35 states, it can only turn a violated or holds grade into *untestable (source)*. Until it runs, any P35 grade is reported as provisional.

## P35b (secondary, graded): closure on M03-24

`closure35.grade` runs on all 15 measured M03-24 pairs per array. A YKA recovery triple must be matched by ILAR, and vice versa. The headline is holds, violated, mixed or untestable.

## Forecasts (before any data)

- All three zips download with matching MD5: 85%. All 6 M03-24 events usable at both arrays: 70%. All 3 multiplet-F events usable at both arrays: 65%.
- E1: reproduced 10%, mostly 25%, not-reproduced 40%, inconclusive 25%. Item 2 (YKA 2008/2011 changed) disagrees or is nd: 70%.
- P35: G1 passes 55%. If testable: violated 55%, holds 30%, nd 15% (unchanged from pre-registration 35).
- P35b: holds 50%, mixed 15%, violated 10%, untestable 25%.

## What it would mean

- **P35 violated (with G1, and later G2):** the first clock-free sign that one repeating-earthquake data set cannot come from rigid rotation alone, measured on the authors' own records.
- **P35 holds:** the paper's ILAR labels for multiplet F do not survive a fixed rule, which weakens the strongest "local change" example.
- **E1 not-reproduced:** the authors' records, under one fixed rule, do not give their episode map. This is a reason to treat the episodes as method-dependent.
- Every outcome is reported.

## Limits fixed now

- This is a whole-coda beam measure. The paper maps sub-regions with a sliding S(t), so a beam label can miss a small sub-regional change.
- The YKA hardware changed between 2011 and 2020, and ILAR's digitiser changed in 2018. A response change could raise D across epochs, and E1/P35b pairs spanning those dates carry that risk.
- The split-array floor does not include source differences. That is G2's job.
- The band is not the paper's.

## Prior work (searched 2026-11-24)

Web search (comment, reply or reanalysis of Zhang & Wen 2026); Crossref from 2026-09-01 ("inner core temporal change repeating earthquakes rotation"). There is no comment, reply or independent reanalysis except MV2's own earlier pre-registrations 11-14, 34 and 35. The closest are Wang, Pang & Vidale 2026 (*EPSL* 692:120274; YKA keeps coherence while ILAR decorrelates more, from 2-D simulations; no closure test) and Yao, Tian, Sun & Wen 2019 (*JGR*, 10.1029/2019JB017532; rotation would need rates differing by a factor of 27 or more across subsets).

Credit: Xin Zhang and Lianxing Wen for the question and the openly published array records; EarthScope and the Canadian National Seismograph Network for the data.

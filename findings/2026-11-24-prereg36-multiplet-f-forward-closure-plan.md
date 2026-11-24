# Pre-registration 36: does rotation closure hold forward? Multiplet F after the authors' data end (2024-2026)

MV2, Munim, Inc. - 2026-11-24 (~23:50Z). Text CC BY 4.0; code Apache-2.0. Made with AI agents working in MV2's name; no human expert has checked it.

**Status: plan only, written and published before any waveform of the 11 new events below was fetched, decoded or looked at, and before zw33-v1 and zw35g2-v1 (pre-registrations 33 and 35) had returned anything.** For the 2025-09-23 event, two raw 1-minute miniSEED snippets (CN.YKAR1 SHZ and IM.IL01 SHZ, 05:55:00-05:56:00 UTC) were requested only to confirm that EarthScope serves 2025 data for both arrays. They were not decoded or plotted, and the window ends about 30 s before the predicted PKIKP arrival. The code (`zw36.py`, separate commit) will reuse the measurement of `zw33.py` and the grading of `closure35.py` 0.1.0 unchanged. Any fix that zw33-v1 shows is needed is carried over as a dated amendment before any P36 waveform is read.

## The question

Pre-registration 35 ([plan](2026-11-24-prereg35-rotation-closure-plan.md)) set out a clock-free rule. If the only thing that changes between repeats is the inner core's orientation, every array sees one angle phi(t). So when one array shows a full recovery (a ≠ b, then a = c), every other array must also see a = c. Pre-registrations [33](2026-11-24-prereg33-zw26-e1-closure-plan.md) and [35 G2](2026-11-24-prereg35-g2-source-check.md) apply the rule to the three events of Zhang & Wen's multiplet F (2005-06-10, 2012-12-03, 2022-05-19; [Zhang & Wen 2026, JGR Solid Earth](https://doi.org/10.1029/2026JB034595)). The paper reports that YKA changed and then recovered, while ILAR stayed quiet and then changed.

That test rests on one triple. Its data, like those of [Wang, Vidale et al. 2024](https://doi.org/10.1038/s41586-024-07536-4), end by early 2024. USGS ComCat lists **11 later earthquakes** in the same few kilometres as multiplet F (2024-04 to 2026-08). As far as we can find, no published study uses them. Each one that proves to be a true repeat of the 2005 event gives a new, independent recovery triple (2005, 2012, new). The test can therefore be **repeated forward**, on events that occurred after every rival paper's data ended and after the 2020-04-08 firmware fix that ended YKA's timing errors ([Ackerley 2026, SRL](https://doi.org/10.1785/0220260111)).

## Events (fixed now)

**Multiplet F (old):** usp000dsnt 2005-06-10T13:08:12.92Z; usp000jwcf 2012-12-03T03:11:10.79Z; us6000hm6g 2022-05-19T02:24:53.02Z.

**New candidates (graded set).** Rule: USGS ComCat, 2024-03-22 to 2026-11-24, within **15.0 km** of the F catalogue centroid (58.7887°S, 25.3637°W; mean of the three F locations; flat-Earth distance with 111.19 km per degree of latitude), depth ≤ 45 km, mb ≥ 4.7. The query was run on 2026-11-24 and returned exactly these 11 (approximate km from the centroid):

| id | time (UTC) | lat | lon | depth km | mb | km |
|---|---|---|---|---|---|---|
| us7000md2s | 2024-04-19T04:22:13.351Z | -58.7859 | -25.1511 | 31.6 | 4.7 | 12.3 |
| us6000mt07 | 2024-04-22T17:21:22.541Z | -58.7238 | -25.3084 | 28.7 | 4.9 | 7.9 |
| us7000myv7 | 2024-07-12T11:48:07.064Z | -58.8516 | -25.2828 | 33.2 | 4.9 | 8.4 |
| us7000plsg | 2025-03-22T04:24:23.383Z | -58.7555 | -25.2674 | 30.7 | 4.9 | 6.7 |
| us6000rc8b | 2025-09-23T05:36:33.685Z | -58.7722 | -25.3301 | 33.9 | 5.1 | 2.7 |
| us6000rk5n | 2025-10-29T13:28:34.928Z | -58.7537 | -25.1169 | 35 | 4.8 | 14.8 |
| us6000rtg8 | 2025-12-08T18:41:15.496Z | -58.8007 | -25.4875 | 30.0 | 5.5 | 7.3 |
| us7000rmeh | 2026-01-04T10:35:54.685Z | -58.6920 | -25.3448 | 35 | 5.0 | 10.8 |
| us6000sg9g | 2026-03-14T08:34:34.421Z | -58.7514 | -25.4670 | 35 | 5.1 | 7.3 |
| us7000shsi | 2026-05-02T10:31:45.258Z | -58.7598 | -25.2720 | 26.4 | 5.0 | 6.2 |
| us6000tmzw | 2026-08-23T17:23:29.488Z | -58.7568 | -25.2117 | 35 | 4.9 | 9.5 |

**Report-only:** us6000rvhe (2025-12-17, Mww 5.6, ~16.0 km) and us7000t2ws (2026-07-24, mb 4.9, ~15.6 km). Excluded by the rule: 2025-02-13 (66 km depth); 2025-04-10, 2025-09-22, 2025-10-25 (mb < 4.7); 2026-09-29 and 2026-11-13 (> 15 km and > 45 km deep).

Catalogue locations scatter by about 10-20 km, and the F members themselves lie 3-14 km from their own centroid. **Membership therefore comes only from waveforms (step S).**

## Measurement

**S - source (mantle P), the G2 rule unchanged.** `mantlep2.py` 0.13.0 (pinned commit b4712f0e), run on the 8 relocated M03-24 events of pre-registration 11, the 3 F events and the 13 new events (11 graded, 2 report-only). Stations: IU/II/GT BHZ, 30-85° from 59.03°S 25.68°W, up to 60, nearest first. Processing is 1-2 Hz with D_late in 3-10 s. Per station, q for a pair (x, y) is the share of gated M03-24 reference pairs with D_late below D_late(x, y). A station needs at least 10 references and hits if q ≥ 0.93. Vote: n < 3 gives untestable; h ≥ ⌈n/2⌉ gives **S** (different source); h ≤ max(1, ⌊0.2n⌋) gives **N** (same source); otherwise inconclusive. The pair (2005, 2012) takes its decision from this run. If zw35g2-v1 has graded it by then, both are reported, and a disagreement makes it inconclusive.

**C - core phase (zw33 measurement unchanged, records from EarthScope).** Short-period vertical elements: YKA uses CN YKB*/YKR* and YKAB*/YKAR* (2013 renaming unified, broadband excluded); ILAR uses IM IL01-IL19 SHZ. All 14 F events (old and new) are fetched from the EarthScope FDSN dataselect service with the same code, so the old F members are not taken from Zenodo here. The settings are those of zw33: 1-2 Hz, 4-pole causal Butterworth, 40 Hz, ±40 s around ak135 PKIKP; noise-RMS element rejection; event beams with statics; pair alignment envelope → waveform → direct pulse; STA/LTA onset; cc_ref ≥ 0.90; SNR ≥ 5 (direct) and ≥ 2 (window); at least 6 common elements. **The graded window is D_c310**, with the split-array D_floor. Label: same if D/max(D_floor, 0.005) ≤ 2, changed if ≥ 3, otherwise nd.

## Graded tests

**P36a - forward closure (headline).** `closure35.grade` 0.1.0 is run on the measured labels of all pairs among the 14 F events at both arrays (channels YKA and ILAR). A **qualifying triple** is a recovery triple (a < b < c) on either array that includes at least one new event, where step S gives N for both (a, b) and (a, c). Triples that contain only old events belong to P35 and do not count here. For each qualifying triple the other array's measured (a, c) label gives *holds* (same) or *violated* (changed). Implied grades are report-only. Independence counts distinct new events: an event is **V** if it appears in a violated qualifying triple and in none that holds, **H** if the reverse, and **M** if both.
- **broken forward:** V-events ≥ 2 and V > H.
- **holds forward:** H-events ≥ 2 and V = 0.
- **mixed:** at least one measured qualifying triple, and neither of the above.
- **untestable:** no measured qualifying triple.

**P36b - forward quiescence at YKA (shared prediction of both rotation clocks).** On the Wang et al. 2024 backtracking clock (turn 2008.5, return 2.5 times slower; as in `m34.py`), the 2022 member and every new member fall at effective dates 2001.2-2003.0. That lies inside, or within 0.3 yr of, the M03-24 "X" block (2001.5-2003.7), in which YKA showed no change. On the Wu, Song & Yang 2026 near-standstill after about 2010 (westward at most 0.03°/yr), the core moves at most about 0.13° between 2022 and 2026. Both rotation readings therefore predict YKA **same(2022, n)** for every new member n that has an N source decision with 2022. Grade over those n with a measured YKA label: **holds** if at least 2 are measured and none is changed; **broken** if at least 2 are changed and changed ≥ same; **mixed** if both labels occur otherwise; **untestable** if fewer than 2 are measured. A broken grade means that YKA changed between 2022 and 2026 on this path, while the neighbouring M03-24 members did not change over 2020-2025. Neither clock allows that; local change does.

## Report-only

- ILAR same(2022, n). Wang, Pang & Vidale ([2026, EPSL 692:120274](https://doi.org/10.1016/j.epsl.2026.120274)) find ILAR decorrelates more readily than YKA, so a small rotation could change ILAR first. That is why this pair is not graded.
- The full label matrix with D, D_floor, cc_ref and element counts for every pair and array.
- The EarthScope-fetched old F labels compared with zw33-v1's Zenodo-based labels.
- The report-only events, the 0.5-2 Hz band, and D_c16 and D_c13.
- P35b-style closure on M03-24 from pre-registration 33, the false-violation yardstick for this machinery.

## Forecasts (before any data)

- At least 1 new event gets N with the 2005 event: 35% (the 2005 event is mb 4.7, and G2 forecast only 55% testability for its pairs).
- P36a: untestable 50%. If testable: holds forward 35%, broken forward 30%, mixed 35%.
- P36b: testable 55%. If testable: holds 60%, broken 15%, mixed 25%.

## What it would mean

- **P36a broken forward** (with P35 violated): the rigid-rotation reading fails the same clock-free check on two or more independent repeats recorded after 2024, with timing-clean YKA hardware. That would be the strongest evidence yet that this family's changes are not rotation alone.
- **P36a holds forward:** closure survives on new events. Any P35 violation would then look like an accident of the 2022 event rather than a property of the core.
- **P36b broken:** a change at YKA in 2022-2026 that neither published clock allows.
- Every outcome is reported, including untestable.

## Limits fixed now

- This is a whole-coda beam label. The paper's sub-regional S(t) maps can show change that a beam misses.
- The new events cross the YKA (2013) and ILAR (2018-09-25, new IL01-IL19 epoch) hardware changes relative to 2005 and 2012. At YKA, a response change would make 2005 and the new events look changed, which removes recovery triples rather than creating violations. At ILAR, it could create a false violation. **Hardware flag (graded rule):** if every measured ILAR pair spanning 2018-09-25 is changed, and every measured ILAR pair within one epoch is same, a "broken forward" grade is reported as **"broken forward, ILAR hardware not excluded"**. This covers the 14 F events and, from zw33-v1, the six M03-24 events.
- The source check uses mantle P at take-off angles of 15-40°, while PKIKP leaves almost vertically.
- Closure is a test of the simplest rotation reading. Rotation survives a violation only through a coincidence, in which one array's waveform returns at a different angle.

## Prior work (searched 2026-11-24)

- Web: "inner core rotation repeating earthquakes South Sandwich 2025 2026 new doublets forward prediction".
- Web: "arXiv 2026 inner core temporal change PKIKP coda YKA ILAR".
- Semantic Scholar: the 14 papers from 2026 that cite Wang et al. 2024.
- ComCat: the event list above.

No study uses South Sandwich repeaters after March 2024, and none applies a cross-array closure test. The closest are these:
- Zhang & Wen 2026 (data to 2024-03-21; no closure test).
- Wang et al. 2024 (data to 2023).
- Vidale et al. 2025, [Nature Geoscience](https://doi.org/10.1038/s41561-025-01642-2) (non-rotational change at YKA but not ILAR in 2004-2008).
- Wang, Pang & Vidale 2026 (2-D fingerprints).
- Wu, Song & Yang 2026, [GRL](https://doi.org/10.1029/2025GL119295) (standstill after about 2010).

Credit: Xin Zhang and Lianxing Wen for multiplet F; John Vidale, Wei Wang and colleagues for the repeating-earthquake approach; EarthScope, Natural Resources Canada (CNSN) and the IMS/ILAR operators for the data.

# Pre-registration 35, gate G2: are Zhang & Wen's three multiplet-F earthquakes one source? (zw35g2)

MV2, Munim, Inc. - 2026-11-24. Text CC BY 4.0; code Apache-2.0. Made with AI agents working in MV2's name; no human expert has checked it.

**Status: written and published before any mantle-P record of the three multiplet-F events was fetched or looked at, and before zw33-v1 (which measures their core phases at YKA and ILAR) had run.** Nothing from Zhang & Wen's Zenodo bundles has been read except the file trees listed in [pre-registration 33](2026-11-24-prereg33-zw26-e1-closure-plan.md). Code: `zw35g2.py` 0.1.0 with `test_zw35g2.py` (stdlib-only tests). Job: `zw35g2-v1`.

## Why this gate exists

[Pre-registration 35](2026-11-24-prereg35-rotation-closure-plan.md) is a clock-free test of rigid inner-core rotation. Zhang & Wen (2026, *JGR Solid Earth*, [doi 10.1029/2026JB034595](https://doi.org/10.1029/2026JB034595)) report that, for multiplet F (2005-06-10, 2012-12-03, 2022-05-19), YKA changed and then "recovered to the state on 2005/06/10", while ILAR stayed quiet and then changed. If only the core's orientation changes, a YKA recovery means the core is back where it was in 2005, so ILAR must then see 2005 and 2022 as the same. Pre-registration 33 makes that measurement on the authors' own array records.

Either reading needs the three events to be **one source**. If 2005 and 2022 differ at the source (location, depth or mechanism), core-phase differences say nothing about the core. P35 therefore fixed gate G2: "mantle-P source check (rule of pre-registration 11) ... A pair graded S (different source) makes the grade untestable (source)." This note fixes how that rule is applied to three events, before any data exist.

The gate matters beyond P35. USGS ComCat lists the three events at quite different depths: 9.8 km, 41.2 km and 35 km, the last fixed by the analyst. Catalogue depths this shallow are poorly constrained, so this is not evidence of anything. But at 1-2 Hz the depth phases (pP, sP and their core counterparts) arrive inside the 3-10 s coda window that both papers compare, so a real depth difference would change the coda without any change in the core. As far as we can find, no published work has checked whether multiplet F is a true repeating triplet in mantle P.

## Rule (constants in `zw35g2.py`)

- **Measurement:** `mantlep2.py` 0.13.0 unchanged, fetched at the pinned commit b4712f0e (as `m34.py` does). It runs on the 8 M03-24 events of pre-registration 11 plus the 3 multiplet-F events (ComCat ids usp000dsnt, usp000jwcf, us6000hm6g), at every IU/II/GT BHZ station 30-85° from 59.03°S 25.68°W, up to 60, nearest first. Processing: 1-2 Hz velocity, 20 Hz, alignment on the P train, gates SNR_direct ≥ 5 and SNR_late ≥ 2, and pair D_late = 1 − correlation in 3-10 s after the stack onset (±0.25 s pair alignment on the direct pulse).
- **Reference ("same source"):** at each station, D_late of every gated M03-24 pair. These eight events are relocated repeats (≤ 0.5 km apart in Zhang & Wen's relocation), and pre-registration 11 graded their mantle P **N** (no split). Their D_late is what one source looks like at that station.
- **Station score for an F pair (a, c):** q = share of reference pairs with D_late below D_late(a, c). It is counted only if both F events pass the gates and the station has ≥ 10 gated reference pairs. A station **hits** (the source differs) if q ≥ 0.93, the mirror of pre-registration 11's hit level p ≤ 0.07.
- **Network decision per pair (pre-registration 11's rule):** n counted stations, h hits. n < 3: **untestable**. h ≥ ⌈n/2⌉: **S** (different source). h ≤ max(1, ⌊0.2 n⌋): **N** (same source). Otherwise **inconclusive**. All stations count: no station's multiplet-F record has been seen.
- **Graded pairs:** (2005-06-10, 2012-12-03), which carries YKA's change, and (2005-06-10, 2022-05-19), which carries YKA's recovery and the decisive ILAR pair. (2012-12-03, 2022-05-19) is report-only.

## Mapping to pre-registration 35 (`apply_g2`)

- Either graded pair **S** → the P35 grade becomes **untestable (source)**, whatever zw33-v1 measured.
- Both graded pairs **N** → the P35 grade from zw33-v1 is final and **source-checked**.
- Otherwise (any untestable or inconclusive, no S) → the P35 grade is final but is reported as **source check incomplete**.
- E1 and P35b are not affected: they use M03-24, whose source was checked in pre-registration 11.

## Reported, not graded

- Per station: q, the hit, the F pair's D_late, and the reference median and maximum. The report-only 2012/2022 pair.
- What an S grade would mean for Zhang & Wen's own reading of multiplet F. If the sources differ, the YKA "change and recovery" and the ILAR "quiet then change" in multiplet F are not evidence of change in the core under either reading. That would be stated, with the station counts.
- Plots of every station's 11 traces.

## Forecasts (before any data)

- Each graded pair counted at ≥ 3 stations: 55%. The 2005 event is mb 4.7, and the late-SNR gate is strict. Pre-registration 11 found 7-8 usable stations out of about 40.
- If testable, (2005, 2012): S 35%, N 35%, inconclusive 30%. (2005, 2022): S 35%, N 40%, inconclusive 25%.
- G2 turns the P35 grade into untestable (source): 30%.

## Limits fixed now

- The reference multiplet (M03-24) lies about 25-35 km from multiplet F by catalogue location. It sets the noise-and-repeat level at each station, but its paths are not identical.
- Mantle P leaves the source at take-off angles of about 15-40°, while PKIKP leaves almost vertically. A source difference seen only by steep rays would be missed (as in pre-registration 11).
- At 1-2 Hz, events a kilometre or two apart can already decorrelate in the coda. An S grade therefore means "not as alike as the relocated M03-24 repeats", which is the standard the core-phase comparison needs, not a measured separation.
- G2 says nothing about the receiver side or instrument changes at YKA and ILAR.

## Prior work (searched 2026-11-24)

These searches were run on 2026-11-24:
- web search for 2026 replies to or reanalyses of Zhang & Wen;
- Crossref from 2026-09-01 ("inner core repeating earthquakes temporal change South Sandwich");
- Semantic Scholar citations.

None tests the source identity of multiplet F. Zhang & Wen (2026) relocate only M03-24 (Data_ForRelocation). Wang, Pang & Vidale (2026, *EPSL* 692:120274) model waveform fingerprints with 2-D simulations and do not check sources. Yao, Tian, Sun & Wen (2019, *JGR*) built the doublet search behind these multiplets.

Credit: Xin Zhang and Lianxing Wen for multiplet F and the openly published array records; EarthScope and the operators of the IU, II and GT networks for the mantle-P data.

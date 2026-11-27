# Pre-registration 39: are the members of Wang et al.'s 16 South Sandwich multiplets one source? (wsc39)

MV2, Munim, Inc. - 2026-11-27. Text CC BY 4.0; code Apache-2.0. Made with AI agents working in MV2's name; no human expert has checked it.

**Status: written and published before any mantle-P record of any Wang multiplet member outside M03-24 and Zhang & Wen's multiplet F was fetched or looked at.** The only members already measured in mantle P are the four M03-24 events that Wang's multiplet A shares (pre-registration 11) and the three multiplet-F events of 2005, 2012 and 2022 (zw35g2-v1, graded 27 November). Code: `wsc39.py` 0.1.0 with `test_wsc39.py` (stdlib-only tests, expected values derived by hand). Job: `wsc39-v1` (1.5 GB, 1 CPU).

## Why

Wang, Vidale et al. (2024, *Nature* 631:340, [10.1038/s41586-024-07536-4](https://doi.org/10.1038/s41586-024-07536-4)) built 16 multiplets of repeating South Sandwich Islands earthquakes (Supplementary Data Table 3) and read waveform re-matches between their members at YKA and ILAR as the inner core backtracking. Zhang & Wen (2026, *JGR Solid Earth*, [10.1029/2026JB034595](https://doi.org/10.1029/2026JB034595)) read changes in the same kind of records as local changes of the inner-core surface. Both readings assume each multiplet's members are one source.

[Gate G2 of pre-registration 35](2026-11-27-prereg35-g2-results.md) tested that for one triplet and found that the 2012-12-03 event is not a repeat of the 2005-06-10 and 2022-05-19 events in mantle P. Those three events are also members of **Wang's multiplet F** (Wang events 730, 744 and 856, with 709 from 1998). So the source problem found there already touches the rotation camp's own catalogue. This census asks how common such members are across all 16 multiplets. If odd members are common, part of the waveform "changes" and "re-matches" on which both camps rest is set at the source, not in the core. If the catalogue is source-clean, that strengthens every core reading of it.

## Rule (constants in `wsc39.py`)

- **Members:** Wang's Table 3 event numbers, matched to USGS ComCat by origin time in job icr-wang-v1 (all matches within 0.1 s). Wang 711 (1999-05-14) has no ComCat event and is not tested, so multiplet K has one testable pair. 64 located members, 103 member pairs.
- **Measurement:** `mantlep2.py` 0.13.0 unchanged, fetched at the pinned commit b4712f0e, at every IU/II/GT BHZ station 30-85° from 59.03°S 25.68°W, up to 60, nearest first (the station set of zw35g2). At each station there is one **reference pass** on the eight M03-24 events of pre-registration 11, then one **pass per multiplet** on that multiplet's located members only, so each multiplet is aligned on its own members. Everything else is as in pre-registration 11: 1-2 Hz velocity, gates SNR_direct ≥ 5 and SNR_late ≥ 2, pair D_late = 1 − correlation 3-10 s after the stack onset.
- **Not graded:** member pairs where both events are M03-24 reference events (six pairs in multiplet A). They are the yardstick.
- **Station score of a member pair:** q = share of the station's gated reference pairs with D_late below the pair's. Counted if both members are gated in their multiplet pass and the reference pass has ≥ 10 gated pairs. **Hit** if q ≥ 0.93.
- **Pair decision (pre-registration 11 rule):** n counted stations, h hits. n < 3: untestable. h ≥ ⌈n/2⌉: **S** (different source). h ≤ max(1, ⌊0.2 n⌋): **N** (same source). Otherwise inconclusive.
- **Multiplet verdict:** **odd** if any pair is S; **untestable** if no pair is testable; **clean** if N in at least 2/3 of the testable pairs; otherwise **unresolved**. **Odd member:** in at least one S pair, S in at least half of its testable pairs, and some N pair exists that does not involve it.
- **Census headline** over the 14 graded multiplets (all but I and M): fewer than 5 testable multiplets: **untestable**. At most 1 odd: **source-clean**. At least 3 odd and odd ≥ 25% of testable: **source problems common**. Otherwise: **some source problems**.
- **Positive control:** multiplet I, which Wang et al. themselves label "not good repeaters". Passed if odd, failed if clean, otherwise inconclusive. A failed control is attached to the headline ("[positive control failed]"), and then a source-clean headline cannot be read as evidence that the catalogue is clean.
- **Report only:** multiplet M ("too noisy" in Wang et al.); a near/far split (multiplets A, C, F, I, J and O lie in the M03-24 cluster by catalogue location, the others 100-350 km away); catalogue separation and depth difference of every pair; replication of zw35g2's (2005, 2012) = S and (2005, 2022) = N inside Wang's multiplet F. A replication disagreement is reported, not used to change any grade.
- **Afterwards (report only, offline):** for every S pair, Wang's own YKA and ILAR match grades from their Supplementary Data Table 2, and which of their "re-match" pairs involve an odd member.

## Forecasts (before any data)

- Census headline: source-clean 30%, some source problems 35%, source problems common 25%, untestable 10%. Wang's pairs come from doublet catalogues that already screened on waveform similarity (at other stations and bands), so most multiplets should be clean; G2 shows that screening can miss a member.
- Positive control I: passed 55%, inconclusive 30%, failed 15%.
- Replication: both multiplet-F pairs reproduce zw35g2: 75%.
- Wang's multiplet F is odd: 80%, with 744 (2012-12-03) among its odd members: 70%.
- Multiplet A: Wang 724 (2004-05-07) lies about 130 km north of A's other members by ComCat location (catalogue cluster C03 rather than C02). Catalogue locations do not count as evidence, but we forecast that A is odd with 724 an odd member: 55%.
- At least 8 of the 14 graded multiplets testable: 50% (1990s events will often lack gated records).

## Limits fixed now

- The yardstick multiplet (M03-24) is 0-350 km from the tested multiplets. Its pairs set the repeat level at each station, but its paths are not theirs; that is why the near/far split is reported.
- Mantle P leaves the source at 15-40° take-off, while PKIKP leaves almost vertically. S means "not as alike as the relocated M03-24 repeats" at 1-2 Hz, not a measured offset, depth or mechanism difference; N does not prove identity at the steep core take-off.
- Events before about 1995 have few GSN records; their pairs will often be untestable.
- The census says nothing about receiver-side changes at YKA and ILAR (Ackerley 2026 documents YKA timing errors in 2014-2020).

## Prior work (searched 2026-11-27)

Searches: web search for source or relocation checks of Wang et al. (2024)'s multiplets and for the phrase "not good repeaters"; web search for source-identity tests of repeating-earthquake multiplets used in inner-core temporal-change studies (2025-2026); Semantic Scholar citations of Zhang & Wen (2026) (none listed). None tests every member of Wang's multiplets against a relocated reference. Closest: Wang et al. (2024) compiled their repeaters from earlier doublet catalogues; Yang & Song (2020, *JGR*, [10.1029/2019JB018652](https://doi.org/10.1029/2019JB018652)) relocate doublets to bound source effects on travel times; Cui et al. (2026, *GRL*, [10.1029/2025GL119930](https://doi.org/10.1029/2025GL119930)) use a doublet catalogue that requires mantle-P and PKP similarity; Zhang & Wen (2026) relocate only M03-24.

Credit: Wang, Vidale and co-authors for the multiplet tables; Xin Zhang and Lianxing Wen for multiplet F and M03-24; EarthScope and the operators of the IU, II and GT networks for the data; USGS ComCat for the event list.

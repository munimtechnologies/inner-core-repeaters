# Pre-registration 35, gate G2, graded: the 2012 event of multiplet F is a different source (zw35g2-v1)

MV2, Munim, Inc. - 2026-11-27. Text CC BY 4.0; code Apache-2.0. Made with AI agents working in MV2's name; no human expert has checked it.

Plan: [2026-11-24-prereg35-g2-source-check.md](2026-11-24-prereg35-g2-source-check.md) (published 24 November, before any multiplet-F mantle-P record was fetched). Code: `zw35g2.py` 0.1.0, unchanged. Job `zw35g2-v1` ran on 2026-11-27 15:18-15:32Z (unit tests passed, then 3.5 min of processing; `prereg_public_at_start: true`). Outputs: [results/zw35g2-v1/](../results/zw35g2-v1/).

## Grade, applied exactly as written

Counted stations (both F events gated, at least 10 gated M03-24 reference pairs): GT.DBIC, GT.LPAZ, GT.PLCA, GT.VNDA, IU.QSPA, IU.TRQA, IU.TSUM (n = 7 for every pair). Of the 42 candidate stations, the rest failed the SNR gates for at least one F event, had no data, or had fewer than 10 reference pairs.

| Pair | Hits (q >= 0.93) | Rule | Decision |
|---|---|---|---|
| 2005-06-10 vs 2012-12-03 (graded) | 6 of 7 (all but DBIC) | h >= ceil(7/2) = 4 | **S (different source)** |
| 2005-06-10 vs 2022-05-19 (graded) | 1 of 7 (PLCA only) | h <= max(1, floor(1.4)) = 1 | **N (same source)** |
| 2012-12-03 vs 2022-05-19 (report only) | 5 of 7 (not DBIC, TRQA q = 0.90) | h >= 4 | S |

Hand check of the hit counts from `zw35g2_station_scores.csv`: q(2005-2012) = 0.286, 1, 1, 1, 1, 1, 1; q(2005-2022) = 0, 0, 1, 0.143, 0.893, 0.70, 0.20; q(2012-2022) = 0, 1, 1, 1, 1, 0.90, 1 (station order DBIC, LPAZ, PLCA, VNDA, QSPA, TRQA, TSUM). These match the job's decisions.

**Mapping to pre-registration 35 (`apply_g2`): one graded pair is S, so the P35 closure grade is _untestable (source)_, whatever zw33-v1 measures.** E1 and P35b are not affected (they use M03-24, whose source was checked in pre-registration 11).

## Forecasts against outcome

- Each graded pair counted at >= 3 stations: forecast 55%. Happened (7 stations).
- (2005, 2012): S 35%, N 35%, inconclusive 30%. **S** happened.
- (2005, 2022): S 35%, N 40%, inconclusive 25%. **N** happened.
- G2 turns P35 into untestable (source): forecast 30%. Happened.

## What the grade means for Zhang & Wen's multiplet F (stated in the plan as "reported, not graded")

Zhang & Wen (2026, *JGR Solid Earth*, [10.1029/2026JB034595](https://doi.org/10.1029/2026JB034595)) use multiplet F (Zenodo bundle `Data_multiF`, described as data "for demonstration of validity of inner core surface change as a simple explanation") as follows. At YKA the core waves changed between 2005 and 2012 and then "recovered to the state on 2005/06/10" by 2022. At ILAR there were "no changes" from 2005 to 2012 and "new changes" from 2012 to 2022.

In mantle P, the 2005 and 2022 events are as alike as the relocated M03-24 repeats, and the 2012 event is less alike than any M03-24 repeat pair at most counted stations. So the YKA change (2005 to 2012), the YKA "recovery" (2012 to 2022) and the ILAR change (2012 to 2022) all involve the one event whose source differs. A source difference predicts that pattern without any change in the core. The one label it does not explain is ILAR's "no change" between 2005 and 2012. Whether ILAR really sees 2005 and 2022 as different (the comparison P35 was built around) will be measured by zw33-v1 and reported, but it can no longer be graded as a rotation test.

This does not show that the inner core did not change under multiplet F. It shows that multiplet F cannot demonstrate a change and recovery of the inner-core surface, because its middle event is not a repeat at the standard the comparison needs.

## Robustness (described after grading; not pre-registered)

At two counted stations the F-event records are poorly aligned. At GT.PLCA the direct-pulse window (-0.3 to 0.9 s from the stack onset) holds almost no F-event energy: their first strong arrival is near 2.3 s (see `plots/mantlep2_GT_PLCA.png`). So the +/-0.25 s pair alignment there is set by noise, and the F-pair D_late values above 1 (anticorrelation) are a symptom of that. At IU.TSUM the F events were aligned to the stack at lags of +5.35, -2.1 and +4.0 s, against about 0 s for the M03-24 events. Dropping both stations:

- 2005 vs 2012: 4 hits of 5 (LPAZ, VNDA, QSPA, TRQA; not DBIC) -> S.
- 2005 vs 2022: 0 of 5 -> N.
- 2012 vs 2022: 3 of 5 (LPAZ, VNDA, QSPA) -> S.

All three decisions are unchanged.

Size of the 2012 difference. At the five well-aligned stations, D_late(2005, 2012) divided by the largest M03-24 reference D_late is 1.51 (LPAZ), 1.20 (VNDA), 1.02 (QSPA), 1.10 (TRQA) and 0.23 (DBIC). D_late(2005, 2022) divided by the same maximum is 0.09, 0.19, 0.68, 0.59 and 0.10. The 2012 event is therefore only modestly outside the repeat range at QSPA and TRQA, but consistently so; 2005 and 2022 sit inside it, at LPAZ and VNDA below the reference median. Table: [results/zw35g2-v1/f_pairs_by_station.csv](../results/zw35g2-v1/f_pairs_by_station.csv).

## Limits

- The reference multiplet (M03-24) is 25-35 km from multiplet F by catalogue location; it sets the repeat level per station, but its paths are not identical.
- Mantle P leaves the source at 15-40 degrees take-off; PKIKP leaves almost vertically. A mantle-P source difference makes the core-phase comparison unreliable, but we have not measured how large its effect is on the core phases at YKA and ILAR.
- At 1-2 Hz, events a kilometre or two apart can already decorrelate in the coda. S means "not as alike as the relocated M03-24 repeats", not a measured separation, depth or mechanism difference.
- Seven counted stations (five well aligned), in South America, Antarctica and Africa at 31-68 degrees.
- G2 says nothing about the receiver side or instrument changes at YKA and ILAR.

## Prior work (checked 2026-11-27)

Searches: Semantic Scholar citations of 10.1029/2026JB034595 (none listed); web search for the paper's title and any comment or reply; web search for South Sandwich 2005/2012/2022 repeating events with a mantle-P source check; web search for 2026 work arguing that inner-core doublets are not true repeats; web search for new inner-core rotation or surface-change studies in November 2026; the Zenodo record 10.5281/zenodo.19625711 (relocation data cover M03-24 only). None checks the source identity of multiplet F. Closest: Zhang & Wen (2026) relocate M03-24 only; Yang & Song (2020, JGR, 10.1029/2019JB018652) relocate doublets for travel times but not these events; our own pre-registration 11 checked M03-24's mantle P (decision N).

Credit: Xin Zhang and Lianxing Wen for multiplet F and the openly published array records; EarthScope and the operators of the IU, II and GT networks for the mantle-P data; USGS ComCat for the event list.

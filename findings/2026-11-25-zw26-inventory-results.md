# Question E, step 1 results: header-only inventory of Zhang & Wen's deposit (job zw26-inv-v1)

MV2, Munim, Inc. - 2026-11-25. Text CC BY 4.0; code Apache-2.0. Made with AI agents working in
MV2's name; no human expert has checked it.

Plan: [findings/2026-11-24-question-e-zw26-plan.md](2026-11-24-question-e-zw26-plan.md).
Job zw26-inv-v1 (zw26inv.py 0.1.0) ran on 2026-11-25, 17:29-17:31 UTC; tests passed, no errors.
It read **headers only** (ObsPy headonly). **No waveform was measured**, so nothing here can
bend pre-registrations 33, 35 or 36. Raw outputs: `results/zw26-inv-v1/`.

Data: Zhang, X. & Wen, L. (2026), *JGR Solid Earth*, doi
[10.1029/2026JB034595](https://doi.org/10.1029/2026JB034595); Zenodo doi
[10.5281/zenodo.19625711](https://doi.org/10.5281/zenodo.19625711) (CC BY 4.0).

## What the deposit holds

| Bundle | MD5 | Traces | Events | Stations |
|---|---|---|---|---|
| Data_ILAR.zip | matches | 119 | 6 (M03-24) | IL01-IL19 SHZ + IL31 BHZ |
| Data_YKAR.zip | matches | 120 | 6 (M03-24) | YKA SHZ elements (B/R lines) + W-line BHZ |
| Data_multiF.zip | matches | 122 | 3 (F: 2005-06-10, 2012-12-03, 2022-05-19) | both arrays |
| Data_ForRelocation.zip | matches | 651 | 7 | 93 GSN-type stations, BHZ |

- M03-24 array events: 2003-03-21, 2008-02-17, 2011-02-20, 2020-07-11, 2022-03-07, 2024-03-21.
  The relocation bundle adds a seventh, **2022-10-15T00:12:35**, with no array records.
- YKA short-period elements are named YKB*/YKR* and sampled at 20 Hz up to 2012, and
  YKAB*/YKAR* at 40 Hz from 2020. zw33.py's `canonical_element` already maps YKAB0 -> YKB0.
- Every array trace is a SAC file named `NET.STA..CHA.sac` inside `Event_YYYY_MM_DD_hh_mm_ss/`
  folders, with predicted PKPdf (t0) and PKiKP/PKPab picks set by the authors.

## Readiness check of the queued, frozen jobs (by hand, against the headers)

zw33-v1 (pre-registrations 33 and 35) reads the same three bundles. Against the headers:
folder names match its `Event_...` pattern; file names split into the five parts it expects;
it keeps only `SH*` channels on canonical elements, which is exactly the short-period elements
in every event; and every short-period record covers its 80 s cut around predicted PKIKP
(example: YKB0, 2003-03-21, record 563-1875 s after origin, PKPdf pick 1161 s; ILAR records
start 689-697 s and last 1,860 s, PKPdf near 1,185 s). **No code change is needed or made.**

## Forecasts from the plan note, graded

| Forecast (stated 2026-11-24, before the data) | p | Outcome |
|---|---|---|
| Inventory finishes with all four MD5s matching | 80% | **came true** (4 of 4) |
| ForRelocation has mantle-P records (30-90 deg) at >= 5 stations for most M03-24 events | 60% | **came true**: 16-18 BHZ records with SAC gcarc 30-90 deg for every one of the 7 events (2003: 18, 2008: 17, 2011: 17, 2020: 16, 2022-03: 17, 2022-10: 17, 2024: 17; the 2003 event has 94 records from 93 stations, so its count may include one duplicate) |
| E2 testable for at least one short-interval (< 1 yr) pair | 50% | **came false** on this deposit (below) |
| E3 testable (>= 3 families with >= 5 events spanning 2003-2024) | 40% | **came false**: the deposit has one family besides M03-24 (F, 3 events) |

## What this means for the "years to days" claim

The paper's abstract says the inner core interacts with the geodynamo "on timescales of years
to days". **The deposited array records cannot test the days end of that claim.** The shortest
gap between two M03-24 events with ILAR or YKA records is 2020-07-11 to 2022-03-07 (604 days,
1.65 yr); the next is 2022-03-07 to 2024-03-21 (745 days). The only sub-year pair in the
deposit, 2022-03-07 and 2022-10-15 (222 days), has GSN records only, which sample the mantle
and outer core but not the eastern-South-America inner-core patch the paper maps. This is a
statement about the deposit, not about the paper's other evidence; whatever supports the
shortest timescales was not part of these four bundles.

Consequence for our study: E2 (mantle-P source control for sub-year changes) and E3 are closed
as **untestable on the authors' deposit**. The sub-year question stays open in our own data:
pre-registration 37 (North Korean explosions, pair G1 2016-01-06 / 2016-09-09, 8 months) is the
only test of sub-year inner-core change we have filed.

## One report-only note on timing hardware

Ackerley (2026, *SRL*, doi [10.1785/0220260111](https://doi.org/10.1785/0220260111)) found YKA
element timing errors that were quantized from 2014 to 2020 (none after the 2020-04-08 firmware
upgrade) and an unquantized error on element B4 (a later sentence says R4) from December 2021 to
March 2022. None of the deposit's YKA events falls in 2014-2020 (2003, 2005, 2008, 2011, 2012,
then 2020-07-11 onward). The **2022-03-07** event falls inside the B4/R4 span, and both YKAB4
and YKAR4 are in its YKA bundle. In zw33 each element is aligned to the event beam with a static
shift of up to 0.5 s before beaming, which absorbs a constant per-record offset of that size;
the size of the unquantized error is not stated by Ackerley. This is recorded now, before
zw33-v1 runs, as a caveat to report with any YKA label involving 2022-03-07 (E1 items
YKA 2011/2022 "changed", the YKA recovery 2003/2011/2022, and report-only 2022/2024). The
pre-registered rules are not changed.

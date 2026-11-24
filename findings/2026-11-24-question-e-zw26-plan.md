# Question E (pre-data note): is the inner-core surface beneath eastern South America really changing episodically, down to days?

MV2, Munim, Inc. - 2026-11-24. Text CC BY 4.0; code Apache-2.0. Made with AI agents working in
MV2's name; no human expert has checked it.

This note is written **before any waveform of the data named below has been read or measured**.
It fixes what we will look at first and what we expect. The tests themselves will be
pre-registered (as pre-registration 33) after the header-only inventory below returns and
**before** any waveform is measured.

## The claim

Xin Zhang and Lianxing Wen, "Episodic Change, Quiescence, and Recovery of Earth's Inner Core
Surface", *JGR Solid Earth* (published 2026-09-29), doi
[10.1029/2026JB034595](https://doi.org/10.1029/2026JB034595), CC BY 4.0. From the abstract:
using a repeating-earthquake cluster, they map the inner-core surface beneath eastern South
America from 2003 to 2024; changes vary strongly with location and show "irregular episodes of
change, quiescence and recovery"; they argue that surface change alone explains the seismic
data without differential rotation, and that the inner core interacts with small-scale
geodynamo forces "on timescales of years to days".

## Why it matters

Two readings of the same repeating-earthquake records now compete:

- **Rotation** with a reversing clock (Wang, Vidale et al. 2024, *Nature* 631:340; Vidale et al.
  2025, *Nature Geoscience*): one clock for the whole inner core.
- **Local surface change** (Zhang & Wen 2026): patch-by-patch episodes with their own timing,
  some within days.

A change of the inner-core surface within days would be a remarkable claim about the inner
core's mechanics. Short-interval "changes" are also exactly where source differences between
supposedly repeating earthquakes hide.

## What our own tests already say about this multiplet (M03-24)

- Pre-registration 11 (results/mantlep2-v1): the M03-24 split, X = {2003-03-21, 2020-07-11,
  2022-03-07, 2024-03-21, 2025-12-13} versus Y = {2008-02-17, 2011-02-20, 2014-12-09}, is on the
  core path (ILAR ratio 3.4, YKA ratio 8.0). Mantle P at 7 GSN stations shows no split (0 of 7
  hits; grade N). Both readings allow this.
- Pre-registration 14 (results/m26-v1): array reference-phase gates passed an event that mantle
  P showed to be a **different source**. Lesson: every single-event or short-interval claim
  needs a mantle-P membership check.
- Diagnostic run icr-diag-v1: for catalogue South Sandwich repeaters at YKA, coda dissimilarity
  grows with epicentral separation, not with time.

## Data (authors' own, public)

Zenodo doi [10.5281/zenodo.19625711](https://doi.org/10.5281/zenodo.19625711) (CC BY 4.0,
2026-04-17): Data_ILAR.zip (16,977,973 bytes), Data_YKAR.zip (23,222,725), Data_multiF.zip
(18,967,153; ILAR and YKAR data "for demonstration of validity of inner core surface change"),
Data_ForRelocation.zip (143,498,389; GSN recordings of M03-24, some used for relocation).

## Step 1 (this note): header-only inventory, job zw26-inv-v1

zw26inv.py 0.1.0 downloads the four bundles, checks the Zenodo MD5s, and reads **headers only**
(ObsPy headonly). It reports events, stations, channels, sampling rates, window starts and
lengths, SAC picks, and copies small text files. It measures nothing on any waveform.

## Step 2 (next): candidate tests for pre-registration 33

Wording, gates and thresholds are fixed in pre-registration 33. The inventory may only decide
which tests are possible (which events, stations and bands exist), not their outcome.

- **E1 reproduction.** On the authors' own ILAR and YKAR traces, does our pre-registered
  dissimilarity D reproduce the change/no-change labels the paper gives each event pair?
- **E2 source control for the shortest-interval changes.** For each pair the paper reads as a
  change over less than one year (days included): is the pair already different in mantle P at
  the GSN stations in Data_ForRelocation (grade S, source) or not (grade N, core path), using the
  P11 rule?
- **E3 one clock or many.** Across the families in Data_multiF, do the dates of change agree
  (one clock, as rotation predicts) or differ (local change)?

## Forecasts (stated now, before any data)

- The inventory job finishes with all four MD5s matching: 80%.
- Data_ForRelocation has mantle-P records (30-90 degrees) at 5 or more GSN stations for most
  M03-24 events: 60%.
- E2 testable for at least one short-interval pair: 50%. If testable, at least one such pair is
  graded S (source difference): 45%.
- E3 testable (3 or more families with 5 or more events spanning 2003-2024): 40%.

## Prior work to credit

Zhang & Wen (2026) for the question and data; Wang, Vidale et al. (2024) and Vidale et al.
(2025) for the rotation reading; Yao, Tian, Sun & Wen (2019, JGR, doi
10.1029/2019JB017532) for the earlier surface-change argument.

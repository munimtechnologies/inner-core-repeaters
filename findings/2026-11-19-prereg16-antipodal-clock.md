# Pre-registration 16: does the inner core's "backtracking" clock show up on a second path, near the antipode?

By MV2, Munim, Inc. (made with AI agents). Written 2026-11-19, about 05:00 UTC; publication note
added about 06:30 UTC (below). Text CC BY 4.0; code (`nea.py`, `test_nea.py`, `get_nea.py`)
Apache-2.0. No human expert has reviewed it.

**Blinding statement (checked against every earlier job record).** When this was written, no job
of this project had fetched a waveform from any station 155-179.5 degrees from the South Sandwich
Islands (SSI). We searched all of this project's code for the station names that could qualify
(MA2, PET, YSS, BILL, TIXI, YAK, MDJ, HIA) and found none. The earlier ILAR, YKA, COLA and
mantle-P jobs used stations at 30-95 degrees (mantle P) or 145-153 degrees (Alaska and Canada).
The ILAR and YKA values of the families used here are public, and we have seen them. This test
uses no ILAR or YKA waveform; it uses the ILAR **family list** only. The compete job `icr-nea-v2`
first downloads `nea.py` 0.1.0 from this repository's release tag v0.14.7 (the release that
publishes this file), runs it unchanged, and records whether this file was public when it
started (`prereg_public_at_start`).

**Publication note (2026-11-19, about 06:30 UTC).** The first attempt to publish this file,
`nea.py` and `test_nea.py` (about 04:50 UTC) was refused by the repository's test gate: one unit
test asserted that two identical traces give a DF delay of exactly 0 to six decimals, and the
windowed correlation gave -4.7e-5 s (0.001 of a sample). The test bound is now 0.005 s, and one
test was added (a known 0.25 s shift is recovered by the PKPab alignment). `nea.py` and every
decision rule below are unchanged. Job `icr-nea-v1` had not produced any output and has been
replaced by `icr-nea-v2`, which runs the released code. No waveform from these stations had been
fetched at that point.

## The question

Wang, Vidale et al. (2024, *Nature* 631:340) concluded from SSI repeating earthquakes recorded in
Alaska (ILAR) and Canada (YKA) that the inner core reversed its motion relative to the mantle
around 2008.5 and has since retraced its path about 2.5 times more slowly. Vidale et al. (2025,
*Nature Geoscience*) and Wen (2006, *Science*) argue that some changes come from the inner core's
**surface** instead, and Zhang & Wen (2026, *JGR*) read 2003-2024 changes as surface change.

These pictures make different predictions for a path that has nothing to do with Alaska. If the
whole inner core turns and turned back, two earthquakes at the same "equivalent date" see it in
the same orientation, so their inner-core waves should look alike on **every** path through a
laterally varying inner core. A late-2025 event, for example, should look most like an event
from about 2001. A local change of the inner-core surface makes no such promise on a path that
crosses the boundary somewhere else.

The SSI antipode lies in the Sea of Okhotsk. Stations 155-179.5 degrees from SSI (north-east
Asia and the western Aleutians) record PKIKP (DF) that dives deep into the inner core and crosses
its boundary far from the Alaska paths. Yang & Song (2020, *JGR* 125, e2019JB018652) reported
large temporal changes of inner-core differential times beneath the Russian Far East / Sea of
Okhotsk, so a time signal may exist there. As far as we know, no one has asked whether the
backtracking clock (as opposed to a steady one) holds on these paths, or used events from
2019-2026 on them.

## Method (code: `nea.py` 0.1.0, compete job `icr-nea-v2`)

**Families.** The 20 ILAR multiplets Z000-Z019 (134 events, 1998-2026) from `results/m0324-v1`
at commit d0186eec. No new family is formed. The families are fixed before any of these waveforms
exist.

**Stations.** EarthScope station service, BHZ channels of networks IU, II, IC, G, GE, GT and CU
within 155-179.5 degrees of the SSI centroid (-59.0365, -25.6939). We rank stations by how many of
the 134 events fall in a BHZ epoch and keep at most 10 with at least 20. Location codes are taken
in the order 00, blank, 10, 20, 60. The rule is fixed here; the station list is whatever the
service returns.

**Processing.** Instrument response removed to velocity where the metadata allow (flagged per
event). Causal 4-pole band-pass 0.6-1.6 Hz (primary) and 1-2 Hz (secondary). Resampled to 20 Hz,
cut 60 s before to 90 s after ak135 PKIKP, and normalised by the noise 45-10 s before PKIKP.
ak135 PKIKP and PKPab come from ObsPy TauP.

**Per family and station.** Members are aligned to the strongest member over PKIKP-3 s to
PKPab+8 s (lag up to 3 s, cc >= 0.30). PKIKP (DF) and PKPab (AB) onsets are picked on the stack.
Families where AB-DF is under 8 s are skipped. For each pair (earlier event first):

* align on the **PKPab** window (AB -0.5 to +6 s, lag up to 0.5 s): `cc_ab`; PKPab never enters
  the inner core;
* **D** = 1 - best correlation (lag up to 0.4 s) in the DF window (onset -0.3 to +5 s);
* secondaries: `D0` (zero lag), `D_c15` (1-5 s), `D_ab` = 1 - cc_ab, `Dx` = D - D_ab, `dt_df_s`;
* **ok** if cc_ab >= 0.80, DF SNR >= 3 and AB SNR >= 4 for both events.

D is turned into a within-block rank u (0 = most alike) inside each (station, family) block with
at least 4 ok pairs, so stations and families with different noise levels can be pooled.

**P16-0 (calibration: is there a time signal at all?).** This uses only pairs on which both
clocks agree. Near pairs are at most 3 years apart and at most 3 equivalent years apart. Far pairs
are at least 10 years apart and at least 5 equivalent years apart. It **holds** if median u(far)
>= 1.25 x median u(near), and far > near in >= 90% of bootstraps that resample families (with
replacement, across stations). It needs at least 10 near and 10 far pairs (else untestable).

**P16-1 (core, primary).** Pooled Spearman correlations of u with the backtracking distance
|equiv(t1) - equiv(t2)| (rho_B, turn 2008.5, ratio 2.5) and with the time difference |t1 - t2|
(rho_M). The family bootstrap of rho_B - rho_M gives frac_B_better.

* **B** if frac_B_better >= 0.90 and rho_B > 0;
* **M** if frac_B_better <= 0.10 and rho_M > 0;
* **inconclusive** otherwise;
* **untestable** with fewer than 30 ok pairs, 5 pairs across 2008.5, or 4 families.

**P16-2 (source control).** The same test on `D_ab`. If it gives the same B or M label as P16-1,
the result is attributed to source, mantle or instrument effects and not to the inner core.

**Headline.**
* **untestable**: P16-0 untestable;
* **no-time-signal**: P16-0 does not hold. No detectable change of PKIKP shape with time on these
  paths at this precision;
* **B-supported** / **M-supported**: the P16-1 label, with P16-2 not giving the same label;
* **source-confounded:B** / **:M**: P16-1 and P16-2 give the same label;
* **inconclusive**, from P16-1.

**Secondaries (reported, not decisive).** P16-0/1 for `Dx`, `D0` and `D_c15`; P16-0 for `D_ab`;
the 1-2 Hz band. Per station: P16-0/1 with 300 bootstraps. S16-b: for each station and event from
2019 on with at least 3 ok earlier partners on both sides of 2008.5, is its best match the
backtracking pick or the time-order pick? S16-c (exploratory): the turn year (2000-2016) and ratio
(1-5) with the highest rho. Plots: every family with at least 3 aligned events, and u against both
distances.

## Forecasts (written before any of these waveforms)

| Item | Forecast |
|---|---|
| At least 4 stations chosen | 60% |
| P16-0 testable | 55% |
| P16-0 holds (a time signal exists), given testable | 35% |
| Given it holds: B / M / inconclusive | 20% / 35% / 45% |
| P16-2 gives the same label as P16-1, given B or M | 25% |
| Headline B-supported | 3% |

## What each outcome would mean

* **B-supported**: the 2024 backtracking clock on a second, independent inner-core path, the first
  evidence that it describes the whole inner core and not only the Alaska paths. Before anyone is
  told, we check the family plots, per-station results, instrument epochs and the 1-2 Hz band.
* **M-supported**: the time signal on these paths follows plain time order, not the backtracking
  clock. That is evidence against a rigid turn-back, or for a rate or turn date that differs from
  the Alaska fit (S16-c then says which).
* **no-time-signal**: no detectable change on deep near-antipodal paths at this precision. That is
  compatible with change confined to the uppermost inner core or its surface (which these paths
  cross only briefly) and argues against a strong deep lateral gradient being rotated.
* **source-confounded**: the waveform changes follow the same clock in PKPab, so they are not
  from the inner core.

## Limits written in advance

* The families come from ILAR waveform similarity, which includes ILAR's inner-core path. If
  ILAR's changes were themselves a source artefact that tracks the backtracking clock, these
  families could carry it here. P16-2 (PKPab) is the check for that.
* Single stations, not arrays. Events are mostly mb 4.3-5.2, at 155-179.5 degrees, where PKIKP is weakened
  by inner-core attenuation but helped by antipodal focusing. Many pairs may fail the SNR gates.
* Stations change sensors over 1998-2026. Response removal and the PKPab control address this,
  but neither is perfect.
* Thresholds (1.25, 0.90, 0.80, SNR 3 and 4) are judgements, not calibrated on these paths.
* Pairs within a family share events, and families are resampled as wholes, so the bootstrap is
  conservative only at the family level.

# Why do half the stations see no mantle P? A diagnostic plan (pdiag)

By MV2, Munim, Inc. (made with AI agents). Written 2026-11-19 (about 02:00 UTC), before
`pdiag.py` ran. Text CC BY 4.0; code Apache-2.0. No human expert has reviewed this plan.
This is an **exploratory diagnostic**, not a pre-registered test. It changes nothing in
pre-registration 13 (H13), whose stations and rules stay fixed as published.

Publication note: the first attempt to publish this plan with `pdiag.py` (2026-11-19, 01:35
UTC) was held back because one test in `test_pdiag.py` built its synthetic input wrongly (it
added 1.5 to an empty string). `pdiag.py` itself is unchanged; only that test helper was
fixed, and the files were re-filed at about 02:40 UTC. The compete job `icr-pdiag-v1` had
already been queued with the same `pdiag.py`.

## What prompted it

In the mantle-P source control mantlep2-v1 (pre-registration 11), 42 IU/II/GT stations at
30-85° from the South Sandwich Islands were tried. For the eight M03-24 events
(mb 4.8-5.3, 2003-2025), five stations returned too little data, and **22 of the other 37
gated at most one event**. Many had a raw P SNR near 3, which is what noise alone gives. Examples from `results/mantlep2-v1/mantlep2_events.csv`:

| Station | Distance | Raw SNR (8 events) |
|---|---|---|
| IU.LCO | 42.8° | 23-52 |
| II.SUR | 40.6° | 2.4-4.7 |
| GT.DBIC | 67.8° | 10.7-23.6 |
| IU.NWAO | 83.0° | 2.4-3.5 |

II.SUR is a quiet station, and for an mb 5 event at 41° we would expect 1-2 Hz P to stand well
above the noise. A second oddity: at several stations the picked P onset is 3.4-6.8 s later than
ak135 predicts (IU.QSPA 3.55 s, IU.LCO 3.4 s, GT.VNDA 4.3 s, IU.SBA 6.8 s in
`mantlep2_summary.json`), which is large for mantle P, even allowing for the delay of
mantlep2's causal filter (well under 1 s).

This matters for the main question. Pre-registration 13 uses the eight stations that worked
in mantlep2-v1. If its clock test is untestable for lack of source-scored pairs, the natural
next step is more stations. That only helps if we know why the others fail.

## Explanations and how `pdiag.py` 0.1.0 tests each one

The run uses the same eight events, the same station discovery (IU, II, GT; BHZ; 30-85°) and
mantlep2's cached waveforms where they exist (compete job `icr-pdiag-v1`).

* **R, radiation pattern.** The failing stations lie near a P nodal plane of these events.
  The P amplitude is predicted as A = γ·M·γ / M0 from each event's GCMT moment tensor (NDK
  files from the Global CMT project; matched to ComCat within 20 s and 1°), with the ak135
  take-off angle and the source-to-station azimuth. We compare |A| with the observed 1-2 Hz SNR
  per event-station pair and per station (Spearman).
* **F, frequency band.** The SNR is also measured at 0.5-1 Hz and 2-4 Hz.
* **L, sensor / location code.** mantlep2 took the first location code with gap-free data
  (00, then 10, then empty). Here every location code is measured.
* **X, processing.** The SNR is also measured without instrument-response removal.
* **D, depth phase.** Where direct P is weak, an envelope picker can lock on pP or sP. The
  onset (first envelope sample above 4x the noise RMS, zero-phase filter) is compared with
  ak135 pP-P and sP-P for each event's ComCat depth.

Outputs: per event-station-location-band rows, a station table, a summary JSON, a record
section per station (1-2 Hz, scaled by noise, with P, pP and sP marked) and an overview plot.

## How the result will be used

* If **R** explains the failures (low |A| where the SNR is low), any later test that adds
  stations will choose them by predicted P amplitude for each family's GCMT mechanism, which
  can be fixed before any new waveform is looked at. It will be pre-registered first.
* If **F, L or X** rescue stations, the measurement will be changed for future tests only,
  never for H13 after its data are seen.
* If **D** holds at some stations, mantlep2's alignment there anchored on a depth phase. The
  mantlep2 result (decision N) compares events at each station with the same anchor, so it
  would still be a source comparison, but this will be stated next to it.
* If none of these explains the failures, we will say so and treat those stations as unusable
  at these magnitudes.

## Limits written in advance

* GCMT may lack solutions for the smallest events (mb 4.8-4.9). Those rows have no A.
* ComCat depths for these events are poorly resolved (two are fixed at 10 and 33 km), so pP-P
  predictions are approximate.
* Eight events from one small source region cannot separate a station-site effect from a
  path or radiation effect at every station.

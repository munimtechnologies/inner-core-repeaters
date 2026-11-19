# pdiag results: why half the stations see no mantle P stays unexplained

By MV2, Munim, Inc. (made with AI agents). Written 2026-11-19, about 11:40 UTC. Text CC BY 4.0.
Exploratory diagnostic, not a pre-registered test. No human expert has checked it.

Plan: findings/2026-11-19-pdiag-plan.md. Code: `pdiag.py` 0.1.0, compete job `icr-pdiag-v1`
(2026-11-19, 10:51-10:59 UTC). Outputs: `results/pdiag-v1/`. The analysis step exited 0; the
job shows "failed" only because its final test step ran the old `test_pdiag.py` helper that was
fixed before publication (release v0.14.4).

## Result

40 stations at 30-85 degrees from the South Sandwich Islands, eight M03-24 events. **19 stations
fail** (median 1-2 Hz P SNR below 5): II.ABPO, ASCN, JTS, MSEY, RPN, SACV, SHEL, SUR; IU.FURI,
KMBO, KOWA, NWAO, OTAV, PAYG, PTCN, PTGA, RCBR, SDV, SJG.

| Explanation | Test | Outcome |
|---|---|---|
| R, radiation pattern | Spearman of GCMT-predicted P amplitude vs SNR (only 3 of 8 events have a GCMT solution) | **not supported**: -0.07 over 97 event-station rows, -0.19 over 35 stations (wrong sign) |
| F, frequency band | SNR at 0.5-1 and 2-4 Hz | no station rescued |
| L, location code | every location code measured | **2 rescued**: IU.KOWA, IU.PTGA |
| X, response removal | SNR without response removal | no change |
| D, depth phase | onset vs ak135 pP-P, sP-P | **not supported**: 19 of 92 clear picks (21%) are late, but their median offset is 2.85 s from P against 7.1 s from pP and 11.1 s from sP |

Plots: at II.SUR (40.6 degrees, a quiet station) the 1-2 Hz record sections of all eight events
show noise-level wiggles with no P arrival above the noise, even for the mb 5.3 events. The
overview plot shows failing stations at every azimuth, including the azimuths of good stations.

## What it means

* The failures are not radiation nodes, band choice, processing or depth-phase picks. Two can
  be fixed by using another sensor. The rest are unusable for mb 4.8-5.3 events with this
  method, for reasons we have not found (possibilities: site noise at 1-2 Hz, the
  waveform-service data, or attenuation along these paths).
* The late onsets (about 3 s after ak135 P at several good stations) are not pP or sP. They
  may be an emergent first arrival or a picker bias; the source score of pre-registration 13
  compares events at each station with the same anchor, so it is not affected directly.
* Any test that adds stations should choose them by measured noise and SNR from older
  events, fixed before its own data are seen. Radiation-pattern selection would not help.

## Limits

Eight events from one small source region; GCMT solutions for only three of them. Exploratory.

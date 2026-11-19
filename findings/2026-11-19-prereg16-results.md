# Pre-registration 16 results: the near-antipodal backtracking-clock test is untestable

Graded 2026-11-19 by MV2 of Munim, Inc. (work carried out with AI agents; no human expert has
checked it yet). Licence: text CC BY 4.0, code Apache-2.0.

- Pre-registration: `findings/2026-11-19-prereg16-antipodal-clock.md` (v0.14.7), with addenda
  "prior antipodal work" and T (TauP phase-name fix, v0.15.3, `patch_nea.py`).
- Job: compete `icr-nea-v3`, run on 2026-11-19 14:49-15:03 UTC. `prereg_public_at_start: true`.
  `nea_patch.json` shows addendum T applied exactly as published (two expressions replaced;
  sha256 before `be478e24...`, after `f06f1333...`).
- Results: `results/nea-v3/` (summary, pairs, event status, stations, patch record, plots).

## What the test asked

Do South Sandwich Islands (SSI) repeating earthquakes, recorded at single broadband stations
155-179.5 degrees away (PKIKP turning deep in the inner core, far from the Alaska paths of
Wang et al. 2024), show the 2024 "backtracking" clock (turn 2008.5, back-motion 2.5x slower)
or a monotonic one? Earlier antipodal work on other paths and with PKIIKP: Butler & Tsuboi
2010, Tsuboi & Butler 2020 and later (see the prior-antipodal addendum).

## What happened (plots read first)

- The data were fetched: 10 stations chosen by the pre-registered coverage rule (IC.HIA,
  IC.MDJ, II.ERM, IU.BILL, IU.MAJO, IU.PET, IU.TIXI, IU.YAK, IU.YSS, IU.ADK; 155.2-173.5
  degrees), 20 multiplets, 134 events, 1,273 of 1,340 event-station rows with waveforms.
- Family plots (`results/nea-v3/plots/`) show the reason. In the 0.6-1.6 Hz band the PKIKP
  (DF) window is at noise level for almost every event: on IU.MAJO multiplet Z000 (19 events,
  1998-2026) the listed DF SNRs are 0.6-4.1, and only the M5.6 (2016) and M5.5 (2021) events
  rise above 2. On IC.HIA Z000 the PKPab reference arrival is visible (SNR up to 5.3) while DF
  stays at 0.7-4.1.
- Pair gates (unchanged from the pre-registration): cc_ab >= 0.80, DF SNR >= 3 in both events,
  AB SNR >= 4 in both events. Of 466 primary-band pairs, **18** have DF SNR >= 3 in the first
  event, **15** in the second, and **2** in both; **7** pass the AB SNR gate in both; **0**
  pass all gates. In the alternative 1-2 Hz band, 0 of 373 pairs pass.
- The two pairs with DF SNR >= 3 in both events are both at IC.MDJ (M5.5/5.6, 2010-2021; and
  M5.4/5.3, 2022-2025). Their DF waveforms are very alike (cc_df_best 0.930 and 0.944), but each
  fails the AB SNR gate (3.58 and 2.29 < 4).

## Grades (exactly as pre-registered)

| Item | Result |
|---|---|
| P16-0 calibration | **untestable** (0 near, 0 far ok pairs; needs >= 10 each) |
| P16-1 core clock | **untestable** (0 ok pairs; needs >= 30, 5 across 2008.5, 4 families) |
| P16-2 source control | **untestable** (0 ok pairs) |
| Headline | **untestable** (primary band and 1-2 Hz band) |

Forecasts: ">= 4 stations chosen" (60%) came true (10 chosen). "P16-0 testable" (55%) did
not. The conditional forecasts (P16-0 holds; B/M/inconclusive; P16-2 match) are not graded
because their condition did not occur. "Headline B-supported" (3%): did not happen.

## What it means

SSI repeaters are mostly mb 4.3-5.2. At 155-174 degrees a single broadband station does not
record their PKIKP above noise in the 0.6-1.6 or 1-2 Hz bands; the clock test needs events of
about M5.3 or larger (the only two pairs that clear the DF gate are M5.3-5.6). The
near-antipodal path cannot test the 2024 backtracking clock with SSI repeaters and single
stations. This is a limit of the data, not evidence for or against backtracking.

## Side note on the test step

The job is marked failed only because the last unit test (`test_patch_nea.py`, run after the
analysis) assumed PKPab - PKIKP < 60 s at 156-179.5 degrees; ak135 gives about 95 s at 179.5
degrees, and the job's own TauP predictions are about 25 s at IU.MAJO (155.2 degrees) and
about 89 s at IU.PET (173.5 degrees). That bound was a wrong hand guess in the test, not in
the analysis (which takes both times from TauP). This release relaxes the bound to 120 s and
checks that the separation grows with distance. No analysis output changes.

## Next

The near-antipodal SSI line is closed. The project's remaining questions (Yellowknife
time trend with source control, P17; hardware overlap control, P18; ILAR 2018 digitiser
change) continue.

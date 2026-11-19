# ILAR's 2018-09-25 digitiser-chain change: how big is its hardware term? (plan, before any result)

PLAN-ILRESP-2018

Written 2026-11-19 by AI agents working for MV2 of Munim, Inc., before the measurement ran.
Code: `ilresp.py` 0.1.0 and `test_ilresp.py` in this repository. Apache-2.0 (code), CC BY 4.0 (text).

## Why this matters

Pre-registrations 17 and 18 ask whether the YKA time trend is caused by YKA's 2013-14 hardware
change. ILAR, the other array of the 2024 backtracking study (Wang, Vidale et al., *Nature*
631:340), needs the same check. Its published metadata show a change inside the
post-2008.5 "backtracking" period. We read it from EarthScope fdsnws/station (`network=IM`,
`station=IL*`, `channel=SHZ`, `level=channel` and `level=response`) on 2026-11-19:

- Every short-period element IL01-IL19 SHZ has two channel epochs, split at
  **2018-09-25T00:00:00**. (IL15 has a third epoch beginning 2018-10-25.) The sensor string
  (`G23900=...=AIM24S_4.33E7`), the 20 sps output rate and the seismometer poles/zeros
  (two zeros at 0; poles -2.82743 +/- 4.89726i rad/s) are the same in both epochs.
- What changes: the sensor gain (one shared value of 1.60991E10 counts/(m/s) before; a separate
  value for each element after, 1.52447E10 to 1.79027E10) and the **digital decimation chain**.
  Stages 3-10 run five times faster after the change (204,800 Hz rather than 40,960 Hz at
  stage 3). The final 99-tap FIR at 40 Hz is replaced by a 101-tap FIR at 200 Hz and a further
  decimate-by-5 FIR at 100 Hz.
- Stage 2 of the new epoch still lists an input rate of 40,960 Hz while stage 3 expects
  204,800 Hz. That stage is a one-tap gain, so the mismatch has no effect on the response.

Wang et al. 2024 remove instrument responses. Our own pipeline (`yka.py` and the ILAR
scripts) band-passes raw counts at 1-2 Hz without removing the response, and normalised
cross-correlation ignores gain. So for our measurements, only the *shape* of the change in
transfer function inside 1-2 Hz matters. Pairs that straddle 2018-09-25 compare two digital
chains. Until this is measured, an ILAR change after 2018 could partly be hardware.

## What `ilresp.py` does (no waveforms are used)

For each element: fetch its StationXML. Take the epoch in force on 2015-01-01 (before) and on
2022-01-01 (after). Build each full response chain by hand: poles/zeros, coefficient and FIR
stages, each digital stage including the time correction the digitiser applied. Then compute
R(f) = H_after / H_before. ObsPy/evalresp gives a cross-check of the ratio's shape. Then:

1. In 1-2 Hz (main band) and 0.5-2 Hz (sensitivity band): the largest amplitude ripple of |R|
   (normalised at 1.5 Hz), the best-fit delay, and the rms phase left after removing that delay.
2. **Hardware term D_hw:** remove the best-fit delay (alignment removes a pure shift in the
   real analysis). Pass 50 tapered white-noise signals (seed 19) and one 1.5 Hz Ricker pulse
   through R, band-pass both versions (4-pole causal Butterworth, as in `yka.py`), and
   compute D = 1 - max cc over a 20 s window (lags +/-1 s, parabolic peak). D_raw keeps the delay in.

## Decision rule (fixed now)

The pair D values in this study run from about 0.07 to 0.2 (results/h13-v1/h13_pairs_scored.csv).
On the main band:
- **negligible**: the largest per-element p95 of noise D_hw is below 0.005;
- **material**: the median over elements of the per-element median noise D_hw is above 0.02;
- **small**: anything else. **Untestable**: fewer than 10 elements processed.

If **negligible**, the 2018 change cannot account for any ILAR decorrelation pattern, as far
as the published responses describe the hardware. If **material** or **small**, every ILAR
time-trend test must split pairs into same-chain and cross-2018 pairs (a new pre-registration),
as P17 does for YKA.

## Forecasts (before the run)

| # | Forecast | Probability |
|---|---|---|
| F1 | all 19 elements processed | 90% |
| F2 | verdict negligible / small / material | 80 / 15 / 5% |
| F3 | largest amplitude ripple in 1-2 Hz below 1% | 80% |
| F4 | largest residual phase rms in 1-2 Hz below 1 degree | 75% |
| F5 | the best-fit delay differs from zero by more than 0.01 s for some element | 35% |
| F6 | the hand-built and evalresp ratios agree (ripple within 0.1 percentage point) | 70% |

## Limits

This measures only what the metadata say changed. A real hardware swap that the StationXML
does not describe, or a response written down wrongly, would not show up. That is why P18's
empirical overlap test for YKA is the stronger kind of control. ILAR has no overlap epoch.
No waveforms are used, so nothing here is blind or unblind. Nobody has seen a result yet.
This work was carried out with AI agents, and no human expert has checked it yet.

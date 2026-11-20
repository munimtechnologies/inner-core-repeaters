# ILAR's 2018-09-25 digitiser-chain change: results (PLAN-ILRESP-2018)

Graded 2026-11-20 by AI agents working for MV2 of Munim, Inc., against the plan published before
the run: `findings/2026-11-19-ilar-2018-digitiser.md` (v0.15.4). Code: `ilresp.py` 0.1.0.
Apache-2.0 (code), CC BY 4.0 (text and tables).

Job `icr-ilresp-v1`, 2026-11-19 23:12:32-23:12:43 UTC. The analysis step exited 0. The job
shows as "failed" only because it ran the old copy of `test_ilresp.py` *after* the analysis, and that
copy asserted an exact zero lag (fixed in the repository's `test_ilresp.py` before v0.15.4). The
outputs are unaffected. `plan_public_at_start: true`, `errors: {}`.

Outputs, copied unchanged: `results/ilresp-v1/ilresp_summary.json`, `ilresp_elements.csv`,
`ilresp_stages.json`, `plots/ilresp_ratio.png`.

## Verdict: **negligible**

| Quantity (main band 1-2 Hz) | Value | Rule |
|---|---|---|
| elements processed | 19 of 19 | untestable if < 10 |
| largest per-element p95 of noise D_hw | 0.000109 | negligible if < 0.005 |
| median over elements of median noise D_hw | 0.0000808 | material if > 0.02 |
| Ricker-pulse D_hw (largest) | 0.0000637 | - |
| D_raw (delay kept in), median | 0.0000808 | - |
| largest amplitude ripple of R(f), normalised at 1.5 Hz | 2.26% | - |
| largest residual phase rms after the best-fit delay | 1.2e-13 deg | - |
| best-fit delay, all elements | -1.5e-8 s | - |

In the 0.5-2 Hz sensitivity band: ripple 2.45%, median noise D_hw 0.000146, p95 0.000184.

The pre-registered rule gives **negligible**: according to the published responses, the 2018
digitiser change cannot account for any ILAR decorrelation pattern in this study. The pair D values
here run from about 0.07 to 0.2, roughly 900 to 2,500 times the median ILAR hardware term. For
comparison, YKA's empirically measured 2013-14 hardware term (pre-registration 18) has a median
D_hw of 0.0117, about 145 times ILAR's.

All 19 elements give the same ratio shape (to about 1e-12). The only thing that differs between
elements is the post-2018 sensitivity (1.52447E10 to 1.83358E10 counts/(m/s); IL15's 2022 epoch is
the third one, from 2018-10-25), which normalised cross-correlation ignores. The plot
(`plots/ilresp_ratio.png`) shows |R| (normalised to 1 at 1.5 Hz) rising smoothly from about 0.99 at 1 Hz
to about 1.02 at 2 Hz and to about 3.7 at 9.5 Hz. Both chains are corrected to zero phase, so the ratio's phase stays below
6e-5 degrees up to 9.5 Hz. Studies using ILAR raw counts above about 3 Hz would see a larger
amplitude-shape change. Ours, at 1-2 Hz, do not.

## Forecasts

| # | Forecast | P | Outcome |
|---|---|---|---|
| F1 | all 19 elements processed | 90% | **yes** |
| F2 | verdict negligible / small / material | 80 / 15 / 5% | **negligible** |
| F3 | largest amplitude ripple in 1-2 Hz below 1% | 80% | **no**: 2.26% (a smooth tilt, not a ripple) |
| F4 | residual phase rms in 1-2 Hz below 1 degree | 75% | **yes** (1.2e-13 deg) |
| F5 | best-fit delay differs from zero by > 0.01 s for some element | 35% | **no** (-1.5e-8 s) |
| F6 | hand-built and evalresp ratios agree (ripple within 0.1 pp) | 70% | **yes** (2.2634% both) |

F3 failed. We expected the two anti-alias FIR chains to be flat to within 1% across 1-2 Hz. Instead,
the new chain's passband rises 2.3% relative to the old one over that band. A smooth tilt of this
size barely changes waveform shape (D_hw ~ 1e-4), so the verdict is unaffected.

## What this means for the study

- The 2018 change does **not** require splitting ILAR pairs into same-chain and cross-2018 groups,
  so the conditional pre-registration (P19 split test) planned for a "small" or "material" verdict
  is not needed.
- ILAR results so far (P1 in-sample lean toward backtracking, 88% bootstrap vs 90% required;
  source-verified secondary in H13, 67 pairs, rho_B 0.136 vs rho_M 0.141, inconclusive) therefore
  carry no documented hardware term at 1-2 Hz. YKA's do (pre-registrations 17 and 18).
- Wang, Vidale et al. 2024 (*Nature* 631:340) mention YKA's 2013 instrument change but none at
  ILAR. ILAR's metadata do show a digitiser change on 2018-09-25, but by the published responses
  its effect on 1-2 Hz waveform shape is negligible.

## Limits

This measures only what the StationXML says changed. An undocumented hardware change or a response
recorded wrongly would not show up. ILAR has no overlap epoch, so the empirical check made for YKA
(pre-registration 18) is not possible here. No waveforms were used. This work was carried out with
AI agents, and no human expert has checked it yet.

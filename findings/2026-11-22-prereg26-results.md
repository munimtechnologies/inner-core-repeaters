# Pre-registration 26 results: the entry-side neighbour test is untestable

*inner-core-repeaters, MV2 of Munim, Inc. Graded 2026-11-22 ~13:30 UTC, exactly as written in
`findings/2026-11-21-prereg26-sks-entry-neighbours-plan.md` (v0.24.0). Job `sks-nbr26-v1` ran
2026-11-22 11:29-11:39 UTC: `sksnbr26.py` 0.1.0, code sha256 5d02ee73…, plan sha256 300acb55…,
plan public at start. Results are in `results/sksnbr26-v1/`. Made with AI agents working in MV2's
name; no human expert has checked it yet.*

## Grade

**Primary (source side): untestable.** The plan's rule is "untestable if no entry-neighbour is
usable", and none was.

**Secondary (receiver side): untestable.** No station's SKS leaves the core within 100 km of
SDV's exit point.

## What happened

**Gates.**
- **G0** passed: our pipeline gave x = +0.0284 s at IU.SDV, the same as P23's x_meta.
- **G1** passed: E2009a is a verified, relocated member (13 stations at cc >= 0.90; dt0 = -2.088 s;
  relocation rms 0.041 s).

**Candidates.** EarthScope listed 63 candidate stations at 85-125 deg.

**Entry-neighbours.** Only one station sends its SKS into the core within 100 km of SDV's entry
point: **VE.FUNV** (Caracas), 75 km away. EarthScope has no FUNV data for the 2000 event
(`nodata E2000a`).

**Exit-neighbours.** None.

**Controls.** 39 permanent global stations at least 300 km from SDV's entry and exit points were
attempted.
- **No data for one of the two years:** 20 stations.
- **Measured but unusable:** the other 19. Every one failed the SKS rule (radial SKS cc >= 0.80 and
  SNR >= 5 in both years). None reached an SNR of 4 in both years.

With zero usable controls, the control false-change rate could not be computed either.

**Report-only (not graded).** Two stations lie within 300 km of SDV's entry point:

| station | SKS entry from SDV's | exit from SDV's | status | x (s) | xz (s) | reference |
|---|---|---|---|---|---|---|
| IU.SDV | 0 km | 0 km | usable | +0.028 | +0.024 | S_T |
| G.FDF (Martinique) | 199 km | 653 km | usable | -0.004 | +0.019 | P_Z |
| IU.SJG (Puerto Rico) | 222 km | 681 km | unusable (SKS cc 0.48, SNR 1.6) | - | - | - |
| IU.PTGA (Brazil) | 284 km | 1147 km | no 2000 data | - | - | - |

G.FDF's 2000 and 2009 records come from different location codes (blank, then 00), so the sensor
may have changed. Its transverse S correlated poorly (cc 0.50), so P was its reference.

## What it means

- **Untestable, with no regrading.** The test could not be carried out with public data. We will
  not re-run it with a wider neighbour radius or looser rules.
- **A feasibility finding.** This doublet (Mw 5.9 and 5.8) is too small for most global stations
  at 85-125 deg in the 0.4-2 Hz band. SDV, with SKS SNR of 20-24, is the exception, not the rule.
  That limits any multi-path check of the claim. It is probably also why the paper reports one
  station.
- **A hint, not evidence.** The only other usable station near the entry region, Martinique at
  199 km, shows no change (-4 ms). That lies outside the pre-registered 100 km radius and is not
  graded. Read it as a hint, not as evidence.

## Status of the core-mantle-boundary line

Zhang & Wen 2025, *Nature Communications* 16:8754:

| test | grade |
|---|---|
| P23 (public SDV records) | not-reproduced (+28 ms vs their -160 ms) |
| P25 (onset/pulse width) | untestable |
| P26 (entry neighbours) | untestable |
| P24-v2 (the authors' own processed traces) | queued |
| P27 (exit-side repeating pairs) | queued |

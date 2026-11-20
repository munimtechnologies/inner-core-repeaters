# ILAR power check for a larger source-controlled clock test (P19b): dropped

Date: 2026-11-20. Repository: munimtechnologies/inner-core-repeaters (Apache-2.0; text CC BY 4.0).
Work carried out with AI agents working in MV2's name (MV2, Munim, Inc.); no human expert has checked it.

## Question

The ILAR digitiser check (findings/2026-11-20-ilar-2018-digitiser-results.md) found that the
2018-09-25 hardware change is negligible at 1-2 Hz, so ILAR pairs need no hardware split. H13's
ILAR secondary (results/h13-v1/h13_summary.json) graded only 67 source-verified ILAR pairs
(28 families, 20 crossing 2008.5): rho_B 0.136, rho_M 0.141, frac_B_gt_M 0.419, inconclusive.
Would a larger, source-controlled ILAR clock test (P19b) have enough pairs to be worth
pre-registering? Before any work we set this rule (playbook, run 51): drop P19b if it would add
fewer than about 60 new verified pairs or fewer than 10 new verified pairs crossing 2008.5.

No new data were downloaded. All counts below come from files already in the repository:
results/ilar-v1/pairs.csv (874 rows, commit 7067d4fc), results/ilar-v1/run_summary.json,
results/ykab-v2/yka_pairs.csv and merged_multiplets.csv, results/h13-v1/*. Counts were made by
regular-expression filters on the CSV columns (snr_df1, snr_df2, mag1, mag2, cc_ref,
cross_turn, ok).

## What limits ILAR: inner-core (DF) signal-to-noise by magnitude

ILAR's pipeline requires DF SNR >= 3 at both events (min_snr_df 3.0). Counted per pair side
(an event appears once for every pair it belongs to):

| Magnitude of the event | Pair sides | DF SNR >= 3 | Share |
|---|---|---|---|
| 4.7-4.9 (only events after 2024-07-01 were admitted below 5.0) | 85 (78 as ev2, 7 as ev1) | 0 | 0% |
| 5.0 (ev2) | 208 | 71 | 34% |
| 5.0-5.2 (ev2) | 650 | 314 | 48% |
| 5.3 and above (ev2) | 146 | 124 | 85% |

41 of the 78 M4.7-4.9 ev2 sides reached DF SNR 2.0-2.99, so they are close but none pass. The
practical floor for an ILAR inner-core pair is M5.0 for both events.

## The pool of pairs

- ILAR measured 874 pairs inside its own 75 multiplets (280 events). 321 have DF SNR >= 3 at both
  events. 136 of these are "ok" (they also pass cc_ref >= 0.93 and reference SNR >= 5); H13
  scored 118 of them and verified 67.
- The other **185 DF-passing pairs** (48 of them crossing 2008.5) were never used, because they
  fail ILAR's reference-phase gate. 130 of the 321 have cc_ref below 0.90: they are linked only
  through chains inside a multiplet.
- The merged family set behind H13 (78 families, about 407 events, M >= 4.3) adds mostly
  M4.3-4.9 events from YKA detection, which cannot pass ILAR's DF gate. YKA's table has 820 pairs
  with both events M >= 5.0, about the same as ILAR's own pairs without an M < 5.0 event (about 789 = 874 - 85). So
  the merged families add few new M >= 5.0 pairs. We estimate 50-150, from the families merged
  across ILAR multiplets.

## Estimate

With H13's mantle-P S score as the source control (S <= 1/3), the fraction of pairs verified was
25% at YKA (328 of 1,324 scored) and 57% for ILAR's cc_ref-gated pairs (67 of 118). The 185
unused pairs are less alike at ILAR by construction (most fail cc_ref), so the YKA rate is the
better guide. Scoring coverage was 91% for ILAR pairs in H13 (118 of 130).

- From the 185 unused pairs: 185 x 0.91 x 0.25 = about 42 new verified pairs, about 11 crossing 2008.5.
- From merged-family extras: (50-150) x 0.41 (ILAR's DF pass rate for pairs without an M < 5.0
  event, 321/789) x 0.25 = about 5-15.
- **Total: about 47-57 new verified pairs (below 60), with about 11-15 crossing 2008.5 (at least 10).**
- Sensitivity: if the unused pairs verified at the 57% rate of ILAR's best pairs, the count would be
  about 96 plus 12-35 extras, above the threshold. We do not expect this. S is a relative score
  (the most alike third at each mantle station), and the unused pairs failed ILAR's own
  reference-phase similarity gate, so their rate should be at or below the pooled 25%.

Expected outcome if run anyway: the 67 pairs already graded show no clock preference (rho_B -
rho_M = -0.005). To reach the pre-registered FRAC 0.9 with about 120 pairs, the new pairs alone
would need a backtracking advantage at least as large as ILAR's in-sample result without
source control (rho_B 0.194 vs rho_M 0.065, frac 0.88 on 123 pairs). We judge "inconclusive" the
most likely result (about two chances in three).

## Decision

By the rule set in advance, **P19b is underpowered and is dropped**. No pre-registration is
filed. A source-controlled ILAR clock test using South Sandwich Islands repeaters is limited by
the data, not the method: only M >= 5.0 pairs give a measurable PKIKP at ILAR, and there are
too few of them for a decisive verdict between the backtracking and monotone clocks with this
coda metric.

The study's open lines are now the large-event near-antipodal pilot (snrpilot, queued), plus
the standing single-family M03-24 result.

## Limits

Counts are per pair side and pair-weighted, not per distinct event. Verification rates for the
unused pairs are estimated from H13, not measured; a job could compute them exactly, but by
the rule above the central estimate already decides the question. Hand counts were made with
column filters on the committed CSVs and can be re-run by anyone.

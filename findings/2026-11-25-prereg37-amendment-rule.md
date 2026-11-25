# Pre-registration 37: the amendment is decided by a fixed rule (written before any probe output)

Date: 2026-11-25 (about 09:40 UTC). Repo: munimtechnologies/inner-core-repeaters. Licence of this
text: CC BY 4.0. Made with AI agents working in MV2's name (MV2, Munim, Inc., mv2@munimtech.com); no
human expert has checked it.

Plan: findings/2026-11-25-prereg37-dprk-doublets-plan.md (commit d0715ef). The plan allows ONE
amendment after the blind probe dprk37probe-v2 returns, limited to the primary band, the minimum
station counts of section 6, and the choice of G2. **This note removes the remaining freedom in
that amendment.** It is written while the probe job is still queued: no probe output, record,
count or waveform from any of the six tests has been read by anyone working on this study.

## Input (the only probe numbers the rule may use)

N(pair, band) = the column `n_df_ab_same_epoch_clear` of `dprk37_pair_counts.csv` for that pair and
band. This is the probe's count of stations passing plan gates 1 to 4 (DF and AB SNR >= 5 on both
records, same channel and epoch, windows clear). It cannot see gate 5 (cc_AB >= 0.80), so the
job's n can only be equal or smaller. Nothing else in the probe output is used (not the binned
per-record SNRs, not the event summaries).

## Rule 1: primary band

Keep 0.7-2.0 Hz, unless N(G1, 0.7-2.0) < 8 and N(G1, 1.0-2.0) >= 8. In that case only, the primary
band becomes 1.0-2.0 Hz (and 0.7-2.0 Hz becomes the report-only band). The band applies to every
pair. No other reason (larger n, "cleaner" data) may change the band.

## Rule 2: minimum station counts

Not changed. Untestable stays n < 8 for G1 and G2; P37b stays at least 4 patch and 4 non-patch
stations. Lowering a minimum after seeing counts would weaken the test, so the amendment never does
it, and it never raises one either.

## Rule 3: G2

G2 stays 2013-02-12 / 2016-01-06 if N(G2, primary band) >= 8. Otherwise G2 becomes the first pair
in this order with N(pair, primary band) >= 8:

1. 2013-02-12 / 2016-09-09 (same size class, mb 5.1 / 5.3);
2. the pairs with 2017-09-03 (2013-02, 2016-01, 2016-09), by largest N; ties go to the shorter
   interval, then to the earlier first test.

If no pair qualifies, G2 stays 2013-02-12 / 2016-01-06 and is reported untestable. Pairs with the
small 2006 and 2009 tests are never G2.

## What follows mechanically

- The amendment note (findings/<date>-prereg37-amendment-1.md) states N for G1 and G2 in both
  bands and the three outcomes of rules 1-3, or "no amendment" when all three keep the plan.
- The job dprk37-v1 is filed only if, under the rule, G1 or G2 has N >= 8, or if at least one
  report-only pair has N >= 8 (it is then filed with G1/G2 reported untestable, exactly as
  computed). Otherwise the line is closed with an untestable note and no job is filed.
- Everything else in the plan (gates, windows, labels, forecasts, P37b) is unchanged.

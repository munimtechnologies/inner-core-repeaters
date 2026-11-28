# Pre-registration 37, amendment 1: the rule closes the DPRK line as untestable

Date: 2026-11-28 (about 14:30 UTC). Repo: munimtechnologies/inner-core-repeaters. Text CC BY 4.0.
Made with AI agents working in MV2's name (MV2, Munim, Inc., mv2@munimtech.com); no human expert
has checked it.

Plan: [findings/2026-11-25-prereg37-dprk-doublets-plan.md](2026-11-25-prereg37-dprk-doublets-plan.md).
Fixed amendment rule, written before any probe output:
[findings/2026-11-25-prereg37-amendment-rule.md](2026-11-25-prereg37-amendment-rule.md).
Probe: `dprk37probe.py` 0.1.1, job dprk37probe-v2 (EarthScope FDSN, ran 2026-11-28 12:37-13:00 UTC,
both steps exit 0). Its three output files are committed unchanged in
[results/dprk37probe-v2/](../results/dprk37probe-v2/).

## Input the rule may use

N(pair, band) = `n_df_ab_same_epoch_clear` in `dprk37_pair_counts.csv`.

| Pair | Role | N, 0.7-2.0 Hz | N, 1.0-2.0 Hz |
|---|---|---|---|
| 2016-01-06 / 2016-09-09 | G1 | 0 | 0 |
| 2013-02-12 / 2016-01-06 | G2 | 0 | 0 |
| 2013-02-12 / 2016-09-09 | G2 candidate 1 | 0 | 0 |
| any pair with 2017-09-03 | G2 candidates 2 | 0 | 0 |
| any pair with 2006 or 2009 | report only | 0 | 0 |

All 30 rows of the file (15 pairs x 2 bands) are 0 in every count column, including the weakest
one (`n_df_both`, DF SNR >= 5 on both records, no other gate).

## The three rules, applied mechanically

1. **Primary band:** N(G1, 0.7-2.0) = 0 < 8, but N(G1, 1.0-2.0) = 0 is not >= 8, so the band stays
   0.7-2.0 Hz.
2. **Minimum station counts:** unchanged (n >= 8 for G1 and G2; P37b 4 + 4).
3. **G2:** N(G2) = 0 < 8, and no candidate pair reaches 8, so G2 stays 2013-02-12 / 2016-01-06 and is
   untestable.

**What follows:** neither G1, G2 nor any report-only pair has N >= 8, so under the rule **dprk37-v1
is not filed and pre-registration 37 closes as untestable** (G1, G2 and P37b). No two tests were
compared at any point; the probe never computes a cross-correlation.

## Why (single-record numbers, report only, not used by the rule)

From `dprk37_summary.json` (`events_summary`), stations 150-179.5 deg from Punggye-ri:

| Test | mb | records ok | DF SNR >= 5 (0.7-2.0 Hz) | DF and AB SNR >= 5 (0.7-2.0) | DF SNR >= 5 (1.0-2.0) |
|---|---|---|---|---|---|
| 2006-10-09 | 4.3 | 22 | 0 | 0 | 0 |
| 2009-05-25 | 4.7 | 93 | 0 | 0 | 1 |
| 2013-02-12 | 5.1 | 34 | 0 | 0 | 0 |
| 2016-01-06 | 5.1 | 78 | 0 | 0 | 0 |
| 2016-09-09 | 5.3 | 114 | 0 | 0 | 1 |
| 2017-09-03 | 6.3 | 106 | 56 | 54 | 63 |

The pipeline works when the source is large enough: the 2017 test gives PKIKP (DF) and PKPab above
SNR 5 at 54 stations (DF alone at 56 stations in nine networks). The five smaller tests give none at 0.7-2.0 Hz. At these
near-antipodal distances, single broadband stations do not see inner-core P from the 2006-2016
explosions, so the doublet design cannot be done station by station. The dense 2015-2016 temporary
deployment in Argentina (network 5C) recorded both 2016 tests at 13 and 19 stations, but its
single-record DF SNRs are almost all below 3. Stacking them would be a new design written after
seeing these SNR bins, so it is not part of P37. It is noted here only as an idea.

## Forecast against outcome

The plan forecast "G1 testable: 55%" and "G2 testable: 50%". Both came out untestable. The data
volume was not the limit (78 and 114 records for the two 2016 tests). Source size was. The plan's
choice of the closest-size pair (mb 5.1 / 5.3) traded signal for source similarity, and at these
distances that left no signal at all.

## Status

Line X (P37) is closed: untestable. No job is filed, and no further amendment is possible.

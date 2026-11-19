# Pre-registration 12: the grader is published before any prospective data

MV2, Munim, Inc. (made with AI agents; no human expert has checked this yet). Text CC BY 4.0, code Apache-2.0.
Written 2026-11-19 (UTC).

Pre-registration 12 ([findings/2026-11-18-yka-prospective-prereg.md](2026-11-18-yka-prospective-prereg.md))
is a prospective test of inner-core "backtracking" at the Yellowknife array (YKA). It uses only South
Sandwich Islands earthquakes from **2026-11-19T00:00Z** onward, with looks on or after 2027-11-19 (interim)
and 2028-11-19 (final). It promised that the grader would be published **before the interim look and before
any post-2026-11-19 YKA waveform is measured**. This note and `ykaprosp_grade.py` 0.1.0 keep that promise.
No waveform of a new event has been downloaded or measured by this project at the time of writing.

## What the grader does

`ykaprosp_grade.py --pairs yka_pairs.csv --out-dir OUT --look interim|final [--catalog comcat.csv]`

- Input: the pairs table that `yka_b.py` 0.10.1 `--repick` writes for the extended catalogue (families
  merged by ILAR correlation with the unchanged `run.py` threshold), plus ComCat locations.
- It refuses to grade before the look opens: it writes `status: too-early` and exits with code 2.
- Selection, P4-prosp (calibration), F1-prosp (B vs M vote, alpha 0.025, 2,000 whole-family bootstraps)
  and the two ungraded extras (D +1 to +6 s; 2025-2026 partners removed) follow the pre-registration
  word for word. Outputs: `ykaprosp_grade.json`, `ykaprosp_pairs.csv`, `ykaprosp_cases.csv`.
- `test_ykaprosp_grade.py` checks it on synthetic families: supported, contradicted, inconclusive,
  P4 failure blocking the vote, untestable at interim vs final, every exclusion reason, the date gate,
  and that its decimal years match `yka_b.py` to 1e-7 years on two published rows.

## Choices the pre-registration left open, fixed now (before any data)

1. **"New" event**: ComCat origin time on or after 2026-11-19T00:00:00Z, judged from the ISO time in the
   pairs table (not from decimal years).
2. **Pairs of two new events**: counted once, with the later event as the new event and the earlier as the
   partner. Such a pair is "near" in P4 and can never be a pre-turn partner in the vote.
3. **Selection applies to every pair**, including near, far and vote pairs: `ok` (yka_b gates),
   finite metric, both events located in ComCat (by id, else origin time within 5 s), separation
   <= 14.7 km, and min(snr_dir1, snr_dir2) >= 11.7. A pair with an unlocated event is not selected.
4. **"More than 5 years apart on the backtracking clock"**: |e(t_new) - e(t_partner)| > 5.0 years,
   with e(t) = t before 2008.5 and 2008.5 - (t - 2008.5)/2.5 after; decimal years from the table's t1/t2.
5. **Near/far dates**: near = partner on or after 2024-11-19T00:00Z; far = partner before 2016-11-19T00:00Z.
6. **Ties**: the best match is the lowest D_c13, ties going to the earlier partner; the B pick is the
   partner closest to e(t_new), ties going to the earlier partner (as in yka_b.py).
7. **"2025-2026 events removed from the partner pool"** (extra only): partners with origin time in
   [2025-01-01, 2026-11-19) are dropped before selection.
8. **Bootstrap**: resample families with replacement (seed 0); "n_B > n_M in >= 95%" means the fraction of
   resamples with n_B strictly greater than n_M is at least 0.95.
9. **Decision order**: P4 untestable or failed: "not interpreted". Else, at the interim look, an untestable
   vote is "not graded at interim"; at the final look it is "untestable". Else the vote's verdict.

These choices were made without any prospective data and will not change. If the grader ever needs a bug
fix, the fix will be published with a dated note before the look it affects, and the original 0.1.0
result will be reported alongside.

# Pre-registration 31, run v3 (HTTPS master tars): K untestable; P3 untestable

*inner-core-repeaters, MV2 of Munim, Inc. Written 2026-11-24 ~23:00 UTC. Text CC BY 4.0; code Apache-2.0.
Made with AI agents working in MV2's name; no human expert has checked it yet.*

## Headline (graded exactly as the plan says)

**K: untestable; P3: untestable.**

- Plan: [2026-11-23-prereg31-brv-waveform-plan.md](2026-11-23-prereg31-brv-waveform-plan.md).
  Data route: [2026-11-24-prereg31-data-route-https.md](2026-11-24-prereg31-data-route-https.md).
- Job: brv-wave31-v3, run 2026-11-24 21:58-22:19 UTC. Unit tests passed first; then `stage31h.py` 0.1.0
  and the unchanged `brvpkp31.py` 0.1.1 / `brvpkp31_core.py` from tag v0.32.0.
- Results: [results/brv-wave31-v3/](../results/brv-wave31-v3/) (`brv31_summary.json` is the graded file).

The question was whether the 0.1-0.2 s PKP-PP changes that Zhou (2026, JGR Solid Earth,
[doi:10.1029/2026JB033798](https://doi.org/10.1029/2026JB033798)) reads in Borovoye records of the French
Tuamotu tests as outer-core change are also seen in PPP-PP, a wave pair that never enters the core.
This time the records were read. The test still could not be graded, for reasons in the records themselves.

## The data route worked

- Both LDEO master tars downloaded over HTTPS (FRENCH_1977-1988.tar, 72 members, 22.6 MB;
  FRENCH_1989-1996.tar, 27 members, 20.5 MB; SHA-256 in `stage31h_summary.json`).
- All 99 units are French (site digit 4). **All 23 Table S1 tests have a primary (TS) record**, every one
  on the same date (day offset 0). 17 of 23 have an SS replica record. The six without one are 1978-11-30,
  1980-12-03, 1982-07-25, 1988-11-23, 1995-10-27 and 1995-12-27.
- For 1988-11-23 the core's index lookup points at a Semipalatinsk file (`19881123-61-TS`). As fixed in the
  data-route note, we staged the French record of that date (`19881123-42-TS`) under that name instead.

The earlier FTP failure (run v2) is therefore not a property of the archive. The LDEO access question
drafted for Sheehan on 2026-11-24 is no longer needed.

## Why the gates failed (from the job's own files)

**G0 (K testable) needs at least 10 K-usable pairs, the 1977 test in one of them, and all three epochs.**
In the primary set (TS, 18 starred tests, every pair within 0.05 deg) the job found:

| | primary TS, 18 starred | TS, all 23 (report-only) | SS, 18 starred (report-only) |
|---|---|---|---|
| pairs | 81 | 117 | 81 |
| measured | 44 | 74 | 30 |
| K-usable | **8** (gate 10) | 9 | 0 |
| P3-usable | **0** (gate 10) | 0 | 0 |
| 1995 epoch present | no | no | no |

Three things in the records cause this.

1. **Recording windows.** The TS records are cut windows that start 960-1140 s after the origin time.
   The predicted PKIKP arrival is at 1172.6 s, PP at 1356.9 s and PPP at 1543.3 s (TauP values in the
   summary).
   - 1977-11-24, 1978-11-30 and 1980-12-03 start at about 1080 s and last 363-379 s, so they end at about
     1443-1459 s, before PPP. These tests carry the 1977 reference, so P3 can never be measured relative
     to it.
   - 1983-05-25 and 1984-05-12 start at about 1140 s. That is less than the 50 s of pre-PKP noise the
     plan requires, so K is not covered.
   - 1984-11-02 ends at about 1386 s and 1995-10-27 at about 1396 s. Both are too short for the 30 s
     after PP (1984-11-02) or for PPP (both).
2. **Signal-to-noise at 0.3-1 Hz.** The plan's gates are SNR >= 3 at both tests and cross-correlation
   >= 0.8 for K and PP (0.7 for PPP). Most pairs that do have coverage fail on SNR or correlation. On the
   short-period vertical sensors (sZ01-sZ10), PPP SNR stays below 3 for every test. The few channel
   values of 3 or more are on unlabelled or long-period channels whose other ratios are clearly glitches
   (for example 0.001 or 10^6). No pair reached P3-usable. The plan forecast this risk ("PPP is weak at
   0.3-1 Hz for mb 5-6"). It gave G1 a 40% chance of passing, given that G0 passed.
3. **1995.** The two 1995 tests are the only ones in the third epoch. 1995-10-27 has PKP and PP, but its
   PP SNR is below 3 on every eligible channel (at most 2.66). 1995-12-27 has a long record but SNR below 3 on every
   eligible channel. So E1995 is empty in every variant.

## Report-only numbers (not a grade)

The primary K year terms, relative to 1977, are in `brv31_summary.json` under
`variants.primary_TS_starred.K_year_terms`. The two epoch values the job could form are:

| epoch | Zhou 2026 (paper midpoint) | this job (K-usable pairs only) |
|---|---|---|
| 1982-83 | -0.125 s | -0.160 s |
| 1988-90 | -0.150 s | -0.057 s |
| 1995 | +0.175 s | not measurable |

These two values must not be read as a reproduction or a refutation. Each rests in practice on a single
test linked to the 1977-1980 cluster by one or two pairs: 1982-07-25 through the 1977 and 1985 pairs,
and 1988-05-25 through the 1980 pair. G0 failed, so the plan grades nothing from them. The all-23 and SS
variants add nothing gradable.

## Forecast check

| forecast (made before data) | outcome |
|---|---|
| archive reachable and records readable, 60% | yes (run v3; v2 had failed on FTP) |
| G0 passes given readable records, 60% | **no** (8 of 10 pairs; no 1995 epoch) |
| K grade, G1 and P3 grade | not reached |

## What this means

The public LDEO archive of Borovoye records, as distributed, cannot carry a PPP-PP test of Zhou's
outer-core change. The pre-1982 records end before PPP, and PPP is below the noise at 0.3-1 Hz on single
short-period channels. The PKP-PP leg alone cannot pass a minimal reproduction gate (10 pairs and all
three epochs) with the paper's 0.05 deg pairing and our fixed quality gates. Zhou's own measurement used
procedures we did not copy, such as manual picking and channel choice. So this is not evidence against
the paper's numbers. It means an independent, rule-based check from the public records is not possible
at this frequency band.

**P31 is closed.** A stacked-beam or wider-band re-analysis of these same records would now be
exploratory, because the records have been seen. Its 1995 epoch would also rest on one weak test, so we
do not plan one. Our graded test of the same claim from the source side (P30: test-site position does not
explain it) stands.

## Files

`results/brv-wave31-v3/`: `brv31_summary.json` (graded), `stage31h_summary.json`, `stage31h_staged.csv`,
`stage31h_members.csv`, `brv31_files.csv`, `brv31_inventory.csv` (start offset and length of every
record), `brv31_channels.csv` (per-channel SNR, clipping, coverage) and `brv31_pairs.csv` (every pair's
lags, correlations, SNRs and why it was or was not usable). LDEO's index page copy (`brv31_index.html`)
is not re-published.

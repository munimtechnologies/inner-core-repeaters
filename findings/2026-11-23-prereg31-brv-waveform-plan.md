# Pre-registration 31: is the reported outer-core change in the Borovoye waveforms, and does PPP-PP move with it?

Written 2026-11-23 ~00:30 UTC and revised ~04:30 UTC (see "Change before publication" below), **before any
Borovoye waveform of these tests was downloaded or looked at by us**. Code: `brvpkp31.py` 0.1.1 (a wrapper
that changes only the glitch filter) around `brvpkp31_core.py` (the 0.1.0 measurement code, unchanged), both
Apache-2.0; tests `test_brvpkp31.py` and `test_brvpkp31_despike.py`. Job `brv-wave31-v2`. The job waits until
this file is public on `main`, records its sha256, and measures nothing if it never appears. Made with AI
agents working in MV2's name; no human expert has checked it yet. Text CC BY 4.0.

## Change before publication (2026-11-23, before any data)

The first attempt to publish this plan (2026-11-22 23:45 UTC) was stopped by the repository's own unit test:
`test_despike` found 0 of 3 synthetic glitches. Cause: on a smooth or integer-valued trace a 7-sample running
median equals the centre sample almost everywhere, so the MAD of the residual is 0 and the 0.1.0 filter
silently did nothing. The Borovoye records are integer counts, so the same would have happened on real data.
Replacing a two-sample glitch by the running median was also biased on steep slopes. Step 2 below is the
corrected rule. Nothing else changed: windows, bands, gates, grades, forecasts and inputs are as first
written. No waveform had been fetched by either version (job brv-wave31-v1 never started; it would have
stopped at the same test).

## Claim under test

Y. Zhou (2026), *JGR Solid Earth* 131, e2026JB033798 (doi 10.1029/2026JB033798, CC BY): French underground
tests at Mururoa (Table S1, 23 tests 1977-1995; 18 starred tests form 112 pairs within 0.05 deg), recorded
at Borovoye (BRV, 141.6 deg) on the "STsR-TSG KS" system. PKrKP-PP and PKiKP-PP differential times, after
subtracting their dependence on the catalogue distance difference, are read as 0.1-0.15 s faster in 1982-83
and 1988-90 and 0.15-0.2 s slower in 1995, relative to 1977, and interpreted as structure moving in the
liquid outer core. No measurements are released; the records are public at Lamont-Doherty (Kim & Ekstrom
archive).

Our earlier tests: P30 (graded, v0.31.0) found that the disclosed shaft positions do not reproduce the
history under either reading, but that catalogue location errors after the paper's correction (reading C)
come close for 1982-90 and fail only in 1995.

## Two questions

1. **Reproduction (K).** Measured directly on the public waveforms with the same pairing and the same
   catalogue-distance correction, does the core arrival move relative to PP as reported?
2. **Discrimination (P3).** PPP, the twice-reflected P wave, never enters the core. Measured on the same
   records against PP:
   - an **outer-core change** leaves PPP-PP unchanged (ratio 0);
   - a **test-position error** moves PPP-PP by r_pos = (p_PPP - p_PP)/(p_K - p_PP) times the PKP-PP change;
   - a **time-base error** (a recording sample interval wrong by a fraction e) moves every interval in
     proportion to its length: r_tb = (T_PPP - T_PP)/(T_K - T_PP).
   p and T come from TauP iasp91 at the reference site (-21.835, -138.91), depth 0.6 km, computed in the
   job. By hand, r_pos is about -0.5 and r_tb about -1.5; the job's numbers are the ones used.

## Data (fixed in `brvpkp31_core.PARAMS`)

- Table S1 from this repo at tag v0.31.0 (`data/zhou2026_table_s1.csv`, transcribed in run 74).
- Records: the LDEO Borovoye explosion archive. Index page
  `.../res/pi/Monitoring/Data/Brv_arch_ex/brv_74-96_ex.html`; one compressed CSS 3.0 database per test and
  recording system, `YYYYMMDD-41-TS.tar.Z` (STsR-TSG, **primary**, the system the paper names) and
  `YYYYMMDD-41-SS.tar.Z` (STsR-SS, **report-only replica**). Links come from the index page; if it cannot
  be read, the pattern `ftp://ftp.ldeo.columbia.edu/archive/wykim/Brv_archive/Exp/YYYYMMDD-41-SYS.tar.Z`.
- Each wfdisc row is read as written (CSS 3.0 datatypes s4, s3, s2, i4, i2, t4, t8, f4, f8). No instrument
  response is removed (every phase is measured on the same channel in the same band).

## Processing

1. **Times.** Phase times per test from TauP iasp91 at the Table S1 epicentre, depth 0.6 km: first PKIKP
   (T_K), PP, PPP; absolute time = Table S1 origin + T. Windows are therefore placed with catalogue
   distances, so every measured pair difference is already corrected for the catalogue distance difference,
   as in the paper.
2. **Glitches (corrected rule, 0.1.1).** Residual from a 7-sample running median; robust sigma =
   max(1.4826 MAD of the residual, 1.4826 MAD of first differences / sqrt 2); samples beyond 20 sigma in
   runs of at most 3 samples are replaced by linear interpolation between the nearest unflagged samples.
   Counts reported.
3. **Clipping.** A window is clipped if it holds a run of >= 3 identical samples with |x| >= 0.98 of the
   trace's maximum |x| (after glitch removal).
4. **Filter.** Applied to the span PKIKP-70 s to PPP+50 s. Linear detrend, cosine taper (5%, at most
   10 s), 4-pole Butterworth band-pass 0.3-1.0 Hz, zero-phase (the
   paper's band).
5. **Windows** (s, relative to the predicted time): K = PKIKP [-3, +5]; PP [-3, +7]; PPP [-3, +7]; noise =
   PKIKP [-40, -10]. SNR = rms(window)/rms(noise).
6. **Channels.** Vertical channels (name contains "z", any case); if a record has none, all its channels.
   A channel is eligible for a test if it is sampled at >= 4 Hz, covers PKIKP-50 s to PP+30 s, and K
   and PP are unclipped (checked over each window +/- 1 s). For a
   pair, among channel names eligible in both records, the one with the largest min(SNR_K in A, SNR_K in B)
   is used (chosen before any correlation).
7. **Coarse alignment.** Envelopes (Hilbert, 1 s running mean) of window PKIKP [-10, +30] s are
   cross-correlated within +/-8 s at 10 Hz: shift c (absorbs clock offsets). Peak at the edge: pair unusable.
8. **Fine lags.** For K, PP and PPP, normalized cross-correlation of record B against A within c +/- 0.6 s
   on a 100 Hz grid (linear interpolation of the filtered traces), parabolic peak. tau_X = lag - c.
   **dK = tau_K - tau_PP**, **dP3 = tau_PPP - tau_PP** (B minus A; positive = later = slower).
9. **Usable.** K-usable: cc_K >= 0.80 and cc_PP >= 0.80, neither at the search edge, SNR_K and SNR_PP >= 3
   in both records. P3-usable: K-usable, record covers PPP+30 s, PPP unclipped, cc_PPP >= 0.70, not at edge,
   SNR_PPP >= 3 in both.
10. **Pairs and year terms.** Primary pair set: the 18 starred tests, every pair within 0.05 deg (Table S1
    epicentres), as the paper. Least-squares year terms from the pair differences, 1977 = 0; pairs within one
    year are skipped; years not linked to 1977 are dropped. Epoch values (1982-83, 1988-90, 1995) = means of
    year terms weighted by the number of tests that year that are in a usable pair. Bootstrap over usable
    pairs (2,000 resamples, seed 31) gives epoch errors.

## Gates

- **G0 (K testable):** >= 10 K-usable pairs, the 1977 test in at least one, and all three epochs present.
- **G1 (P3 testable):** >= 10 P3-usable pairs, all three epochs present for both K and P3 in that pair set,
  >= 20 valid bootstrap resamples.

## Grades

**K reproduction** (epochs against the paper's midpoints -0.125 / -0.150 / +0.175 s; same arithmetic as P30):
- *reproduced*: all three signs right (|value| >= 0.02 s) and rms(paper - ours) <= 0.05 s;
- *partly*: >= 2 signs right and rms < 0.75 x rms(paper) (= 0.114 s);
- *not-reproduced*: anything else; *untestable* if G0 fails.

**P3 discrimination** (P3-usable pairs; K and P3 epochs from the same pairs). For each hypothesis H with
ratio r_H (core 0, position r_pos, time-base r_tb): chi2_H = sum over the three epochs of
(P3_e - r_H K_e)^2 / var_e, where var_e is the bootstrap variance of P3_e - r_H K_e (floor (5 ms)^2).
H "fits" if chi2_H <= 7.815 (95%, 3 degrees of freedom).
- grade = the single hypothesis that fits (*core*, *position* or *time-base*); *unresolved* if more than
  one fits; *none-fits* if none; *untestable* if G1 fails.

Headline: "K: <grade>; P3: <grade>". If K is not reproduced, all three hypotheses predict P3 near 0 and
the P3 grade is expected to be *unresolved*; it is reported as it comes.

## Report-only (no grade changes)

- All 23 tests instead of the 18 starred.
- The SS system as an independent replica (same rules) where SS records exist.
- Per-pair slope of dP3 on dK (least squares through 0 and total least squares).
- Per test: system, channels, sample rates, start time relative to origin, record length, glitch counts.
- Every pair's lags, correlations, SNRs and the channel used.

## Forecasts (before data)

- Archive reachable from the job and records readable: 60%.
- G0 passes, given readable records: 60%.
- K: reproduced 15%, partly 30%, not-reproduced 55% (given G0).
- G1 passes, given G0: 40% (PPP is weak at 0.3-1 Hz for mb 5-6).
- P3, given G1: core 25%, position 25%, time-base 15%, unresolved or none-fits 35%.

## What each outcome would mean

- **K reproduced + P3 core:** independent support, from the same public data, for a change inside the
  outer core: remarkable, and the author would be told first.
- **K reproduced + P3 position or time-base:** the reported change is real in the records but comes from
  where the tests were or from the recording clock, not from the core.
- **K not reproduced:** the public waveforms, measured with the paper's pairing and correction, do not
  show the reported history.
- **untestable:** reported as such; a later version may change only what the plan names as failing, by
  public amendment before results.

## Limits

- Broadband cross-correlation in 0.3-1 Hz, not the paper's single-frequency multitaper at 0.8 Hz.
- K is the whole first core-arrival group (PKIKP window), not separated into PKrKP and PKiKP.
- 1-D iasp91 times and ray parameters; Table S1 epicentres; no relocation.
- Glitch and clip rules are simple and fixed; a record with long glitches may still pass, and a glitch in
  the first or last three samples of the processed span is not detected (far from every phase window).
- One station, one or two systems; no instrument response removed.

Credit: Ying Zhou set up the question and the event selection; the records exist thanks to the Borovoye
archive rescue at Lamont-Doherty (Kim, Richards, Ekstrom and colleagues).

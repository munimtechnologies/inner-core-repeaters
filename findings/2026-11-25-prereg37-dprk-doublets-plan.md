# Pre-registration 37: do inner-core waves from the North Korean explosions change between tests?

MV2, Munim, Inc. (mv2@munimtech.com). Written 2026-11-25 by AI agents working in MV2's name; no human
expert has reviewed it. Text CC BY 4.0; code Apache-2.0. Repository:
https://github.com/munimtechnologies/inner-core-repeaters

**Blindness.** This plan was written before any output of the blind probe `dprk37probe.py` 0.1.1
(commit 0fdf6cb, job dprk37probe-v2, still queued when this was written) existed. No DPRK waveform,
SNR, station count or pair count had been seen by anyone working on this repository when the
measures, gates, pairs and decision rules below were fixed. This is stricter than the scoping note
(findings/2026-11-25-line-x-dprk-explosion-doublets-scoping.md), which planned to write P37 after
the probe.

**One amendment window.** When the probe returns (binned single-record SNRs and pair counts only,
never a comparison of two tests), ONE dated amendment may change only: the primary band (to the
probe's other band, 1.0-2.0 Hz), the minimum station counts in section 6, and which pair is G2
(only to another pair of the four tests 2013-02, 2016-01, 2016-09, 2017-09). Nothing else changes,
and no amendment is allowed after the P37 job has started.

## 1. Question and why it matters

Recent inner-core studies disagree about what happened after about 2008: a slow backward rotation
(Wang, Vidale et al. 2024, Nature 631:340), near standstill after an eastward phase (Wu, Song & Yang
2026, GRL 53:e2025GL119295), or local, episodic change of the inner-core surface (Zhang & Wen 2026,
JGR Solid Earth, 10.1029/2026JB034595; both: Vidale et al. 2025, Nat. Geosci.). Every one of these
uses earthquake doublets, and critics (Yang et al. 2026, Commun. Earth Environ.,
10.1038/s43247-026-03961-7) argue that doublet differences can come from the sources themselves.

The six underground tests at Punggye-ri (2006-2017) are explosions at one site with relative
positions published to about 100 m. Stations in South America and Antarctica lie 150-179.5 deg away,
where PKIKP (DF, through the inner core) and PKPab (AB, outer core only) arrive on the same record
seconds apart. AB is a built-in control for the source, the receiver, clocks and the mantle. A web
and Semantic Scholar search on 2026-11-25 found no study that has used these tests for temporal
change of Earth's core (searches listed in section 9).

Wow, if it works: the first inner-core time-lapse test with 21st-century explosions, in exactly the
years the rival readings disagree, with an 8-month pair (2016-01 / 2016-09) that directly tests
"episodic" change on a timescale no rotation model predicts.

Honest limit fixed now: near-antipodal DF rays are weakly sensitive to rigid rotation (they cross
the inner core almost through its centre), so a null result says little about rotation. They are
sensitive to change where they cross the inner-core boundary (ICB).

## 2. Data (fixed)

- Sources: USGS ComCat ids usc000f5t0 (2013-02-12 02:57:51.49, mb 5.1), us10004bnm (2016-01-06
  01:30:01.48, mb 5.1), us10006n8a (2016-09-09 00:30:01.44, mb 5.3), us2000aert (2017-09-03
  03:30:01.76, mb 6.3), plus usp000eurb (2006) and usp000gxgc (2009), coordinates as in
  dprk37probe.py 0.1.1. The 2017-09-03 03:38 collapse is never used.
- Stations and records: exactly the probe's selection and cache (EarthScope FDSN, BHZ before HHZ,
  location 00, '', 10; 150-179.5 deg; records DF-90 s to AB+120 s; cache layout
  /data/dprk37/<eventid>/NET.STA.LOC.CHA.mseed, data_key dprk37). The P37 job loads
  dprk37probe.py at commit 0fdf6cb and fetches any missing record the same way.
- Predicted times: ak135, TauP, surface source, as in the probe (`branches`).

## 3. Pairs

- **G1 (primary): 2016-01-06 / 2016-09-09** (8 months; mb 5.1 / 5.3).
- **G2 (secondary): 2013-02-12 / 2016-01-06** (2.9 years; mb 5.1 / 5.1, the closest sizes).
- Report-only: 2013-02 / 2016-09, every pair with 2017-09 (yield about ten times larger, so the
  source spectrum differs most), and every pair with 2006 or 2009 (small events).

## 4. Station gates for a pair (all must hold)

1. Both records "ok" in the probe sense (one trace, no gap, covers the window).
2. The same location code, channel and channel-epoch start for both tests (no instrument change).
3. windows_clear (probe definition): AB-1 s and BC-1 s fall after DF+4 s.
4. SNR (probe definition: rms of [phase-1, phase+4] s over rms of [DF-60, DF-5] s, 4-pole causal
   band-pass) >= 5 for DF and >= 5 for AB on both records, in the primary band 0.7-2.0 Hz.
5. After alignment (section 5): cc_AB >= 0.80.

## 5. Measurements (per gated station, per pair)

Processing: demean, linear detrend, 5% taper, 4-pole causal band-pass 0.7-2.0 Hz (secondary
1.0-2.0 Hz, report-only), cubic interpolation to 100 Hz. No response removal (gate 2 makes
responses equal). Predicted DF and AB times come from the fixed site (41.30N, 129.05E) for both
tests, so catalogue location noise never enters the measurement.

- **AB alignment:** normalised cross-correlation of [AB-0.5, AB+5] s, lags within +-1.0 s, parabolic
  sub-sample peak (lag and value): tau_AB, cc_AB.
- **DF:** window [DF-0.5, DF+4] s, lags within tau_AB +- 0.4 s, parabolic peak: tau_DF, cc_DF.
- A station is dropped for the pair if either best lag lies on the edge of its search range, or a
  window falls outside a record.
- **dt = tau_DF - tau_AB** (s, sign: positive means DF later in the second test relative to AB).
- **Excess decorrelation X.** For a record with SNR r, the expected best correlation of two noisy
  copies of the same signal is sqrt(f_a f_b) with f = max(0, 1 - 1/r^2). D_DF = 1 - cc_DF,
  E_DF = 1 - sqrt(f_DF,a f_DF,b), and likewise for AB. **X = (D_DF - E_DF) - (D_AB - E_AB).**
  X > 0 means DF changed more than AB beyond what noise explains.
- Source-position budget per station, reported: |p_AB - p_DF| (s/deg, TauP) x 0.0027 deg (300 m,
  three times the published relative precision).

## 6. Pair statistics and labels (G1 and G2 separately)

n = number of gated stations. Bootstrap over stations: 2000 resamples, Python random.Random(37);
the CI is the 2.5% and 97.5% linear-interpolation quantiles of the resampled statistic.

- **Untestable** if n < 8.
- **W (waveform):** m_X = median X, CI = bootstrap 2.5-97.5% of the median.
  "changed" if m_X >= 0.03 and CI lower > 0.01; "same" if CI upper < 0.03; else "inconclusive".
- **T (timing):** fit dt = c + a cos(baz) + b sin(baz) by least squares over stations (baz = back-
  azimuth from station to source; this soaks up any source-position shift, whose effect varies as a
  cosine of azimuth; with n < 8 the fit is skipped and the pair is untestable anyway). CI of c from
  the same bootstrap (refit per resample).
  "shift" if |c| >= 0.020 s and the CI excludes 0; "no shift" if the CI lies within
  [-0.020, +0.020] s; else "inconclusive".
- **Headline per pair:** "DF changed" if W = changed or T = shift; "no change at this level" if
  W = same and T = no shift; otherwise "inconclusive". The headline of P37 is G1's.

## 7. P37b: the patch Zhang & Wen report as changing

For every station, the job computes the two points where DF crosses the ICB (ak135 ICB depth
5153.5 km; TauP get_pierce_points for PKIKP; positions along the source-station great circle). It
computes, the same way, the ICB crossing points of DF from Zhang & Wen's multiplet F centroid
(58.7887S, 25.3637W, surface source) to YKA (62.493N, 114.605W) and to ILAR (64.7714N, 146.8866W):
four reference points. A station is a **patch station** if either of its DPRK-DF ICB points lies
within 10.0 deg (angular distance) of any reference point; otherwise a non-patch station.

- Untestable for a pair if fewer than 4 patch or fewer than 4 non-patch stations are gated.
- d = median X (patch) - median X (non-patch); p = share of 10,000 label permutations
  (Python random.Random(37)) with a difference >= d.
- "patch excess" if d >= 0.03 and p <= 0.05; "no patch excess" if d <= 0.01; else "inconclusive".
- Reading against Zhang & Wen's timeline uses only their episode dates as printed, quoted in the
  results note; no re-dating.

## 8. Forecasts (written before any data)

- G1 testable: 55%. Given testable: W changed 15%, same 45%, inconclusive 40%; T shift 10%,
  no shift 55%, inconclusive 35%; headline "DF changed" 20%.
- G2 testable: 50%; given testable, headline "DF changed" 25%.
- P37b testable for G1: 30%; given testable, patch excess 20%.

## 9. Prior work (checked 2026-11-25)

- Searches: Google-type web search "North Korea nuclear test PKIKP inner core temporal change
  repeating explosions South America"; web search "DPRK explosions PKP inner core temporal change OR
  time-lapse 2026"; web search "Punggye-ri PKiKP OR PKIKP inner core boundary nuclear test antipodal
  Argentina Chile stations"; Semantic Scholar citations of Wang et al. 2024 (14 papers in 2026; none
  uses explosions) and of Zhang & Wen 2026 (no citations yet). Earlier runs: five more searches.
- Closest prior work: explosion-based inner-core change with Cold War tests only (Vidale & Earle
  2005, Mururoa at NORSAR; Vidale et al. 2019, Novaya Zemlya at LASA; Wang & Vidale 2022, Amchitka
  at LASA); Zhou 2026 (JGR SE, 10.1029/2026JB033798, French tests at Borovoye, outer core); Butler
  2026 (EPSL 683:119963, antipodal waves, static ICB structure, not temporal).
- Simulations of what changes look like: Wang, Pang & Vidale 2026, EPSL 692:120274.

## 10. Code and order of steps

1. This plan (its own commit). 2. dprk37.py 0.1.0 and test_dprk37.py (a separate commit; the unit
tests are the job's first run step). 3. The blind probe returns; the one amendment, if any, is
committed. 4. Only then is the job dprk37-v1 filed (data_key dprk37, records already cached). The
results note reports every label exactly as computed, including untestable.

## 11. Limits fixed now

Near-antipodal DF is weakly sensitive to rigid rotation. The tests differ in yield, so source
spectra differ; DF is lower-frequency than AB after inner-core attenuation, so a source difference
probably decorrelates AB more than DF and biases X toward "same". Single stations are noisier than
the arrays used by the rival studies. The thresholds (0.03 in X, 20 ms in c) are judgements, not
calibrated on these paths.

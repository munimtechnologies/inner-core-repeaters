# Line X: North Korean nuclear tests as inner-core doublets (scoping note, probe only)

Date: 2026-11-25. Repository: munimtechnologies/inner-core-repeaters (Apache-2.0; text CC BY 4.0).
Made with AI agents working in MV2's name (MV2, Munim, Inc., mv2@munimtech.com). No human expert has checked this yet.

This note is **not** a pre-registration. It records, before any waveform is read, a new question and the blind
feasibility probe (`dprk37probe.py`, job `dprk37probe-v1`) that will decide whether a pre-registration (number 37)
is worth writing.

## The question

Every test of whether Earth's inner core changes over years rests on "repeating" sources. Earthquake doublets are
never exactly repeats: their positions, depths and rupture can differ (that is the core of the 2026 artifact
critiques, Yang et al. 2026, https://doi.org/10.1038/s43247-026-03961-7, and of our own gate G2,
findings/2026-11-24-prereg35-g2-source-check.md). Nuclear explosions are cleaner sources, and they have been used
before for the inner core, but only Cold War tests: Novaya Zemlya and Amchitka recorded at LASA (Vidale et al. 2019,
https://doi.org/10.1029/2019GL083774; Wang & Vidale 2022), Mururoa at NORSAR (Vidale & Earle 2005) and Mururoa
at Borovoye (Zhou 2026, https://doi.org/10.1029/2026JB033798, the outer core).

North Korea's six underground tests at Punggye-ri (2006-10-09, 2009-05-25, 2013-02-12, 2016-01-06, 2016-09-09,
2017-09-03; USGS ComCat ids usp000eurb, usp000gxgc, usc000f5t0, us10004bnm, us10006n8a, us2000aert) are
explosions at one site, with relative positions known to about 100 m from published relocations. Seen from them,
Chile, Argentina, Bolivia, southern Brazil and the Antarctic Peninsula lie 150-180 degrees away, where PKIKP (DF,
through the inner core) and PKPab (AB, outer core only) arrive tens of seconds apart on the same record. AB shares
the source, the site and the mantle near both ends with DF, so it is a built-in control. The EarthScope station
service lists 160 broadband stations at these distances for the January 2016 test (51 in Chile's C1 network,
18 in C, 26 in a 2016-2018 temporary array ZR, plus GT, IU, II, G, BL, ON, AI, WA; query of 2026-11-25).

The question for pre-registration 37: **between 2013 and 2017, do these explosions' DF waves change relative to
their AB waves?** These are the years in which the published readings disagree: steady backtracking 2008-2023
(Wang et al. 2024, https://doi.org/10.1038/s41586-024-07536-4), a near standstill after about 2010 (Wu, Song &
Yang 2026, https://doi.org/10.1029/2025GL119295), and episodic change of the inner-core surface (Zhang & Wen 2026,
https://doi.org/10.1029/2026JB034595). The 2016-01 / 2016-09 pair, 8 months apart with nearly equal magnitudes
(mb 5.1 and 5.3), is the cleanest short-interval test of "episodic" change with sources that cannot be relocation
artifacts.

Prior-work check (2026-11-25): four web searches ("North Korean nuclear tests PKIKP inner core temporal change",
"North Korea nuclear explosions core phases PKP inner core rotation doublet South America", "Wen Zhang DPRK
nuclear tests PKiKP PKIKP inner core boundary change", "nuclear explosion doublet inner core Punggye-ri DPRK
PKP-DF PKP-AB differential travel time") found no study that uses the North Korean tests for inner-core or
outer-core temporal change. To be repeated, with Semantic Scholar citation checks, before any article.

## Honest expectations

- These paths are nearly antipodal, so DF turns deep in the inner core, close to Earth's centre. A rigid rotation
  of a few tenths of a degree moves structure there very little; the paths are weakly sensitive to rotation and
  more sensitive to change near where DF crosses the inner-core boundary (under north-east Asia and under southern
  South America). A null result would therefore say little about rotation and more about local surface change.
- The 2017 test (mb 6.3) was 1.0-1.2 magnitude units larger than 2013 and 2016, so its source spectrum differs at 1-2 Hz;
  DF relative to AB cancels most of that, but the 2017 pairs will be report-only unless the probe shows otherwise.
- The 2006 (mb 4.3) and 2009 (mb 4.7) tests are probably too small for clear DF at these distances.

## The probe (blind)

`dprk37probe.py` 0.1.0 fetches, for each test, every BHZ/HHZ record at 150-179.5 degrees from EarthScope FDSN
(BHZ before HHZ; location 00, empty, 10), and computes on each single record the rms signal-to-noise ratio of the
DF and AB windows (-1 to +4 s around the ak135 prediction) against -60 to -5 s before DF, in 0.7-2 Hz and 1-2 Hz.
It writes per-record SNR, channel epoch and whether the DF, BC and AB windows overlap, and for each pair of tests
only counts of stations where both records pass SNR >= 5. **It never compares two tests**: no correlation, time
difference, amplitude ratio or overlay of two records is computed, written or plotted. Records are cached in the
job's data volume (key `dprk37`) for the later test.

Pre-registration 37 will be written only from these counts, and will fix the measure (DF-minus-AB differential
time change and DF coda dissimilarity relative to AB coda dissimilarity), the station set, the graded pairs, the
decision rules and forecasts, before any two records are compared.

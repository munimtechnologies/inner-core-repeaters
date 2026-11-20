# Second path pilot (chilepilot): can Chile -> Central Asia arrays test the inner-core clock?

*inner-core-repeaters, MV2 of Munim, Inc. Written 2026-11-20, before any waveform for this
pilot was downloaded. Code: `chilepilot.py` (tests: `test_chilepilot.py`). Licence: code
Apache-2.0, text CC BY 4.0. Made with AI agents working in MV2's name; no human expert has
checked it yet.*

## Why a second path

The 2024 *Nature* claim that Earth's inner core turned around in about 2008.5 and has been
"backtracking" since (Wang, Vidale et al. 2024, 631:340) rests on repeating South Sandwich Islands
(SSI) earthquakes recorded at two Alaskan/Canadian arrays, ILAR and YKA. Our own pre-registered
tests on that path have now reached its data limit:

- YKA's multi-family time trend is explained by the array's 2013 hardware change
  (findings/2026-11-19-prereg17-results.md, prereg18-results.md).
- ILAR's 2018 digitiser change is negligible (findings/2026-11-20-ilar-2018-digitiser-results.md),
  but a source-controlled ILAR clock test has too few measurable pairs
  (findings/2026-11-20-ilar-p19b-power-check.md).
- The single-family M03-24 result (core waves split the family; mantle waves do not) stands.

If the inner core turns as a rigid body, its clock (when it reoccupies a past position) is the
same on every path. If the changes are local deformation of the inner-core surface, different
paths need not agree. A source-controlled test on an independent path is the cleanest way to tell
these apart.

## Prior work on other paths (searched 2026-11-20)

- Zhang et al. 2008 (*GJI* 174:993): one Aleutian doublet at BOSA (South Africa) and one Kuril
  doublet at BDFB (Brazil), about 0.1 s change over 6-7 years.
- Yang & Song 2020 (*JGR* 125, e2019JB018652): systematic global doublet search; large changes
  under four Northern-Hemisphere regions; a Chile -> Kazakhstan path (station CHKZ) with steady
  changes and a sharp lateral velocity gradient in the top of the inner core.
- Yang & Song 2023 (*Nat. Geosci.* 16:182): paths that changed before about 2009 show little
  change since; a multidecadal oscillation.
- Yu 2016 (*GJI* 204:1204): Tonga-Kermadec-Vanuatu repeaters recorded in Europe; mostly null.
  Deep Tonga repeaters are rare (Yu & Wen 2012, *BSSA*).
- Wu, Song & Yang 2026 (*GRL*, 10.1029/2025GL119295): event-based sequences on SSI -> COLA and
  Kuril -> PLCA (Argentina), 1982-2024; an oscillation with no abrupt changes.
- Tkalčić et al. 2024 (*GRL*): data selection of doublets can change the inferred period.

What we did not find: a source-controlled (mantle-wave verified), multi-family, array-beam test
of the backtracking clock on the Chile -> Central Asia path that includes 2009-2026 earthquakes,
with each array's hardware changes measured. That is the gap this line aims at.

## Why Chile -> Central Asia

- Same distances as SSI -> ILAR/YKA (about 143-165 deg), so PKIKP, its coda and PKPbc can be
  measured exactly as in our ILAR/YKA code (1-2 Hz, coda 1-6 s after PKIKP).
- Three arrays with open data at EarthScope: Makanchi (KZ MK01-MK09, short-period, from 2000),
  Karatau (KZ KK01-KK09, from 2001) and KNET (KN, broadband, from 1991). No restricted flag in
  the station metadata. MKAR's broadband channel changed sensor on 2011-08-08 (KS-54000 to
  CMG3-TB); the short-period elements show one epoch from 2000.
- Many sources: ComCat lists 638 M >= 5.0 events (depth <= 300 km) between 18S and 42S before
  2008.5 since 1995, and 1,047 after 2008.5 between 26S and 42S alone (counts fetched
  2026-11-20; the 2010 Maule and 2015 Illapel sequences dominate the later ones).

## What the code does

1. ComCat: 18-42 S, 66-76 W, depth <= 300 km, M >= 4.5, 1995-01-01 to 2026-11-20.
2. Four magnitude bins [4.5, 5.0), [5.0, 5.5), [5.5, 6.0), [6.0, 9.9). In each bin, up to 6 events
   before 2008.5 and 6 after, evenly spaced through the time-ordered list (deterministic).
3. For each event and array: elements with a BHZ/SHZ/HHZ epoch covering the event; array centre
   distance must lie in 143-165 deg. ak135 PKIKP predicted per element (TauP). One bulk
   waveform request per event and array.
4. Each element: response removed to velocity, causal 4-pole band-pass (primary 1.0-2.0 Hz as
   ILAR/YKA; alt 0.6-1.6 Hz), 20 Hz, aligned on its own predicted PKIKP, normalised by the rms of
   -45 to -10 s. Beam = mean of at least 3 elements, renormalised the same way.
5. Read-outs per beam: rms SNR in -3 to +7 s (PKIKP and PKPbc at these distances), in the coda
   window 1-6 s, and in a noise-only control window -58 to -48 s; median element SNR.
6. Catalogue proxy for repeaters: pairs with both M >= 5.0, epicentres within 5 km, |dM| <= 0.2
   and |dz| <= 10 km, counted overall and straddling 2008.5 (and with both M >= 5.5). ComCat
   locations are uncertain by 5-20 km, so this is a loose proxy, not a repeater count.
7. Outputs: `chile_rows.csv`, `chile_summary.json`, `plots/chile_snr_vs_mag.png`, the catalogue.

## Decision rule (fixed now)

- **Array verdict** (primary band, PKP window): the lowest bin with at least 4 beams whose median
  beam SNR >= 3: `testable-from-M<lo>` or `none`.
- **Availability**: per array, the fraction of in-range picked events with a beam, before and
  after 2008.5. An array counts only if both fractions are >= 0.6.
- `build-catalogue`: the best counting array is testable from M4.5 or M5.0. Next: a Chile
  multiplet catalogue (PKPbc and mantle-P similarity, the source-score method of
  pre-registration 13), the array's hardware epochs measured, then a pre-registered Chile clock
  test (pre-registration 20) before any PKIKP pair is measured.
- `build-catalogue-m5.5`: testable from M5.5 and at least 30 straddling M >= 5.5 close pairs.
- Otherwise `close-chile-path`, reported as such.

## Forecasts (made before any data)

| # | Forecast | Probability |
|---|---|---|
| F1 | At least one array has beams for >= 70% of in-range events both before and after 2008.5 | 60% |
| F2 | Best array verdict: M4.5 / M5.0 / M5.5 / M6.0 / none | 15 / 35 / 25 / 10 / 15% |
| F3 | Over beams with >= 5 elements, the beam SNR beats the median element SNR in >= 75% | 75% |
| F4 | Noise-only control: median beam value in [0.8, 1.3] | 85% |
| F5 | At least 300 straddling close pairs (M >= 5.0) in the catalogue | 55% |

## Limits known in advance

- Alignment uses ak135 predictions only; 3-D structure can misalign elements by a few tenths of a
  second, which lowers beam SNR at 1-2 Hz (most for the 200 km KNET aperture).
- 6 events per bin and side is a small sample; ComCat magnitudes mix mb and Mw.
- PKPbc falls inside the PKP window at these distances, so the window SNR is not DF alone.
- Northern Chile events reach KKAR and KNET at 139-146 deg, near the PKP caustic; the distance
  gate removes the nearest of them.
- This pilot makes no claim about the inner core.

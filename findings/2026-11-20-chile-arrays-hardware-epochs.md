# Chile -> Central Asia arrays: documented hardware epochs (prep for pre-registration 20)

*inner-core-repeaters, MV2 of Munim, Inc. Written 2026-11-20, while `icr-chilepilot-v1` waits for
compute and before any Chile-path waveform was examined. Text CC BY 4.0. Made with AI agents
working in MV2's name; no human expert has checked it yet.*

## Why this comes first

On the SSI -> YKA path an apparent inner-core time trend turned out to be the array's 2013
recording-system change (findings/2026-11-19-prereg18-results.md, median hardware-only coda
difference D_hw 0.0117). ILAR's 2018 digitiser change was measured as negligible (D_hw 8e-5,
findings/2026-11-20-ilar-2018-digitiser-results.md). Our rule since then: measure the instrument
before believing a trend. For the Chile -> Central Asia path (findings/2026-11-20-chile-path-plan.md)
this note maps every documented hardware epoch of the three candidate arrays, so that
pre-registration 20 can be written around them.

## Source

EarthScope FDSN station service, `level=channel`, `format=text`, fetched 2026-11-20:
`net=KZ&sta=MK*,KK*&cha=?HZ` and `net=KN&cha=?HZ`. Table: `results/chile-hw-epochs/chile_array_epochs.csv`
(vertical channels used or considered for beams; KN HHZ and LHZ epochs follow the same switch
dates as BHZ and are not listed; KKL, a 1997-98 temporary KZ station, and ERPT, ended 2000, are
omitted).

## What the metadata shows

| Array | Beam elements | Documented changes 1995-2026 on those elements |
|---|---|---|
| Makanchi (KZ) | MK01-MK09 SHZ, Geotech GS-21, 40 sps | none: one epoch each from 2000-08-28, open |
| Karatau (KZ) | KK01-KK09 SHZ, Geotech GS-21, 40 sps | none: one epoch each from 2001-12-11, open |
| KNET (KN) | 11 stations, STS-2 BHZ, 40 sps | one, on 1995-01-05/06 (digitiser RT72A-02 -> RT72A-08, 20 -> 40 sps), then one open epoch |

Changes exist only on channels the pilot does not beam:

- MKAR BHZ (the array's broadband site): KS-54000 -> CMG3-TB on 2011-08-08; MK32 (GS-13) ends
  2011-08-07.
- KKAR BHZ: KS-54000 -> CMG3-TB on 2018-09-20.

`chilepilot.py` beams only MK01-09, KK01-09 and the KN stations, so neither change enters its beams.

## What it means for pre-registration 20

1. On the metadata, the short-period Makanchi and Karatau beams have **no documented hardware
   boundary anywhere in 2000/2001-2026**, so a reoccupation-clock test across the 2008.5 turn on
   these arrays does not have to cross one. That is the property YKA lacked (2013) and ILAR
   only partly had (2018, negligible).
2. KNET is metadata-flat from 1995-01-06. P20 should start KN pairs on 1995-01-07.
3. MKAR and KKAR broadband channels, and MK32, stay out of every P20 beam.
4. **Metadata-flat is not hardware-flat.** An undocumented digitiser, gain or timing change would
   not show here (the YKA lesson was a documented change; the open risk is an undocumented one).
   P20 will therefore include a data-side check fixed in advance: on each array, the median
   element-to-beam coherence and the beam noise spectrum (-45 to -10 s before PKIKP, 1-2 Hz) as a
   function of year, with a step of more than 3 robust standard deviations treated as a candidate
   hardware boundary that splits pairs, exactly as the YKA 2013 boundary did.
5. No `ilresp`-style response comparison is needed for MK/KK (one response each); for KN the
   only switch predates the study window.

## Limits

This is a metadata map, not a measurement. It makes no claim about the inner core. Whether the
arrays have the signal and the data to test the clock at all is decided by `icr-chilepilot-v1`
under its own pre-registered rule.

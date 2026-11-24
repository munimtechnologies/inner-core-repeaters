# Rival-paper check, 2026-11-24 (evening): two September 2026 papers we had missed

MV2, Munim, Inc. Text CC BY 4.0. This was written with AI agents working in MV2's name, and no
human expert has checked it yet.

We ran this check under our own rule: look for new rival papers every few runs. It found two
papers from September 2026 that bear on work already in this repo. Earlier runs had not seen
either one. We record them here so that our own claims are read in their light.

## Searches run (2026-11-24)

- Crossref: `query.bibliographic=inner core temporal change repeating earthquakes`, from 2026-10-01.
- arXiv export API: `all:"inner core" AND (rotation OR repeating)`, newest first.
- Semantic Scholar: the citations of Zhang & Wen 2026 (JGR SE, 10.1029/2026JB034595). It has none yet.
- Semantic Scholar: the citations of Wang et al. 2024 (Nature, 10.1038/s41586-024-07536-4). This search found both papers.
- Semantic Scholar: the citations of Zhang & Wen 2025 (Nat Commun, 10.1038/s41467-025-63814-3) and of Zhou 2026 (JGR SE, 10.1029/2026JB033798). The second has none yet.
- Web search for the exact title of Yang et al. 2026.

## 1. Yang, Li, Ringler, Jiang, Wu & Song (2026)

"Evaluation of seismological evidence for decadal changes in Earth's lowermost mantle and outer
core", *Communications Earth & Environment*, published 2026-09-21,
[doi 10.1038/s43247-026-03961-7](https://doi.org/10.1038/s43247-026-03961-7) (open access, CC BY).

- **What the abstract says:** the authors re-examined the earthquake-pair waveforms from two
  recent studies. They tested source repeatability independently and checked station quality
  against co-located sensors. They conclude that the reported changes "can be largely explained by
  artifacts".
- **The two studies:** their reference list includes Zhang & Wen 2025 (Nat Commun, the ULVZ
  change at IU.SDV) and Zhou 2022 (Commun Earth Environ, Kermadec SKS). These are the only
  earthquake-pair studies of deep-mantle or outer-core change that it cites, so they are almost
  certainly the two it re-examined. We have not yet read the full text, and we will confirm this
  before relying on it.
- **Zhou 2026 is not covered:** the list does not cite Zhou 2026 (Borovoye / Mururoa PKP).

**What this means for our work:**

- **Our article [ulvz-decadal-change-sdv-not-reproduced](https://www.munimtech.com/research/ulvz-decadal-change-sdv-not-reproduced)**
  (tests 23-24, SDV) does not cite Yang et al. 2026, and it should have. Yang et al. published
  first that the SDV change can be explained by artifacts. The part of our article that may still
  be new is the specific mechanism we measured. The authors' published 2009 traces sit in a frame
  rotated 3.6 degrees from their 2000 frame, while the sensor turned only 0.55 degrees. Both
  public and published traces give shifts of the opposite sign (+21 to +28 ms against the
  reported -160 ms). Once the full text has been read, we will update the article (same slug). It
  will cite Yang et al. 2026 as the closest prior work and claim only what they did not show.
- **Pre-registration 32** (Zhou 2022 Kermadec SKS, job ks-sks32-v1, still queued) is no longer
  the first independent re-examination of that claim. If it runs, its grade will be reported as a
  second, pre-registered check, compared point by point with Yang et al. 2026.
- **Our Borovoye article** ([outer-core-change-test-site-position](https://www.munimtech.com/research/outer-core-change-test-site-position))
  tests Zhou 2026, which Yang et al. do not treat. Its main point, about source position and
  repeatability, is in the same spirit as Yang et al., so it will cite them as related work.

## 2. Ackerley (2026)

"Detection and Modeling of Timing Errors at Seismic Arrays", *Seismological Research Letters*,
2026-09-01, [doi 10.1785/0220260111](https://doi.org/10.1785/0220260111) (CC BY).

- **The errors:** some elements of the Yellowknife array (YKA) had relative timing errors of up to
  0.4 s. From 2014 to 2020 the errors came in steps of 0.125 s (five samples at 40 sps). They grew
  until the datalogger was reset.
- **The fix:** none have appeared since a digitizer firmware upgrade on 8 April 2020.

**What this means for our work:**

- **Our YKA beams are immune by construction.** `yka.py` 0.2.0 (run by `m0324_check.py` for the M03-24 labels) estimates per-element,
  per-event statics within +/-0.5 s against the event's own PKIKP-group beam before it stacks,
  and `zw33.py` aligns elements to the event beam the same way. A timing error that is constant
  within one event's window (a few minutes) is removed. Ackerley's errors only change at resets,
  so this applies.
- **The one member inside the window:** among the M03-24 members, only the 2014-12-09 event (our
  "Y" label) falls in Ackerley's 2014 to April 2020 window. That label sets the lower edge of the
  one-turn window in amendment 1b (`wsy1b`). The per-element statics stored for that event are an
  independent check: on the affected elements they should show multiples of 0.125 s.
- **What Ackerley's table would add:** we will compare the table, once read, with those statics
  before we rely on the 2014-12-09 label again.
- **Zhang & Wen's data are outside the window.** Their multiplets M03-24 (2003, 2008, 2011,
  2020-07-11, 2022, 2024) and F (2005, 2012, 2022) have no event between 2014 and 2020-04-08.

## Also noted

The Nature paper "Gravitational torque drives multidecadal variations in length of day"
(2026-09-23, [doi 10.1038/s41586-026-10999-2](https://doi.org/10.1038/s41586-026-10999-2)) cites
Wang et al. 2024. It is relevant background for the rotation reading and does not change any of
our tests.

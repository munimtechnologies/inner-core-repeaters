# Pre-registration 27: does SKS at IU.SDV change over a decade on other earthquake pairs that leave the core near the same exit point?

Written 2026-11-21 ~23:40 UTC, before any catalogue query or waveform for this test was downloaded. Code:
`sdvexit27.py` 0.1.0 (Apache-2.0). The job waits until this file is public on `main` and records its
sha256. Made with AI agents working in MV2's name; no human expert has checked it yet. Text CC BY 4.0.

## Claim under test

Zhang & Wen 2025, *Nature Communications* 16:8754 (doi 10.1038/s41467-025-63814-3) report a decade-scale
change of an ultra-low velocity zone (ULVZ) at the core-mantle boundary (CMB). It rests on one Aleutian
doublet (2000-07-07 Mw 5.9 and 2009-12-17 Mw 5.8, their AI-0009) recorded at one station, IU.SDV in
Venezuela (94.9 deg): radial SKS 160 ms earlier in 2009. Their models put the change near the SKS
**entry** point (source side). They also say a change near the SKS **exit** point (receiver side) "could
explain the seismic data equally well".

- Pre-registration 23 re-measured SDV from the public records: +28 ms, graded *not-reproduced*.
- Pre-registration 24 (the authors' own traces) and 25 (onset and pulse width) are queued.
- Pre-registration 26 (queued) tests the entry side with other stations whose SKS enters the core near
  SDV's entry point.

This test covers the **exit side** and does not depend on the outcome of 24-26.

## Idea

An SKS leg leaves the core about 10 deg of arc from the receiver, along the back-azimuth. Other earthquakes
at a similar back-azimuth and distance from SDV, but hundreds of km along the arc from the doublet, send
SKS out of the core within 100 km of the doublet's exit point while entering it far away. A receiver-side
ULVZ change big enough to move the doublet's SKS by 0.1-0.16 s, ongoing in time as the paper argues
("vigorous" flow), should also show in SDV records of **other repeating pairs** on those paths spanning
five years or more. If none of them changes, the exit-side reading loses its support, independently of how
the doublet itself is measured. Pairs at other back-azimuths are controls. Together they also show whether
SDV itself (sensor, clock, the 2009-09-18 metadata step) produces apparent SKS changes.

## Data and processing (fixed in `sdvexit27.PARAMS`)

- **Catalogue:** USGS ComCat CSV, earthquakes M >= 5.3, depth <= 100 km, 85-125 deg from IU.SDV.
  - Earlier events: 1999-11-06 to 2009-09-17 (IU.SDV location 00, before the metadata step).
  - Later events: 2009-12-17 to 2026-11-01.
  - The query files' sha256 are recorded.
- **Candidate pairs** (A earlier, B later):
  - epicentres <= 15 km apart, catalogue depths within 25 km, magnitudes within 0.8;
  - B - A >= 5 years;
  - the doublet itself is excluded here (it is the G0 check).
- **Pre-selection** (fast spherical geometry from SDV):
  - wedge pairs: A within +/-20 deg of the doublet's azimuth seen from SDV, at 85-110 deg;
  - control candidates: A at least 25 deg off that azimuth, at 88-110 deg. Ranked as below, each event used
    at most once, top 40 kept.
- **Geometry:** TauP iasp91 SKS pierce points from event A (catalogue depth, clamped to 0-700 km) to SDV.
  Entry and exit are the first and last CMB pierce points. d_entry and d_exit are the arc distances (CMB
  radius 3,480 km) to the doublet's entry and exit points (from E2000a). Geometry is computed for at most
  400 wedge events, in rank order.
- **Classes:**
  - **near:** d_entry <= 100 km and d_exit <= 100 km (shares both regions; secondary).
  - **exit:** d_exit <= 100 km and d_entry > 100 km (**primary**).
  - **control:** a control candidate with d_entry > 300 km and d_exit > 300 km.
  - **other:** not used.
- **Separation gate:** |tS - tSKS| >= 10 s, and tSKKS - tSKS >= 6 s (or no SKKS).
- **Ranking within a class:** epicentral separation (smaller first), then the smaller magnitude of the
  pair (larger first). Each event is used at most once per class.
- **Verification:** `sdv23.membership`, unchanged except for the number of stations: mantle P, 0.8-2 Hz, at
  the 12 nearest IU/II stations at 30-90 deg. A pair is a member if >= 3 stations have cc >= 0.90, and it is
  relocated if >= 4 do; both are required.
- **Caps** (pairs checked / pairs verified), in this order:
  - exit: 25 / 12;
  - control: 15 / 8;
  - near: 10 / 6.
  Time budget 240 min.
- **Measurement at IU.SDV** (identical to the SDV branch of pre-registration 26):
  - `sdv23.process`: response removed, displacement, 0.4-2 Hz zero-phase, upsampled to 200 Hz.
  - Windows on record A: SKS_R and SKS_Z [-2, +3] s, S_T [-3, +10] s, P_Z [-1, +4] s around the
    largest peaks near the iasp91 times.
  - Record B is shifted by the relocation dt0, then fine lags are taken within +/-0.6 s.
  - x = SKS_R lag - reference lag - location-correction difference (reference: S_T if cc >= 0.90, else
    P_Z, else abs). xz is the same for SKS_Z.
- **Frames.** Every primary pair spans the 2009-09-18 metadata step, which turns SDV's horizontals by 3 deg
  in the station records. Pre-registration 23 measured the real relative rotation between 2000 and 2009 from
  long-period waves: only 0.55 deg. So the **primary frame is "nom"**: both records rotated with the
  pre-step azimuths, BHN 0 / BHE 90 deg. The "meta" frame (each record's own metadata, as in 26) is reported
  with the same grades.
- **Usable pair:** SKS_R cc >= 0.80 and not at the search edge, SKS_R SNR >= 5 in both records, x finite.

## Gates

- **G0:** the pipeline reproduces P23 on the doublet: the doublet's meta-frame x against S_T is within
  0.010 s of P23's x_meta (+0.0284 s, read at run time from results/sdv23-v2/sdv23_summary.json).
- **G1:** E2009a is a verified, relocated member (sdv23 rules, 30 stations).

## Pair grade

- **stable:** |x| <= 0.04 s.
- **changed:** |x| >= 0.08 s, and xz has the same sign with |xz| >= 0.04 s.
- **unclear:** any other usable pair.
- **unusable:** a pair that is not usable.

## Primary grade (exit class, nom frame)

- **untestable:** G0 or G1 fails, or no exit pair is usable.
- **change:** all three of these hold:
  - at least half of the usable exit pairs (and at least one) are *changed*;
  - all changed pairs have the same sign;
  - the controls are quiet: >= 3 usable controls, with <= 20% of them at |x| >= 0.08 s.
- **stable:** no usable exit pair is *changed*, and at least 2/3 of them are *stable*.
- **unclear:** anything else.

The same rule grades the **near** class (secondary) and the meta frame (reported). Also reported:
- the median x of all usable pairs (a station-side check);
- every pair's x, xz, references, SNR and correlation;
- the list of candidates with their geometry.

## Forecasts (before data)

- G0 passes: 80% (same pipeline as 26).
- At least one usable exit pair: 45%. Usable pairs need M >= 5.3 events with enough SKS signal at 95 deg
  and cc >= 0.90 mantle-P repeats 5+ years apart. Both are uncertain.
- Usable controls >= 3: 40%.
- If testable: stable 55%, unclear 35%, change 10%.
- At least one usable near pair: 25%.

## What each outcome would mean

- **stable:** SKS paths leaving the core within 100 km of the doublet's exit point did not change by more
  than 0.04 s over 5-25 years. That removes support for the receiver-side reading of the claimed ULVZ
  change. With 23 (not reproduced) and 26 (entry side), the public data would then show no decadal change
  in this region from any side.
- **change:** the first evidence, on independent earthquake pairs, of a decade-scale change under the
  Caribbean side of this path. That would be remarkable. Because the doublet itself showed +28 ms in P23,
  any claim would need the geometry explained, and the authors told first.
- **unclear / untestable:** reported as such, with no re-grading under other thresholds.

## Limits

- Exit-side sampling only. A change confined to the doublet's exact exit point, narrower than the pairs'
  spread, would escape.
- Pairs span different years; a change that happened only between 2000 and 2009 is covered fully only by
  pairs with A early and B after 2009.
- Catalogue locations and depths have errors; mantle-P membership and relocation correct the
  differential location only.
- Beyond about 105 deg, SPdKS/SKPdS can merge with SKS. Zhang & Wen's own group showed this for a different
  SKS doublet claim (ESS Open Archive, doi 10.22541/essoar.167276417.78056130). Pairs are graded as written
  and reported with their distance.
- One station; iasp91 ray geometry; the Fresnel zone is argued, not modelled.

Credit: Zhang & Wen 2025 set up the measurement and published their processed traces. Cui, Wang, Vidale
et al. 2026 (GRL e2025GL119930) tested CMB stability with 231 repeating-earthquake pairs.

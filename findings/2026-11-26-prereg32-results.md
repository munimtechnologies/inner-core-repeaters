# Pre-registration 32 results: Zhou's 1 s Kermadec SKS advance is reproduced in her band, but it is not in the direct SKS pulse

Graded 2026-11-26 from job `ks-sks32-v1` (ks32.py 0.1.0, finished 2026-11-25 23:21:49Z, 0 errors,
unit tests passed). Plan: `findings/2026-11-24-prereg32-kermadec-sks-plan.md` (public before any
waveform was downloaded; sha256 71f0bcff..., checked by the job at start). Text CC BY 4.0, code
Apache-2.0. Made with AI agents working in MV2's name (MV2, Munim, Inc.); no human expert has
checked it yet.

Files: `results/ks32-v1/ks32_summary.json`, `ks32_phases.csv`, `ks32_stations.csv`, `stdout.txt`
(job output, unchanged), and `ks32_grade.json`, `ks32_bands.csv`, `ks32_band_shift.csv` (values
copied from those files, rounded to 4 decimals; medians and differences computed from them).

## Grades (exactly as pre-registered)

Primary pair 1997-05-03 / 2018-09-10, four stations, all four testable in every band.
x = median over passing references of d(SKS) - d(ref), after the iasp91 location correction;
negative = SKS earlier in 2018.

| Station | Zhou 2022 (Table S2) | x long, 20-100 s (T1) | x mid, 5-20 s (reported) | x short, 2-10 s (T2) |
|---|---|---|---|---|
| COLA | -1.05 | -0.944 | -0.301 | -0.001 |
| ULN | -1.0 | -1.064 | -0.330 | +0.415 |
| INK | -0.8 | -0.907 | -0.373 | -0.482 |
| YAK | -1.1 | -0.753 | -0.563 | +0.356 |
| median | | -0.926 | -0.351 | +0.178 |

- **T1: reproduced.** k = 4 stations at or below -0.5 s (rule: k >= 3). Our values are within
  0.11 s of Zhou's at COLA, ULN and INK, and 0.35 s smaller at YAK.
- **T2: absent.** k = 0 stations at or below -0.5 s (rule: k <= 1). INK is the closest (-0.48 s).
- **Joint: window-content.** As the plan states, this means the long-window advance is not in the
  direct SKS pulse, which supports Zhang & Wen's reading over the outer-core reading.

## Reported, not graded

- **The advance shrinks with frequency at every station** (median -0.93, -0.35, +0.18 s from
  long to short band). Measured against Pdiff alone, which passes every gate in every band at
  all four stations, SKS - Pdiff is -0.42/-0.51/-0.53/-0.72 s (long), +0.09/-0.33/-0.11/-0.24 s
  (mid) and +0.27/+0.41/-0.23/+0.80 s (short) at COLA/ULN/INK/YAK.
- **The frequency dependence sits in the SKS window.** From long to short band, the SKS residual
  moves by +0.92, +0.93, +0.73 and +1.61 s; the Pdiff residual in the same records moves by only
  +0.24, +0.01, +0.42 and +0.09 s (`ks32_band_shift.csv`). A relative mislocation, an origin-time
  error or a clock error moves every band by the same amount, so none of them makes this
  difference. (Post-hoc description, not a pre-registered test.)
- **Uncorrected x_long** (catalogue location ignored): -1.11, -1.74, -1.35, -1.40 s.
- **SKKS minus references (long, ULN only):** +1.08 s. SKKS failed the SNR gate at the other
  three stations.
- **Other pairs:** all 20 station-pairs of 1997-1999, 1997-2006, 1999-2006, 1999-2018 and
  2006-2018 were untestable in the long band (14 SKS failed a gate, 6 had no passing reference),
  so the planned time series says nothing about when a change appeared.

## Forecasts scored

| Forecast (before data) | p | Came true? |
|---|---|---|
| T1 reproduced | 55% | yes |
| T2 absent, given T1 reproduced | 45% | yes |
| Joint window-content | 25% | yes |
| iasp91 correction moves x_long by < 0.2 s at every station | 70% | **no** (0.17, 0.68, 0.44, 0.65 s) |
| Per-reference x_long agree within 0.3 s at each station | 50% | **no** (spreads 0.76, 1.27, 1.59, 0.56 s) |
| Amendment 1: at least one station's per-reference spread > 0.5 s | 55% | yes (all four) |

## What this does and does not show

The ~1 s advance Zhou reported is in the public data and is reproduced by independent code in her
band and windows. It is not in the direct SKS pulse at 2-10 s period, and it is about a third as
large at 5-20 s. A change of P-wave speed along the SKS ray in the outer core would advance the
direct pulse too. So the long-period advance comes from something else in the long window.

The test does not say what that is. Candidates, none tested here:
(1) SPdKS/SKPdS, which separate from SKS at 100-106 degrees and change between the events (Zhang &
Wen, ESSOAr preprint 2023); (2) differences between the two earthquakes themselves (depth phases,
source time function, radiation), as Yang et al. (2026, Communications Earth & Environment) argue
for source repeatability; (3) a band-dependent response difference: at COLA, ULN and YAK the 1997
and 2018 records come from different sensors (blank location code vs 00), and Yang et al. report
co-located-sensor differences above 1 s at ULN and YAK in 2018. An outer-core anomaly lying off the
ray path, inside the broad long-period sensitivity zone only, is also not excluded; it would no
longer be a change along SKS's path.

## Limits

One pair, four stations; Zhou's other six pairs are not re-measured. USGS catalogue locations and
iasp91 are used for the correction, and the per-reference values disagree by 0.56-1.59 s, a sign
that the relative location is imperfect (SS, the most slowness-sensitive reference, differs most).
This does not affect the band comparison, which uses the same correction in every band. The
second pre-registered reading, F (amendment 1: a joint fit of an SKS-only term and the relative
source position), is filed as job ks32f-v1 and will be reported separately; it does not change
T1, T2 or the joint grade. The 0.1-0.5 Hz windows may still hold some SPdKS energy. COLA's sensor
azimuth changed 11 days before the 2018 event; each event is rotated with its own metadata.

## Prior work

Zhou 2022, Commun. Earth Environ. 3:97 (doi 10.1038/s43247-022-00432-7): the claim.
Zhang & Wen 2023, ESSOAr (doi 10.22541/essoar.167276417.78056130/v1): no SKS difference at
>= 0.2 Hz after their own relocation; SPdKS-SKPdS changes. Our T2 agrees with this, using
catalogue locations, open code and rules fixed in advance. Zhou's reply (doi
10.22541/essoar.167458062.24504379/v1). Yang, Li, Ringler, Jiang, Wu & Song 2026, Commun. Earth
Environ. (doi 10.1038/s43247-026-03961-7): the reported changes "can be largely explained by
artifacts" of source repeatability and station issues. Searches on 2026-11-26 (web, Crossref from
2026-01-01, Semantic Scholar citations of Zhou 2022: 6 papers) found no other re-measurement.
What is new here: an independent, pre-registered, open-code measurement that reproduces Zhou's
number in her own band and windows at all four stations and, in the same pipeline, shows that it
shrinks with frequency and is gone from the direct pulse, with the frequency dependence confined
to the SKS window.

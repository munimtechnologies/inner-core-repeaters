# Pre-registration 16, addendum G: how a P16 result will be read, given where the paths cross the inner-core boundary

By MV2, Munim, Inc. (made with AI agents). Written 2026-11-19, about 08:00 UTC. Text CC BY 4.0;
code (`geom.py`, `test_geom.py`) Apache-2.0. No human expert has reviewed it.

**Status when written.** Job `icr-nea-v2` (pre-registration 16, released code `nea.py` 0.1.0 at
tag v0.14.7) is queued and has produced no output. No near-antipodal waveform has been fetched or
seen. The geometry described here has not been computed yet either: `geom.py` runs as its own
job (`icr-geom-v1`) after this file is public. Nothing in pre-registration 16 changes: its
stations, processing, gates, calibration, clock test, source control and headline stay exactly as
published. This addendum only fixes, in advance, what each headline will be **taken to mean**.

## Why an addendum is needed

Pre-registration 16 says that near-antipodal PKIKP "crosses [the inner-core] boundary far from the
Alaska paths' entry and exit points". We had not computed that. It matters for the reading:

* A rigid turn-back of the whole inner core (Wang, Vidale et al. 2024, *Nature* 631:340) makes the
  "equivalent date" clock the same on every path.
* A change of the inner-core **surface** (Wen 2006, *Science*; Zhang & Wen 2026, *JGR Solid
  Earth*, doi:10.1029/2026JB034595, whose key points say the surface change "varies greatly in
  geographical location") gives the clock only where a path crosses the changing patch.

So a backtracking clock on the near-antipodal paths speaks against a single local surface patch
only if those paths cross the inner-core boundary (ICB) **away from** where the ILAR and YKA paths
cross it. Every PKIKP path crosses the ICB twice: once on the source side (entry) and once on the
receiver side (exit). Near-antipodal paths enter close to the source (both they and the Alaska
paths start at the same South Sandwich Islands source), so a shared source-side patch is possible
and must be checked, not assumed.

## What `geom.py` computes (ObsPy TauP, ak135)

From the multiplet centroid (-59.0365, -25.6939), at the median catalogue depth of the ILAR
multiplet events (results/m0324-v1 at commit d0186eec), for the ILAR and YKA array centres and for
every candidate station of nea.py (BHZ of IU, II, IC, G, GE, GT, CU at 155-179.5 degrees; nea.py's
cover and ranking rule reproduced to mark the stations it is expected to choose):

* the source-side ICB **entry** point and receiver-side **exit** point (lat, lon), the bottoming
  point and its depth below the ICB;
* the inner-core path length, the share of it within 300 km of the ICB (`upper_ic_frac`), and its
  length-weighted mean distance from the spin axis (how far rigid rotation moves the material it
  samples);
* for each candidate, the distance **on the ICB** (km, at the ak135 inner-core radius) from its
  entry and from its exit to the nearest of the four ILAR/YKA crossing points;
* a path class: **separate** if both distances are >= 250 km; otherwise **shared-entry**,
  **shared-exit** or **shared-both**. 150 km and 400 km are reported as sensitivity only.

250 km is a judgement: roughly two to three times our rough estimate of the first Fresnel-zone
radius of 1 Hz PKIKP at the ICB (about 100 km). It is not calibrated.

## Pre-registered readings

Weights are the numbers of ok pairs per station in `nea_pairs.csv` of `icr-nea-v2`; `f_sep` is the
weighted share of stations whose path class is **separate** (`geom.py --grade` computes it).

| P16 headline | f_sep >= 0.75 | f_sep < 0.75 (or no weights) |
|---|---|---|
| B-supported | **R1**: the backtracking clock on ICB crossings separate from ILAR/YKA's; favours a motion of the whole inner core over one local surface patch | **R2**: compatible with both pictures; no claim beyond "the clock is seen on a second set of stations" |
| M-supported | **R3-separate**: on separate crossings the time signal follows plain time order; the 2024 clock does not describe these parts of the inner core | **R3-shared**: plain time order on crossings shared with ILAR/YKA; in conflict with the backtracking clock seen there |
| no-time-signal | **R4**: read with `upper_ic_frac` against ILAR's (a much smaller upper-inner-core share supports "change confined near the surface") | **R4** |
| untestable, inconclusive, source-confounded | no reading beyond the headline | same |

Even R1 is not proof of rotation: a surface change that is coherent over thousands of kilometres
would also give it. R1 only rules out the simplest local-patch version.

## Forecasts (before any geometry or P16 output)

| Item | Forecast |
|---|---|
| Every expected P16 station has its exit >= 250 km from every ILAR/YKA crossing | 90% |
| Every expected P16 station has its entry >= 250 km from ILAR's and YKA's entries | 70% |
| Unweighted share of expected stations classed **separate** >= 0.75 | 70% |
| Median `upper_ic_frac` of expected stations below ILAR's | 85% |

Reasoning: near 180 degrees PKIKP is nearly vertical in the mantle and outer core, so it enters
the inner core within roughly 5-15 degrees of the source and dives towards the centre. At
ILAR's ~150 degrees the ray stays in the upper few hundred kilometres and enters much further from
the source. These are rough numbers we have not computed; the job will.

## Files

* `geom.py` 0.1.0, `test_geom.py` (unit tests without network).
* Results will be published under `results/geom-v1/` (`geom_paths.csv`, `geom_summary.json`,
  `plots/geom_icb_crossings.png`), and the reading of P16 under `results/nea-v2/geom_reading.json`.

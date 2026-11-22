# Pre-registration 24 results (job sdv-zframe24-v1): graded "other" / "untestable" because of a set-pairing defect in our code

*inner-core-repeaters, MV2 of Munim, Inc. Graded 2026-11-22 from job `sdv-zframe24-v1` (sdvz24.py
0.1.0 with sdv23.py 0.1.1, release v0.21.0; started 2026-11-22 00:20 UTC on the Dell, finished
00:31 UTC, all steps exit 0, 20 unit tests passed). Plan:
findings/2026-11-21-prereg24-sdv-zenodo-frame-plan.md (public before the job ran:
`plan_public_at_start` true, sha256 dbca3571...). Results: results/sdvz24-v1/. Licence: code
Apache-2.0, text CC BY 4.0. Made with AI agents working in MV2's name; no human expert has
checked it yet.*

## Short version

The job downloaded the authors' published IU.SDV traces (Zhang & Wen 2025, Zenodo
10.5281/zenodo.16421554; MD5 matched, gate G0) and ran without errors. But the graded pair
compares **unlike files**: the authors' *plain* 2000 traces (`2000.R.sac`, `2000.T.sac`) against
their *band-passed* 2009 traces (`2009.R.sac.bp`, `2009.T.sac.bp`). The literal pre-registered
grades are below. They describe that mismatch, not the authors' processing or the Earth. A
corrected run that pairs like with like (amendment 1) is queued as `sdv-zframe24-v2`.

## What the zip contains (results/sdvz24-v1/sdvz24_manifest.csv)

| Event | Plain files | Band-passed files |
|---|---|---|
| 2000-07-07 | 2000.R/T/Z.sac, 100 Hz | 2000.R/T/Z.sac.bp, **1000 Hz** |
| 2009-12-17 | 2009.R/T/Z.sac, 100 Hz | 2009.R/T.sac.bp **100 Hz**, 2009.Z.sac.bp 1000 Hz |

All R/T files carry SAC headers baz 323.63, cmpaz(R) 143.63, cmpaz(T) 233.63 (the standard
great-circle rotation), and identical phase markers (SKSac 1430.84 s, S 1470.40 s).

## The defect

The plan groups traces into sets by "same file name apart from the component code". For file
names without a channel code, sdvz24.set_key 0.1.0 instead fell back to folder|rate|start. For
2000 the plain (100 Hz) and band-passed (1000 Hz) files differ in rate and stayed separate. For
2009 both are 100 Hz, so they fell into one set, and the `.bp` R and T (read after the plain
files) overwrote the plain ones. The primary set of each event was the first sorted key:
plain for 2000, the merged set (in effect `.bp`) for 2009. The unit tests did not cover file
names like these.

## Literal grades (pre-registered rules)

| Part | Value | Grade |
|---|---|---|
| Frame gates G0-G3 | MD5 ok; R/T sets for both; raw location 00 for both; fit c 0.99997 (2000) and 0.978 (2009), both "proper" | pass |
| Frame: rho = A(2009) - A(2000) | 178.02 - 0.17 = **177.85 deg** | **other** |
| Combined, with pre-registration 23's delta (+0.55 deg, its G5 passed) | eps = 177.30 deg | misrotated (by the rule) |
| Combined, as coded (Rayleigh-Z delta; c 0.78 / 0.66 < 0.8) | - | unresolved |
| Measurement gate G4 | S_T cc 0.795 (< 0.90), SKS_R cc 0.786 (< 0.80) | fail |
| Measurement | x0 = -13.6 ms (not usable) | **untestable** |

The code applied the combined rule with the Rayleigh-Z delta only. The plan says to use
pre-registration 23's long-period delta when its gate G5 passed, which it did. Both readings
are listed; neither is meaningful here.

## What the frame fits do show (results/sdvz24-v1/sdvz24_fits.csv; reported, not graded)

| Authors' set | Assumed channel-1 azimuth A | Lag vs raw | c |
|---|---|---|---|
| 2000 plain | 0.17 deg | 0.00 s | 0.99997 |
| 2000 band-passed | 180.24 deg | -0.64 s | 0.979 |
| 2009 band-passed (merged set) | 178.02 deg | -0.62 s | 0.978 |

1. **The authors' plain 2000 traces are exactly our raw 2000 record** rotated with nominal
   north/east (A = 0.17 deg, c = 0.99997, no time shift).
2. **Both band-passed sets look like our raw records reversed in polarity and shifted by about
   0.63 s** (A near 180 deg). The fit compares 0.4-2 Hz zero-phase versions of both, so this is
   consistent with the authors' band-pass being causal (one-pass) rather than zero-phase. A
   causal filter delays and reshapes every pulse. It cannot by itself make one phase earlier
   than another in only one event, unless the two events' spectra differ (the paper reports a
   broader 2009 SKS wavelet). The corrected run measures this directly.
3. Comparing the two band-passed sets with each other (post hoc, from the same table): rho =
   178.02 - 180.24 = **-2.22 deg**. That is within 1 deg of the 3-deg metadata step (frame
   grade "metadata-frame" by the plan's rule), so each event was probably rotated with its own
   station metadata. With pre-registration 23's measured relative rotation of +0.55 deg, the
   2009 band-passed R/T would be misrotated by about 2.8 deg. Pre-registration 23 found that a
   3-deg relative rotation moves the SKS-minus-S time by about 19 ms. This is a post-hoc
   reading of a secondary set, **not a grade**.

## Forecasts

| # | Forecast | P | Outcome |
|---|---|---|---|
| H1 | Usable R/T sets for both events | 85% | yes |
| H2 | Fits with c >= 0.95 for both | 70% | yes |
| H3 | metadata / same / other-or-untestable | 50 / 35 / 15 | other (code defect) |
| H4 | consistent / misrotated / unclear-unresolved | 40 / 35 / 25 | not meaningful (code defect) |
| H5 | reproduced, given testable | 65% | untestable |

## What happens next

Amendment 1 (findings/2026-11-22-prereg24-amendment-1.md) fixes the grouping so it does what
the plan says. It changes no threshold, window or rule. It re-runs the pinned sdvz24.py
0.1.0 on the plain pair (primary, graded) and on the band-passed pair (secondary). It also
measures how the authors' band-passed traces are delayed relative to their plain traces in each
event and window. Job `sdv-zframe24-v2` (sdvz24b.py 0.1.0) reuses the same downloaded data.

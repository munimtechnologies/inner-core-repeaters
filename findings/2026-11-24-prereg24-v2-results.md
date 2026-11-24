# Pre-registration 24 results (job sdv-zframe24-v2): the authors' own IU.SDV traces do not show the 160 ms SKS advance

*inner-core-repeaters, MV2 of Munim, Inc. Graded 2026-11-24 from job `sdv-zframe24-v2`
(sdvz24b.py 0.1.0 running the pinned sdvz24.py 0.1.0 and sdv23.py 0.1.1, code commit dbbd91d;
started 2026-11-23 21:02 UTC on the Dell, finished 21:10 UTC, all steps exit 0, 10 unit
tests passed). Plan: findings/2026-11-21-prereg24-sdv-zenodo-frame-plan.md (public at start,
sha256 dbca3571...). Amendment 1: findings/2026-11-22-prereg24-amendment-1.md (public at
start, sha256 03ad21a3...). Results: results/sdvz24-v2/. Licence: code Apache-2.0, text
CC BY 4.0. Made with AI agents working in MV2's name; no human expert has checked it yet.*

## Short version

Zhang & Wen (2025, *Nature Communications* 16:8754) report that radial SKS at IU.SDV arrived
**160 ms earlier** for the 2009 Aleutian earthquake than for the 2000 one, and read it as a
decade-scale change of an ultra-low velocity zone. They published their processed traces
(Zenodo 10.5281/zenodo.16421554). We measured those traces with the method fixed in advance.

- **Frame: metadata-frame.** Their plain traces are the raw IU.SDV records rotated with each
  event's own station metadata: channel-1 azimuth 0.17 deg for 2000 and 357.14 deg for 2009
  (metadata 0 and 357 deg), both with fit correlation 0.99997 and no time shift. rho = -3.03 deg.
- **Combined: misrotated.** Our test 23 measured that the sensor physically turned by only
  +0.55 deg between the events. So their 2009 radial/transverse pair is rotated 3.58 deg away
  from their 2000 pair (eps = -3.58 deg; rule: misrotated if |eps| >= 2 deg).
- **Measurement: not-reproduced.** On their own traces, radial SKS minus transverse S changed
  by **+20.7 ms** (later in 2009, not earlier). The rule needed <= -100 ms. All gates passed:
  SKS correlation 0.987, transverse S correlation 0.991, none at the search edge.

The band-passed traces they also published give the same grades: metadata-frame (rho -2.23
deg), misrotated (eps -2.78 deg), not-reproduced (+27.4 ms).

## Grades (pre-registered rules, amendment 1)

| Part | Primary (plain .sac, graded) | Secondary (.sac.bp, reported) |
|---|---|---|
| Gate A1 (like with like) | pass | - |
| Frame gates G0-G3 | pass (MD5, sets, raw 00, c 0.99997 / 0.99997, same convention) | pass (c 0.979 / 0.978) |
| rho = A(2009) - A(2000) | -3.03 deg | -2.23 deg |
| Frame grade | **metadata-frame** | metadata-frame |
| eps = rho - 0.55 deg | -3.58 deg | -2.78 deg |
| Combined grade | **misrotated** | misrotated |
| Measurement gates G1, G4 | pass | pass |
| x0 = lag(SKS_R) - lag(S_T) | **+0.0207 s** | +0.0274 s |
| Measurement grade | **not-reproduced** | not-reproduced |

## Every window on their traces (primary, relative to transverse S)

From results/sdvz24-v2/sdvz24b_summary.json (`lags0`; each lag minus the S_T lag, 2009 minus
2000). The last column is our test 23 on the raw public records (metadata rotation,
results/sdv23-v2/sdv23_pairs.csv).

| Window | Their traces (P24-v2) | Public records (P23) |
|---|---|---|
| radial SKS | +20.7 ms | +28.4 ms |
| radial after SKS | +31.9 ms | +39.0 ms |
| radial S-Scd-ScS start | -38.3 ms | -38.6 ms |
| radial S-Scd-ScS middle | +15.3 ms | +15.0 ms |

Their traces give the same numbers as the raw records to within 8 ms. Nothing in their
processing adds a 160 ms SKS advance that the raw records lack. As in test 23, only the start
of the radial S-Scd-ScS wave moves early (-38 ms; the paper reports -80 ms there). That window
was not graded.

## Rotation sensitivity (results/sdvz24-v2/plain/sdvz24_sens.csv)

Rotating their 2009 traces by da from -20 to +20 deg moves x from -0.129 s to +0.129 s,
6.1 ms per degree near zero. Between da = -3.5 and +3.5 deg, roughly the size of the
misrotation, x stays between -0.001 and +0.042 s. Reaching -0.100 s takes da = -17 deg.
Getting -0.160 s would take more than 20 deg, beyond the grid. Their 3.6-deg misrotation
therefore moves this measurement by about 20 ms. It does not create or remove a 160 ms change.

## Their band-pass

Compared with their plain traces, their band-passed traces are reversed in polarity. They are
delayed by 0.59-0.65 s on R and T and 0.24-0.34 s on Z, identically in both events. The
2009-minus-2000 difference of that delay is at most 7.8 ms in any window: R SKS -1.95 ms, T
SKS +0.41 ms, Z SKS +3.75 ms. So the filter is not zero-phase, but it treats both events the
same and cannot make SKS earlier in one of them.

## Forecasts

| Forecast | P | Outcome |
|---|---|---|
| A1 passes | 90% | yes |
| Primary frame metadata / same / other | 55 / 30 / 15 | metadata-frame |
| Combined misrotated, given metadata-frame | 80% | yes |
| Primary measurement not-reproduced / reproduced / untestable | 80 / 10 / 10 | not-reproduced |
| Secondary measurement not-reproduced | 65% | yes |
| Band-pass delay difference in SKS window >= 30 ms | 15% | no (2 ms) |
| Plan H5: reproduced, given testable | 65% | no |
| Plan H6: abs(dx/dda) <= 5 ms/deg | 75% | no (6.1; band-passed 5.0) |
| Plan H7: S_R_early lag changes >= 30 ms between da = -3 and +3 | 30% | no (3.6 ms) |

## What it means

Every line of the test now points the same way. The public records (test 23, +28 ms) and the
authors' own published traces (this test, +21 ms) give no 160 ms SKS advance when SKS is
measured against the transverse S wave that the paper reports unchanged. The authors rotated
each event with its own metadata, which leaves the 2009 traces 3.6 deg off. That matters for
the S-Scd-ScS part of the claim (the radial S wave is sensitive to mixing in transverse
energy) but only by about 20 ms for SKS. The paper does not state the reference against which
its 160 ms was measured. If it was an absolute alignment, for example on P or on a catalogue
origin, the difference from ours would come from that reference, not from the traces. Tests
25, 26 and 27 (an onset reading, entry-side neighbour stations, and exit-side pairs) were all
untestable with public data, so no other route confirms or rules out a change.

## Limits

One station, one earthquake pair. This test grades a transverse-referenced differential
measurement on the authors' traces; it does not test the Earth directly. The authors' exact
measurement procedure for the 160 ms is not given in the paper. Made with AI agents working
in MV2's name; no human expert has checked it yet.

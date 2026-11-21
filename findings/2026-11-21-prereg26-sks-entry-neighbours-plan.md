# Pre-registration 26: do other SKS paths through the same core-mantle-boundary entry region change between 2000 and 2009?

Written 2026-11-21 ~22:00 UTC, before any waveform of this test is downloaded (station lists only were
looked at, to know which stations exist). Code: `sksnbr26.py` 0.1.0 (Apache-2.0); the job waits until this
file is public on `main` and records its sha256. Made with AI agents working in MV2's name; no human
expert has checked it yet.

## Claim under test

Zhang & Wen 2025, *Nature Communications* 16:8754 (doi 10.1038/s41467-025-63814-3) report a decade-scale
change of an ultra-low velocity zone (ULVZ). It rests on one Aleutian doublet (2000-07-07 Mw 5.9 and
2009-12-17 Mw 5.8, their AI-0009) recorded at one station, IU.SDV (Venezuela, 94.9 deg). There the radial SKS
arrives 160 ms earlier in 2009 with a broader wavelet. Their models put the changing ULVZ near the
**SKS entry point** on the core-mantle boundary (CMB), on the source side. They also state that a change near
the SKS exit point (receiver side) "could explain the seismic data equally well". The paper does not
report SKS from any other station for this doublet.

Our pre-registration 23 (results/sdv23-v2) re-measured SDV from the public records and found +28 ms, not
-160 ms (graded *not-reproduced*). Pre-registrations 24 (authors' own traces) and 25 (onset and pulse-width
reading) are queued. This test does not depend on their outcome.

## Idea

SKS legs from the same source to stations at nearly the same azimuth enter the core at nearly the same
place. Seen from the Aleutian source, VE.FUNV (Caracas) is about 4 deg from SDV in azimuth and 1 deg
further away, and IU.PTGA (Brazil) is about 3 deg off in azimuth at 109 deg distance. By rough spherical
estimate, their SKS entry points lie within tens of km of SDV's. At 0.4-2 Hz that is inside the first Fresnel
zone of the CMB conversion (about 100 km half-width). Their exit points, near the receivers, are hundreds of
km from SDV's. So:

- A source-side ULVZ change large enough to move SDV's SKS by 0.1-0.16 s should also show on these
  *entry neighbours*.
- A receiver-side change would not.
- Nothing has changed if neither the neighbours nor SDV show it (as P23 found at SDV).

The paper leaves the side unresolved, and these neighbours can resolve it. If a neighbour changed, it
would be the first evidence for the change on a second, independent path.

## Data and processing (fixed)

- **Events:** E2000a and E2009a as in sdv23.py 0.1.1 (pinned commit dbbd91d, fetched at run time and its
  sha256 recorded).
- **Relocation and origin-time shift (dt0):** mantle P at IU/II stations at 30-90 deg, computed with
  `sdv23.membership`, unchanged.
- **Candidate stations:** every station that EarthScope's FDSN station service lists with a BHZ channel at
  both event times, at 85-125 deg from E2000a.
- **Geometry:** SKS pierce points from TauP, iasp91 model, E2000a source depth. Entry and exit points are the
  first and last pierce points at the CMB depth, placed on the source-station great circle (spherical).
  - d_entry: arc distance between a station's SKS entry point and SDV's, at the CMB radius of 3480 km.
  - d_exit: the same for the exit points.
- **Classes:**
  - **entry-neighbour:** d_entry <= 100 km.
  - **exit-neighbour:** not an entry-neighbour, and d_exit <= 100 km.
  - **control:** d_entry > 300 km and d_exit > 300 km, permanent global networks IU, II, GE, G, GT, MN only,
    at most 40 (sorted by network and station code).
  - **other:** the rest. Others within 300 km of SDV's entry or exit point are measured and reported but not
    graded.
- **Separation gate:** a station is measured only if TauP gives |tS - tSKS| >= 10 s and
  tSKKS - tSKS >= 6 s (or no SKKS).
- **Waveforms:** downloaded from EarthScope. Location preference 00, then blank, then 10, then others; the
  location must have Z and a horizontal pair (N/E or 1/2). IU.SDV uses sdv23's own download windows, so its
  numbers can be compared with P23.
- **Processing:** `sdv23.process`, unchanged: response removed, displacement, 0.4-2 Hz zero-phase, upsampled
  to 200 Hz. Rotation to R/T uses each record's own station metadata.
- **Windows,** set on the 2000 record around the iasp91 prediction (peak search +/-5 s; S search -5/+8 s):
  - SKS_R: radial SKS, -2/+3 s around the peak;
  - SKS_Z: vertical SKS, -2/+3 s;
  - S_T: transverse S or Sdiff, -3/+10 s;
  - P_Z: vertical P or Pdiff, -1/+4 s.
- **Lags:** the 2009 record is first shifted by round(dt0 x 200 Hz) samples. The fine lag is the normalised
  cross-correlation peak within +/-0.6 s, refined with a parabola (as sdv23.window_lag).
- **x (per station, s; negative = SKS earlier in 2009):** the SKS_R lag minus the reference lag, minus the
  location-correction difference loc_pred(dx, az, p_SKS - p_ref).
  - The reference is S_T if its cc >= 0.90 and it is not at the edge; otherwise P_Z on the same terms;
    otherwise "abs" (lag minus dt0 minus loc_pred(dx, az, p_SKS)).
  - xz is the same quantity for SKS_Z, with the same reference.
- **Usable station:**
  - SKS_R cc >= 0.80 and not at the edge;
  - SKS_R SNR >= 5 in both records (peak over RMS of the window P-120 s to P-20 s);
  - x finite.

## Gates

- **G0:** our pipeline reproduces P23 at SDV, i.e. SDV's x against S_T is within 0.010 s of P23's x_meta
  (+0.0284 s, read from results/sdv23-v2/sdv23_summary.json at run time).
- **G1:** E2009a is a verified member and relocated (sdv23 rules).

## Station grade

- **stable:** |x| <= 0.04 s.
- **changed:** |x| >= 0.08 s, and xz has the same sign with |xz| >= 0.04 s.
- **unclear:** anything else that is usable.
- **unusable:** not usable.

## Primary grade (source-side test)

- **untestable:** G0 or G1 fails, or no entry-neighbour is usable.
- **shared-change:** at least one usable entry-neighbour is *changed*, and the controls are quiet. Quiet
  means at least 5 usable controls, with at most 20% of them at |x| >= 0.08 s.
- **no-shared-change:** every usable entry-neighbour is *stable*.
- **unclear:** otherwise.

The secondary grade applies the same rule to exit-neighbours. The control false-change rate and every
station's x are reported in full.

## Forecasts (before data)

- **Testability:** G0 passes 80%. At least one usable entry-neighbour 50% (FUNV's 2000 data and both
  stations' SNR at Mw 5.8 are the risks).
- **If testable:** no-shared-change 60%, unclear 30%, shared-change 10%.

## What each outcome would mean

- **no-shared-change:** paths that enter the core within 100 km of SDV's SKS entry point did not change by
  >= 0.04 s over the same decade. The paper's preferred source-side ULVZ change is not supported. Only a change
  confined to a strip narrower than the neighbour separation would escape, and so would a receiver-side
  change. Together with P23, the public data show no decadal change on this path.
- **shared-change:** an independent path through the same CMB region changed. This would be the first
  multi-path evidence of a decade-scale ULVZ change, which is remarkable. SDV itself showed +28 ms in P23, so
  the geometry of such a change would need explaining before any claim.
- **unclear or untestable:** reported as such. No re-grading with other thresholds.

## Limits

- Only one doublet.
- TauP iasp91 ray geometry. Fresnel-zone overlap is argued, not modelled, and we do not model the ULVZ.
- A neighbour's S or P reference can carry its own path changes; the abs reference depends on station
  clocks.
- Metadata rotation only (P23 found SDV's real relative rotation to be 0.55 deg, not 3).
- Beyond about 105 deg (e.g. IU.PTGA at about 109 deg), SPdKS/SKPdS energy can merge with SKS. Wen's
  group showed this for a different SKS doublet claim (ESS Open Archive, doi 10.22541/essoar.167276417.78056130).
  A change at such a station is still graded as written, but it is reported with this caveat.

Credit: Zhang & Wen 2025 set up the measurement and published their processed traces. Cui, Wang, Vidale et
al. 2026 (GRL e2025GL119930) tested CMB stability with 231 repeating-earthquake pairs.

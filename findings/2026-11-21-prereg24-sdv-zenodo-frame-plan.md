# Pre-registration 24: which sensor orientation do the published IU.SDV waveforms behind the reported decadal ULVZ change assume?

Written 2026-11-21, **before the authors' data file or any waveform for this test was
downloaded**. Code: `sdvz24.py` 0.1.0 (uses `sdv23.py` 0.1.1 from the same commit). Job:
`sdv-zframe24-v1`. Licence: text CC BY 4.0, code Apache-2.0. Made with AI agents working in
MV2's name (MV2, Munim, Inc.). No human expert has checked it yet.

## Background

Zhang & Wen 2025 (*Nature Communications* 16:8754, doi 10.1038/s41467-025-63814-3) report
the only published decadal change of an ultra-low velocity zone at the core-mantle
boundary. It rests on one Aleutian doublet (AI-0009: 2000-07-07 Mw 5.9 and 2009-12-17 Mw 5.8)
recorded at one station, IU.SDV in Venezuela (94.9 deg). On the **radial** component SKS
arrives 160 ms earlier in 2009 and S-Scd-ScS shifts 80 -> 30 -> 0 ms. The **transverse**
component is unchanged. The paper removes the instrument response and rotates to the
great-circle path, but it does not discuss sensor orientation.

Pre-registration 23 (this repository) found that IU.SDV's metadata change the azimuths of
location 00 (STS-1) on 2009-09-18, between the two events: BHN 0 -> 357 deg, BHE 90 -> 87
deg. Location 10 changes 0 -> 343 deg. Pre-registration 23 measures whether the sensor
physically rotated. It re-measures the SKS change under three rotations, but it cannot know
which rotation the authors used.

The authors have published their processed data: Zhang & Wen, "Preprocessed Doublet Data of
AI-0009", Zenodo doi 10.5281/zenodo.16421554 (CC BY 4.0, 2025-07-25). The file
`data_IU.SDV.zip` (26,812,269 bytes, MD5 9ba0853da7a537a7a8f1117407675112) holds "the IU.SDV
station subset plus filtered three-component data" rotated "to great-circle paths". **So the
rotation they applied can be measured directly**: compare their radial and transverse traces
with the raw IU.SDV records.

## Questions

1. **Frame (primary).** For each event, which channel-1 azimuth did their rotation assume?
   Did they rotate 2009 in the same frame as 2000 (as with nominal N/E, or with the old
   metadata for both), or in the new metadata frame (3 deg apart)?
2. **Combined.** Given the physical rotation of the sensor between the events, is their 2009
   radial/transverse pair consistent with their 2000 pair, or misrotated, and by how much?
3. **Measurement on their own traces.** Does the SKS-minus-transverse-S differential time
   on their published traces show the 160 ms advance? How does it change if the 2009 traces
   are rotated by a few degrees?

## Method (all constants in `sdvz24.PARAMS`)

**Data.** Download `data_IU.SDV.zip` from Zenodo and check the MD5 (gate G0). Extract it
safely (members with `..` or absolute paths are skipped; nested zips one level). Read
every file's headers into a manifest (all SAC header fields that bear on orientation and
timing: cmpaz, cmpinc, baz, az, gcarc, o, b, e, khole, kcmpnm, t0-t3, user0-2). Assign
traces to an event when their start lies within -1 h .. +2 h of its USGS origin. Group them
into sets: same file name apart from the component code, or else same folder, rate and
start. The primary set of an event is the first set (sorted by name) that has R and T
covering SKS - 40 s to S + 60 s. Without one, it is the first set with R and T.

**Raw.** IU.SDV BH? locations 00 and 10 from the EarthScope FDSN service. Response removed
to displacement, pre-filter 0.1/0.2/7/9 Hz, 0.4-2 Hz zero-phase band-pass (2 corners),
rotated to geographic N/E with **each event's own station metadata**. The authors' traces get
the same 0.4-2 Hz band-pass.

**F. Frame fit.** On a 50 Hz grid, window SKS - 30 s .. S + 60 s after origin (TauP iasp91,
master location, 31 km), or its part that their traces cover (at least 40 s). For each lag
within +/- 5 s and each convention ("proper", or "reflected" = their T negated), the rotation
a_fit of our metadata N/E that best matches their (R, T) is solved in closed form (maximising
R.R' + T.T'). We keep the lag/convention with the highest normalised vector correlation c.
The channel-1 azimuth their processing assumed is **A = metadata azimuth + a_fit**. Done for
location 00 (graded) and 10 (reported).
- rho = A(2009) - A(2000), wrapped to +/-180 deg, primary sets, location 00.

**RZ (for the combined grade).** Rayleigh-Z azimuth of each event's location-00 sensor
(0.02-0.05 Hz, group velocity 4.2-3.0 km/s, search +/-45 deg around the metadata), as in
pre-registration 23. delta_RZ = a(2009) - a(2000) if both |c| >= 0.8.

**M. Measurement on their traces.** Their primary-set R and T (T negated if the fit says
"reflected") on a 200 Hz grid, SKS - 40 s .. S + 60 s after each event's own origin. The
windows come from their 2000 traces, as in pre-registration 23: SKS_R [-2, +3] s around the
largest radial peak within +/-5 s of SKS; S_T [-3, +10] s around the transverse peak in
S - 5 .. S + 8 s; S_R_early [-2, +1.5] s and S_R_mid [+1.5, +5] s around the radial peak there;
postSKS_R [+5, +11] s. 2009 is aligned to 2000 on S_T within +/-5 s (integer samples), then
each lag is measured within +/-0.6 s (normalised cross-correlation, parabolic refinement).
**x0 = lag(SKS_R) - lag(S_T)** (negative = SKS earlier in 2009 relative to S). Sensitivity:
the 2009 R/T rotated by da = -20 .. +20 deg in 0.5 deg steps (same alignment shift), x(da) and
every window lag recorded.

## Decision rules (pre-registered)

**Frame grade** (`classify_frame`). Gates, all required, else **untestable**: G0 MD5 matches;
G1 both events have an R/T set; G2 location-00 raw data fit for both; G3 c >= 0.95 for
both, with the same convention.
- **same-frame**: |rho| <= 1.0 deg (their 2009 rotation used the same frame as 2000, as with
  nominal N/E or a single set of azimuths);
- **metadata-frame**: |rho - (-3)| <= 1.0 deg (each event rotated with its own metadata);
- **other**: anything else.

**Combined grade** (`grade_combined`), eps = rho - delta, wrapped:
**consistent** if |eps| <= 1.0 deg, **misrotated** if |eps| >= 2.0 deg, else **unclear**.
delta is pre-registration 23's long-period delta if its gate G5 passes when we grade. If not,
it is this job's delta_RZ if both |c| >= 0.8. Otherwise the grade is **unresolved**.

**Measurement grade** (`grade_reproduce`). Gates: G1 both events cover every window plus
5.6 s of search margin; G4 S_T cc >= 0.90 and SKS_R cc >= 0.80, neither at the search
edge. **reproduced** if x0 <= -0.10 s, else **not-reproduced**. This checks whether a
transverse-referenced differential measurement on the authors' own traces shows the
claimed advance. It is not a test of the Earth.

Reported, not graded: a_fit, A, c, lag and convention for every set and both locations;
the SAC baz/cmpaz headers; x(da) and the slope dx/dda at 0; x at da = -17, -3, 0, +3, +17; the
da that minimises |x|; the change of the S_R_early lag between da = -3 and +3; T/R amplitude
ratios in the SKS and S windows.

## Forecasts (before data)

- H1 both events have usable R/T sets covering SKS to S: 85%.
- H2 frame fits with c >= 0.95 for both events: 70%.
- H3 frame grade: metadata-frame 50%, same-frame 35%, other or untestable 15%.
- H4 combined grade (when graded): consistent 40%, misrotated 35%, unclear/unresolved 25%.
- H5 measurement reproduced (x0 <= -0.10 s), given testable: 65%.
- H6 |dx/dda| at 0 <= 5 ms per degree (SKS is nearly insensitive to small rotations): 75%.
- H7 S_R_early lag changes by >= 30 ms between da = -3 and +3 deg: 30%.

Prior physics (as in pre-registration 23): a 3 deg error mixes 5% of the transverse into the
radial. That cannot by itself move a clean radial SKS pulse by 160 ms unless the transverse
energy in the SKS window is ~20 times the radial. It can distort the radial S-Scd-ScS, where
SH is strong. A **misrotated** grade would therefore matter most for the S-Scd-ScS part of the
claim and for the D'' anisotropy reading. A **consistent** grade would remove orientation as
an explanation of the published traces.

## Limits known in advance

One station, one doublet. The band-pass is applied a second time to traces that are already
filtered. This changes waveform shape slightly but not the rotation angle. The fit absorbs
any difference between their back-azimuth and ours into a_fit. That difference is common to
both events (0.6 km apart), so it cancels in rho. Location corrections (about 16 ms at most
for SKS - S, from the paper's 0.6 km relocation) are not applied in M. If the zip lacks R/T
for one event, everything is untestable and we say so.

## Prior work

Zhang & Wen 2025 (above) and their Zenodo data (above); Cui, Wang, Vidale et al. 2026 GRL 53
e2025GL119930 (CMB stable over decades in 231 core-reflection pairs; their one apparent
change, under Siberia, was a local earthquake); Ringler et al. 2013 SRL 84:926 (GSN sensor
orientation errors). Pre-registration 23 and amendment 1 (this repository). We found no
comment, reply or replication of the ULVZ claim (searched 2026-11-21).

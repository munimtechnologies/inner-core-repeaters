# Pre-registration 32, amendment 1: the source-fingerprint reading (F)

Written 2026-11-24, **before any waveform of pre-registration 32 was downloaded**. At the time of
writing, job `ks-sks32-v1` was still queued for compute and had not started. The graded T1, T2 and
joint headline of the plan (`findings/2026-11-24-prereg32-kermadec-sks-plan.md`, code `ks32.py`
0.1.0) is **unchanged**. This amendment adds one separately labelled reading, **F**, computed by
`ks32f.py` 0.1.0 from the job's own output file `ks32_phases.csv`. Text CC BY 4.0, code
Apache-2.0. Made with AI agents working in MV2's name (MV2, Munim, Inc.); no human expert has
checked it yet.

## Why

The plan corrects each phase with iasp91 predictions at the USGS locations of the two
earthquakes. For this pair the USGS depths differ by 6.7 km (108.3 km in 1997, 115.0 km in 2018).
Zhang & Wen's relocation (ESSOAr 2023) puts the 2018 event only 1.45 km deeper, and 14.20 km
from the 1997 event. Their stated direction is ambiguous ("204.33°" and "NW").

Depth matters unequally for the four reference phases. SKS leaves the source as an S wave, so a
deeper source shortens it by about 0.2 s per km. Pdiff and PP leave as P waves and shorten by only
about 0.10-0.12 s per km. SKS minus a **P-type** reference therefore changes by about **0.1 s per
km** of depth difference. SKS minus an **S-type** reference (Sdiff, SS) changes by about 0.03 s
per km or less. These are rough hand estimates from iasp91 source velocities; they are not used
in any grade.

So a few kilometres of relative depth error can make a large part of a 1 s "SKS advance" relative
to Pdiff and PP, and almost none of it relative to Sdiff. Horizontal offsets matter most for SS,
whose slowness (about 14 s/deg) is much larger than SKS's. An outer-core change moves SKS against
**every** reference by the same amount. These three effects leave different fingerprints across
the four references and the four station azimuths, so they can be separated without trusting any
catalogue location.

This also bears on a plan forecast. We forecast a 70% chance that the iasp91 correction moves
x_long by less than 0.2 s at every station. Because of the 6.7 km USGS depth difference, **we now
expect that forecast to fail**: our rough estimate is a correction near 0.6-0.7 s for the P-type
references and under 0.1 s for the S-type ones. The forecast stays as written and will be scored
as written.

## Reading F (rules fixed here)

**Data.** Primary pair E1997/E2018 from `ks32_phases.csv`. For station s and reference k
(P1 = P/Pdiff, PP, S1 = S/Sdiff, SS), the raw change of each phase is d + dp (event B relative to
event A, each after its own origin time). The value used is r(s,k) = raw(SKS) - raw(k). Origin
time and station clock cancel in r. A pair (s, k) is used only if both SKS and k pass ks32's own
gates in that band (not at the search edge; correlation >= 0.90 in the long band and >= 0.80 in
the short and mid bands; peak-to-noise >= 3 in both events).

**Model.** r(s,k) = X + [dT_SKS - dT_k](u_N, u_E, dh), where X is a change that only SKS feels,
(u_N, u_E) is the horizontal position of B relative to A in km, and dh is how much deeper B is
in km. For each phase, dT = -(dT/dDelta)/R * (u_N cos az_s + u_E sin az_s) + (dT/dh) dh.
dT/dDelta and dT/dh are central differences of the iasp91 (TauP) first arrival, taken at E1997's
USGS depth and each station's distance (steps 0.05 deg and 2 km). az_s is the azimuth from the
event to the station. R = 6371 km.

**Fit.** Unweighted least squares for (X, u_N, u_E, dh). The formal sigma of X is
s^2 (G^T G)^-1, with s^2 = RSS / (n - 4). The condition number is that of the design matrix
after scaling each column to unit length.

**Grade F (long band, 0.01-0.05 Hz, Zhou's band):**
- **untestable** if n < 7 pairs, or fewer than 3 stations, or condition number > 100, or
  sigma_X > 0.5 s, or any of these is not finite;
- **core-term** if X <= -0.5 s and X + 2 sigma_X < 0: an SKS-only advance survives after the
  relative source position and depth are fitted from the reference phases themselves;
- **no-core-term** if X > -0.25 s: the advance is accounted for by relative source position;
- **intermediate** otherwise.

The same fit is **reported, not graded**, for the short (0.1-0.5 Hz) and mid (0.05-0.2 Hz)
bands. Also reported: the fitted (u_N, u_E, dh) compared with USGS (B is 5.15 km north, 0.87 km
east and 6.7 km deeper) and with Zhang & Wen (14.20 km, 1.45 km deeper); X_fixed_USGS (the mean
of r minus the USGS location term, close to the plan's x); and every residual.

**How F combines with the headline.** F does not change T1, T2 or the joint grade. If T1 is
reproduced and F is no-core-term, the result will be stated as "the long-window advance is
reproduced, but it is fitted by relative source position", whatever T2 says. If T1 is reproduced
and F is core-term, the advance is not a location artefact in the long band. Whether it is in the
direct pulse is then T2's question.

## Forecasts (before data)

- F: core-term 25%, intermediate 15%, no-core-term 30%, untestable 30%.
- The fitted dh is closer to Zhang & Wen's 1.45 km than to USGS's 6.7 km: 55%.
- At least one station's long-band per-reference x values (plan) spread by more than 0.5 s: 55%.

## Known limits

The fit has four unknowns and at most 16 values (4 stations x 4 references). The four stations
lie within about 60 degrees of azimuth (roughly 319 to 16 degrees, northwest to north-northeast
of the source, by hand calculation), so the east-west offset is the weakest unknown. The partial derivatives are 1-D (iasp91) and linearised.
Lateral structure near the source enters SKS and the references differently, and the fit cannot
separate that from X. All phases share one record, so a timing error of the station clock
cancels, but a band-dependent response phase does not. Metadata: COLA's sensor azimuths in the
station metadata changed on 2018-08-30 (BH1 304 -> 41 deg), eleven days before the 2018 event.
Rotation to R and T uses each event's own metadata, as in the plan. A metadata error would mix
components and change waveforms more than times; we do not test the orientation here. The work is
ours and has not been reviewed.

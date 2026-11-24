# Pre-registration 32: is the reported 1 s SKS advance of the Kermadec doublet in the direct SKS pulse?

Written 2026-11-24, **before any waveform for this test was downloaded**. Code: `ks32.py` 0.1.0
(with `sdv23.py` 0.1.1 downloaded unchanged from tag v0.32.0 for data access, response removal
and rotation). Job: `ks-sks32-v1`. Licence: text CC BY 4.0, code Apache-2.0. Made with AI
agents working in MV2's name (MV2, Munim, Inc.); no human expert has checked it yet.

## The claim and the dispute

Zhou (2022, *Communications Earth & Environment* 3:97, doi 10.1038/s43247-022-00432-7)
reported that SKS from the Kermadec earthquake of 2018-09-10 (Mw 6.9, 115 km) arrived about
**1 s earlier** than SKS from the 1997-05-03 earthquake (Mw 6.9, 108 km) at the same place:
"When any of the reference waves (Pdiff, Sdiff, PP, or SS) from the two earthquakes are
aligned, the SKS wave from the second earthquake arrives about 1 s earlier". It was read as a
transient 1-1.5% P-speed change, a 2-3% density deficit and fast flow (~40 km/yr) in the upper
outer core. The paper's band is 20-100 s ("band-pass filtered between 20 and 100 s period");
its Table S2 gives SKS windows (s after origin) COLA 1411-1463, ULN 1439-1482, INK 1447-1494,
YAK 1429-1471, and changes -1.05, -1.0, -0.8 and -1.1 s.

Zhang & Wen (ESSOAr preprint, 2023, doi 10.22541/essoar.167276417.78056130/v1) studied the
same four stations. They found "no discernable relative travel time differences" in SKS at
higher frequency, and "evident changes of waveform and travel time" in SPdKS-SKPdS, which
fall inside Zhou's long, low-frequency window. Their reading: ULVZ sampling at the base of the
mantle, not the outer core. Zhou's response (doi 10.22541/essoar.167458062.24504379/v1): the
sensitivity kernels of the measurement windows "can only be created by having a structure
anomaly in the Earth's outer core". Neither preprint is in a journal (Crossref, 2026-11-24).
We found no independent re-measurement. Cui et al. (2026, GRL, doi 10.1029/2025GL119930)
call the question debated.

## Questions

**T1 (reproduce).** With public data, Zhou's band and her Table S2 windows, is SKS at the four
stations about 1 s earlier in 2018, relative to the reference phases?

**T2 (decisive).** Is the advance in the **direct SKS pulse**? We measure SKS and the same
reference phases at 0.1-0.5 Hz in 7-s windows around each pulse's peak. An outer-core speed
change along SKS's path moves the direct pulse by the same ~1 s. A change of other phases
inside a 40-50 s window at 20-100 s period does not.

## Method (all constants in `ks32.PARAMS`)

**Data.** BH? of IU.COLA, IU.ULN, CN.INK, IU.YAK from the EarthScope FDSN service, for the
four USGS events at the same place (M >= 5.8, 70-160 km depth, within 40 km): 1997-05-03
(108.3 km), 1999-04-20 (Mw 6.5, 95.7 km), 2006-05-16 (Mw 7.4, 152 km), 2018-09-10 (115 km).
Location 00 if present, else blank. Every trace is resampled to 20 Hz, then the response is
removed to velocity (pre-filter 0.002/0.004/2/4 Hz), band-passed (2 corners, zero phase) and
rotated to R/T with each event's own station metadata. The record runs 600-2600 s after the
origin.

**Predictions.** iasp91 (TauP), each event at its own USGS location, depth and origin time.
For each phase, d = (T_B,obs - T_B,pred) - (T_A,obs - T_A,pred), measured by normalised
cross-correlation of A (earlier event) and B (later event). Phases: SKS (R), SKKS (R,
reported only), P1 = P/Pdiff (Z), PP (Z), S1 = S/Sdiff (T), SS (T). Origin-time errors cancel
in d(SKS) - d(ref). The iasp91 term corrects the 1997/2018 differences in location and depth.
The uncorrected value is reported too.

**Long band (T1), 0.01-0.05 Hz.** SKS window: Zhou's Table S2 window if it contains the
iasp91 SKS time of the 1997 event at that station. It is kept at the same offsets from each
event A's predicted SKS (identical for the 1997-2018 pair). Otherwise the window is SKS -10 to
+40 s. Reference windows: prediction -10 to +40 s. Search +/- 8 s around the predicted
difference.

**Short band (T2, graded), 0.1-0.5 Hz; mid band 0.05-0.2 Hz (reported).** For each phase,
window -3 to +4 s around the largest |amplitude| of event A within +/- 5 s of its prediction.
Event B is first moved by the predicted difference plus c, where c is the station's median
long-band reference change (the alignment on the references). Then it is searched within
+/- 1.5 s.

**Station value.** x = median over passing references of d(SKS) - d(ref); negative means
SKS earlier in the later event. A measurement passes if it is not at the search edge, its
correlation is >= 0.90 (long band) or >= 0.80 (short and mid bands), and its peak-to-noise
ratio (noise window 60-10 s before the prediction) is >= 3 in both events. SKS must pass and
at least one reference must pass, or the station is untestable in that band.

## Decision rules (pre-registered; primary pair 1997-05-03 / 2018-09-10)

- **T1** (`grade_count`, threshold -0.5 s, at least 3 testable stations, else
  **untestable**): k = number of stations with x_long <= -0.5 s. **reproduced** if k >= 3,
  **not-reproduced** if k <= 1, **partial** if k = 2.
- **T2**, same rule on x_short: **present** (k >= 3), **absent** (k <= 1), **mixed** (k = 2),
  or untestable.
- **Joint reading** (`grade_joint`): T1 reproduced and T2 present -> **direct-SKS-advance**:
  the direct SKS pulse moved, as an outer-core change requires, and the window-content
  explanation does not account for it. T1 reproduced and T2 absent -> **window-content**: the
  long-window advance is not in the direct SKS pulse, which supports Zhang & Wen's reading
  over the outer-core reading. T1 not-reproduced -> **not-reproduced**. Anything else ->
  **mixed** or **untestable**.

**Reported, not graded:** every phase lag, correlation and SNR; per-reference x; x without
the iasp91 correction; SKKS minus references; the mid band; the window actually used and its
source; and the long-band x for every other pair: 1997-1999, 1997-2006, 1999-2006, 1999-2018,
2006-2018. Zhou's claim is a transient, so this series asks when any change appeared.

## Forecasts (before data)

- T1: reproduced 55%, partial 15%, not-reproduced 20%, untestable 10%.
- T2 given T1 reproduced: present 35%, absent 45%, mixed or untestable 20%.
- Joint: direct-SKS-advance 20%, window-content 25%, not-reproduced 20%, mixed or
  untestable 35%.
- The iasp91 location/depth correction moves x_long by less than 0.2 s at every station: 70%.
- Per-reference x_long values agree within 0.3 s at each testable station: 50%.

## Known limits

Four stations and one pair are graded; Zhou's other six pairs and the remaining SKS records
are not re-measured here. USGS catalogue locations are used. A relocation error moves
SKS - SS and SKS - PP more than SKS - Pdiff, because their slownesses differ. The per-reference
values show this, but we do not relocate. COLA changed sensor between the events (KS-54000 in
1997, a new location-00 sensor in 2018), and INK changed digitizer in 2007. A relative timing
between phases in the same record cancels the instrument's absolute timing, but not
differences in response phase between frequency bands. The 0.1-0.5 Hz pulses at 100-106 deg
may still hold SPdKS where it is close to SKS. A **present** grade would then still need the
SPdKS separation checked; an **absent** grade would be the stronger statement. Response
removal follows the station metadata. The work is ours, not reviewed.

## Prior work

Zhou 2022 (above); Zhang & Wen 2023 preprint and Zhou's 2023 response (above); Wen &
Helmberger 1998 (SPdKS and ULVZ); Cui, Wang, Vidale et al. 2026 GRL (CMB stable over decades in
repeating-earthquake core reflections); Wang et al. 2024 preprint (doi
10.22541/au.172538614.44654756/v1, no temporal change in high-frequency waves scattered near
the core-mantle boundary). This repository's tests 23-31 examine the same two groups' later
claims (Zhang & Wen 2025 at IU.SDV; Zhou 2026 from nuclear-test PKP).

# Pre-registration 9: the YKA hardware control (same earthquake, two systems)

MV2, Munim, Inc. (made with AI agents). Written 2026-11-17 (about 21:45 UTC), while YKA
stage B (icr-ykab-v1, pre-registration 8) was still queued. **No YKA pair, coda
similarity or stage-B number had been computed. No waveform from the overlap window had
been looked at.** Text CC BY 4.0, code (`ykagen.py`) Apache-2.0.

## Why this control is needed

Every out-of-sample pair at the Yellowknife array (YKA) puts a 1990-2008 event next to a
2025-2026 event, so every such pair crosses YKA's 2013-2014 hardware change. FDSN station
metadata (EarthScope, read 2026-11-17) shows what changed and what did not:

| | legacy YKB*/YKR* | current YKAB*/YKAR* |
|---|---|---|
| Sensor | Geotech S-13 short-period | Geotech S-13 short-period |
| Electronics | CNSN Yellowknife digitiser | YKA pre-amp + Guralp digitiser |
| Sampling | 20 sps | 40 sps |
| Recording | 1989-01-25 to 2014-04-02 | 2013-06-26 onward |
| Piers | same coordinates (elevations differ by about 1 m) | |

`yka.py` applies no instrument correction. It uses the same causal 1-2 Hz filter on
both systems. A difference in pre-amplifier or digitiser phase response therefore enters
every cross-generation pair. Three effects are possible:

* it raises the coda dissimilarity D, the quantity being tested;
* it lowers the direct-pulse correlation, so pairs fail stage B's cc ≥ 0.90 gate;
* it flips polarity.

Pre-registration 2 only covered polarity (the flip test). The two systems overlapped
for about nine months, so the same earthquakes were recorded by both. That allows a
direct measurement of this effect.

## Design (fixed here)

* **Events.** USGS ComCat, 2013-06-27 to 2014-04-01, M ≥ 4.7, depth ≤ 120 km, 128.5-145.5°
  from YKA (the full distance range of `yka.py`'s travel-time grid, not only the South
  Sandwich Islands). Events in the SSI box are flagged and reported as a subset.
* **Data.** Legacy channels (`YKB*,YKR*` SHZ) come from NRCan first, then EarthScope. Current
  channels (`YKAB*,YKAR*` SHZ) come from EarthScope first, then NRCan. Each generation
  is processed **separately** with `yka.py` v0.2.0's `process_event`, unchanged: causal
  1-2 Hz filter, per-element ak135 cut, noise normalisation, per-element statics, beam.
  `yka.py` is fetched at commit 7067d4f.
* **Comparison.** The current beam is the reference. This mirrors stage B exactly:
  `align_member` with the flip test (allowed, because the generations differ), an STA/LTA
  onset on the two-beam stack, then `yka_pair`. That gives the direct-pulse (-0.3 to
  +0.9 s) alignment and correlation cc_ref, and **D_hw = 1 − cc** in the coda window
  (+1 to +6 s), plus the three sensitivity windows. Element-by-element D on the shared
  piers is reported, not graded.
* **Usable event.** Stage B's SNR gates, applied to both generations: direct SNR ≥ 5 and
  coda SNR ≥ 2. The cc ≥ 0.90 gate is not applied, because it is itself under test (H2).
* **Testable** if there are ≥ 5 usable events.

## Hypotheses and forecasts

| | Statement | Rule | Forecast |
|---|---|---|---|
| T | Testable | ≥ 5 usable events | 80% |
| H1 | The hardware floor is small | median D_hw ≤ 0.05 | 65% |
| H2 | The hardware change does not block the stage-B gate | cc_ref ≥ 0.90 in ≥ 80% of usable events | 75% |
| H3 | No polarity reversal | the flip is chosen in ≤ 10% of usable events | 85% |
| H4 | The floor is below the real-pair floor | median D_hw < median YKA D of same-generation ok pairs less than 2 years apart in stage B (icr-ykab-v1, `yka_pairs.csv`) | 75% |

H4 is graded when both jobs have finished. Its stage-B side is fixed now: same-generation
pairs (gen1 = gen2) that pass all stage-B gates and are separated by less than 2.0 years.

**Reported, not graded:**

* the timing offset between the two systems (median and 10-90% range of the beam lag, in seconds);
* D_hw against SNR;
* the SSI subset;
* per-element D;
* overlay plots of the 8 highest-SNR and 6 lowest-SNR usable events.

## How the result is used (fixed now)

* **If H1, H2 and H4 hold:** the hardware change is too small to matter. Stage B's
  cross-generation D values are read as they stand. This is the expected case: same
  sensors, and only the electronics changed.
* **If H1 or H4 fails:**
  * The new-event vote F1-Y is still valid. Both of its candidate partners are legacy
    events and the new event is current-system, so the hardware term enters both sides
    of every comparison.
  * P2-Y and P1-Y are different. They rank pairs across the 2013 change against pairs
    that do not cross it, so they are confounded. They will be reported as **not
    interpretable** for cross-generation pairs, and only same-generation pairs will be
    graded in a follow-up pre-registration.
* **If H2 fails:** stage B's cc ≥ 0.90 gate removes cross-generation pairs for
  instrumental reasons. Every result of stage B that is untestable will say so. A
  follow-up must correct for the response difference measured here, under a new
  pre-registration, before any gate is relaxed.
* **If H3 fails** (most events flip): this confirms that the flip test is needed. It does
  not change any stage-B grade, because stage B already tests for the flip.

## Why it matters beyond this project

YKA is one of the two arrays behind the 2024 *Nature* backtracking result (Wang et al.),
which hand-corrected YKA timing errors in 2013-2020. A same-event, two-system
measurement of the waveform and timing differences across YKA's hardware change
directly bounds an instrumental confound that every waveform-similarity study on this
array has to rule out. It is reproducible from public data with the published code.

The job is `icr-ykagen-v1` (data_key `icr-ykagen`, separate from the stage-B volume).

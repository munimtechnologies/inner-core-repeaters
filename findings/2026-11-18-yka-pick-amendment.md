# Amendment: the YKA onset pick (pre-registration 8 and the M03-24 check)

By MV2, Munim, Inc. (made with AI agents). Written 2026-11-18, **before** `icr-ykab-v1`
(YKA stage B, [pre-registration 8](2026-11-17-yka-stageb-prereg.md)) has returned, and before
the re-pick job `icr-m0324chk-v2` has run. Text CC BY 4.0, code Apache-2.0.

## What went wrong, and what we have already seen

`icr-m0324chk-v1` ran `yka.py` v0.2.0, unchanged, on the nine M03-24 (Z004) events at the
Yellowknife array (YKA), as [pre-registered](2026-11-17-m0324-check-prereg.md). The onset
picker (`xcorr.pick_onset`) takes the **first** STA/LTA ≥ 6 trigger within ±6 s of the ak135
PKIKP time. On the stacked 9-event beam, where incoherent noise is strongly suppressed, it
triggered **4.8 s before ak135** on weak coherent energy. The large impulsive arrival came
4.9-6.5 s after the pick, i.e. at the ak135 PKIKP time.

YKA's own single-event quicklooks from stage A (2011-02-20 and 2025-12-15, both 139°, in
[`results/yka-v1/plots`](../results/yka-v1/plots/)) show PKIKP at YKA as strong and impulsive,
with its onset about 1 s before ak135 and nothing comparable 5 s earlier. We therefore read
the large arrival as PKIKP. We read the weak early energy most likely as scattered PKP
precursors, which are expected on the 120-145° range (it is not proven; see "Open question").
Consequences for v1:

* The "direct pulse" window (−0.3 to 0.9 s after the pick) held only that weak energy, so
  the direct-SNR gate (≥ 5) failed for **all 36 pairs** (direct SNR 0.8-4.5). Y1 and Y0 were
  untestable.
* The pre-registered "coda" window (+1 to +6 s) held mostly pre-PKIKP energy plus the first
  ~1.1 s of PKIKP.

**Not blind.** We have looked at the v1 stack plot and its ungated numbers. The wide windows
of v1 overlap the real PKIKP: ungated D(2003, 2025-12-13) = 0.044 and D(2011, 2025-12-13) =
0.017 (the same direction as the ILAR anomaly), D(2003, 2011) = 0.009. The re-pick test
below is therefore a **corrected, partly unblinded** measurement. It cannot count as an
independent confirmation on its own, and it will not be described as one.

## Rule A: pick validity for YKA stage B (`icr-ykab-v1`, already queued, unchanged)

Stage B uses the same picker on every merged ILAR family seen at YKA, so the same failure
can happen there. Fixed now, before any stage B output is read:

1. A multiplet's stack pick is **valid** if |`pick_minus_pred_s`| ≤ 2.0 s
   (`yka_multiplets.csv`). YKA timing errors are at most 0.2 s (stage A statics), and PKIKP
   onsets sit within about 1 s of ak135.
2. **Primary grading** of P4-Y, F1-Y / F1x-Y, P2-Y, P1-Y and P6 uses only pairs from
   multiplets with a valid pick. Grading on all multiplets is also reported and labelled
   "includes invalid picks".
3. If more than 25% of the multiplets with ≥ 2 YKA members have invalid picks, or if F1-Y is
   untestable and at least one would-be F1-Y family has an invalid pick, stage B is re-run
   as `icr-ykab-v2` with the constrained picker of `yka_repick.py` (the same picker, search
   ak135 ± 2.0 s). That run then becomes primary. Every other rule, gate, window and
   forecast of pre-registration 8 stays as written. The v2 run is described as "after an
   amendment for mis-picks".
4. The pick diagnostics (original ±6 s pick, constrained pick, envelope onset) are
   reported for every multiplet.

## Rule B: the M03-24 re-pick test (`icr-m0324chk-v2`)

* **ILAR (R0, C1, C2): unchanged.** v1 could not run them because the ILAR beams cached by
  `icr-m0324-v1` are on another machine (compete dataset volumes are per host). v2 first
  rebuilds them with the published `m0324.py` (commit 7067d4f), then runs `m0324_check.py
  --no-yka` unchanged. R0 now also checks that the rebuilt beams reproduce the published D
  values (max |ΔD| < 0.001).
* **YKA (Y1R, Y0R):** `m0324_check.run_yka` is run unchanged, except that the onset is
  picked within ak135 ± 2.0 s (`yka_repick.py`). The gates (direct cc ≥ 0.90, direct SNR
  ≥ 5, coda SNR ≥ 2), the windows (−0.3 to 0.9 s alignment, +1 to +6 s D) and the decision
  rules are those of Y1 and Y0. Y1R/Y0R are interpreted only if the pick is valid
  (|pick − ak135| ≤ 2 s).
* How it will be described: if Y1R is testable and A > B, the result is "the 2025-12-13
  anomaly also appears at YKA in a corrected, partly unblinded measurement". It is not a
  refutation of backtracking, and author contact still waits for a blind test (stage B, or
  new 2026 events of the patch). If A < B, the ILAR anomaly is not supported at YKA.

### Forecasts (written before v2 ran)

* Pick valid: 90%.
* Y1R testable: 75%. If testable, A > B: 70%. This is high because the v1 windows already
  overlapped PKIKP.
* Y0R testable: 70%. If testable, holds: 30%. v1 ungated D(2003, 2011) was the smallest of
  the three.
* ILAR R0 holds after the rebuild: 85%. C1 robust: 55%. C2 window-robust: 50% (unchanged).

## Also re-filed: `icr-bcdf-v2`

`bcdf.py` v0.5.0, unchanged, with Q1-Q3 and the forecasts as pre-registered. v1 found no
cached ILAR data on its host. v2 first rebuilds the cache with the published `run.py` (same
settings as `icr-ilar-v1`).

## Open question (not part of any decision rule)

The weak energy about 5 s before PKIKP at YKA repeats across the M03-24 events. Two readings
are possible:

* **(a)** scattered PKP precursors from near the core-mantle boundary, expected at 120-145°;
* **(b)** a weak, very early DF on a polar path, as at ILAR (where DF is about 4.35 s early
  and weak at 1-2 Hz, and PKPbc dominates). The big arrival would then be a later phase.

The quicklooks favour (a). Beam slowness and back-azimuth of both arrivals at YKA would
settle it, and that needs its own pre-registration.

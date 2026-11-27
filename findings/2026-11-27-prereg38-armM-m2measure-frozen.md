# Pre-registration 38, arm M: measurement front end frozen (m2measure 0.1.0)

MV2, Munim, Inc. (mv2@munimtech.com), 2026-11-27. Text CC BY 4.0, code Apache-2.0.
Made by AI agents working in MV2's name; no human expert has checked it.

**Status: FROZEN, before any M2 waveform is fetched.** Code: `m2measure/m2measure.py` 0.1.0 with
`m2measure/test_m2measure.py` (standard library; expected values derived by hand, listed below).
It writes the input file of the decision core frozen in
[findings/2026-11-27-prereg38-armM-plan.md](2026-11-27-prereg38-armM-plan.md) (m2core v0.1.0, commit
063ca49). It changes nothing in that plan. No YKA or ILAR record of any Wang et al. 2024 multiplet
other than M03-24 (~ Wang's A) has been fetched or measured by us. As the plan says, the job
(`m2-v1`) is filed only after zw33-v1, zw36-v1 and dprk37probe-v2 have run. (zw35g2-v1, also named
there, has already run.) Since then wsc39-v2 has been queued ahead of it as well.

## What it does

Each multiplet is measured in its own pass, on its ComCat-matched members only. The multiplets
are those of Wang, Vidale et al. 2024 (Nature 631:340, https://doi.org/10.1038/s41586-024-07536-4)
Supplementary Data Table 3, except I and M. Multiplet A uses 734 (2008-02-17) for the misprinted
724, as in pre-registration 39 amendment 1. Event matches are the same as in wsc39.py 0.2.0.
Wang 711 has no ComCat match, so K has two located members.

- **YKA.** D = D_c310 (1 - cc in the 3-10 s late PKIKP coda) and `ok` from yka.py v0.2.0. They come
  through m0324_check v0.8.0's `run_yka`, with yka_repick 0.9.2's +/-2 s picker. This is exactly
  how m34.py ran M03-24. Only the event list changes.
- **ILAR.** The events are beamed with run.py v0.1.0 (`fetch_event` / `process_event`) at the
  pinned commit 7067d4f. This is the code that ilar-v1 and m0324.py ran. Repeater edges among a
  multiplet's members come from run.py's `detect_pairs`, using ilar-v1's thresholds (beam cc >= 0.90
  in the PKPbc+ab window, lag up to 8 s, snr_group >= 3.0). Only events with snr_group >= 3.0 enter.
  Then `ilar_frame` and `all_D` of m0324_check v0.8.0 run unchanged. `ok` = D finite and, for both
  events, snr_group >= 3.0, reference-window SNR >= 5.0 and DF-window SNR >= 3.0. These SNRs are
  computed in the multiplet's frame, as run.analyze_multiplet does.

## Choices the plan left open, fixed here

1. **No catalogue-distance limit on ILAR edges.** run.py used 60 km and m0324.py 90 km. Arm M's
   membership is Wang's, and the plan forbids catalogue-proximity rules. So pairs inside a
   multiplet are compared whatever their ComCat separation. The nine M03-24 events are all within
   70 km of each other, so their edges are unchanged.
2. **ILAR members that no edge connects to the multiplet's highest-SNR event** are not aligned
   (m0324_check's `ilar_frame`). Their pairs are written with D = nan and ok = false, and m2core
   ignores them.
3. **Gate, exactly.** The same code runs first on the nine M03-24 (Z004) events from
   z004_events.csv, with fresh downloads (data_key m2, ILAR beams in proc_m2_v1). It must
   reproduce our published rows: results/m0324chk-v2/yka_z004_pairs.csv (column D_c310) and
   results/m0324chk-v2/ilar_z004_matrix.csv (column D_reproduced). At each array, the max
   |difference| must be <= 0.001 over the pairs common to both, and at least 90% of the published
   finite rows must be present. m34-mirror-v1 already reproduced the YKA rows exactly (36 of 36,
   difference 0.0).
4. **m2core runs only if the gate passes.** It runs inside the same job (m2core v0.1.0, from next to
   the script or downloaded at commit 063ca49) and writes <out>/m2core/. If the gate fails, m2core
   is not run and nothing is graded until a recount note explains why.
5. **Time budget.** 200 min. Multiplets are measured in letter order, and the outputs are
   rewritten after each one. A multiplet skipped for time gets rows with ok = false. m2core treats
   these as unmeasured, and the summary lists them.

## Hand checks (test_m2measure.py)

- decimal_year("2003-03-21T10:47:47.620Z") = 2003.2176708 (z004_events.csv convention).
- 14 measured multiplets (no I, M). A has 7 members including usp000fz4d (Wang 734). K has 2.
  No event is in two multiplets.
- ILAR gates: (3.0, 5.0, 3.0) passes, 2.9 / 4.9 / 2.9 fail, and a nan D fails. A pair without D,
  or with a member outside the frame, gets nan / false.
- frame_snr on [1,1,1,1,2,2,2,2], noise 0-4: window 4-8 gives 2.0. Shift -1 gives
  sqrt(3.25). Shift +1 runs off the end and gives nan. A zero beam gives nan.
- YKA rows take D_c310 and ok from the pipeline row, in either event order. A missing pair gives
  nan / false.
- Gate: a 0.0005 difference passes and 0.002 fails. No common pairs fails. With 36 published
  rows, 32 present fails and 33 present passes.

## Forecast (before the run)

- Gate passes: 70%. (YKA already reproduced exactly in m34. ILAR beams have never been rebuilt
  from fresh downloads.)

## Job (to file after the queue clears)

`m2-v1`: python:3.12; `pip install numpy scipy matplotlib obspy requests`; run
`python -m unittest -v test_m2measure`, `python fetch_pinned.py --dest v080`, then
`python m2measure.py --mods v080 --data-dir /data --out-dir /out --time-budget-min 200`.
Network on, data_key m2, 1 CPU, 2 GB, timeout 260 min.

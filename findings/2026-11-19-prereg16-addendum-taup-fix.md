# Pre-registration 16, addendum T: a phase-name bug stopped the first run before any waveform

By MV2, Munim, Inc. (made with AI agents). Written 2026-11-19, about 13:20 UTC, before any
waveform for pre-registration 16 has been fetched. Text CC BY 4.0; code (`patch_nea.py`,
`test_patch_nea.py`) Apache-2.0. No human expert has reviewed it. Marker: ADDENDUM-T-TAUP-FIX.

## What happened

Compete job `icr-nea-v2` ran `nea.py` 0.1.0 from tag v0.14.7 on 2026-11-19 (12:49-12:55 UTC).
It confirmed that pre-registration 16 was public when it started (`prereg_public_at_start: true`),
found the station list, and chose 10 stations: IC.HIA, IC.MDJ, II.ERM, IU.BILL, IU.MAJO, IU.PET,
IU.TIXI, IU.YAK, IU.YSS and IU.ADK.

It then stopped working for every event. `nea.py` asked ObsPy's TauP for the phases
`["PKIKP", "PKPab"]`. ObsPy does not accept the branch suffix "ab", so every travel-time call
raised `ValueError: Invalid phase name: ab could not be parsed in PK...`. Of the 1,340
event-station rows, 1,339 have status `taup-error`; the other one is IU.ADK with no channel
epoch for that event. Because the travel time comes before the download, **no waveform was
requested, downloaded or looked at**. Every family has `n_in = 0`, and the job reported the
headline `untestable`. That label comes from the bug, not from the data, so **pre-registration
16 is not graded by this run.** The raw outputs are in `results/nea-v2-failed/`.

The unit tests passed because they use synthetic traces and never call TauP. That gap is what let
the bug through.

## The fix (and nothing else)

`patch_nea.py` downloads nothing new. It rewrites exactly two expressions in the released
`nea.py` 0.1.0 and refuses to run unless each matches exactly once:

| Released (v0.14.7) | Patched |
|---|---|
| `phase_list=["PKIKP", "PKPab"]` | `phase_list=["PKIKP", "PKP"]` |
| `next((a.time for a in arr if a.name == "PKPab"), NAN)` | `max((a.time for a in arr if a.name == "PKP"), default=NAN)` |

In ObsPy, the phase `PKP` returns every P ray that turns in the outer core: the AB branch and,
close to 155 degrees, the BC branch. BC (where it exists) always arrives before AB, so the latest
`PKP` arrival is the AB arrival that pre-registration 16 meant by "ak135 PKPab".

Nothing else changes: not the families, the station rule, the windows, the band, the gates
(cc_ab >= 0.80, SNR 3 and 4, AB-DF >= 8 s), the rank transform, P16-0, P16-1, P16-2, the
headline rules, the secondaries or the forecasts. Addenda G (geometry) and H (prior antipodal
work) stand as written. `VERSION` stays "0.1.0"; the job writes `nea_patch.json` with the SHA-256
of `nea.py` before and after patching.

## The rerun

Compete job `icr-nea-v3`:

1. downloads `nea.py` and `test_nea.py` from tag v0.14.7;
2. runs `patch_nea.py`, which waits up to 120 minutes for this addendum to be public on the main
   branch and exits with code 3 if it is not;
3. runs the patched `nea.py`, exactly as `icr-nea-v2` did;
4. runs `test_patch_nea` (including a TauP check that the AB-DF time is between 0 and 60 s at
   156, 165, 175 and 179.5 degrees) and the released `test_nea`.

## Forecast graded now

One forecast in pre-registration 16 depends only on station metadata, which `icr-nea-v2`
fetched: **"At least 4 stations chosen" (60%) came true, with 10 stations.** All the other
forecasts wait for `icr-nea-v3`.

## Lesson for this project

Any job that calls TauP must have a unit test that calls TauP with the exact phase list it uses,
run inside the job's container. `test_patch_nea.py` adds that test for this job.

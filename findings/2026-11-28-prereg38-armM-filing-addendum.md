# Pre-registration 38, arm M: filing addendum (order only; no rule or code changes)

MV2, Munim, Inc. (mv2@munimtech.com), 2026-11-28 (written ~06:30Z, before `m2-v1` runs and before any
waveform of a Wang et al. multiplet other than M03-24 has been fetched for arm M). Text CC BY 4.0, code Apache-2.0.
Made with AI agents working in MV2's name; no human expert has checked it.

## What changes

The [arm M plan](2026-11-27-prereg38-armM-plan.md) says "the job waits until zw33-v1, zw35g2-v1, zw36-v1 and
dprk37probe-v2 have run", and the [frozen front-end note](2026-11-27-prereg38-armM-m2measure-frozen.md) repeats
that order. zw33-v1 and zw35g2-v1 have run. zw36-v1 and dprk37probe-v2 have waited in the shared compute queue
since 2026-11-25 for free memory. **We file `m2-v1` now instead of after them.** This is the only change.

## Why this does not touch the test

- The wait was a queue-ordering choice, not a decision rule. Neither the plan nor the frozen note makes any arm M
  rule depend on the outputs of zw36-v1 (pre-registration 36, a different question) or dprk37probe-v2
  (pre-registration 37, North Korean explosion records). They use different data and none of their outputs is
  read by m2measure or m2core.
- **The code is the frozen code, byte for byte.** The job does not use a local copy. `get_frozen.py` downloads
  `m2measure/m2measure.py` and `m2measure/test_m2measure.py` from the freezing commit
  (a7b58bdb1d27cfe0fce09fd54ced519d557e5346) and stops unless their git blob SHA-1s equal the ones GitHub lists
  for that commit (eb92f4cf... and 61e7144c...). `fetch_pinned.py` (the helper named in the frozen job draft)
  fetches the v0.8.0 pipeline at commit 7067d4f as before. m2core is fetched at its frozen commit 063ca49 by
  m2measure itself.
- The gate is unchanged: the same code must first reproduce our published M03-24 rows at both arrays (max
  |difference| <= 0.001, >= 90% of rows present), or nothing is graded.
- Every rule, forecast and threshold of arm M (M2a, M2b), and of arm T
  ([plan](2026-11-28-prereg38-armT-plan.md), gradet 0.1.0) which reads m2-v1's output, stays as frozen.

## One thing we now know that bears on the gate (stated before it runs)

zw33-v1 ([results](2026-11-28-prereg33-results.md)) showed our plain STA/LTA onset picker firing at the edge of
its +/-6 s search on Zhang & Wen's 1-2 Hz YKA beams. m2measure does not use that path: at YKA it uses yka_repick
0.9.2's constrained +/-2 s picker on our own beams, the setup that reproduced all 36 published M03-24 YKA rows
exactly in m34-mirror-v1. If it misfires on other multiplets, the per-multiplet `pick_valid` flag and the `ok`
column will show it; we will report it and not re-tune anything.

## Job

`compete/m2/job.json`: name m2-v1, python:3.12, setup `pip install numpy scipy matplotlib obspy requests`;
run `python -m unittest -v test_get_frozen`, `python get_frozen.py --dest .`,
`python -m unittest -v test_m2measure`, `python fetch_pinned.py --dest v080`,
`python m2measure.py --mods v080 --data-dir /data --out-dir /out --time-budget-min 200`; network on;
data_key m2; 1 CPU; 2 GB; 260 min.

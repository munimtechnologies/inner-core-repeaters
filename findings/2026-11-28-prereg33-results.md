# Pre-registration 33 graded: E1 inconclusive, P35 and P35b untestable (zw33-v1)

MV2, Munim, Inc. (mv2@munimtech.com), 2026-11-28. Text and data CC BY 4.0, code Apache-2.0.
Made with AI agents working in MV2's name; no human expert has checked it.

Plan: [2026-11-24-prereg33-zw26-e1-closure-plan.md](2026-11-24-prereg33-zw26-e1-closure-plan.md), written before
any waveform in Zhang & Wen's (2026, JGR Solid Earth, https://doi.org/10.1029/2026JB034595) Zenodo bundles
(https://doi.org/10.5281/zenodo.19625711, CC BY 4.0) was read. Job `zw33-v1` (zw33.py 0.1.0, closure35.py 0.1.0,
unchanged) ran on the Dell on 2026-11-28, 01:16-01:30Z: tests passed, run exit 0, all three zips matched their
Zenodo MD5s, and all six M03-24 and three multiplet-F events were present at both arrays. Outputs are committed
unmodified in [`results/zw33-v1/`](../results/zw33-v1/).

## Grades, exactly as the plan states

| test | grade | why |
|---|---|---|
| E1 (8 graded items) | **inconclusive** | 0 agree, 0 disagree, 8 nd in the graded window D_c310 (1-2 Hz) |
| P35 | **untestable** | gate G1 failed: YKA 2005-06-10 / 2022-05-19 is nd (ratio D / floor = 2.79, between the 2 and 3 cuts), so no YKA recovery is reproduced |
| P35b (closure on M03-24) | **untestable** | 0 usable triples |

Forecasts: E1 inconclusive had 25%; G1 passing had 55% (failed).

## Why almost everything is nd (diagnosis, report-only)

29 of the 30 M03-24 pairs (15 per array) and 3 of the 6 multiplet-F pairs (all at ILAR) in the graded 1-2 Hz
band fail the pre-registered **direct-window gate** (direct-pulse SNR >= 5 in both beams, cc_ref >= 0.90). The
failing direct SNRs are 0.8-4.2, while the same YKA beams have window SNRs of 40-118. The likely cause is our
onset picker, not the data: in the 1-2 Hz band the STA/LTA onset of 14 of the 15 YKA M03-24 pairs lies at -5.5
to -6.0 s from the predicted PKIKP time, at or next to the edge of the +/-6 s search, and all 15 ILAR pairs lie
at -4.0 to -6.0 s. The direct window is then cut on pre-signal noise. The pairs whose onset
lies inside (YKA 2022-03-07 / 2024-03-21 at -2.3 s; all YKA multiplet-F pairs at -1.8 to -1.9 s) pass the gate.
Because the coda window is placed after the onset, a mispicked onset also moves the "3-10 s" window onto the
direct arrival, so the ungated zw33 D values are **not** late-coda values and are not used anywhere. (A possible
physical reason, not checked: South Sandwich to YKA is about 139 degrees, inside the range where PKP precursors
arrive ahead of PKIKP.)

We do not fix the picker and re-run zw33. The paper's labels are now known to us, so any retuned picker would be
chosen after seeing the outcome. The frozen arm M front end (`m2measure` 0.1.0), which uses the YKA pipeline that
reproduced our published M03-24 rows exactly, is the instrument for these records.

## Report-only numbers that pass every gate

- YKA 2022-03-07 / 2024-03-21: **same** in all three windows (D_c310 = 0.0041, floor 0.0255). The paper calls
  this pair quiescent; agree (report-only item).
- YKA multiplet F, D_c310: 2005/2012 changed (D 0.370, floor 0.010), 2012/2022 changed (0.300, floor 0.0096),
  2005/2022 nd (0.035, floor 0.0125). D(2005, 2022) is 0.095 of D(2005, 2012), so the 2022 event is much
  closer to 2005 than 2012 is, as the paper's recovery reading says. Under the frozen ratio rule this still
  does not count as "same". In D_c16 and D_c13 the 2005/2022 pair is "changed".
- ILAR multiplet F: every pair fails the direct gate, so the decisive P35 pair (ILAR 2005 / 2022) was not
  measured.

## What this closes

Line E (E1) is closed **inconclusive** because our own measurement failed its gates. This is a failure of our
pipeline on these records, and it says nothing for or against the paper. P35 stays untestable. It was already
untestable for source reasons (multiplet F's 2012 event fails the mantle-P source check), and now also by G1.

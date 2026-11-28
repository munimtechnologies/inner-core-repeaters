# Pre-registration 40, addendum 2: grader bug fix (array name case), gradeq40 0.1.1

MV2, Munim, Inc. (mv2@munimtech.com), 2026-11-28 ~10:15Z. Text CC BY 4.0, code Apache-2.0.
Made with AI agents working in MV2's name; no human expert has checked it.

## The bug

`gradeq40.py` 0.1.0 (commit c1fc29f) keeps a row only if its `array` field equals `"YKA"`. The pinned input,
`results/h13-v1/h13_pairs_scored.csv` at commit 4c44a95 (git blob 3bc1eacf..., 275,396 bytes, re-checked
2026-11-28 through the GitHub contents API), writes the array in lower case: `yka`. So 0.1.0 would have
selected **no pairs** for either set and graded the test "untestable" for a code reason, not a data reason.
The unit tests passed because their synthetic rows used `"YKA"`.

Found 2026-11-28 ~10:05Z by reading the grader while preparing for its result, before the job `gradeq40-v1`
had run (it was still waiting in the shared compute queue). See addendum 1 for what had been displayed by then.

## The fix (mechanical)

- `select()` and the row count now compare the array name case-insensitively (`is_yka`: strip, upper case,
  equal to `"YKA"`). This is what the plan already says ("array YKA").
- The grader now waits for this addendum, as well as the plan, to be public on `main` before grading.
- Nothing else changed: sets P and C, the source cut, the timing-era exclusion, the thresholds (20 pairs,
  4 families, 5 long pairs, f_pos 0.95), the seeds (0 for P, 1 for C), 2,000 resamples, the headline order,
  the report-only lines and the input commit and blob are identical to 0.1.0.
- New tests: any case of `YKA` is kept and `ILAR` in any case is not; a table in the exact header and
  spelling of the real file (lower-case array, CRLF line ends, synthetic values) gives the hand-derived
  selections; the end-to-end grade is the same for `YKA` and `yka`.

## Job

`gradeq40-v2` replaces `gradeq40-v1` (only the newest version of a job name runs): python:3.12, no setup,
`python -m unittest -v test_gradeq40` then `python gradeq40.py --out-dir /out`, network on, 1 CPU, 0.5 GB,
60 min. If v1 runs anyway, its "untestable" output is a pipeline failure and is not reported as a grade.

## Rule learned

Every grader that reads an existing table gets one test built on that table's exact header and value
spelling, not on fixtures written from memory.

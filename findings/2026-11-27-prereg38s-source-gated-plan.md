# Pre-registration 38, arm WS: does the turn signal in Wang et al.'s grades survive the source gate?

MV2, Munim, Inc. (mv2@munimtech.com), 2026-11-27. Text CC BY 4.0, code Apache-2.0.
Made by AI agents working in MV2's name; no human expert has checked it.

**Status: written after arm W was graded ([results](2026-11-27-prereg38-W-results.md)) and before
any result of the source census wsc39 (pre-registration 39, job wsc39-v2) exists.** It is a labelled
follow-up of a graded result, not a new confirmatory test. Code: `gradeclock38s/gradeclock38s.py`
0.1.0 with hand-derived offline tests `gradeclock38s/test_gradeclock38s.py`. It imports the frozen
`gradeclock38.py` 0.1.0 unchanged (same clock grids, 500 date shuffles, seed 38, W1-W4 rules).
Job `gradeclock38s-v1` is filed only after wsc39-v2's summary is committed to
`results/wsc39-v2/wsc39_summary.json`, which the job reads.

## Why

Arm W found that Wang et al. (2024)'s own by-eye grades at YKA and ILAR favour one turn clock far
beyond chance, at both arrays. A turn clock explains a pair graded "different" by the core having
moved. But a pair can also look different because its two earthquakes are not one source, and
[gate G2 of pre-registration 35](2026-11-27-prereg35-g2-results.md) found exactly that for the 2012
member of a triplet inside Wang's multiplet F. If the grades that carry the turn involve such members,
the turn signal is partly a source signal. This arm removes them, by a rule fixed now.

## Rule

Each Table 2 pair gets a flag from wsc39-v2: **odd-member** (either event is an odd member of any
multiplet), **S** (the pair itself decided not one source), **N**, **inconclusive**, **untestable**,
**yardstick** (both events are M03-24 reference events: Wang 718, 734, 742, 848, 854) or
**untested** (no decision, e.g. events 855 and 858, which Table 3 omits).

- **Gate:** the job reproduces arm W's counts and best test-grid scores (143 rows; 193 and 141
  comparisons; 0.124352 and 0.120567), and wsc39's summary must be version 0.2.0, finished, with 734
  in multiplet A. Otherwise nothing is graded.
- **WS1 (primary):** drop every pair flagged odd-member or S (rows are removed, so their events also
  leave the shuffle pool). Grade W1-W4 exactly as arm W. If nothing is dropped, WS1 is arm W:
  "unchanged".
- **WS2 (strict):** keep only pairs flagged N or yardstick. Grade W1-W4; an array with fewer than 30
  comparisons is untestable, as in arm W.
- **Reading per array:** *survives* if W1 structure and W2 turn-needed; *lost* if W1 no-structure or W2
  not-needed; otherwise *weakened*. **Headline (WS1):** survives at both arrays: "survives"; lost at
  both: "lost"; else "mixed".

## How we will read it

- *survives*: the turn signal in the published grades does not come from earthquakes that fail the
  mantle-P source test. That removes the most direct source objection to the backtracking reading of
  these grades (it does not make the grades independent; they remain by-eye grades).
- *lost* or *mixed*: some or all of the turn signal rests on pairs whose members are not one source;
  wherever arm W is cited we will say so.
- WS2 is the stricter, smaller test; if it is untestable, WS1 carries the conclusion.

## Forecasts (before any census result)

WS1 headline: survives 60%, mixed 25%, lost 10%, unchanged 5%. WS2 testable at YKA 50%, at ILAR 35%;
where testable, survives 55%.

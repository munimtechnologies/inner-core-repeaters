# Pre-registration 30, amendment 1: inputs committed, a report-only starred-18 run

*inner-core-repeaters, MV2 of Munim, Inc. Written 2026-11-22 ~17:45 UTC, in the same commit
series that adds the inputs. No prediction from Table S1 has been computed by us, and the graded
job (pkp-pair30-v1) has not returned. Code: `star30.py` 0.1.0 + `test_star30.py` (Apache-2.0).
Text CC BY 4.0. Made with AI agents working in MV2's name; no human expert has checked it yet.*

## What is added

- `data/zhou2026_table_s1.csv`: Zhou's (2026) Table S1, transcribed verbatim (23 tests, the
  table's own header, stars kept).
- `data/zhou2026_methods.json`: `distance_correction = "not-stated"`.

## Two things the table taught us before any computation

1. **The table has 23 rows, not 18.** Its caption says that stars mark the 18 tests with the best
   data, which formed the 112 pairs with differential epicentral distances below 0.05 deg. The 5
   unstarred tests (1985-05-08, 1989-05-11, 1989-06-03, 1991-05-18, 1991-07-15) were used only in
   the closest-pair (< 2 km) analysis. The plan says "all S1 pairs within 0.05 deg", and
   `pkppair30.py` 0.1.0 reads every row, so **the graded headline stays as written: all 23 rows.**
   We add a **report-only** run on the 18 starred rows (`star30.py`, job `pkp-star30-v1`). It
   calls the unchanged `pkppair30.py` from tag v0.29.0 and only drops unstarred rows before
   parsing. It is not graded and will not replace the headline whichever way it comes out.
2. **The coordinates are not from the CMR database.** The caption gives Marshall et al. (1993)
   before 1990 and the ISC Bulletin after 1990. The plan's text named CMR/Yang et al. (2002); that
   was our error. The code does not use the catalogue's name, only its numbers.

## The methods entry

The browser copy of the paper's Methods sentence on distance correction was cut off, so we do
**not** know whether the paper states one. Under the plan, anything other than an explicit "none"
or "catalogue-1D" is "not-stated", and both readings are the headline (R:<grade> / C:<grade>).
If the sentence is found later, the results note will quote it as context. We will not switch the
headline to one reading after seeing the grades.

## Logged now for a later test (not part of P30)

PKiKP-PP at Borovoye spans roughly two minutes of record. A sample-interval error of a fraction
ε scales that differential by ε, so ε = 0.1% alone gives about 0.1 s, the size of the reported
signal. Clock offsets cancel in a differential time, but the clock rate does not. The Borovoye
systems changed their nominal sampling (e.g. 0.032 s to 0.024 s in mid-1981, from the Zenodo copy
of the archive). A later pre-registration may test this with PPP-PP on the same records. The three
explanations predict different PPP-PP / PKiKP-PP ratios: about zero for a core change, the
ray-parameter ratio for a position error, and the ratio of the two differential lengths for a
time-base error.

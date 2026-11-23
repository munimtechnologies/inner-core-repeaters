# Question C, step 2 results: what public catalogues say about where the French tests were

*inner-core-repeaters, MV2 of Munim, Inc. Written 2026-11-23 ~13:00 UTC from job `pkp-loc29-v1`
(pkploc29.py 0.1.0, ran 2026-11-23 11:52-11:53 UTC on the Blade). Results: `results/pkploc29-v1/`.
This is a **scan, not a graded test**. It grades nothing and changes no grade. Text CC BY 4.0, code
Apache-2.0. The shaft table derives from English Wikipedia (CC BY-SA 4.0). Made with AI agents
working in MV2's name; no human expert has checked it yet.*

The scan was planned in `findings/2026-11-22-question-c-shaft-positions.md` as context for
pre-registration 30. Test 30 was graded first, from its own job (`results/pkppair30-v1/`), so
this note only records what the scan found. It came back about 27 hours after it was queued, because
the compute queue was full.

## What ran

- No errors (`errors: []`), 0.53 min of run time.
- **Wikitext check of the shaft table:** 123 of 142 rows matched the live Wikipedia series pages.
  The 19 rows that did not are all 1989-91 tests. The 1989-91 page came back as 91 bytes (0
  coordinates), so those rows were *not checked*. They were not found wrong.
- **Exact sensitivities** (iasp91, from the Blacknest reference point 21.835 S, 138.910 W;
  141.66 deg to Borovoye at azimuth 331.8 deg), per km of shot movement toward Borovoye:
  PKiKP-PP **+36.3 ms**, PKIKP-PP +38.8 ms, PPP-PP **-15.4 ms**. A position error therefore moves
  PPP-PP by **-0.425 times** its effect on PKiKP-PP. This is the ratio pre-registration 31 uses to
  tell position from core change on the waveforms.
- **Catalogues:** USGS ComCat matched 77 tests and the ISC Bulletin 79.
- **Text sources that gave no positions:** the Yang et al. appendix (180 French Polynesia lines,
  none parsed with a date), the official CIVEN and AVEN test lists (dates and zones, no
  coordinates), and the OCR of Blacknest AWE Report O2/95 (no dated lines parsed). So the scan has
  no Yang-minus-shaft or Blacknest-minus-shaft comparison.

## What it found: modern catalogue epicentres scatter far more than the signal

For every test both in a catalogue and in the shaft table, the job converts the catalogue-minus-
shaft offset into the PKiKP-PP time it would fake at Borovoye. Per epoch (ms, from
`pkp29_summary.json`, epochs.*_minus_shaft_PKiKP_PP_ms):

| epoch | ISC: n, mean, sd | USGS: n, mean, sd |
|---|---|---|
| 1977 | 3, -773, 1055 | 3, -384, 522 |
| 1982-83 | 7, +134, 549 | 7, +87, 350 |
| 1988-90 | 13, +166, 305 | 13, -10, 206 |
| 1995 | 4, +312, 193 | 4, +77, 112 |

The within-epoch scatter is 112 to 1055 ms. That is several times the whole 0.1-0.2 s history
reported by Zhou (2026). At 36 ms per km, it means ISC and USGS put these tests 3 to 30 km from
their disclosed shafts along the path to Borovoye. Restricting to tests of 20 kt or more does not
change the picture (ISC sd 221-329 ms, USGS sd 74-281 ms).

None of this reaches Zhou's analysis directly. Her Table S1 takes pre-1990 epicentres from
Marshall et al. (1993), the Blacknest locations, which are much better than these teleseismic
catalogues. It does matter for the **after-1990 tests**, which Table S1 says come from the ISC
Bulletin. Those include both 1995 tests, the epoch on which test 30's reading-C grade turned.
Two details from `pkp29_tests.csv`:

- For 1995-10-27 (Aepytos), Table S1 gives 21.852 S, 138.952 W. The ISC hypocentre the job
  downloaded is 21.7955 S, 138.8818 W and USGS gives 21.891 S, 138.983 W. Each is roughly 5-10 km
  from Table S1 (our hand arithmetic, not a job output).
- For 1995-12-27 (Themisto), Table S1 gives 21.792 S, 139.042 W. The job has ISC 21.7553 S,
  138.9997 W and USGS 21.881 S, 138.973 W. Again, neither is within 5 km (hand arithmetic).

So the Table S1 position of the 1995 tests does not match today's ISC or USGS epicentres, and
the shaft table flags both shots as rough ("pro-forma") positions. Of the four epochs, the 1995
epoch has the least certain geometry from every source we have.

## What it means for the outer-core question

- It strengthens the conclusion of test 30: no public location set is good to the ~1 km that
  PKP-PP pairing at Borovoye needs. Catalogue position cannot settle the claim either way.
- It makes the waveform test (pre-registration 31, `findings/2026-11-23-prereg31-brv-waveform-plan.md`,
  v0.32.0) the deciding step. PPP-PP on the same records, at -0.425 times the PKiKP-PP position
  effect, separates position from core change without trusting any catalogue. Job
  `brv-wave31-v2` is queued.

## Files

`results/pkploc29-v1/pkp29_summary.json` (summary, sensitivities, per-epoch statistics) and
`results/pkploc29-v1/pkp29_tests.csv` (142 tests: shaft, USGS and ISC positions, along/across-path
km and iasp91 PKIKP-PP, PKiKP-PP and PPP-PP relative to the reference). The job's text copies of
third-party documents (Wikipedia pages, Yang appendix lines, CIVEN/AVEN lists, Blacknest OCR) are
not re-published here.

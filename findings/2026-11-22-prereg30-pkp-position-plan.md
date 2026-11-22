# Pre-registration 30: can test-site position alone make the reported outer-core change?

*inner-core-repeaters, MV2 of Munim, Inc. Written 2026-11-22 ~14:00 UTC. We have not read
Table S1 of the paper, and no position prediction for its events has been computed. Code:
`pkppair30.py` 0.1.0 + `test_pkppair30.py` (Apache-2.0). Text CC BY 4.0. Made with AI agents
working in MV2's name; no human expert has checked it yet.*

## Claim under test

Y. Zhou (2026), *Temporal changes in the Earth's outer core from nuclear-test PKP waves*, JGR Solid
Earth 131, e2026JB033798 (doi:10.1029/2026JB033798, CC BY):

- Data: French underground tests recorded at Borovoye, Kazakhstan (BRVK, 141-142 deg away), from
  the LDEO Borovoye archive (Kim & Ekström 1996).
- Pairs: 112 pairs from 18 tests whose catalogue epicentres agree within 0.05 deg.
- Result: PKrKP-PP and PKiKP-PP double differences are about 0.1-0.15 s faster in 1982-83 and
  1988-90, and 0.15-0.2 s slower in 1995, relative to 1977. The paper reads this as change in a
  region of the shallow outer core more than 700 km across.

If true, this would be the first seismic evidence of outer-core change on a 2-5 year scale.

## Why position matters

At 141.6 deg, PP leaves the source with a ray parameter of about 6 s/deg and PKiKP/PKP with about
2 s/deg. A shot moved 1 km toward BRVK therefore changes PKP-PP by about +35 ms. The 0.05 deg
pairing tolerance allows 5.5 km, which is up to about 0.2 s: as large as the reported signal.

The catalogue epicentres behind the pairing are from the CMR database (Yang et al. 2002). They are
not the disclosed shaft positions. France has published which shaft each test used, and our table
`data/mururoa_shafts_1977_1996.csv` (v0.28.0) gives those positions. The rim shafts sit at named
points on the reef. Most lagoon shots after 1981, and all 1995 Mururoa shots, have positions that
Wikipedia flags as rough ("pro-forma"), and these carry `approx = 1`.

A hand estimate from tests of 20 kt or more (scan note of 2026-11-22) found that raw position
alone gives about -36, +18 and -90 ms for the three epochs. That does **not** match the sign
pattern. But that estimate did not use Zhou's 18 tests and did not use her pairing. This test does
both.

## Inputs (fixed now; transcribed later)

- `data/zhou2026_table_s1.csv`: Zhou's Table S1, transcribed verbatim from the open-access
  supporting information. Columns: date, time, lat, lon, mag, name, and any others. It is
  committed by a later run. That run copies the table and computes nothing from it.
- `data/zhou2026_methods.json`: `{"distance_correction": "none" | "catalogue-1D" | "not-stated",
  "quote": "<the paper's own sentence(s)>"}`.
  - "none": the paper says the pair differential times were not corrected for the pair's
    distance difference.
  - "catalogue-1D": the paper says they were corrected using the catalogue epicentres in a 1-D
    model.
  - "not-stated": anything else.
- `data/mururoa_shafts_1977_1996.csv`, as committed in v0.28.0.

The job polls the repo's main branch for the first two files. Without them it ends "no-input" and
grades nothing.

## Method (fixed)

1. **Parse and match.** Parse S1, then match each S1 test to the shaft table: same UT date, else
   +/-1 day; when several tests share a date, take the nearest origin time.
2. **Pairs.** All S1 pairs whose S1 epicentres are within 0.05 deg (great circle). Report the count
   against Zhou's 112.
3. **Travel times.** T(lat, lon) = t(phase) - t(PP) from TauP, iasp91, source depth 0.6 km, to
   BRVK (53.0581 N, 70.2828 E), for two observables:
   - **PKiKP-PP (primary).**
   - **PKIKP-PP (report-only),** standing in for PKrKP, whose ray parameter iasp91 does not give.
     The PKP-PP sensitivity is set mostly by PP.
4. **Per-test anomaly from position alone, two readings:**
   - **R (raw):** x_i = T(shaft_i). No distance correction, so the true position stays in the data.
   - **C (catalogue-corrected):** x_i = T(shaft_i) - T(S1_i). What is left after a correction made
     with catalogue positions.
5. **Inversion.** One term per calendar year, with 1977 = 0, by least squares on the pair
   differences x_j - x_i. Pairs within one year are skipped, and years not linked to 1977 through
   pairs are left out. Epoch value = the year terms in the epoch, weighted by number of tests:
   - E1982_83 (1982-1983);
   - E1988_90 (1988-1990);
   - E1995.
6. **Report-only:**
   - the same epochs using only rim tests (approx = 0);
   - the median, 90th percentile and maximum |T(S1_j) - T(S1_i)| over the pairs, i.e. how big a
     PKP-PP difference the 0.05 deg tolerance itself allows.

## Gates

- **G0:** at least 14 of the S1 tests match the shaft table.
- **G1:** S1 contains a 1977 test and at least one test in each epoch.

If either gate fails, the grade is *untestable*.

## Grade (per reading, primary observable)

Zhou's reported values, the midpoints of the quoted ranges: E1982_83 = -0.125 s, E1988_90 =
-0.150 s, E1995 = +0.175 s.

- A predicted epoch value has a sign only if |value| >= 0.02 s.
- signs_match = the number of epochs whose predicted sign equals Zhou's.
- rms_res = rms(Zhou - predicted) over the three epochs.
- rms_Zhou = 0.151 s.

Grades:

- **explains:** signs_match = 3 and rms_res <= 0.05 s.
- **partly:** signs_match >= 2 and rms_res < 0.75 x rms_Zhou.
- **does-not-explain:** otherwise.
- **untestable:** a gate fails, or an epoch value is missing.

**Headline:**

- If the paper states "none", the headline is reading R.
- If it states "catalogue-1D", the headline is reading C.
- If it is "not-stated", both grades are the headline, written R:<grade> / C:<grade>. We do not
  pick one afterwards.

## Forecast (before data)

Testable (G0 and G1 pass): 65%. The risk is S1 events that the shaft table cannot match, e.g.
safety trials.

If testable:

| | does-not-explain | partly | explains |
|---|---|---|---|
| reading R | 60% | 25% | 15% |
| reading C | 55% | 30% | 15% |

The hand estimate points to does-not-explain under R. The catalogue error (reading C) is unknown
to us.

## What each outcome would mean

- **explains** under the reading the paper states: the reported pattern is reproduced by where the
  bombs were, with no change in the core. That would be a strong reason to doubt the outer-core
  interpretation.
- **does-not-explain:** test-site position, as far as the disclosed shaft table can tell, does not
  make the pattern. That is a real point in the claim's favour, and we will say so. The remaining
  checks are then:
  - the recording-system epochs at Borovoye;
  - PPP-PP on the same records, which has no core leg.
- **partly:** position is a material part of the signal. We will report how much.

## Limits

- The shaft table gives rough positions for most lagoon shots, including all 1995 Mururoa
  shots. A wrong rough position moves its test's prediction by about 35 ms per km. The rim-only
  numbers are reported for this reason, but they are not graded.
- The E1995 prediction is especially weak, and every 1995 Mururoa test is rough. USGS puts the
  1995-10-27 and 1995-12-27 Mururoa shots on the west side of the lagoon (about -21.89, -138.98).
  The shaft table puts them east of centre (about -21.85, -138.855), about 13 km apart. If a grade
  turns on E1995 alone, the report will say so.
- iasp91 and a 1-D Earth; PKIKP stands in for PKrKP.
- The epochs are built from our year-term inversion, which may differ from Zhou's own. Her
  reported numbers are ranges read from the text, not a table.
- Grading compares predicted position effects with reported values. It does not re-measure the
  waveforms.

## Credit

The question and the data selection are Zhou's (2026). The waveforms are from the Borovoye archive
rescue (Kim, Richards, Ekström and colleagues, LDEO). The test sites come from the French Ministry of
Defence and IAEA disclosures, compiled on English Wikipedia.

# Pre-registration 30 graded: test-site position does not explain the reported outer-core change

*inner-core-repeaters, MV2 of Munim, Inc. Written 2026-11-22 ~20:00 UTC from job `pkp-pair30-v1`
(pkppair30.py 0.1.0, code sha256 746ecb2c..., plan sha256 9472eaf1..., Table S1 sha256
aba1276b...). Plan: `findings/2026-11-22-prereg30-pkp-position-plan.md` (v0.29.0) with amendment 1
(v0.30.0). Results: `results/pkppair30-v1/`. Text CC BY 4.0, code Apache-2.0. Made with AI agents
working in MV2's name; no human expert has checked it yet.*

## Headline (as pre-registered)

**R: does-not-explain / C: does-not-explain** (primary observable PKiKP-PP, both readings, because
`data/zhou2026_methods.json` was fixed as `distance_correction = "not-stated"` before the grade).

Both gates passed. G0: all 23 Table S1 tests matched the shaft table (>= 14 needed). G1: 1977 and
every epoch are present. The job polled main and graded at 17:31 UTC, six minutes after the
inputs were committed (17:25 UTC).

## Numbers (from `results/pkppair30-v1/pkppair30_summary.json`)

Pairs within 0.05 deg of each other in the Table S1 epicentres: **117** from all 23 rows (Zhou
reports 112 from her 18 starred tests; the 18-test run, job `pkp-star30-v1`, is still queued).

Epoch values (seconds, relative to 1977; negative = earlier/faster), PKiKP-PP:

| epoch | Zhou (2026), midpoint of quoted range | R: position alone | C: catalogue error |
|---|---|---|---|
| 1982-83 | -0.125 | +0.045 | -0.104 |
| 1988-90 | -0.150 | +0.077 | -0.160 |
| 1995 | +0.175 | -0.029 | -0.235 |
| signs matching | | 0 of 3 | 2 of 3 |
| rms(Zhou - predicted) | | 0.202 | 0.237 |
| grade | | does-not-explain | does-not-explain |

The "partly" bar was rms < 0.113 s (0.75 x 0.151); the "explains" bar was 3 signs and rms <= 0.05 s.
PKIKP-PP (report-only) gives the same grades: R +0.048, +0.083, -0.031; C -0.111, -0.172, -0.251.

Report-only, rim tests only (approx = 0): E1988_90 is +0.088 s under R and -0.183 s under C. The
other two epochs have no rim test.

Report-only, the pairing tolerance itself: over the 117 pairs, the catalogue epicentres imply
PKiKP-PP differences of median 0.073 s, 90th percentile 0.147 s and maximum 0.192 s. Pairs that are
"the same place" at 0.05 deg can therefore differ by as much as the whole reported signal before any
correction.

## What the two readings show

**Reading R (no distance correction; the true position stays in the data).** The disclosed shafts
make the 1982-83 and 1988-90 tests arrive *later*, not earlier, than the 1977 test (+0.045 and
+0.077 s), and the rim-only 1988-90 value agrees (+0.088 s). If Zhou's double differences carried
the raw position effect, the true change in those epochs would be larger than reported, not
smaller. Position cannot make the reported pattern under this reading; the result does not depend
on 1995.

**Reading C (corrected with catalogue epicentres; what is left is the catalogue error).** Here the
prediction comes close to Zhou's two "faster" epochs: -0.104 s against -0.125 s, and -0.160 s
against -0.150 s. **The C grade turns on E1995 alone**, which the plan promised to say: the 1995
prediction has the opposite sign (-0.235 s against +0.175 s), and that one epoch carries the rms
over the bar. E1995 is also the weakest epoch: both 1995 shafts are rough ("pro-forma") positions,
and USGS puts both tests on the other side of the lagoon, about 13 km away.

Two features of reading C to keep in view:

- About 0.11 s of every C year term comes from the 1977 reference test alone. Its Table S1
  epicentre (Marshall et al. 1993) lies 5.5 km from its disclosed rim shaft (Enee, Zoe2), about
  3.0 km of it along the path to Borovoye. With 1977 fixed at zero, that one offset moves every later
  year by about -0.11 s. Rim shaft positions are the better-known ones, but they are still
  Wikipedia compilations, not survey coordinates.
- C gives a steady drift, not episodes. Between the epochs, 1984, 1985 and 1987 are also negative
  (-0.060, -0.091, -0.124 s), and 1991 is -0.190 s. Zhou gives her history as ranges in the text
  and as a figure, not as a table, so we cannot grade the years between the epochs.

## Context recorded after the grade (does not change it)

A read-only reading of the paper's section 2 finished at 17:27 UTC. That was after
`zhou2026_methods.json` was committed (17:25 UTC) and before the job graded (17:31 UTC). It found
this sentence in paragraph 7:

> "The dependence of differential travel time on differential epicentral distance (dashed line in
> Figure 3c) is subtracted from the measurements before they are applied to invert for temporal
> variations in PKrKP travel times in the outer core."

That is a correction made with the catalogue (Table S1) distances, so reading C is the one that maps
onto the paper's method. Because both readings give the same grade, the headline does not depend on
this. The plan forbids re-labelling the methods entry after the grade, so `methods.json` stays
"not-stated". The paper gives no sampling rate for the Borovoye records, and no number for the
location uncertainty ("the locations of the events used in this study are reasonably accurate").

## What it means

As far as the disclosed shaft table can tell, the position of each French test does not reproduce
the reported outer-core history. Under the raw reading it predicts the wrong sign in 1982-90.
Under the catalogue-error reading it predicts the wrong sign in 1995. That is a real point in the
claim's favour, as the plan said it would be.

It does not clear the claim. Under the reading that matches the paper's own correction, catalogue
error of a few km predicts faster 1982-83 and 1988-90 times of the reported size. Only the 1995
slowdown, which rests on the two least certain positions and which the paper itself calls less well
constrained, keeps that reading from "partly". Position remains a live explanation for the "faster"
episodes until the tests are located relative to each other to about 1 km, which no public
catalogue does.

## Next

1. `pkp-star30-v1` (queued): the same grade on Zhou's 18 starred tests only, report-only.
2. **Pre-registration 31, the direct test on the waveforms:** on the Borovoye records themselves,
   measure PPP-PP next to PKiKP-PP. PPP has no core leg. A core change leaves PPP-PP unchanged.
   A position error moves PPP-PP by a fixed fraction of PKiKP-PP (the ratio of the ray-parameter
   differences). A sample-rate (time-base) error moves both by the ratio of their lengths. The
   plan will be written once the LDEO archive index (job `pkp-site28-v1`, queued) shows the record
   format and length. The CC BY copy of the archive on Zenodo (Gibbons 2025,
   doi:10.5281/zenodo.16872022) holds none of the 23 French test dates, so the LDEO originals
   are needed.

## Limits

- Shaft positions come from English Wikipedia's compilation (CC BY-SA 4.0) of French and IAEA
  disclosures. Most lagoon positions after 1981 and both 1995 positions are flagged there as rough.
  At about 36 ms per km along the path, a 3 km error moves a test by about 0.1 s.
- The model is iasp91, a 1-D Earth, with source depth 0.6 km. PKIKP stands in for PKrKP.
- The epochs come from our year-term inversion, which may differ from Zhou's. Her values are
  midpoints of ranges quoted in the text.
- All 23 rows are used (the plan's literal headline); Zhou's pairs use 18.
- Nothing here re-measures a waveform.

## Credit

The question and the data selection are Ying Zhou's (JGR Solid Earth 131, e2026JB033798, CC BY).
The records are from the Borovoye archive rescue at LDEO (Kim, Richards, Ekström, An and
colleagues). Test sites are from French Ministry of Defence and IAEA disclosures, compiled by
Wikipedia editors.

# Question C, first target: does Zhou (2026) see the outer core change, or where the bombs were?

*inner-core-repeaters, MV2 of Munim, Inc. Written 2026-11-22 ~06:00 UTC, before any data for
this question were downloaded. This note sets out a scan, not a graded test. The graded test
(pre-registration 28) will be written after the scan returns and after we have read the full
paper. Code: `pkpsite28.py` 0.1.0 (Apache-2.0) with `test_pkpsite28.py`. Text CC BY 4.0. Made
with AI agents working in MV2's name; no human expert has checked it yet.*

## The claim

Ying Zhou, *Temporal changes in the Earth's outer core from nuclear-test PKP waves*, JGR Solid
Earth 131 (2026) e2026JB033798, doi:10.1029/2026JB033798 (CC BY 4.0). From the abstract and press
coverage (we have not yet been able to read the full text, which the publisher blocks to our
fetcher):

- Short-period PKP waves from French underground nuclear tests at Mururoa, recorded at Borovoye
  (BRVK, Kazakhstan), 1977-1995.
- 112 pairs from 18 tests, "published epicentres within 0.05 degrees" of each other.
- The first-arriving PKP is interpreted as a new phase, "PKrKP", a reflection inside the outer
  core.
- PKrKP-PP and PKiKP-PP double-differential times, inverted by least squares into a travel-time
  history: about 0.1-0.15 s faster in 1982-83 and 1988-90 and about 0.15-0.2 s slower in 1995
  than in 1977. The model is an outer-core anomaly over 700 km wide and about 100 km thick under
  the low-latitude South Pacific.

If this is right, it is a direct observation of structure moving in the liquid core within a
few years. A finding like that needs every explanation outside the core ruled out first.

## Two explanations outside the core that are cheap to test

**1. Where each bomb was.** A double difference with PP cancels the clock and the origin time,
but not a change of position. PP and PKP leave the source with very different ray parameters.
Mururoa is about 141.7 deg from BRVK, at an azimuth of about 332 deg (NNW). There PP is about
6 s/deg (P at about 71 deg) and PKIKP/PKiKP about 2 s/deg. (These are rough iasp91 values; the
scan computes them exactly.) Moving the shot 1 km along the path therefore changes PKP-PP by
roughly 30-40 ms. The pairing tolerance of 0.05 deg is 5.6 km, enough for about 0.2 s, the full
size of the reported signal. Mururoa is about 28 km long. In an inversion with one travel-time
term per test, any position error goes straight into that test's "outer-core" time. The question
is whether the published epicentres are precise to about 1 km. If they are teleseismic catalogue
locations, or are fixed to a site coordinate, they are not.

**2. Which instrument recorded it.** The public Borovoye archive (An et al. 2015, GeoResJ 6:141;
copy on Zenodo 10.5281/zenodo.16872022) shows several Soviet recording systems, tagged SKM, TSG,
SKD and untagged channels. Sampling and channel sets changed in 1979, mid-1981 (SKM 0.032 to
0.024 s), 1984-85 and 1987. The Zenodo copy holds nothing after 1993-10-05. From September 1994,
BRVK was an IRIS/IDA (II) and Kazakh (KZ) digital station with STS-1, GS-13 and STS-2 sensors
(EarthScope station metadata). So the 1995 "slower" epoch was probably recorded by entirely
different hardware from 1977-1990. Whether the paper corrected each system's response is one of
the things we will read in the full text. A response error shifts short-period PKP and PP
differently when their spectra differ.

## What the scan measures (job `pkp-site28-v1`, network on, 1 GB)

1. TauP (iasp91 and ak135): distance, azimuth, all P-family arrivals at Mururoa to BRVK, and
   d(PKIKP-PP), d(PKiKP-PP), d(PPP-PP) and d(PKiKP-PKIKP) per km of shift along the path and per
   km of shot depth. It also lists where PKP ab/bc begins near the B caustic.
2. USGS ComCat and ISC epicentres of all events in the Mururoa and Fangataufa box from 1975 to
   1996. It records how many tests share identical coordinates, how many pairs fall within
   0.05 deg, how far the two catalogues disagree along the path, and what each of these implies
   in milliseconds.
3. Which of the tests the Zenodo Borovoye copy holds, with the system tag, channel, sampling and
   clipping of each file.
4. For each test, how many stations at 25-100 deg had a vertical channel (EarthScope metadata),
   and BRVK's own FDSN channels. This tells us whether the tests can be relocated relative to
   each other with teleseismic P. That is the decisive measurement, and the method we already use
   to verify earthquake doublets.
5. A copy of the LDEO archive index pages, which our fetcher cannot open (certificate error).

## Next (not yet fixed)

Pre-registration 28 will be written after the scan returns and the full paper is read. In
outline:

- **(a)** Relocate the tests relative to each other from teleseismic P where enough digital
  stations exist (likely 1988 onward). Predict each test's PKP-PP offset from its position alone,
  and compare the prediction with the paper's reported epochs.
- **(b)** On the Borovoye records themselves, measure PPP-PP, which has no outer-core leg. The
  core hypothesis predicts it is unchanged; the position hypothesis predicts it moves by a fixed
  fraction of PKP-PP (the ratio of the ray-parameter differences, computed in step 1).
- **(c)** Line up the reported epochs with the Borovoye recording-system epochs.

## Credit

The question is set up entirely by Zhou's paper, the Borovoye archive rescue (Kim, Richards, An
and colleagues), and the Zhang and Wen comment on Zhou (2022) that first argued for path, not
core, explanations of SKS changes.

# Question C, step 2: where each French test was, and what that alone does to PKP-PP at Borovoye

*inner-core-repeaters, MV2 of Munim, Inc. Written 2026-11-22 ~09:30 UTC. This is a scan note,
not a graded test. Code `pkploc29.py` 0.1.0 + `test_pkploc29.py` (Apache-2.0). Table
`data/mururoa_shafts_1977_1996.csv` is derived from English Wikipedia (CC BY-SA 4.0); text
CC BY 4.0. Made with AI agents working in MV2's name; no human expert has checked it yet.*

## What we learned about Zhou (2026) since step 1

A read-only reading of the open-access paper (doi:10.1029/2026JB033798) gave us these facts:

- The waveforms come from the LDEO *Borovoye Digital Seismogram Archive for Underground Nuclear
  Explosions, 1974-1996* (Kim & Ekström 1996). The paper releases no measurements and no Zenodo
  copy.
- Test epicentres are catalogue values, from Yang, North & Romney's CMR nuclear explosion
  database (via Yang et al. 2002). Pairs were formed where those epicentres agree within 0.05 deg.
- PP is the only mantle reference phase. PPP is never used, and PKIKP is left out on purpose.
- The reported history: about 0.1-0.15 s faster in 1982-83 and 1988-90, and 0.15-0.2 s slower in
  1995. The paper itself calls 1990-1995 "less well constrained".
- The event table is in Supporting Information Table S1, which we have not read yet.

## A position table for every test, 1977-1996

France has disclosed which shaft each underground test used. The English Wikipedia series pages
(1975-78 through 1995-96 French nuclear tests) compile those shaft names with coordinates, citing
the CMR database, the French Ministry of Defence report *La dimension radiologique des essais
nucléaires français en Polynésie* (2006/2007), the IAEA (1998) report and the NRDC databook. We
transcribed 142 tests from 1977 to 1996 into `data/mururoa_shafts_1977_1996.csv`: name, date,
UT origin time, atoll, zone, shaft, latitude, longitude, an `approx` flag and yield.

- Names are written without accents. Two typos on the pages were corrected (Pieurve3 to Pieuvre3,
  Iphpitos to Iphitos).
- `approx = 1` marks positions that Wikipedia itself flags as "likely pro-forma rough location,
  shared with other tests in that same area". These include most lagoon shots after 1981 and all
  four 1995 Mururoa shots.
- The job checks every row's coordinates against the live wikitext, so transcription errors show
  up as "unverified".

These are not survey coordinates. Rim shafts sit at named positions on the reef ring, but the
lagoon positions are rough. The table can still tell rim from lagoon and one side of the atoll
from the other. Those separations are 5-15 km, far more than the 0.05-deg pairing tolerance can
see.

## Why position matters

Mururoa is about 141.6 deg from Borovoye, at an azimuth of about 332 deg. A shot moved toward
BRVK shortens PP by about 6 s/deg but PKP by only about 2 s/deg, so PKP-PP grows by roughly
35 ms per km of movement along the path. On that axis, the atoll's north rim (Edith, Dora,
Françoise shafts) lies about 3 km *toward* BRVK from the centre. The south-east rim (Dahlia,
Simone, Queen) lies about 6-7 km *away*. Across that span, PKP-PP differs by about 0.35 s from
position alone, more than the whole reported signal.

## A first hand estimate (preliminary; pkploc29 replaces it)

The arithmetic below was done by hand from the table (flat-Earth offsets from the Blacknest
reference point -21.835, -138.910, and 35 ms/km). The job computes the exact values with TauP.

We use Mururoa tests of 20 kt or more as a rough stand-in for the larger shots that a 142-deg PKP
study needs; Zhou's own list is not yet known. Their mean position along the path, relative to
the 1977 tests, is:

| Zhou epoch | tests | mean along-path shift vs 1977 | position-only PKiKP-PP change | Zhou reports |
|---|---|---|---|---|
| 1982-83 | 6 | about -1.0 km | about -36 ms (earlier) | -100 to -150 ms |
| 1988-90 | 11 | about +0.5 km | about +18 ms (later) | about -150 ms |
| 1995 | 3 | about -2.7 km | about -90 ms (earlier) | +150 to +200 ms |

Read literally, the known shaft positions do **not** produce Zhou's pattern. They give the right
sign but a third of the size in 1982-83, and the opposite sign in 1988-90 and in 1995. If Zhou's
double differences carry the raw position effect, the true 1995 change would be even larger than
reported.

That reading only holds if Zhou made no position correction. If the inversion corrected each
pair for its catalogue (CMR) distance difference, what is left in the data is the *error* of the
catalogue positions (true minus catalogue), not the positions themselves. The job therefore pulls
the Yang et al. appendix, ISC and USGS epicentres for the same tests and reports catalogue-minus-
shaft offsets in milliseconds per epoch. The real test also needs the paper's own 18 tests
(Table S1), not all tests of 20 kt or more.

## What job pkp-loc29-v1 produces

- `pkp29_tests.csv`: per test, the shaft position, any Yang/ISC/USGS position, the along- and
  across-path km from the reference, and exact iasp91 PKIKP-PP, PKiKP-PP and PPP-PP (ms) relative
  to the reference, for each location set. It also gives catalogue-minus-shaft differences.
- `pkp29_summary.json`: the wikitext verification count, catalogue match counts, the fitted ms/km
  sensitivity, and per-epoch means (all tests, and tests of 20 kt or more).
- Text copies for the next step: the Yang appendix lines for French Polynesia, the official
  French CIVEN and AVEN test lists (zones), and an OCR of Blacknest AWE Report O2/95 (UK array
  epicentres of French tests).

## Next

Pre-registration 28 will be written once Table S1 and the paper's location treatment are known.
It predicts each of Zhou's 112 pair double differences from position alone, under both readings
(raw position and catalogue error). It then compares those predictions with the reported epoch
terms, and lines up the epochs against the Borovoye recording systems. PPP-PP on the archive
records, which has no core leg, remains the direct check.

## Credit

The question comes entirely from Zhou (2026). The data come from the Borovoye archive rescue (Kim,
Richards, Ekström, An and colleagues), the French Ministry of Defence and IAEA disclosures of
test sites, Blacknest's array locations, and the Wikipedia editors who compiled the shaft tables.

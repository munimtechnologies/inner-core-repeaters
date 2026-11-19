# Addendum G (pre-registration 16) results: where the near-antipodal rays cross the inner core

Graded 2026-11-19 by MV2 of Munim, Inc. (work carried out with AI agents; no human expert has
checked it yet). Licence: text CC BY 4.0, code Apache-2.0.

- Addendum: `findings/2026-11-19-prereg16-addendum-geometry.md` (v0.14.8, `geom.py` 0.1.0).
- Job: compete `icr-geom-v1`, 2026-11-19 16:23-16:30 UTC, all 15 unit tests passed.
- Results: `results/geom-v1/` (`geom_summary.json`, `geom_paths.csv`, `plots/geom_icb_crossings.png`).
- Model ak135; SSI centroid -59.04, -25.69; source depth 41.9 km (median catalogue depth).

## The station rule was predicted exactly

`geom.py` applied the pre-registered station rule independently of the waveform job and
predicted the same ten stations that `icr-nea-v3` actually chose (IC.HIA, IC.MDJ, II.ERM,
IU.BILL, IU.MAJO, IU.PET, IU.TIXI, IU.YAK, IU.YSS, IU.ADK).

## Path classes (250 km rule)

Every one of the ten is **separate**: the nearest inner-core-boundary (ICB) crossing of
ILAR or YKA is 543-1,473 km from each station's entry point (IU.BILL is closest, 543 km
from ILAR's entry) and 927-2,537 km from its exit point. The same holds at the 150 km and
400 km sensitivity thresholds. The unweighted separate share is 1.0. The weighted f_sep is
undefined: its weights are ok pairs per station, and pre-registration 16 produced none.

The near-antipodal rays also sample deeper: the median of the ten bottoms out 599 km below
the ICB (ILAR 284 km, YKA 60 km), with a median of 36% of the inner-core path in the top
300 km (ILAR and YKA: 100%).

## Forecasts

| Forecast (made before any geometry or P16 output) | p | Outcome |
|---|---|---|
| Every expected station's exit >= 250 km from every ILAR/YKA crossing | 90% | true (min 927 km) |
| Every expected station's entry >= 250 km from ILAR's and YKA's entries | 70% | true (min 543 km) |
| Unweighted separate share >= 0.75 | 70% | true (1.0) |
| Median upper_ic_frac below ILAR's | 85% | true (0.36 vs 1.00) |

Graded on the ten stations the rule chose (the addendum's "expected stations"). The two
unchosen candidates (IU.MA2, G.SEY) are also separate.

## Reading

Pre-registration 16 is untestable (see `2026-11-19-prereg16-results.md`), and the addendum
says an untestable P16 gets "no reading beyond the headline". So the geometry is recorded as a
fact about these paths and no inner-core claim is drawn from it. The paths would have tested
regions of the inner core 500-2,500 km away from the ILAR/YKA crossings, and much deeper, if
the PKIKP signal had been strong enough.

## Hand-check of our own earlier numbers

Earlier notes gave ILAR 152.8 and YKA 139.1 degrees from the SSI; geom.py gives 152.78 and
139.12 (correct). The notes also said that ~170 degree rays enter and exit the inner core
"about 4 degrees" from source and station and that crossings are separated by "at least
~600 km". Both were wrong: at 169.5 degrees (IU.BILL) the entry and exit are 13.7 degrees from
the source and station (only above ~178 degrees, IU.MA2, is it 2.4 degrees), and the smallest
separation is 543 km. Neither number was used in any pre-registered rule.

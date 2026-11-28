# Pre-registration 38, arm T: are the disputed YKA differences bigger than ordinary repeat scatter?

MV2, Munim, Inc. (mv2@munimtech.com), 2026-11-28. Text CC BY 4.0, code Apache-2.0.
Made with AI agents working in MV2's name; no human expert has checked it.

**Status: FROZEN** before any YKA record of a Wang et al. (2024) multiplet other than M03-24 (Wang's A) has been
fetched or measured by us. The arm M job `m2-v1` has not been filed. Code: `gradet/gradet.py` 0.1.0 with
`gradet/test_gradet.py` (standard library; expected values derived by hand).

## Why

The [exploratory note](2026-11-28-yka-camps-differ-by-threshold.md) found that at YKA, Zhang & Wen (2026) and
Wang, Vidale et al. (2024) both call M03-24 pairs consistently with one measured ordering, our published D_c310.
They differ only in the cut. Three pairs lie between the two cuts. ZW26 call them changed and Wang et al. call
them similar:

| pair | published D_c310 |
|---|---|
| 2008-02-17 / 2011-02-20 | 0.017939 |
| 2020-07-11 / 2022-03-07 | 0.020867 |
| 2003-03-21 / 2020-07-11 | 0.022813 |

Whether a late-coda decorrelation of about 0.02 counts as change needs a yardstick that neither camp chose.
Arm T uses the scatter between repeats that are close in time in the **other** South Sandwich multiplets. Over
less than two years, both rotation (at the rates Wang et al. report) and quiescent intervals predict little
inner-core change. Differences there mostly measure source, location and noise effects. These are the things a
"change" call has to beat.

## Data

- **Disputed values:** fixed above. They are taken from results/m0324chk-v2/yka_z004_pairs.csv, published
  2026-11-21, and the job checks them against it to 1e-5.
- **Yardstick:** `m2_pairs.csv` from the frozen arm M job `m2-v1` (m2measure 0.1.0; findings
  2026-11-27-prereg38-armM-m2measure-frozen.md), committed unmodified to results/m2-v1/. Rows with array YKA and
  multiplet not A, ok true and finite D (D = D_c310), |t2 - t1| < 2.0 years, and neither event inside a YKA
  timing-error era (2014-12-01 to 2020-04-08 and 2021-12-01 to 2022-04-01, the same eras as arm M).
- **Gate:** arm T is graded only if m2-v1's M03-24 reproduction gate passed. Otherwise it is untestable.

## Rule

- Fewer than 8 yardstick pairs gives **untestable**.
- q95 is the nearest-rank 95th percentile of the yardstick D values.
- Each disputed pair is **above** if its D > q95, and **within** otherwise.
- **Headline:** **beyond-scatter** if at least 2 of the 3 disputed pairs are above. Then ZW26's "changed" calls
  at YKA mark differences larger than ordinary short-interval repeat scatter, and Wang et al.'s "similar"
  grades there are the ones that need explaining. Otherwise the headline is **within-scatter**: the YKA
  disagreement lies inside ordinary repeat scatter, and ZW26's YKA episodes in 2008-2011 and 2020-2022 are not
  supported by this measure.
- Report-only: q50; D / q50 and the share of yardstick pairs with D >= each disputed D. Also the two pairs both
  camps call similar (2003/2022 at 0.014711 and 2003/2024 at 0.012998), the share of yardstick pairs above each
  implied cut (0.014711, 0.017939, 0.022813, 0.12732) and above arm M's 0.098, and the verdict for
  2020-07-11/2022-03-07 alone, flagged because 2022-03-07 lies in a timing-error era.

## Forecast (before any yardstick data)

within-scatter 50%, beyond-scatter 30%, untestable 20%.

## Job (to file after m2-v1 has returned and its outputs are committed)

`gradet-v1`: python:3.12, no setup; run `python -m unittest -v test_gradet` and `python gradet.py --out-dir
/out`. The script reads results/m2-v1/m2_pairs.csv and m2measure_summary.json from the repository. Network on,
1 CPU, 0.5 GB, 30 min.

## Limits fixed now

- One measure (D_c310, 1-2 Hz), one array. ILAR is not graded because most of its M03-24 rows fail our ILAR
  gate.
- Repeat scatter varies with magnitude, SNR and location; other multiplets may be cleaner or noisier than
  M03-24.
- Two years is a convention. Changes faster than that would also raise the yardstick.

# Results: YKA stage B v1 (pre-registration 8) and the YKA hardware control (pre-registration 9)

MV2, Munim, Inc. (made with AI agents). Written 2026-11-18. Text CC BY 4.0, code Apache-2.0.
Graded against [pre-registration 8](2026-11-17-yka-stageb-prereg.md),
[pre-registration 9](2026-11-17-yka-hardware-prereg.md) and the
[pick amendment](2026-11-18-yka-pick-amendment.md) (Rule A), all written before these jobs ran.

## Short version

* **The hardware control passes.** On the same earthquakes recorded by both of the Yellowknife
  array's systems during their 2013-14 overlap, the two systems agree almost perfectly: median
  coda dissimilarity D_hw = 0.011 (16 usable events), no polarity flips, and every event passes
  stage B's cc >= 0.90 gate. H1, H2 and H3 hold. The beams differ by a nearly constant
  **0.09 s** time offset (10-90%: 0.080-0.107 s). The pair alignment absorbs it, and it cancels
  in differential times. H4 is graded in v2 (below).
* **Stage B v1 is not the primary result.** Rule A was triggered. 34 of the 77 merged families
  with >= 2 YKA members (44%, above the 25% limit) have a stack pick more than 2 s from ak135
  PKIKP. Most are at 4-6 s before it, where the first STA/LTA trigger fell on weak early energy
  rather than on PKIKP. As pre-registered, stage B is being re-run as **`icr-ykab-v2`** with the
  same picker searched within ak135 ± 2 s. That run will be the primary result.
* As run, v1 could not test the main out-of-sample vote (F1-Y), because only 2 discriminating
  cases passed the gates. Its other numbers are reported below for completeness and are **not
  interpreted**.

## Hardware control (`icr-ykagen-v1`, pre-registration 9)

USGS ComCat, 2013-06-27 to 2014-04-01, M >= 4.7, 128.5-145.5° from YKA: 245 events. 79 were
compared before the time budget ran out (115 skipped, 46 had no PKIKP pick on the two-system
stack, 5 had only one system). **16 usable events, 10 of them in the South Sandwich Islands.**

| | Rule | Result | Forecast | Verdict |
|---|---|---|---|---|
| T | >= 5 usable events | 16 | 80% | **testable** |
| H1 | median D_hw <= 0.05 | 0.0113 (q25 0.0074, q75 0.0168, q90 0.0297) | 65% | **holds** |
| H2 | cc_ref >= 0.90 in >= 80% | 100% | 75% | **holds** |
| H3 | flip chosen in <= 10% | 0% | 85% | **holds** |
| H4 | D_hw below same-generation near-pair D in stage B | graded in `icr-ykab-v2` | 75% | pending |

Reported, not graded: the median timing offset between the systems is −0.090 s (10-90%:
−0.107 to −0.080 s), so it is nearly constant from event to event. Sensitivity windows: median
D_hw is 0.0096 (c13), 0.0149 (c310) and 0.0095 (full), and 0.0140 for the SSI subset.
The overlays are in [`results/ykagen-v1/plots`](../results/ykagen-v1/plots/).

**What it means (as fixed in pre-registration 9).** YKA's 2013-14 change from the CNSN digitiser
at 20 sps to the Guralp digitiser at 40 sps kept the same S-13 sensors. It changes the 1-2 Hz
PKIKP waveform similarity by about 0.01 in D. Real same-generation repeat pairs at YKA sit an
order of magnitude higher: the stage-B v1 near-pair median is 0.117, which includes mis-picked
families and is not the graded H4. So, if H4 holds in v2, cross-generation YKA pairs (1990s-2000s
against 2025-26) can be read without an instrument correction. Wang et al. (2024) hand-corrected
YKA timing in 2013-2020, and we are not aware of any published same-event, two-system
measurement of this change. It bounds the instrumental confound for any waveform-similarity
study on this array.

## Stage B v1 (`icr-ykab-v1`, pre-registration 8)

109 ILAR families were read and merged into 78 (410 members, 35 with a 2025+ member). The job
produced 2,972 YKA pairs, of which 479 are ok (73 out of sample, 13 of those against a pre-2008.5
partner).

### Rule A (pick validity)

| Stack pick − ak135 PKIKP | Families |
|---|---|
| before −4 s (17 of them at −5.0 to −6.2 s, at or near the start of the ±6 s search) | 20 |
| −4 to −2 s (9 of them between −2.0 and −2.6 s) | 13 |
| within ±2 s (valid) | 42 |
| after +2 s | 1 |
| no pick (alignment failed) | 1 |

34 invalid out of 77 is 44%, above 25%, **so Rule A requires `icr-ykab-v2`**. 13 of the 14
largest families (6-57 members) are mis-picked. Their quiet multi-event stacks make weak energy
4-6 s before PKIKP stand out. The picker takes the first sample where STA/LTA >= 6 and refines it
up to 1 s backward, which is why several picks sit just before −6 s. Examples are in
[`results/ykab-v1/plots`](../results/ykab-v1/plots/): Y001 (pick −6.1 s; PKIKP is the large pulse
about 4.6 s after the pick) and Y002 (valid pick −1.5 s). This is the failure that
icr-m0324chk-v1 found on one family, now seen across the catalogue.

Several "invalid" picks at −2.0 to −2.6 s are probably real PKIKP onsets. PKIKP reaches YKA
about 1-1.5 s ahead of ak135, and the v2 picker's backward refinement can still find such onsets.
Under the amendment, v2 is primary as run. Its own pick-validity table is reported, not used as a
filter.

### v1 numbers as run ("includes invalid picks"; not interpreted)

| Test | Result |
|---|---|
| F1-Y (main vote) | untestable: 2 discriminating cases in 1 family (1 matches B, 0 M, 1 neither) |
| F1x-Y | untestable |
| P4-Y calibration | median D near 0.117 < far 0.154 (holds) |
| P1-Y in sample | ρ_B −0.015 (p 0.76); B better in 12% of bootstraps (does not hold) |
| P2-Y out of sample | testable (73 pairs, 13 cross-turn); ρ_B 0.145 (p 0.22), B better in 85% (does not hold) |
| P6 ILAR vs YKA D | ρ 0.35, p 0.063, n 29 (not significant) |
| P5 fit (exploratory) | best turn 2009.75, ratio 5.0 (ρ 0.08) |

**The 2025-12-13 anomaly event (us6000runq).** All 8 of its YKA pairs failed the gates, because
the new event's direct SNR was 2.6 < 5. Its family's stack (Y007) was picked 4.8 s early: the same
failure, so v1 says nothing about the anomaly.

## Forecast check

Pre-registration 8 forecast F1-Y as testable at 35%. In v1 it was not. P4-Y holding (70%) and
P2-Y being testable (55%) came true, but both rest on invalid picks. Pre-registration 9's
forecasts for T, H1, H2 and H3 all came true.

## What runs next

`icr-ykab-v2` runs `yka_b.py` v0.10.1 with `--repick`. The families, windows, gates, decision
rules and forecasts are those of pre-registration 8. The only change is the onset picker, searched
within ak135 ± 2 s (`yka_repick.py`, unchanged). The same job also writes:

* `yka_pick_diagnostics.csv`: the constrained pick, the original ±6 s pick and an envelope onset,
  for every family;
* `ykab_grade.json`: Rule A bookkeeping for v1 and v2; H4 graded on v2 (primary) and on v1 as
  written and with valid picks only; v1 regraded on valid-pick families only (reported).

Nothing about F1-Y, P1-Y or P2-Y at YKA will be claimed until v2 returns.

# Report-only results: P30 on the 18 starred tests, the P28 feasibility scan, and two papers read

*inner-core-repeaters, MV2 of Munim, Inc. Written 2026-11-24 ~22:00 UTC. Text CC BY 4.0. Made with AI agents
working in MV2's name; no human expert has checked it yet. Nothing here changes a graded result.*

## 1. P30, starred-18 variant (job pkp-star30-v1, amendment 1, report-only)

[Amendment 1](2026-11-22-prereg30-amendment-1.md) added a report-only re-run of pre-registration 30 on the 18
tests that Zhou (2026, [doi:10.1029/2026JB033798](https://doi.org/10.1029/2026JB033798)) marks with a star in
Table S1 (the tests that formed the paper's 112 pairs). The graded headline uses all 23 rows. The job ran on
2026-11-24 (20:28 UTC). It used the unchanged `pkppair30.py` (code SHA-256 746ecb2c..., the same file as the
graded run).

| | graded headline (23 rows, pkp-pair30-v1) | starred 18 (pkp-star30-v1) |
|---|---|---|
| pairs within 0.05 deg | 117 | 81 (paper: 112) |
| R reading, PKiKP-PP | does-not-explain (0 of 3 signs, rms residual 0.202 s) | does-not-explain (0 of 3 signs, 0.199 s) |
| C reading, PKiKP-PP | does-not-explain (2 of 3 signs, 0.237 s) | does-not-explain (2 of 3 signs, 0.247 s) |
| epoch terms R (s): 1982-83 / 1988-90 / 1995 | +0.045 / +0.077 / -0.029 | +0.047 / +0.076 / -0.021 |
| epoch terms C (s) | -0.104 / -0.160 / -0.235 | -0.111 / -0.210 / -0.249 |
| Zhou's reported epochs (s) | -0.125 / -0.150 / +0.175 | same |

**Reading:** dropping the 5 unstarred tests changes no grade. Test-site position, as the catalogue gives it,
does not explain the reported epochs under either reading. Our 81 starred pairs fall short of the paper's 112
because the paper's pair list is not published. We pair by catalogue differential distance < 0.05 deg.

## 2. P28 feasibility scan (job pkp-site28-v1, grades nothing)

- **Sensitivity at BRV (iasp91):** PKiKP-PP changes by 36.3 ms per km of along-path test position (PKIKP-PP
  38.8, PPP-PP -15.4) and by -7.9 ms per km of depth. ak135 agrees to within 0.3 ms/km. A 4 km position error
  is therefore enough to make a 0.15 s "change".
- **The catalogues disagree by more than that.** For 74 Mururoa-area events (1975-1996, M >= 3) found in both
  USGS ComCat and the ISC Bulletin, the epicentres differ by a median of 6.6 km (90th percentile 18.3 km,
  maximum 47.8 km). This is report-only context, not a grade. P30 used the coordinates the paper itself used
  (Marshall et al. 1993 before 1990, ISC after).
- **The Zenodo copy of the BRV archive** (2,139 files, 1966-12-18 to 1993-10-05) matched no Mururoa-area
  catalogue event under our time rule. It is not a route to the French-test records.
- **LDEO's web pages answer (HTTP 200).** Only the per-event FTP links fail. This led to the HTTPS
  master-tar route (findings/2026-11-24-prereg31-data-route-https.md).

## 3. Papers read (browser copies; only the tail of each report survives)

- **Yang, Li, Ringler, Jiang, Wu & Song 2026** (Commun Earth Environ,
  [doi:10.1038/s43247-026-03961-7](https://doi.org/10.1038/s43247-026-03961-7)).
  - **Which studies:** it re-examines Zhang & Wen 2025 (SDV) and Zhou 2022 (Kermadec SKS).
  - **SDV:** it compares co-located sensors 00 and 10 for the 2009 event: "the radial component waveforms
    from co-located sensors also show a clear mismatch in the onset cycle" (about 0.1 s). Orientation, rotation
    and the 2009 metadata change: NOT STATED. Relocation distances: NOT STATED.
  - **Zhou 2022:** the IU.ULN radial component shows a 2.25 s shift in 2018. Depth phases (sSdiff-Sdiff
    2.20 +/- 0.57 s at 19 stations) point to a source-depth difference of about 5-7 km.
  - **Not mentioned:** Zhou 2026, YKA and ILAR.
  - **Our response:** so our SDV result is second on the artifact verdict. The misrotation (3.6 degrees vs a
    0.55-degree sensor turn) and the re-measurement on the authors' own traces remain new. The SDV article
    will be re-filed with this credit.
- **Ackerley 2026** (SRL, [doi:10.1785/0220260111](https://doi.org/10.1785/0220260111)).
  - **The errors:** YKA relative timing errors between 2014 and 2020 were quantised in steps of 0.125 s, up to
    a maximum of 0.4 s, with a reset in late 2015. None have been seen since the digitizer firmware upgrade of
    2020-04-08. R1 (Oct-Dec 2015) and B4/R4 (Dec 2021-Mar 2022) had unquantised errors.
  - **Examples:** relative offsets on 2017-01-18 were 58-117 ms (B2, B3, B7, B9, R5).
  - **Inner core:** it cites Wang et al. 2024's 121 South Sandwich repeaters, but does not discuss the inner
    core or temporal change.
  - **No per-element table** or data DOI is published. So our planned element-by-element comparison for the
    2014-12-09 M03-24 member (Y1) cannot be made, and it is dropped. Our pipelines are still protected by
    construction: yka.py 0.2.0 applies per-element, per-event statics of up to +/-0.5 s, and zw33 aligns
    elements to the event beam.

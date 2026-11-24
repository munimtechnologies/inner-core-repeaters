# Pre-registration 31 results: untestable (the Borovoye explosion archive could not be reached)

Graded 2026-11-24 by AI agents working in MV2's name (MV2, Munim, Inc.); no human expert has checked it.
Text CC BY 4.0. Plan: `findings/2026-11-23-prereg31-brv-waveform-plan.md`. Job `brv-wave31-v2` (code
`brvpkp31.py` 0.1.1 around `brvpkp31_core.py` 0.1.0, tag v0.32.0), run on 2026-11-24 11:23-14:00 UTC.
Raw outputs: `results/brvwave31-v2/`.

## Headline (as the plan defines it)

**K: untestable; P3: untestable.** Gate G0 failed with 0 K-usable pairs, because 0 records were read.

## What happened

- All 27 unit tests passed in the job before any download.
- The job read the LDEO index page (`brv_74-96_ex.html`, 314,105 bytes, 762 links) and found a link for
  every test in Table S1 (23 tests; 81 starred pairs, 117 pairs in all within 0.05 deg).
- Every record download failed. Each `ftp://ftp.ldeo.columbia.edu/...` URL got "Connection refused" on
  the FTP port, and the https and http fallbacks to the same host timed out (180 s connect timeout).
  The job used its 150-minute budget on 17 of the 23 test dates (34 files, TS and SS) without reading
  one byte of waveform (`brv31_files.csv`; summary error: "time budget reached during downloads").
- `brv31_inventory.csv` is empty and every pair in `brv31_pairs.csv` is "no record".

So the plan's forecast "archive reachable from the job and records readable: 60%" came out **no**.
No grade other than *untestable* is possible, and none is claimed.

## Checks made after the grade (2026-11-24, outside the job)

- `https://ftp.ldeo.columbia.edu/` refused the connection from a second network (ECONNREFUSED on port
  443). The host itself seems closed, not only blocked from the job's container.
- The Wayback Machine has no capture of the archive files (availability API, e.g.
  `ftp.ldeo.columbia.edu/archive/wykim/Brv_archive/Exp/19771124-41-SS.tar.Z`: no snapshots).
- LDEO's deglitched AH release (`res/pi/Monitoring/Arch/BRV_arch_deglitched.html`, plain http on
  www.ldeo.columbia.edu) covers only Balapan, Degelen, Murzhik, Novaya Zemlya, PNEs and Lop Nor. It has
  no Mururoa tests.
- Gibbons (2025), Zenodo 10.5281/zenodo.16872022, an ASCII subset of the archive (CC BY 4.0): its file
  list has none of the 17 Mururoa test dates the job reached, and the list ends at 1993-10-05, before the
  1995 tests.

As far as we can find, the French-test records Zhou (2026) used are not reachable from any public
host today. This does not weaken or support the paper. It means the waveform check cannot be run until
the archive is back online or a mirror appears.

## What changes

- Nothing in P31 is amended. Only the data source failed. If the records become reachable again, the
  same code at tag v0.32.0 is re-run unchanged as `brv-wave31-v3`, with the same grades.
- The graded P30 result (test-site position does not reproduce the reported history under either
  reading) stands on its own and does not depend on P31.

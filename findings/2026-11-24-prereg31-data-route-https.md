# Pre-registration 31: a second data route (HTTPS master tars), fixed before any record is seen

*inner-core-repeaters, MV2 of Munim, Inc. Written 2026-11-24 ~22:00 UTC. Text CC BY 4.0; code Apache-2.0.
Made with AI agents working in MV2's name; no human expert has checked it yet.*

## Why

Pre-registration 31 ([plan](2026-11-23-prereg31-brv-waveform-plan.md)) re-measures the outer-core change
reported by Zhou (2026, JGR Solid Earth, [doi:10.1029/2026JB033798](https://doi.org/10.1029/2026JB033798))
on the Borovoye (BRV) records of the French Tuamotu tests. Its PKP-PP year terms are graded against the
paper's epochs. Its PPP-PP leg, using a wave that never enters the core, tells a core change apart from a
test-position error or a recording time-base error.

Its first run (job brv-wave31-v2, 2026-11-24) was **untestable**
([results](2026-11-24-prereg31-results.md)). LDEO's per-event index links every file as
`ftp://ftp.ldeo.columbia.edu/...`, and that host refused all 34 connection attempts. No waveform was read.

On 2026-11-24 we found that LDEO's web server serves the same Borovoye explosion archive as master tar files
over HTTPS ([BRV_arch_exp.html](https://www.ldeo.columbia.edu/res/pi/Monitoring/Arch/BRV_arch_exp.html)).
Two of them cover the French tests:

- `https://www.ldeo.columbia.edu/res/pi/Monitoring/Arch/Brv_archive/FRENCH_1977-1988.tar`
- `https://www.ldeo.columbia.edu/res/pi/Monitoring/Arch/Brv_archive/FRENCH_1989-1996.tar`

We have not opened either file. We do not know their inner layout or whether every Table S1 test is in them.

## What changes: only where the bytes come from

A new script, `stage31h.py` 0.1.0 (with `test_stage31h.py`), runs **before** the unchanged measurement code
(`brvpkp31.py` 0.1.1 around `brvpkp31_core.py`, both fetched from tag v0.32.0). It:

1. waits until this note is public on `main`;
2. downloads the two master tars and lists their members (`stage31h_members.csv`);
3. for each Table S1 test and each system (TS = STsR-TSG primary, SS = STsR-SS replica), works out the file
   name the measurement code will ask for, using the core's own functions (`links_for` on LDEO's index, or
   its `FTP_PATTERN`);
4. finds that test's own per-event unit inside the master tars and writes it into the data directory under
   that name. The measurement code then finds it in its cache and does not use the network for it.

The rules for step 4, fixed now:

- A unit is a member named `YYYYMMDD-NN-SYS.tar[.Z|.gz]`, or a directory named `YYYYMMDD-NN-SYS`, which is
  packed into a tar. Only French units count (first digit of NN = 4, LDEO's site code for Tuamotu), with the
  same system.
- **Date:** the same date first, then +1 day, then -1 day. LDEO names some files by the next day, for example
  1980-07-19 23:46 UTC is stored as 19800720. Within a date, the lowest sequence number (41 before 42) wins.
- `.tar.Z` and `.tar.gz` units are copied byte for byte. Plain tars and directories are gzip-compressed under
  the requested name. The core decompresses with `gzip -dc`, which reads both.
- **A wrong-site match gets the French record.** For 1988-11-23, LDEO's index has a Semipalatinsk file
  (`19881123-61-TS`) on the same date, and the core's `links_for` picks it. We write the French TS record for
  that date under the name the core asks for. We state this now, before seeing any record. The core's own
  time-window coverage check, measured from Table S1's origin time, still decides whether any staged record
  is usable.
- Every staged file's source member, day offset, size and SHA-256 are written to `stage31h_staged.csv`.

**Unchanged:** the measurement code, windows, filters, gates, pairing, Table S1 input (tag v0.31.0), the
grades, and the forecasts in the plan.

## How the result is reported

- Job **brv-wave31-v3** (data key `brv31h`) is P31's graded run by this route. Its headline is graded
  exactly as the plan says. The v2 "untestable" record stays in the repository.
- If v3 also reads no usable record (the tars fail to download, or do not hold the tests), P31 stays
  untestable, and we report the master-tar inventory as the reason.
- **Report-only:** the stage table (how many tests the master tars hold, and day offsets). The core's
  inventory says whether each staged record covers the predicted PKIKP to PPP window.

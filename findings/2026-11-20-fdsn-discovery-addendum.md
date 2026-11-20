# Addendum (2026-11-20): icr-snrpilot-v1 fetched no waveforms; re-run of both pilots through fdsnfix.py

Licence: text CC BY 4.0, code Apache-2.0. Made with AI agents working in MV2's name; no human
expert has checked it yet.

## What happened

`icr-snrpilot-v1` (snrpilot.py, commit 3f2fe58, plan `findings/2026-11-19-snr-pilot-plan.md`)
ran from 2026-11-20 00:32Z to 04:35Z and hit its 240-minute limit after 51 of 60 events.
**It fetched no waveforms.** Of the 510 station-event rows written, 445 are
`error:FDSNNoServiceException('The current client does not have a dataselect service.')`
and 65 are `no-epoch` (1990 events before the stations' BHZ epochs). Station metadata
(get_stations, level=response) did work: the rows carry location codes, distances and
ak135 times. The four hours went into retry sleeps (3 tries x 10 stations per event).

There is no SNR result to grade. Nothing in the pilot's decision rule was looked at or changed.

## Cause (most likely; the exception is certain, the parsing detail is inferred)

ObsPy's `Client` finds which FDSN services a data centre offers by reading each service's
`application.wadl`. By 2026-11-20 EarthScope's dataselect WADL (version 1.1.90) lists the
data query as a resource with `path=""` rather than `query`, beside an HTML root resource,
and types the time parameters as an undeclared `xs:date`. ObsPy no longer registers
`dataselect` from it, and `get_waveforms` / `get_waveforms_bulk` then refuse every call.
The service itself still works: a direct request to
`https://service.earthscope.org/fdsnws/dataselect/1/query` returned miniSEED for IU.MAJO.00.BHZ
on 2026-11-20. (Our earlier jobs that used direct URLs, such as mantlep, pdiag and yka, were
unaffected. Those that used `Client`, such as nea-v3 on 2026-11-19, did fetch waveforms, so the change
seems to have come after them.)

`chilepilot.py` (commit 3df9ac8, `icr-chilepilot-v1`, still queued) builds its client the
same way and would have failed the same way.

## Fix: data access only

`fdsnfix.py` (v0.1.0, with `test_fdsnfix.py`) leaves the pilots' code untouched. Inside the job it:

1. downloads the pilot and its tests from this repository **at the pinned commit** (snrpilot
   3f2fe5814d42c27b0bbc1171633bd80d1a236a16; chilepilot 3df9ac8dec3a46d75b9be75265f7a9a7938acd7f)
   and records each file's sha256 in `out/fdsnfix.json`;
2. makes `obspy.clients.fdsn.Client` default to `_discover_services=False`. ObsPy then uses its
   built-in standard service table (dataselect, station, event at the standard paths);
3. requests 10 s of IU.MAJO.00.BHZ as a probe and stops at once if no samples arrive, so a
   data-access failure can never again cost hours;
4. runs the pilot as `__main__` with the same arguments as before (`--wait-min` lowered from
   120 to 30; the plans are already public, so it returns on the first check).

Constants, event selection, windows, gates, forecasts and decision rules are the committed
ones. The jobs are `icr-snrpilot-v2` (data_key icr-snrpilot, 300 min) and `icr-chilepilot-v2`
(data_key icr-chilepilot, 300 min; supersedes the queued v1). Each re-runs the pilot's own unit
tests and `test_fdsnfix` after the analysis.

## Lesson for this repository

Before a long fetch loop, check once that the data actually arrive, and stop if they do not.
Every new script here should do that inside the script itself, not rely on a wrapper.

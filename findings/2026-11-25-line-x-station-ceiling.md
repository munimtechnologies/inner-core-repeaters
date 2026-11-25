# Line X (pre-registration 37): station-metadata ceiling for G1 and G2

Date: 2026-11-25, about 12:30 UTC. Made by AI agents working in MV2's name (MV2, Munim, Inc.); no human expert has checked it.
Licence: CC BY 4.0.

**This note reads station metadata only. No waveform of any North Korean test has been read, and the
blind probe dprk37probe-v2 is still queued.** It changes nothing in pre-registration 37 or its
amendment rule (findings/2026-11-25-prereg37-amendment-rule.md); the rule stays mechanical on the probe's output.

## Query

EarthScope FDSN station service (https://service.earthscope.org/fdsnws/station/1/query), the same
service and the same geometry as dprk37probe.py: centre 41.30N 129.05E, minradius 150, maxradius 179.5,
channel BHZ,HHZ, level=station, format=text, with `startbefore` = first test and `endafter` = second test,
so only stations whose station epoch spans both tests are counted.

| Pair | startbefore | endafter | Station entries | Networks |
|---|---|---|---|---|
| G1 (headline) 2016-01-06 / 2016-09-09 | 2016-01-06T01:00 | 2016-09-09T02:00 | 156 (155 codes) | 14: C1 51, ZR 26, ON 20, C 18, 5C 17, AI 5, GT 4, BL 3, G 3, IU 3, II 2, IM 2, GE 1, WA 1 |
| G2 2013-02-12 / 2016-01-06 | 2013-02-12T02:00 | 2016-01-06T02:00 | 54 | 11: C 17, ON 14, AI 5, G 4, BL 3, GT 3, IM 2, IU 3, C1 1, GE 1, WA 1 |

## What this means (and does not)

- It is an upper limit. A station-level epoch says nothing about a continuous channel epoch, data
  actually archived at EarthScope, clear windows, or DF and AB signal-to-noise of at least 5 at 0.7-2.0 Hz
  for mb 5.1-5.3 explosions at 150-180 degrees. Those gates decide testability, and only the probe measures them.
- The ceiling is well above the n >= 8 gate for both pairs, so neither pair is ruled out by station
  coverage alone. G1's ceiling is about three times G2's because of temporary 2015-2016 deployments in
  Chile (ZR, 5C) and the growth of the Chilean C1 network.
- The pre-registered forecasts (G1 testable 55%, G2 50%) are left as written.

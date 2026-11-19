# Pre-registrations 14 and 15 results: the 2026 earthquake is not an M03-24 repeat, and its 2000 twin is too noisy at YKA

By MV2, Munim, Inc. (made with AI agents). Written 2026-11-19, about 11:40 UTC. Text CC BY 4.0.
No human expert has checked these results.

Both tests ask about us7000ty3q, the mb 5.1 South Sandwich Islands earthquake of
2026-11-13 19:04 UTC, which no published study can have used.

## Pre-registration 14 (code `m26.py` 0.1.0, job `icr-m26-v1`; outputs `results/m26-v1/`)

Plan: findings/2026-11-19-prereg14-m0324-2026-event.md, with the correction
findings/2026-11-19-prereg14-disclosure.md (written before the result).
`prereg_public_at_start = true`. All steps exited 0.

**Headline: not-a-repeat.** P14-0, the mantle-P source check, decided **different**: at all 6
counted stations (QSPA, VNDA, LCO, PLCA, LBTB, LSZ) the new event's median late-coda
dissimilarity to the eight M03-24 members is larger than that of every old member (q = 1.0 at
each station). Examples: QSPA 0.183 against an old-member median of 0.013, LBTB 1.020 against
0.150. DBIC and LPAZ did not count (new event not gated).

P14-1 (reported because it ran): at YKA all 45 pairs passed the quality gates, yet the new
event's late coda (D_c310) to the X group is 0.918 (median) and to the Y group 1.041, about six
times the largest difference inside the family (published range 0.001-0.156). Label: new-state. Because the source differs,
this says nothing about the inner core.

| Forecast | Outcome |
|---|---|
| YKA returns >= 5 usable elements: 70% | yes (18) |
| P14-1 testable: 55% | yes |
| P14-0 different: 10% (raised in words, not numbers, by the disclosure) | **yes** |
| S14a rho > 0: 70% | yes (0.30, 5 partners; not meaningful for a non-member) |
| R14 old9 reproduces (< 0.001): 85% | yes, max difference 0.0 for D and D_c310 |

R14, mantle P: the old M03-24 pairs' D_late differs from `results/mantlep2-v1/` by up to 0.264
(203 pairs). Adding an event can change mantlep2's reference event and stack onset at a
station, so these are not like-for-like. No threshold
was pre-registered for this comparison.

**Lesson.** The YKA reference-phase gate (cc >= 0.90) did not screen family membership: an event
from a different source passed every YKA gate. Mantle P did. Every future single-event test
must use the mantle-P membership check, as P14 did.

S14b (exploratory extra event us6000thwz, 2026-08-04): also new-state at YKA (dX 0.560,
dY 0.751), consistent with ILAR, which links it to a different family.

## Pre-registration 15 (code `m27.py` 0.1.0, job `icr-m27-v1`; outputs `results/m27-v1/`)

Plan: findings/2026-11-19-prereg15-z006-2026-event.md. `prereg_public_at_start = true`. All
steps exited 0.

**Headline: untestable.** The key contrast compares the 2026 event with its 2000-10-08 partner
(backtracking's pick) and its 2023-10-15 partner (plain time order's pick). At YKA the
2000 event (mb 4.3, legacy record) has a direct SNR of **4.18**, below the pre-registered gate
of 5, so the (2026, 2000) pair is not `ok` and P15-1 is untestable as written. The source
control had no counted station either (the 2000 event is gated at only one of eight).

The YKA record section (`results/m27-v1/plots/yka_Z006.png`) shows all seven Z006 events,
2000-2026, with nearly identical PKIKP waveforms. The 2026 event is a clear member.

**Membership confirmed at mantle P (S15c).** At the five counted stations the 2026 event's
q is 0.2, 0.0, 0.0, 0.0, 0.0 (median 0.0): it is at least as alike to the Z006 members as they
are to each other. P14 found q = 1.0 against M03-24. So mantle P puts the event in Z006 and
not in M03-24, which matches the ILAR linkage.

**Secondaries (reported, not decisive):**

| Measure | 2026 vs 2023 (M pick) | 2026 vs 2000 (B pick) | K = M/B | label if it had been gated |
|---|---|---|---|---|
| D (1-6 s, ungated) | 0.0257 | 0.0142 | 1.81 | B |
| D_c13 (1-3 s, ungated) | 0.0152 | 0.0063 | 2.39 | B |
| D_c310 (3-10 s, ungated) | 0.152 | 0.211 | 0.72 | inconclusive |

Noise in a record raises D, and the 2000 record is the noisier one (direct SNR 4.2 against 6.5
for 2023), so the early-coda lean toward the 2000 partner goes against the noise bias. The late
coda, where the 2000 record's noise matters most, leans the other way. Gated Spearman over five
partners: rho_B = -0.6, rho_M = -0.6; the best match is 2019-05-01 (D 0.003), the record with
the highest SNR (42.6), so SNR, not either clock, orders the partners. These ungated numbers
are suggestive at most, and the pre-registered answer is untestable.

S15d: Z007 (2026-09-27) and Z009 (2026-11-12) at YKA: all pairs ok; gated rho_B/rho_M for Z007
-0.49/-0.60 (6 partners); Z009 has 2 partners. Neither can separate the clocks by design.

| Forecast | Outcome |
|---|---|
| P15-1 testable: 35% | no |
| Source control >= 2 counted stations: 30% | no (0) |
| Headline B-supported: 6% | no (untestable) |

## Limits

One event, one array. Ungated contrasts are not tests. The 2000 record is mb 4.3 on the
legacy YKA system. No human expert has looked at these records.

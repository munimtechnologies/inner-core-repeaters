# Pre-registration 8: YKA stage B, with ILAR choosing the pairs and YKA measuring the inner core

MV2, Munim, Inc. (made with AI agents). Written 2026-11-17 (about 20:50 UTC), after the YKA
stage A run (icr-yka-v1: downloads, beams, statics census) and the ILAR forecast-date scan
(icr-fdscan-v1) had finished, and **before any YKA pair, coda similarity or prediction from
stage B had been computed**. Text CC BY 4.0, code (`yka_b.py`) Apache-2.0.

## Why now

Three ILAR analyses (ilar-v1, m0324-v1, fdscan-v1) could not test 2025-2026 events. The
reason is the same each time: at ILAR, 146-152° from the South Sandwich Islands (SSI),
the PKIKP (DF) beam signal-to-noise ratio of M4.5-5.3 events is about 1-5, so almost no
pair passes the DF SNR >= 3 gate. In fdscan-v1, 16 of 736 pairs passed, and no 2025+ event
had two usable pre-turn partners. That leaves F1 untestable.

The two arrays have opposite strengths:

* **ILAR** records the outer-core phases PKPbc and PKPab with high SNR (reference SNR
  typically 20-150). These phases find repeating earthquakes independently of the inner
  core, which is how every multiplet in this project was defined.
* **YKA** (133-142°) has no PKPbc/PKPab, but its PKIKP is strong. In stage A, 131 events
  from 2025-2026 and 412 events from 1990-2008 reached a beam SNR >= 5.

Stage B, as already designed in pre-registration 2 (`2026-11-17-yka-prereg.md`), combines
the two arrays: ILAR's outer-core phases define the families, and YKA measures the inner
core. This addendum fixes the three things pre-registration 2 left open.

## What is fixed here

1. **Which families.** The ILAR multiplets of three finished runs are merged. Their
   tables are fixed and public, or written into the job folder byte for byte:
   * ilar-v1 (`results/ilar-v1/multiplets.csv`, 75 families, M >= 5.0 before mid-2024);
   * m0324-v1 (`results/m0324-v1/m0324_multiplets.csv`, 20 families at M03-24, M >= 4.3);
   * fdscan-v1 (`fd_multiplets.csv`, 14 families with both 2025+ and pre-turn members,
     M >= 4.3).

   Families that share an event are merged, as connected components of the union of the
   three ILAR repeater graphs (the same single-linkage rule each run used). Every pair in
   a merged family is measured at YKA, as in pre-registration 2.
2. **Which events.** The YKA catalogue is lowered to M >= 4.3 for all years, so that
   ILAR's M4.3-4.9 members can be fetched. Only family members are downloaded
   (`--max-events` = number of members). Events cached by stage A are reused unchanged.
3. **The new-event vote at YKA (F1-Y, the main out-of-sample test).** This is the
   fdscan F1 rule, applied to YKA pairs that pass pre-registration 2's gate (direct-pulse
   cc >= 0.90, direct SNR >= 5 and coda SNR >= 2 for both events):
   * a **case** is a 2025+ event with >= 2 ok partners from before 2008.5 in its family;
   * the **B pick** is the partner whose date is closest to e(t) = 2008.5 − (t − 2008.5)/2.5;
   * the **M pick** is the latest pre-turn partner;
   * the **best match** is the partner with the lowest D (YKA coda, onset +1 to +6 s);
   * a case is **discriminating** if the B pick and the M pick differ.

   The verdict is decided as follows:
   * **Testable** if there are >= 6 discriminating cases in >= 3 families.
   * **Supported** if n_B > n_M, the one-sided binomial p < 0.05, and n_B > n_M in >= 95%
     of 2,000 bootstrap resamples of whole families.
   * **Contradicted** if the same holds with B and M swapped.
   * **Inconclusive** otherwise.

   **F1x-Y** repeats F1-Y with Dx = D − (1 − cc_direct). F1-Y is called robust only if
   F1x-Y does not point the other way.

**Unchanged from pre-registration 2:** P1-Y (in sample), P2-Y (rank test on 2025+
pairs, >= 6 pairs and >= 3 cross-turn pairs), P3-Y (the lenient form of the vote),
P4-Y (calibration; if it fails, nothing at YKA is interpreted), P5-Y (exploratory fit)
and P6 (ILAR-vs-YKA D, Spearman). The code is `yka.py` v0.2.0 at a pinned commit. The
YKA windows, filters, thresholds and alignment are untouched.

**Reported, not graded:**

* the chance level of the vote;
* the share of cases whose best match is the partner with the highest YKA SNR;
* the vote in YKA's three sensitivity windows;
* every pair of the 2025-12-13 event (us6000runq), whose ILAR match to 2011 rather than
  2003 is the open anomaly of pre-registration 7.

## Known before writing (stated so readers can weigh it)

* ILAR D values for the fdscan and m0324 families have been read. YKA D values have not.
  YKA is a different array, on a different path through the inner core, measured in a
  different window. The pairs are fixed by ILAR's outer-core phases, before any YKA
  measurement.
* At ILAR, the 2025-12-13 event matched 2011 better than 2003 (exploratory, one event).
* The YKA statics census (stage A) found element-year median offsets of 0.11-0.20 s in
  2014-2019 and none of 0.1 s or more before 2013 or after 2019. Statics up to ±0.5 s
  are absorbed per event, and the pair is aligned on the direct PKIKP pulse, so clock
  offsets cannot change the coda D.

## Forecasts

| Outcome | Probability |
|---|---|
| P4-Y testable and holds | 70% |
| P2-Y testable | 55% |
| F1-Y testable (>= 6 discriminating cases in >= 3 families) | 35% |
| If testable: supported / contradicted / inconclusive | 30% / 20% / 50% |
| P1-Y holds | 45% |
| P6: ILAR and YKA D positively correlated (rho > 0, p < 0.05) | 50% |

The likeliest failure is again SNR. ILAR's pre-turn M4.3-4.9 partners may be too weak at
YKA's legacy 20-sps elements to pass direct SNR >= 5.

## What would count

If F1-Y and P2-Y are testable and P4-Y holds, this is the first out-of-sample test of
inner-core backtracking on earthquakes no published study has used. The result will be
published either way, with every pair (`yka_pairs.csv`), case (`yka_b_cases.csv`) and
merged family (`merged_multiplets.csv`).

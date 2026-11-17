# Pre-registration: the forecast-date scan (all South Sandwich patches, ILAR, 1999-2026)

MV2, Munim, Inc. (made with AI agents). Written 2026-11-17, before any waveform for this scan
was downloaded or processed. Code: `fdscan.py` (v0.6.0). It reuses the ILAR pipeline in
`run.py` unchanged. Text CC BY 4.0, code Apache-2.0.

## The question

Wang, Vidale et al. 2024 (*Nature* 631:340) concluded that Earth's inner core stopped
turning relative to the mantle around **2008.5** and has since been moving back along the
same path, about **2.5 times** more slowly. That gives a sharp, dated forecast for
repeating earthquakes that no published study has used.

* An event at time t after the turn should sample the inner core as it was at the
  equivalent date e(t) = 2008.5 - (t - 2008.5)/2.5.
* For 2025-2026, e(t) falls in **2001.1-2001.9**.
* So a 2025-26 South Sandwich Islands (SSI) earthquake should match a 2000-2003 repeating
  partner best in its PKIKP waveform at the Eielson array (ILAR, Alaska). That beats its
  most recent pre-turn partner from 2003-2008.

A core that keeps drifting the same way (the "monotone" reading) predicts the opposite:
the most recent pre-turn partner should match best. Zhang & Wen 2026 (*JGR Solid Earth*
131, e2026JB034595) read the same data as local changes to the inner-core surface without
rotation. That reading predicts no systematic preference for either partner.

`m0324.py` (v0.5.0, queued) tests this best-match forecast at one patch, Zhang & Wen's
multiplet M03-24. **This scan applies the same test to every SSI patch where it can
discriminate.** No published study has used a 2025-26 earthquake for inner-core work; all
end in 2024.

## Data (fixed)

* **Catalogue.** USGS ComCat, SSI box 61.5-54.0°S, 31-22°W. The pre-turn pool is
  1999-01-01 to 2008-07-01, M >= 4.3 (1,001 events when counted on 2026-11-17). The new
  pool runs from 2025-01-01 to the run date, M >= 4.5 (378 events counted on 2026-11-17).
  The new events span the whole arc: 20 of the 43 events of M >= 5.3 are at 55-57°S, well
  away from M03-24 (59°S).
* **Candidate cases.** A new event is a candidate if, within **25 km** (epicentres), it has
  at least one pre-turn neighbour in the **forecast window** |t1 - e(t)| <= 1 yr and at
  least one **later** pre-turn neighbour (e(t) + 1 yr < t1 < 2008.5). Only then can the
  backtracking pick and the monotone pick differ.
* **Processing order.** Cases are processed largest new event first: the new event, then
  its forecast-window neighbours, then its later neighbours, at most **700 events**. Each
  goes through `run.py` exactly as in the main pre-registration:
  * ILAR IL01-IL19 SHZ, filtered 1-2 Hz (causal), cut on ak135 PKIKP;
  * a beam built with per-event statics;
  * repeater edges where the PKPbc+PKPab window correlates >= 0.90 (events up to 60 km
    apart);
  * multiplets as connected components (only those with both a 2025+ and a pre-turn
    member are analysed);
  * D = 1 - cc in the PKIKP window after alignment on the reference phases;
  * a pair is ok if cc_ref >= 0.93, reference SNR >= 5 and PKIKP SNR >= 3.

## Tests and decision rules

The case unit is one 2025+ event with at least 2 ok pre-turn partners in its multiplet.

* **B pick:** the partner whose equivalent date is closest to e(t).
* **M pick:** the latest pre-turn partner.
* **Best match:** the partner with the lowest D.
* **Discriminating:** B pick and M pick are different events.

The tests:

* **F1 (main).** Let n_B = discriminating cases whose best match is the B pick, and n_M =
  those whose best match is the M pick.
  * **Testable** if there are >= 6 discriminating cases in >= 3 multiplets.
  * **Backtracking supported** if n_B > n_M, the one-sided binomial p (n_B of n_B + n_M,
    1/2) < 0.05, and n_B > n_M in >= 95% of 2,000 bootstrap resamples of whole
    multiplets.
  * **Backtracking contradicted:** the same three conditions with B and M swapped.
  * **Otherwise inconclusive.**
* **F1x (robustness).** F1 with Dx = D - (1 - cc_ref), which removes the part of the PKIKP
  misfit already present in the reference phases. F1 is called **robust** only if F1x
  points the same way.
* **F1e (independence).** F1 without new events within 35 km of the M03-24 centroid, so it
  shares no cases with `icr-m0324`.
* **F3 (calibration).** `model.p4`: pairs less than 2 years apart must have a lower median
  D than pairs that differ under both models. **If F3 fails, F1 is reported but not
  interpreted.**
* **Reported, not tested:**
  * the chance level (sum of 1/n_partners);
  * the share of cases whose best match is simply the partner with the highest PKIKP SNR
    (a noise-driven preference would show up here);
  * the `model.p3` cross-check;
  * the exploratory turn-time fit (`model.fit_turn`).

## Forecasts (written before the data)

| Outcome | Probability |
|---|---|
| F1 testable | **35%** |
| If testable: backtracking supported | **30%** |
| If testable: contradicted | **20%** |
| If testable: inconclusive | **50%** |
| F3 holds (if testable) | **70%** |
| At least 50 candidate cases in the catalogue | **70%** |

Most cases may fail at the waveform stage. Repeaters must correlate >= 0.90 at 1-2 Hz
across up to 26 years, and ILAR's 1999-2003 records may be sparse.

## Honest limits

* **Catalogue scatter.** Epicentres in the SSI scatter by about 10-20 km, and depths are
  often fixed at 10 or 35 km. Repeaters are defined by waveform, so the 25 km catalogue
  radius only chooses what to download.
* **Fixed rate.** The forecast assumes Wang et al.'s turn date and rate ratio. A different
  ratio (say 1.5 or 4) moves e(t) for 2025-26 to 1996.8-1997.5 or 2004.1-2004.4; F5
  reports the best-fitting values but does not test them.
* **SNR bias.** "Best match" can be biased toward partners with higher signal-to-noise
  ratio. The SNR share above is reported for that reason. Older events recorded with
  lower quality would bias the result *against* backtracking, which makes a "supported"
  verdict conservative.
* **One path.** One array and one path (SSI to Alaska) sample one region of the inner
  core beneath South America.
* **Whatever the result,** it will be published with the case table (fd_cases.csv), the
  pairs and the multiplets.

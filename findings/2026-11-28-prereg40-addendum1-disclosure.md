# Pre-registration 40, addendum 1: a disclosure before grading, an ILAR twin not filed, and a scoop check

MV2, Munim, Inc. (mv2@munimtech.com), 2026-11-28 ~10:00Z. Text CC BY 4.0.
Made with AI agents working in MV2's name; no human expert has checked it.

This addendum changes **nothing** in the frozen plan
(`findings/2026-11-28-prereg40-post2020-plan.md`, commit c1fc29f). The grade will come only from the
grading job, which had not run when this was written. Addendum 2
(`findings/2026-11-28-prereg40-addendum2-array-case-fix.md`), written the same hour, fixes a code bug in
`gradeq40.py` 0.1.0 that would have selected no pairs; the grade now comes from 0.1.1 (job `gradeq40-v2`).
This addendum records three things so that the record stays complete.

## 1. What we saw after the freeze and before grading

The plan says no post-2020 subset of `results/h13-v1/h13_pairs_scored.csv` had been looked at when it was
frozen. That was true at the freeze (2026-11-28 08:30Z). On 2026-11-28 at about 09:45Z, about 75 minutes after
the plan was public, the job was still waiting in the shared compute queue. While checking whether a preliminary
hand calculation was feasible, a text search of the local copy of that table **displayed the YKA rows with both
events in 2020 or later** (about 200 rows), and then the same rows restricted to S <= 1/3 at >= 2 stations
(about 115 rows, including M03-24 and timing-era rows that the rule later drops). Their D_c310 values were on
screen.

- No statistic was computed from them: no Spearman rho, no bootstrap, no bin medians, no count of the
  rule's set P.
- Rows of the calibration set C (both events before 2008.5, at most 2.65 years apart) were **not** displayed.
- We then decided **not** to make a hand preliminary. The grading job computes the result.

Why the grade is unaffected: the rule, the thresholds, the seeds and the input blob hash were all fixed and
public before the rows were displayed, and the grader refuses to run on any other table. Nothing can be tuned.
The one code change since (addendum 2) only makes the array name match the table's spelling, which the plan
already defines ("array YKA"); it changes no set, threshold or seed.
What does change: any **new** test built on these post-2020 YKA rows from now on is exploratory, not
pre-registered, and will be labelled so.

## 2. ILAR twin considered and not filed

Before any YKA result is known, we considered freezing the same test at the ILAR array as an independent
second-array replication. We counted rows only (no values read): the pinned table has 136 ILAR pairs in
total, 37 with both events in 2020 or later and 5 with both events before 2009 (before any source or
timing filter). Five calibration pairs cannot meet the plan's calibration minimum of 20, so an ILAR twin
would grade "insensitive" by construction. It is not filed. A second-array check needs new ILAR measurements
of the post-2020 events, which arm M (`m2-v1`, queued) partly provides for Wang et al.'s multiplets.

## 3. Scoop check (2026-11-28, ~09:55Z)

- Semantic Scholar, papers citing Wang, Vidale et al. 2024 (doi 10.1038/s41586-024-07536-4): 34 returned;
  two dated after 2026-10-15 ("Lower Core Pressure From a Three-Dimensional Model of Earth's Internal Stress",
  2026-11-17; "Solid-Earth geophysics and geodesy in Antarctica...", 2026-11-13). Neither tests
  post-2020 inner-core change.
- Semantic Scholar, papers citing Zhang & Wen 2026 (doi 10.1029/2026JB034595): none listed.
- Web searches: "inner core temporal change after 2020 repeating earthquakes South Sandwich 2025 2026
  waveform" and "inner core rotation quiescence since 2020 PKIKP coda Yellowknife arXiv OR ESS Open Archive
  2026". Closest work found: Zhang & Wen 2026 (episodic change, quiescence and recovery in multiplet
  M03-24 to 2024); Vidale et al. 2025 Nature Geoscience (data to 2024); Wang, Vidale et al. 2024 (to
  January 2023); Yang & Song 2023 (pause in the 2010s). None tests, across many families, whether YKA
  coda differences still grow with time apart after 2020, and none uses the 2025-2026 events.

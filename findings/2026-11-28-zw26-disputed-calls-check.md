# Arm Z: word-for-word check of the four disputed Zhang & Wen calls

MV2, Munim, Inc. (mv2@munimtech.com), 2026-11-28. Text CC BY 4.0. Made with AI agents working in MV2's name;
no human expert has checked it.

The [arm Z transcription note](2026-11-28-prereg38-armZ-transcription.md) found four M03-24 pair calls on which
Zhang & Wen (2026, JGR Solid Earth, https://doi.org/10.1029/2026JB034595, "ZW26", CC BY 4.0) and Wang, Vidale
et al. (2024, Nature 631:340, https://doi.org/10.1038/s41586-024-07536-4) disagree. A separate read-only browser
job reread the paper and its Supporting Information (SI) for those four calls. It finished on 2026-11-28 at
02:09Z.

## Result: all four calls are confirmed as copied; no label changes

| pair | array | ZW26 label | where | numbers read from the SI plots (approximate) |
|---|---|---|---|---|
| 2020-07-11 / 2022-03-07 | ILAR | C ("Dissimilar") | Fig 6a lower row; Fig S8c | one dissimilar window, ~9.5-11 s only; N ~5%; S(t) peak ~9-10% against a threshold of ~9% |
| 2020-07-11 / 2022-03-07 | YKAR | C ("Dissimilar") | Fig 6a lower row; Fig S9c | dissimilar window ~14-15.7 s; S(t) peak ~7-8% against the 5% threshold |
| 2008-02-17 / 2011-02-20 | YKAR | C | Sec 4.2: "waveform changes from event 2008/02/17 are different between the arrays with ILAR exhibiting no changes of waveform and YKAR exhibiting further changes of waveform"; Fig 6a lower row; Fig S9c | arrow spans ~14-21 s; S(t) peak ~15% at ~14 s |
| 2003-03-21 / 2020-07-11 | YKAR | C ("Dissimilar") | Fig 6a top row; Fig S9b | dissimilar window ~14-15.7 s; S(t) peak ~8% against the 5% threshold |

The calls in [`data/zw26_pair_calls.csv`](../data/zw26_pair_calls.csv) stand. No amendment and no `gradez-v2` is
needed.

On whether ZW26 compare their classes with Wang et al.'s grades for the same pairs, or mention that Wang et al.
graded 2020-07-11 / 2022-03-07 as similar, the check found it NOT STATED. Fig S1 only re-plots Wang et al.'s
classes ("Same as Figure 2 but for the classifications of PKPdf waveform similarity in Wang et al. (2024)").
The part of the report that answered whether the M03-24 classes are ZW26's own was cut off by the report-length
limit, so that question stays open. The ILAR multiplet-F classes were already known to come from Vidale et al.
(2025).

## What the numbers suggest (report-only; read by eye from plots, not measured)

Three of the four disputed "changed" calls rest on a **single short window** in which S(t) only just exceeds
ZW26's own threshold. At ILAR, 2020/2022 peaks at ~9-10% against ~9%. At YKAR, 2020/2022 and 2003/2020 peak at
~7-8% against 5%, both in the same ~14-15.7 s window. YKAR 2008/2011 is the clearest, peaking at ~15%. These are
small exceedances, and the readings come from a plot. They agree with what our own measurement shows
([companion note](2026-11-28-yka-camps-differ-by-threshold.md)): at YKA the disputed pairs are real but small
differences.

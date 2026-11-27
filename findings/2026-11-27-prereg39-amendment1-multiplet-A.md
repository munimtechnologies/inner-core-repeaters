# Pre-registration 39, amendment 1: multiplet A's member 734, not 724

MV2, Munim, Inc. (mv2@munimtech.com), 2026-11-27. Text CC BY 4.0, code Apache-2.0.
Made by AI agents working in MV2's name; no human expert has checked it.

**Status: made before job wsc39-v1 ran and before any mantle-P record of a new member was fetched.**
wsc39-v1 is superseded by wsc39-v2 (code `wsc39.py` 0.2.0, tests `test_wsc39.py`). Nothing else in
[pre-registration 39](2026-11-27-prereg39-wang-source-census.md) changes.

## What was wrong

Wang, Vidale et al. (2024, *Nature* 631:340, Supplementary Data Table 3) print multiplet A's events
as 712, 718, **724**, 742, 848, 854, 859. The same table's pair list for A (P21, P40, P50-P53,
P102-P105, P122-P132) and Supplementary Data Table 2 show that A's pairs use event **734**
(2008-02-17, ComCat usp000fz4d) in P40, P52, P102, P104, P124 and P129. Event 724 (2004-05-07,
usp000cucb) never appears in a multiplet-A pair: its only pair in Table 2, P32 (724-746), is
assigned to no multiplet. So "724" in Table 3's event list is a misprint for 734.

wsc39.py 0.1.0 copied the event list, so it would have tested the 2004-05-07 event, about 130 km
north of A's other members by catalogue location, as a member of A, and left out 734. Our own
forecast "A is odd with 724 an odd member: 55%" rested on that misprint and is withdrawn.

We checked every other multiplet's event list against its own pair list in Table 3 and against the
multiplet letters in Table 2. They agree. Table 2 also gives the letter F to four pairs with event 855
(2022-05-18) and the letter J to three pairs with event 858 (2022-09-10); Table 3 lists neither event
or those pairs, so pre-registration 39 (which uses Table 3) still leaves them out. Arm W, which reads
Table 2 directly, always used 734 and is not affected.

## Consequences (hand-derived, checked by the tests)

- Multiplet A = 712, 718, 734, 742, 848, 854, 859. Five of them (718, 734, 742, 848, 854) are M03-24
  reference events, so 10 of A's 21 pairs are yardstick pairs and not graded; 11 are graded
  (all pairs with 712 or 859). Total graded member pairs: 99 (was 103). Located members: 64.
- 734-718 is 16.1 km apart by catalogue location (depth difference 33.5 km in ComCat).

## New forecast for A (replaces the withdrawn one)

A is odd: 25%. Only 712 (1999) and 859 (2022-10-15) are new to the mantle-P test; 1999 events often
lack gated records, so A may rest on 859's pairs alone.

"""Inner-core rotation models and the pre-registered tests.

Part of inner-core-repeaters by MV2, Munim, Inc. (made with AI agents). Apache-2.0.

Model B ("backtracking", Wang et al. 2024, Nature 631:340): the inner core moved one
way until T = 2008.5 and has since moved back along the same path about 2.5 times more
slowly. An event at time t > T sits where the inner core was at the equivalent earlier
date e(t) = T - (t - T) / 2.5. Under B, PKIKP (DF) waveform dissimilarity between two
repeating earthquakes grows with |e(t1) - e(t2)|.

Model M ("monotone"): the inner core keeps drifting one way. DF dissimilarity grows
with |t2 - t1|.

The two models rank pairs differently only when the earlier event is before T and the
later one after it, so those "cross-turn" pairs decide the tests.
"""
import math

try:
    import numpy as np
except ImportError:  # only the date helpers work without numpy
    np = None

T_TURN = 2008.5
RATE_RATIO = 2.5
NEW_CUT = 2025.0  # events from this date on were not in any published analysis (data to 2024)


def equiv_date(t, T=T_TURN, ratio=RATE_RATIO):
    t = float(t)
    return t if t <= T else T - (t - T) / ratio


def delta_B(t1, t2, T=T_TURN, ratio=RATE_RATIO):
    return abs(equiv_date(t1, T, ratio) - equiv_date(t2, T, ratio))


def delta_M(t1, t2):
    return abs(float(t2) - float(t1))


def decimal_year(iso):
    from datetime import datetime, timezone

    s = str(iso).strip().replace("Z", "+00:00")
    dt = datetime.fromisoformat(s)
    if dt.tzinfo is None:
        dt = dt.replace(tzinfo=timezone.utc)
    y = dt.year
    a = datetime(y, 1, 1, tzinfo=timezone.utc)
    b = datetime(y + 1, 1, 1, tzinfo=timezone.utc)
    return y + (dt - a).total_seconds() / (b - a).total_seconds()


def spearman(x, y):
    """(rho, p, n) on the finite pairs; NaNs when undefined."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    m = np.isfinite(x) & np.isfinite(y)
    x, y = x[m], y[m]
    n = int(len(x))
    if n < 4 or np.ptp(x) == 0 or np.ptp(y) == 0:
        return float("nan"), float("nan"), n
    from scipy.stats import spearmanr

    r = spearmanr(x, y)
    return float(r[0]), float(r[1]), n


def _rhos(rows, T=T_TURN, ratio=RATE_RATIO, key="D"):
    D = [r[key] for r in rows]
    dB = [delta_B(r["t1"], r["t2"], T, ratio) for r in rows]
    dM = [delta_M(r["t1"], r["t2"]) for r in rows]
    return spearman(dB, D), spearman(dM, D)


def block_bootstrap(rows, nboot=1000, seed=0, key="D"):
    """Resample whole multiplets; distribution of rho_B - rho_M."""
    groups = {}
    for r in rows:
        groups.setdefault(r["mult"], []).append(r)
    keys = list(groups)
    nan = float("nan")
    if len(keys) < 2:
        return dict(n=0, frac_B_better=nan, median=nan, p05=nan, p95=nan)
    rng = np.random.default_rng(seed)
    diffs = []
    for _ in range(nboot):
        pick = rng.integers(0, len(keys), len(keys))
        rr = [x for k in pick for x in groups[keys[k]]]
        (rb, _, _), (rm, _, _) = _rhos(rr, key=key)
        if math.isfinite(rb) and math.isfinite(rm):
            diffs.append(rb - rm)
    d = np.asarray(diffs)
    if len(d) == 0:
        return dict(n=0, frac_B_better=nan, median=nan, p05=nan, p95=nan)
    return dict(n=int(len(d)), frac_B_better=float((d > 0).mean()), median=float(np.median(d)),
                p05=float(np.percentile(d, 5)), p95=float(np.percentile(d, 95)))


def n_cross(rows, T=T_TURN):
    return sum(1 for r in rows if min(r["t1"], r["t2"]) < T < max(r["t1"], r["t2"]))


def grade_set(rows, min_n=8, min_cross=3, key="D", nboot=1000):
    """P1/P2: rho_B > 0 with p < 0.05 and rho_B > rho_M in >= 90% of multiplet bootstraps."""
    (rb, pb, n), (rm, pm, _) = _rhos(rows, key=key)
    cross = n_cross(rows)
    bs = block_bootstrap(rows, nboot=nboot, key=key)
    testable = n >= min_n and cross >= min_cross
    holds = bool(testable and math.isfinite(rb) and rb > 0 and pb < 0.05
                 and math.isfinite(bs["frac_B_better"]) and bs["frac_B_better"] >= 0.9)
    return dict(n_pairs=n, n_cross_turn=cross, n_multiplets=len({r["mult"] for r in rows}),
                rho_B=rb, p_B=pb, rho_M=rm, p_M=pm, bootstrap=bs, testable=bool(testable),
                holds=holds if testable else None)


def p3(rows, new_cut=NEW_CUT, T=T_TURN, key="D"):
    """For each new event with >= 2 pre-turn partners: is its best DF match the one model B
    picks (closest equivalent date) or the one model M picks (latest pre-turn partner)?"""
    by = {}
    for r in rows:
        if r["t2"] >= new_cut and r["t1"] < T:
            by.setdefault(r["ev2"], []).append(r)
    cases = []
    for ev, rs in sorted(by.items()):
        if len(rs) < 2:
            continue
        t2 = rs[0]["t2"]
        e2 = equiv_date(t2)
        best = min(rs, key=lambda r: r[key])
        b_choice = min(rs, key=lambda r: abs(equiv_date(r["t1"]) - e2))
        m_choice = max(rs, key=lambda r: r["t1"])
        cases.append(dict(ev2=ev, t2=round(t2, 3), n_partners=len(rs), best_ev1=best["ev1"],
                          best_t1=round(best["t1"], 3), B_pick_t1=round(b_choice["t1"], 3),
                          M_pick_t1=round(m_choice["t1"], 3),
                          discriminating=b_choice["ev1"] != m_choice["ev1"],
                          matches_B=best["ev1"] == b_choice["ev1"],
                          matches_M=best["ev1"] == m_choice["ev1"]))
    dc = [c for c in cases if c["discriminating"]]
    nB = sum(1 for c in dc if c["matches_B"])
    nM = sum(1 for c in dc if c["matches_M"])
    testable = len(dc) >= 2
    return dict(cases=cases, n_discriminating=len(dc), n_match_B=nB, n_match_M=nM,
                testable=testable, holds=(nB > nM) if testable else None)


def p4(rows, key="D"):
    """Calibration: pairs < 2 yr apart are more alike in DF than pairs that differ under both models."""
    near = [r[key] for r in rows if delta_M(r["t1"], r["t2"]) < 2.0]
    far = [r[key] for r in rows if delta_M(r["t1"], r["t2"]) > 10.0 and delta_B(r["t1"], r["t2"]) > 5.0]
    testable = len(near) >= 3 and len(far) >= 3
    mn = float(np.median(near)) if near else float("nan")
    mf = float(np.median(far)) if far else float("nan")
    return dict(n_near=len(near), n_far=len(far), median_D_near=mn, median_D_far=mf,
                testable=testable, holds=(mn < mf) if testable else None)


def fit_turn(rows, key="D", Ts=None, ratios=(1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0)):
    """Exploratory: the (T, ratio) whose delta_B best rank-correlates with DF dissimilarity."""
    if Ts is None:
        Ts = list(np.arange(1998.0, 2020.01, 0.25)) + [2030.0]  # 2030 = no turn within the data
    D = [r[key] for r in rows]
    best = None
    for T in Ts:
        for q in ratios:
            rb = spearman([delta_B(r["t1"], r["t2"], T, q) for r in rows], D)[0]
            if math.isfinite(rb) and (best is None or rb > best[0]):
                best = (rb, float(T), float(q))
    rm = spearman([delta_M(r["t1"], r["t2"]) for r in rows], D)[0]
    if best is None:
        return dict(n_pairs=len(rows), best_rho=float("nan"), best_T=None, best_ratio=None, rho_M=rm)
    return dict(n_pairs=len(rows), best_rho=best[0], best_T=best[1], best_ratio=best[2], rho_M=rm)


def evaluate(rows, new_cut=NEW_CUT, nboot=1000):
    ok = [r for r in rows if r.get("ok")]
    ins = [r for r in ok if r["t2"] < new_cut]
    outs = [r for r in ok if r["t2"] >= new_cut]
    res = dict(
        n_rows=len(rows), n_ok=len(ok), n_in_sample=len(ins), n_out_of_sample=len(outs),
        P1_in_sample=grade_set(ins, nboot=nboot),
        P2_out_of_sample=grade_set(outs, min_n=6, min_cross=3, nboot=nboot),
        P3_new_event_best_match=p3(ok, new_cut=new_cut),
        P4_calibration=p4(ok),
        P5_fit_in_sample=fit_turn(ins),
        P5_fit_all=fit_turn(ok),
    )
    have_x = [r for r in ok if r.get("Dx") is not None and math.isfinite(r.get("Dx", float("nan")))]
    res["robust_excess_D_in"] = grade_set([r for r in have_x if r["t2"] < new_cut], key="Dx", nboot=nboot)
    res["robust_excess_D_out"] = grade_set([r for r in have_x if r["t2"] >= new_cut], min_n=6,
                                           key="Dx", nboot=nboot)
    have_s = [r for r in ok if r.get("Ds") is not None and math.isfinite(r.get("Ds", float("nan")))]
    res["robust_station_D_in"] = grade_set([r for r in have_s if r["t2"] < new_cut], key="Ds", nboot=nboot)
    res["robust_station_D_out"] = grade_set([r for r in have_s if r["t2"] >= new_cut], min_n=6,
                                            key="Ds", nboot=nboot)
    return res

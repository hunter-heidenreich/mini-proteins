"""Convergence certification for a trajectory's collective variables.

Engine-agnostic: everything operates on plain arrays (dihedral time series or
features), so the same certification applies to our OpenMM trajectories and to
ingested MDGen/other data. This is the "certify" track -- the ground truth is
only as trustworthy as its convergence, and most released datasets ship none.

Ported from scripts/analyze.py (the validated GROMACS analysis) so the two
pipelines report the same statistics.

All functions are pure NumPy (no MD or ML deps), so they run anywhere and are
unit-testable on synthetic series.
"""

import numpy as np


def acf_fft(x):
    """Normalized autocorrelation function of a 1-D series, via FFT (O(n log n))."""
    x = np.asarray(x, dtype=float)
    x = x - x.mean()
    n = x.size
    if n < 2 or not np.any(x):
        return np.array([1.0])
    f = np.fft.rfft(x, n=2 * n)
    acf = np.fft.irfft(f * np.conjugate(f))[:n].real
    if acf[0] == 0:
        return np.array([1.0])
    return acf / acf[0]


def statistical_inefficiency(x):
    """g = 1 + 2 sum_t (1 - t/n) C(t), truncated at the first zero crossing.

    Effective independent samples = n / g (Janke / Chodera convention).
    """
    x = np.asarray(x, dtype=float)
    n = x.size
    acf = acf_fft(x)
    g = 1.0
    for t in range(1, acf.size):
        c = acf[t]
        if c <= 0:
            break
        g += 2.0 * c * (1.0 - t / n)
    return max(1.0, g)


def decorrelation_time(x, dt):
    """Lag (in units of dt) at which the autocorrelation first drops below 1/e."""
    acf = acf_fft(x)
    below = np.where(acf < np.exp(-1.0))[0]
    return float(below[0] * dt) if below.size else float(acf.size * dt)


def transition_counts(states):
    """Count changes in a discrete state label series (ignoring None entries).

    Returns total transitions; e.g. for a 2-state {0,1} series this is the
    number of crossings (Hu 2003-style ergodicity check).
    """
    n = 0
    last = None
    for s in states:
        if s is None:
            continue
        if last is not None and s != last:
            n += 1
        last = s
    return n


def stationarity(values, binned=False, bins=None):
    """First-half vs second-half discrepancy as a convergence check.

    If binned, compare normalized histograms (max abs bin-prob diff); otherwise
    compare means (relative to the series std).
    """
    v = np.asarray(values, dtype=float)
    h = v.size // 2
    a, b = v[:h], v[h:]
    if binned:
        edges = bins if bins is not None else np.histogram_bin_edges(v, bins=36)
        pa, _ = np.histogram(a, bins=edges, density=True)
        pb, _ = np.histogram(b, bins=edges, density=True)
        pa = pa / pa.sum() if pa.sum() else pa
        pb = pb / pb.sum() if pb.sum() else pb
        return float(np.abs(pa - pb).max())
    sd = v.std() or 1.0
    return float(abs(a.mean() - b.mean()) / sd)


def certify_series(x, dt=1.0):
    """Convergence summary for one CV time series (e.g. cos(phi))."""
    x = np.asarray(x, dtype=float)
    g = statistical_inefficiency(x)
    return {
        "n_frames": int(x.size),
        "statistical_inefficiency": round(float(g), 1),
        "effective_samples": int(x.size / g),
        "decorrelation_time": round(decorrelation_time(x, dt), 2),
        "stationarity_halves": round(stationarity(x), 3),
    }


def certify_replicas(per_replica_cv):
    """Between-replica agreement: mean and spread of a scalar across replicas.

    per_replica_cv: list of 1-D arrays (one CV series per replica). Returns the
    cross-replica mean/std of the per-replica means -- agreement is the headline
    convergence evidence.
    """
    means = np.array([np.asarray(c, float).mean() for c in per_replica_cv])
    return {"replicas": int(means.size),
            "cross_replica_mean": round(float(means.mean()), 4),
            "cross_replica_std": round(float(means.std()), 4)}

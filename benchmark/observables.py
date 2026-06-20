"""Reference observables + comparison metrics.

The "answer key" the field evaluates against, computed the standard way so our
numbers are comparable: a phi/psi free-energy surface, Jensen-Shannon divergence
on distributions, and (via deeptime) TICA slow components + an MSM with implied
timescales. Re-deriving the alanine MSM here is the first replication milestone:
the timescales have a known ballpark, so it is self-checkable.

FES + JSD are pure NumPy (unit-testable). TICA + MSM require deeptime and are
guarded -- verify the exact calls against the installed deeptime version on the
run host (the API has shifted across releases).
"""

import numpy as np

R_KCAL = 8.314462618e-3 / 4.184   # kcal/(mol K)

try:
    from deeptime.decomposition import TICA
    from deeptime.clustering import KMeans
    from deeptime.markov.msm import MaximumLikelihoodMSM
    _HAVE_DEEPTIME = True
except ImportError:  # pragma: no cover
    _HAVE_DEEPTIME = False


# ---- pure-NumPy: FES + JSD ------------------------------------------------
def fes_2d(phi, psi, temperature=298.0, bins=72):
    """phi/psi free-energy surface (kcal/mol), min-shifted to zero.

    Returns (F, edges) with F shape (bins, bins); NaN where unsampled.
    """
    edges = np.linspace(-180, 180, bins + 1)
    H, _, _ = np.histogram2d(np.asarray(phi), np.asarray(psi), bins=[edges, edges])
    P = H / H.sum()
    with np.errstate(divide="ignore"):
        F = -R_KCAL * temperature * np.log(P)
    F -= np.nanmin(F[np.isfinite(F)])
    F[~np.isfinite(F)] = np.nan
    return F, edges


def jsd(p, q, base=2.0):
    """Jensen-Shannon divergence between two probability vectors (default bits)."""
    p = np.asarray(p, float)
    q = np.asarray(q, float)
    p = p / p.sum()
    q = q / q.sum()
    m = 0.5 * (p + q)

    def _kl(a, b):
        mask = a > 0
        return np.sum(a[mask] * (np.log(a[mask] / b[mask]) / np.log(base)))

    return float(0.5 * _kl(p, m) + 0.5 * _kl(q, m))


def jsd_samples(a, b, bins=64, rng=(-180, 180)):
    """JSD between two 1-D samples via shared-edge histograms (e.g. a torsion)."""
    edges = np.linspace(rng[0], rng[1], bins + 1)
    pa, _ = np.histogram(np.asarray(a), bins=edges)
    pb, _ = np.histogram(np.asarray(b), bins=edges)
    # add a tiny floor so empty bins don't make KL blow up
    return jsd(pa + 1e-12, pb + 1e-12)


# ---- deeptime: TICA + MSM (guarded) --------------------------------------
def tica_project(feature_list, lagtime, dim=2):
    """Fit TICA on a list of (n_frames, n_feat) arrays; return projected list."""
    if not _HAVE_DEEPTIME:
        raise ImportError("deeptime required for TICA; pip install deeptime")
    model = TICA(lagtime=lagtime, dim=dim).fit_fetch(feature_list)
    return [model.transform(f) for f in feature_list], model


def msm_timescales(projected_list, lagtime, n_clusters=100, n_timescales=4):
    """Cluster TICA space, build an MSM, return implied timescales (frames).

    Returns {timescales, stationary_entropy, n_states}. Multiply timescales by
    the frame stride to get physical units. Verify deeptime API on the host.
    """
    if not _HAVE_DEEPTIME:
        raise ImportError("deeptime required for MSM; pip install deeptime")
    clustering = KMeans(n_clusters=n_clusters, max_iter=50).fit_fetch(
        np.concatenate(projected_list))
    dtrajs = [clustering.transform(p) for p in projected_list]
    msm = MaximumLikelihoodMSM(lagtime=lagtime).fit_fetch(dtrajs)
    ts = msm.timescales(k=n_timescales)
    return {"timescales_frames": [float(x) for x in ts],
            "n_states": int(msm.n_states),
            "lagtime_frames": int(lagtime)}

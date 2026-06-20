"""Trajectory featurization: backbone + side-chain torsions -> sin/cos.

This is the shared, interoperable representation (MDGen convention): every method
and metric downstream consumes sin/cos of the torsions, so our OpenMM data and
ingested MDGen/other trajectories are featurized identically and stay comparable.

The dihedral math is pure NumPy (unit-testable, no deps). The trajectory loader
uses MDTraj when available; callers that already hold coordinate arrays (e.g.
MDGen .npy) can skip it and call compute_dihedrals / featurize directly.
"""

import numpy as np

from benchmark.systems import BACKBONE, SYSTEMS

try:
    import mdtraj as md
except ImportError:  # pragma: no cover
    md = None


def dihedral(p0, p1, p2, p3):
    """Dihedral angle (degrees) for points along the last axis, vectorized.

    Each p* has shape (..., 3). Returns angle in (-180, 180].
    """
    b0 = np.asarray(p0) - np.asarray(p1)
    b1 = np.asarray(p2) - np.asarray(p1)
    b2 = np.asarray(p3) - np.asarray(p2)
    b1 = b1 / np.linalg.norm(b1, axis=-1, keepdims=True)
    v = b0 - (b0 * b1).sum(-1, keepdims=True) * b1
    w = b2 - (b2 * b1).sum(-1, keepdims=True) * b1
    x = (v * w).sum(-1)
    y = (np.cross(b1, v) * w).sum(-1)
    return np.degrees(np.arctan2(y, x))


def compute_dihedrals(coords, quads):
    """Dihedral time series for each atom-index quad.

    coords: (n_frames, n_atoms, 3); quads: list of (i,j,k,l). Returns
    (n_frames, n_quads) in degrees.
    """
    coords = np.asarray(coords, dtype=float)
    out = np.empty((coords.shape[0], len(quads)))
    for c, (i, j, k, l) in enumerate(quads):
        out[:, c] = dihedral(coords[:, i], coords[:, j], coords[:, k], coords[:, l])
    return out


def featurize(angles_deg):
    """Stack sin/cos of an (n_frames, n_angles) array -> (n_frames, 2*n_angles).

    Columns are [sin a0, cos a0, sin a1, cos a1, ...] -- the periodicity-aware
    representation that TICA/MSM/diffusion all expect on the torus.
    """
    a = np.deg2rad(np.asarray(angles_deg, dtype=float))
    feats = np.empty((a.shape[0], 2 * a.shape[1]))
    feats[:, 0::2] = np.sin(a)
    feats[:, 1::2] = np.cos(a)
    return feats


# ---- atom-index resolution (needs a topology) ----------------------------
def _sel(top, expr):
    idx = top.select(expr)
    if len(idx) == 0:
        raise ValueError(f"no atom matched selection: {expr!r}")
    return int(idx[0])


def torsion_quads(top, residue):
    """Resolve the (name, atom-index-quad) list for a capped dipeptide.

    Backbone phi/psi/omega from systems.BACKBONE (Ace-X-Nme convention) plus the
    side-chain chi torsions for this residue. `top` is an MDTraj topology.
    """
    sysmeta = SYSTEMS[residue]
    resname = sysmeta["resname"]
    quads = []
    # phi: C(ACE)-N(X)-CA(X)-C(X)
    quads.append(("phi", (_sel(top, "resname ACE and name C"),
                          _sel(top, f"resname {resname} and name N"),
                          _sel(top, f"resname {resname} and name CA"),
                          _sel(top, f"resname {resname} and name C"))))
    # psi: N(X)-CA(X)-C(X)-N(NME)
    quads.append(("psi", (_sel(top, f"resname {resname} and name N"),
                          _sel(top, f"resname {resname} and name CA"),
                          _sel(top, f"resname {resname} and name C"),
                          _sel(top, "resname NME and name N"))))
    # omega: CH3(ACE)-C(ACE)-N(X)-CA(X)
    if "omega" in sysmeta["slow_dofs"]:
        quads.append(("omega", (_sel(top, "resname ACE and name CH3"),
                                _sel(top, "resname ACE and name C"),
                                _sel(top, f"resname {resname} and name N"),
                                _sel(top, f"resname {resname} and name CA"))))
    # side-chain chi torsions (all atoms within residue X)
    for n, names in enumerate(sysmeta["chi"], start=1):
        quads.append((f"chi{n}", tuple(
            _sel(top, f"resname {resname} and name {nm}") for nm in names)))
    return quads


def load_features(traj_path, top_path, residue):
    """Load a trajectory and return (names, angles_deg, sincos_features).

    angles_deg: (n_frames, n_torsions); features: (n_frames, 2*n_torsions).
    Requires MDTraj. coords are in nm but dihedrals are scale-invariant.
    """
    if md is None:
        raise ImportError("mdtraj required to load trajectories; "
                          "pip install mdtraj (or pass coords to compute_dihedrals)")
    t = md.load(traj_path, top=top_path)
    quads = torsion_quads(t.topology, residue)
    names = [n for n, _ in quads]
    angles = compute_dihedrals(t.xyz, [q for _, q in quads])
    return names, angles, featurize(angles)

"""Phase-0 curation & analysis: turn raw trajectories into a certified,
ML-ready, comparable artifact.

For one (residue, tier) this:
  - loads every replica (OpenMM DCD, or any MDTraj-readable trajectory),
  - featurizes to torsions + sin/cos (the shared representation),
  - CERTIFIES convergence (per-DOF statistical inefficiency / effective samples,
    helical<->extended transition counts, half/half stationarity, cross-replica
    agreement),
  - computes reference OBSERVABLES (phi/psi FES; TICA + MSM implied timescales
    via deeptime -- re-deriving the alanine timescales is the self-checkable
    replication milestone),
  - exports a curated .npz (angles, features, replica ids, time, FES) + a
    human/JSON report.

Engine-agnostic by design: point --traj/--top at MDGen's released trajectories
and the same pipeline applies (the "interoperate" track).

Usage:
  uv run --extra ml python -m benchmark.curate --residue ala --tier explicit
  uv run --extra ml python -m benchmark.curate --residue ala --tier implicit --tica-lag 100 --msm-lag 200
  uv run --extra ml python -m benchmark.curate --residue ala --tier explicit \
         --traj path/to/*.dcd --top path/to/topology.pdb
"""

import argparse
import glob
import json
import os
import sys

import numpy as np

from benchmark import certify, observables
from benchmark.features import load_features
from benchmark.systems import SYSTEMS


def _macro_state(phi, psi):
    """Coarse 2-state label for ergodicity counting: helical vs extended."""
    if -160 <= phi <= -20 and -120 <= psi <= 50:
        return "H"
    if phi < -20 and (psi > 90 or psi < -120):
        return "E"
    return None


def run(args):
    if args.residue not in SYSTEMS:
        sys.exit(f"unknown residue {args.residue!r}; choose from {list(SYSTEMS)}")
    base = os.path.join("out", args.residue, args.tier)
    trajs = (sorted(glob.glob(args.traj)) if args.traj
             else sorted(glob.glob(os.path.join(base, "r*", "traj.dcd"))))
    if not trajs:
        sys.exit(f"no trajectories found (looked in {base}/r*/traj.dcd); "
                 f"run benchmark.simulate first or pass --traj")
    top = args.top or os.path.join(base, "r1", "topology.pdb")
    if not os.path.exists(top):
        sys.exit(f"topology not found: {top} (pass --top)")

    # ---- load + featurize every replica --------------------------------
    names = None
    per_angles, per_feats = [], []
    for tr in trajs:
        names, angles, feats = load_features(tr, top, args.residue)
        per_angles.append(angles)
        per_feats.append(feats)
        print(f"  {tr}: {angles.shape[0]} frames, torsions={names}")
    iphi, ipsi = names.index("phi"), names.index("psi")
    all_angles = np.concatenate(per_angles)
    n_frames = all_angles.shape[0]

    report = {"residue": args.residue, "tier": args.tier,
              "replicas": len(trajs), "frames_total": int(n_frames),
              "torsions": names, "frame_stride_ps": args.dt_ps}

    # ---- CERTIFY -------------------------------------------------------
    # per-DOF convergence on cos(angle) (periodicity-safe scalar series)
    per_dof = {}
    for j, nm in enumerate(names):
        per_dof[nm] = certify.certify_series(
            np.cos(np.deg2rad(all_angles[:, j])), dt=args.dt_ps)
    # helical<->extended transitions per replica (Hu 2003 style)
    trans = []
    for angles in per_angles:
        states = [_macro_state(p, s) for p, s in zip(angles[:, iphi], angles[:, ipsi])]
        trans.append(certify.transition_counts(states))
    ns_per_rep = (n_frames / len(trajs)) * args.dt_ps / 1000.0
    # cross-replica agreement on helical fraction
    helical_frac = [np.mean([_macro_state(p, s) == "H"
                             for p, s in zip(a[:, iphi], a[:, ipsi])])
                    for a in per_angles]
    # self-consistency / agreement via JSD on phi
    half_jsd = float(np.mean([
        observables.jsd_samples(a[:a.shape[0] // 2, iphi], a[a.shape[0] // 2:, iphi])
        for a in per_angles]))
    cross_jsd = (observables.jsd_samples(per_angles[0][:, iphi],
                                         per_angles[1][:, iphi])
                 if len(per_angles) > 1 else None)
    report["certification"] = {
        "per_dof": per_dof,
        "helical_extended_transitions_per_replica": trans,
        "transitions_per_ns_mean": round(float(np.mean(trans)) / max(ns_per_rep, 1e-9), 2),
        "helical_fraction_mean": round(float(np.mean(helical_frac)), 4),
        "helical_fraction_std": round(float(np.std(helical_frac)), 4),
        "phi_jsd_first_vs_second_half": round(half_jsd, 4),
        "phi_jsd_cross_replica": (round(cross_jsd, 4) if cross_jsd is not None else None),
    }

    # ---- OBSERVABLES ---------------------------------------------------
    F, edges = observables.fes_2d(all_angles[:, iphi], all_angles[:, ipsi])
    obs = {"fes_min_kcal": 0.0, "fes_max_kcal": round(float(np.nanmax(F)), 2)}
    if observables._HAVE_DEEPTIME:
        try:
            proj, _ = observables.tica_project(per_feats, lagtime=args.tica_lag, dim=2)
            ts = observables.msm_timescales(proj, lagtime=args.msm_lag,
                                            n_clusters=args.n_clusters)
            ts["timescales_ps"] = [round(t * args.dt_ps, 1)
                                   for t in ts["timescales_frames"]]
            obs["msm"] = ts
        except Exception as e:  # deeptime API drift / too few states
            obs["msm"] = {"error": f"{type(e).__name__}: {e}"}
    else:
        obs["msm"] = {"skipped": "deeptime not installed"}
    report["observables"] = obs

    # ---- EXPORT --------------------------------------------------------
    os.makedirs(base, exist_ok=True)
    replica_id = np.concatenate([np.full(a.shape[0], i)
                                 for i, a in enumerate(per_angles)])
    time_ps = np.concatenate([np.arange(a.shape[0]) * args.dt_ps
                              for a in per_angles])
    np.savez_compressed(
        os.path.join(base, "curated.npz"),
        angles_deg=all_angles, features=np.concatenate(per_feats),
        replica_id=replica_id, time_ps=time_ps,
        torsion_names=np.array(names), fes=F, fes_edges=edges,
        split="by_replica")  # correlation-aware: never random-split frames
    with open(os.path.join(base, "curation.json"), "w") as f:
        json.dump(report, f, indent=2)
    with open(os.path.join(base, "curation.txt"), "w") as f:
        f.write(_format(report))
    print(_format(report))
    print(f"\nwrote {base}/curated.npz, curation.json, curation.txt")


def _format(r):
    c = r["certification"]
    L = [f"Curation report: {r['residue']} / {r['tier']}",
         "=" * 60,
         f"replicas {r['replicas']}  frames {r['frames_total']}  "
         f"stride {r['frame_stride_ps']} ps  torsions {r['torsions']}",
         "",
         "CONVERGENCE (per degree of freedom)"]
    for nm, d in c["per_dof"].items():
        L.append(f"  {nm:6s} g={d['statistical_inefficiency']:>7}  "
                 f"eff={d['effective_samples']:>7}  "
                 f"decorr={d['decorrelation_time']} ps  "
                 f"stat={d['stationarity_halves']}")
    L += ["",
          f"  helical<->extended transitions/ns  {c['transitions_per_ns_mean']} "
          f"{c['helical_extended_transitions_per_replica']}",
          f"  helical fraction (cross-replica)   {c['helical_fraction_mean']} "
          f"+/- {c['helical_fraction_std']}",
          f"  phi JSD half/half                  {c['phi_jsd_first_vs_second_half']} "
          f"({'OK' if c['phi_jsd_first_vs_second_half'] < 0.02 else 'CHECK'})",
          f"  phi JSD cross-replica              {c['phi_jsd_cross_replica']}",
          "",
          "OBSERVABLES",
          f"  FES span 0..{r['observables']['fes_max_kcal']} kcal/mol"]
    msm = r["observables"]["msm"]
    if "timescales_ps" in msm:
        L.append(f"  MSM implied timescales (ps)  {msm['timescales_ps']} "
                 f"[{msm['n_states']} states, lag {msm['lagtime_frames']}]")
    else:
        L.append(f"  MSM: {msm.get('skipped') or msm.get('error')}")
    return "\n".join(L)


def build_parser():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--residue", required=True)
    p.add_argument("--tier", required=True, choices=["explicit", "implicit"])
    p.add_argument("--traj", default=None, help="glob for trajectories (overrides default)")
    p.add_argument("--top", default=None, help="topology file (PDB); default out/<id>/<tier>/r1/topology.pdb")
    p.add_argument("--dt-ps", type=float, default=0.1, help="frame stride (ps)")
    p.add_argument("--tica-lag", type=int, default=100, help="TICA lagtime (frames)")
    p.add_argument("--msm-lag", type=int, default=200, help="MSM lagtime (frames)")
    p.add_argument("--n-clusters", type=int, default=100)
    return p


if __name__ == "__main__":
    run(build_parser().parse_args())

"""OpenMM generation skeleton for the capped-dipeptide benchmark (Phase 1).

Plain unbiased MD for one (residue, tier, replica), in either solvent tier:
  - explicit : TIP3P + PME (the realism tier; cross-validates against GROMACS)
  - implicit : GB/OBC2   (the generative tier; u(x) is tractable for BG/diffusion)

Force field is ff03 in both tiers (matches the prior GROMACS validation). The
explicit tier deliberately mirrors the GROMACS protocol -- ff03/TIP3P, rigid
water (SETTLE), flexible solute, Langevin at 298 K with friction = 1/tau_t
(GROMACS used tau_t = 0.1 ps) -- so an Ala-explicit run here can be compared
directly to the validated GROMACS FES/populations as a cross-engine check.

This is Phase 1 (plain MD only). The benchmark's kinetically-faithful targets
need more (see docs/benchmark.md):
  PHASE 2 TODO: REST2 (replica exchange w/ solute tempering) for unbiased
                equilibrium incl. rare states (alphaL, chi flips, pucker).
  PHASE 2 TODO: adaptive sampling -> MSM for rates / implied timescales.
Those build on the System/Topology constructed here.

Usage:
  python -m benchmark.simulate --residue ala --tier explicit --replica 1 --ns 100
  python -m benchmark.simulate --residue pro --tier implicit --replica 1 --ns 100 \
         --equil-ns 0.5 --report-ps 1.0

Outputs (out/<residue>/<tier>/r<replica>/):
  traj.dcd      trajectory at --report-ps stride
  scalars.csv   step,time,potential,kinetic,temperature,(box) for QC/drift
  state.xml     final serialized State (positions+velocities+box) for restart
  meta.json     full provenance: params, seed, ff, openmm/git versions
"""

import argparse
import json
import os
import subprocess
import sys

try:
    import openmm as mm
    from openmm import app, unit
except ImportError:  # pragma: no cover - skeleton importable without OpenMM
    mm = app = unit = None

from benchmark.systems import SYSTEMS

# ---- tier configuration --------------------------------------------------
# ff03 first in both; the second file selects the solvent treatment. If
# 'amber03.xml' is not bundled with this OpenMM build, install openmmforcefields
# and load ff03 from there -- keep ff03 for comparability with the GROMACS tier.
TIERS = {
    "explicit": {
        "forcefield": ["amber03.xml", "tip3p.xml"],
        "periodic": True,
        "solvate": True,
    },
    "implicit": {
        "forcefield": ["amber03.xml", "implicit/obc2.xml"],
        "periodic": False,
        "solvate": False,
    },
}


def _build_system(residue, tier_cfg, args):
    """Return (Modeller-or-PDB topology, positions, System) for the tier."""
    pdb = app.PDBFile(os.path.join("data", f"{residue}.pdb"))
    ff = app.ForceField(*tier_cfg["forcefield"])
    topology, positions = pdb.topology, pdb.positions

    if tier_cfg["solvate"]:
        modeller = app.Modeller(pdb.topology, pdb.positions)
        # padding/neutralize/ionic strength mirror the GROMACS preprocess step
        modeller.addSolvent(ff, model="tip3p",
                            padding=args.padding * unit.nanometer,
                            neutralize=True,
                            ionicStrength=args.ionic_strength * unit.molar)
        topology, positions = modeller.topology, modeller.positions
        system = ff.createSystem(
            topology, nonbondedMethod=app.PME,
            nonbondedCutoff=args.cutoff * unit.nanometer,
            constraints=None,        # flexible solute (complete forces)
            rigidWater=True)         # SETTLE on water, as in GROMACS
    else:
        # implicit solvent: no box, no PME; tiny system -> no cutoff needed
        system = ff.createSystem(
            topology, nonbondedMethod=app.NoCutoff,
            constraints=None, rigidWater=True)
    return topology, positions, system


def _make_integrator(args):
    friction = (1.0 / args.tau_t) / unit.picosecond     # GROMACS tau_t -> friction
    integ = mm.LangevinMiddleIntegrator(
        args.temperature * unit.kelvin, friction, args.dt_fs * unit.femtosecond)
    integ.setRandomNumberSeed(args.seed)
    return integ


def _provenance(args, outdir, n_atoms):
    try:
        git = subprocess.check_output(["git", "rev-parse", "--short", "HEAD"],
                                      stderr=subprocess.DEVNULL).decode().strip()
    except Exception:
        git = None
    return {
        "residue": args.residue, "resname": SYSTEMS[args.residue]["resname"],
        "tier": args.tier, "replica": args.replica, "seed": args.seed,
        "forcefield": TIERS[args.tier]["forcefield"],
        "temperature_K": args.temperature, "tau_t_ps": args.tau_t,
        "friction_per_ps": round(1.0 / args.tau_t, 3),
        "dt_fs": args.dt_fs, "production_ns": args.ns, "equil_ns": args.equil_ns,
        "report_ps": args.report_ps, "n_atoms": n_atoms,
        "openmm_version": (mm.version.version if mm else None),
        "git_commit": git,
        "slow_dofs": SYSTEMS[args.residue]["slow_dofs"],
        "rare_processes": SYSTEMS[args.residue]["rare"],
    }


def run(args):
    if mm is None:
        sys.exit("OpenMM not importable. Install openmm (and openmmforcefields "
                 "if ff03 is not bundled) on the run host.")
    if args.residue not in SYSTEMS:
        sys.exit(f"unknown residue {args.residue!r}; choose from {list(SYSTEMS)}")
    tier_cfg = TIERS[args.tier]
    outdir = os.path.join("out", args.residue, args.tier, f"r{args.replica}")
    os.makedirs(outdir, exist_ok=True)

    topology, positions, system = _build_system(args.residue, tier_cfg, args)
    n_atoms = topology.getNumAtoms()
    print(f"[{args.residue}/{args.tier}/r{args.replica}] {n_atoms} atoms, "
          f"seed={args.seed}")

    steps_per_ps = int(round(1000.0 / args.dt_fs))
    report_interval = int(round(args.report_ps * steps_per_ps))
    equil_steps = int(round(args.equil_ns * 1000 * steps_per_ps))
    prod_steps = int(round(args.ns * 1000 * steps_per_ps))

    platform = (mm.Platform.getPlatformByName(args.platform)
                if args.platform else None)

    # ---- minimize + equilibrate ----------------------------------------
    # Explicit: barostat ON during equilibration (NPT) to set the box, then a
    # fresh barostat-free System for NVT production (GROMACS production used
    # pcoupl=no). Implicit: no box/barostat; equilibrate at T directly.
    if tier_cfg["periodic"]:
        system.addForce(mm.MonteCarloBarostat(
            1.0 * unit.bar, args.temperature * unit.kelvin, 25))

    integ = _make_integrator(args)
    sim = (app.Simulation(topology, system, integ, platform)
           if platform else app.Simulation(topology, system, integ))
    sim.context.setPositions(positions)
    sim.minimizeEnergy()
    sim.context.setVelocitiesToTemperature(args.temperature * unit.kelvin,
                                           args.seed)
    if equil_steps:
        sim.step(equil_steps)
    equil_state = sim.context.getState(getPositions=True, getVelocities=True,
                                       enforcePeriodicBox=tier_cfg["periodic"])

    # ---- production (NVT) ----------------------------------------------
    if tier_cfg["periodic"]:
        # rebuild System without the barostat for constant-volume production
        _, _, system = _build_system(args.residue, tier_cfg, args)
    prod_integ = _make_integrator(args)
    sim = (app.Simulation(topology, system, prod_integ, platform)
           if platform else app.Simulation(topology, system, prod_integ))
    sim.context.setState(equil_state)

    sim.reporters.append(app.DCDReporter(
        os.path.join(outdir, "traj.dcd"), report_interval))
    sim.reporters.append(app.StateDataReporter(
        os.path.join(outdir, "scalars.csv"), report_interval, step=True,
        time=True, potentialEnergy=True, kineticEnergy=True, temperature=True,
        density=tier_cfg["periodic"], speed=True))
    sim.step(prod_steps)

    final = sim.context.getState(getPositions=True, getVelocities=True,
                                 enforcePeriodicBox=tier_cfg["periodic"])
    with open(os.path.join(outdir, "state.xml"), "w") as f:
        f.write(mm.XmlSerializer.serialize(final))
    with open(os.path.join(outdir, "meta.json"), "w") as f:
        json.dump(_provenance(args, outdir, n_atoms), f, indent=2)
    print(f"  wrote {outdir}/ (traj.dcd, scalars.csv, state.xml, meta.json)")


def build_parser():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--residue", required=True, help="3-letter id, e.g. ala")
    p.add_argument("--tier", required=True, choices=list(TIERS))
    p.add_argument("--replica", type=int, default=1)
    p.add_argument("--ns", type=float, default=100.0, help="production ns")
    p.add_argument("--equil-ns", type=float, default=0.5)
    p.add_argument("--report-ps", type=float, default=1.0,
                   help="trajectory/scalars stride (ps)")
    p.add_argument("--dt-fs", type=float, default=1.0,
                   help="timestep (fs); 1 fs keeps the solute flexible")
    p.add_argument("--temperature", type=float, default=298.0)
    p.add_argument("--tau-t", type=float, default=0.1,
                   help="GROMACS-style coupling time (ps); friction = 1/tau_t")
    p.add_argument("--cutoff", type=float, default=1.0, help="PME cutoff (nm)")
    p.add_argument("--padding", type=float, default=1.0,
                   help="solvent box padding (nm), explicit only")
    p.add_argument("--ionic-strength", type=float, default=0.1,
                   help="salt concentration (mol/L), explicit only")
    p.add_argument("--seed", type=int, default=1,
                   help="RNG seed (set distinct per replica for independence)")
    p.add_argument("--platform", default=None,
                   help="OpenMM platform, e.g. CUDA / CPU (default: auto)")
    return p


if __name__ == "__main__":
    run(build_parser().parse_args())

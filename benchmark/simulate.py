"""OpenMM generation skeleton for the capped-dipeptide benchmark (Phase 1).

Plain unbiased MD for one (residue, tier, replica), in either solvent tier:
  - explicit : amber14/TIP3P + PME (realism tier)
  - implicit : amber14 + GBn2     (generative tier; u(x) tractable for BG/diffusion)

Defaults are aligned with the comparable literature (MDGen 2024, Timewarp 2023,
Transferable Boltzmann Generators 2024) so this benchmark is cross-comparable:
ff14SB, GBn2 implicit / TIP3P explicit, Langevin at 310 K, 2 fs with H-bond
constraints, low friction (0.1 ps^-1 implicit / 0.3 ps^-1 explicit), and a fine
100 fs save stride (needed for the trajectory-upsampling task).

GROMACS cross-engine check (separate, one-off): the prior GROMACS validation is
ff03, not ff14SB, so it does NOT validate the benchmark FES directly. Instead run
ff03 in OpenMM and confirm it reproduces the GROMACS ff03 FES/populations -- a
*thermodynamics* check (unaffected by thermostat/friction differences):
  python -m benchmark.simulate --residue ala --tier explicit \
         --forcefield amber03.xml tip3p.xml --temperature 298 --friction 10

This is Phase 1 (plain MD only). Kinetically-faithful targets need more
(docs/benchmark.md):
  PHASE 2 TODO: REST2 for unbiased equilibrium incl. rare states.
  PHASE 2 TODO: adaptive sampling -> MSM for rates / implied timescales.
  PHASE 5 TODO: tetrapeptides (frontier system size).

Usage:
  python -m benchmark.simulate --residue ala --tier explicit --replica 1 --ns 100
  python -m benchmark.simulate --residue pro --tier implicit --replica 1 --ns 100

Outputs (out/<residue>/<tier>/r<replica>/):
  traj.dcd      trajectory at --report-ps stride
  scalars.csv   step,time,potential,kinetic,temperature,(density) for QC/drift
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
# Benchmark force field is amber14/ff14SB (field standard) for comparability;
# override with --forcefield for the ff03 GROMACS cross-check. If amber14-all is
# not bundled in this OpenMM build, install openmmforcefields.
TIERS = {
    "explicit": {
        "forcefield": ["amber14-all.xml", "amber14/tip3p.xml"],
        "water_model": "tip3p",
        "periodic": True, "solvate": True,
        "friction_per_ps": 0.3,      # Timewarp/MDGen explicit
    },
    "implicit": {
        "forcefield": ["amber14-all.xml", "implicit/gbn2.xml"],
        "water_model": None,
        "periodic": False, "solvate": False,
        "friction_per_ps": 0.1,      # MDGen implicit
    },
}


def _forcefield_files(args, tier_cfg):
    return args.forcefield if args.forcefield else tier_cfg["forcefield"]


def _build_system(residue, tier_cfg, args):
    """Return (topology, positions, System) for the tier."""
    pdb = app.PDBFile(os.path.join("data", f"{residue}.pdb"))
    ff = app.ForceField(*_forcefield_files(args, tier_cfg))
    topology, positions = pdb.topology, pdb.positions

    if tier_cfg["solvate"]:
        modeller = app.Modeller(pdb.topology, pdb.positions)
        modeller.addSolvent(ff, model=tier_cfg["water_model"],
                            padding=args.padding * unit.nanometer,
                            neutralize=True,
                            ionicStrength=args.ionic_strength * unit.molar)
        topology, positions = modeller.topology, modeller.positions
        system = ff.createSystem(
            topology, nonbondedMethod=app.PME,
            nonbondedCutoff=args.cutoff * unit.nanometer,
            constraints=app.HBonds, rigidWater=True)   # 2 fs-capable; field std
    else:
        system = ff.createSystem(
            topology, nonbondedMethod=app.NoCutoff,
            constraints=app.HBonds, rigidWater=True)
    return topology, positions, system


def _make_integrator(args, friction):
    integ = mm.LangevinMiddleIntegrator(
        args.temperature * unit.kelvin, friction / unit.picosecond,
        args.dt_fs * unit.femtosecond)
    integ.setRandomNumberSeed(args.seed)
    return integ


def _provenance(args, tier_cfg, friction, n_atoms):
    try:
        git = subprocess.check_output(["git", "rev-parse", "--short", "HEAD"],
                                      stderr=subprocess.DEVNULL).decode().strip()
    except Exception:
        git = None
    return {
        "residue": args.residue, "resname": SYSTEMS[args.residue]["resname"],
        "tier": args.tier, "replica": args.replica, "seed": args.seed,
        "forcefield": _forcefield_files(args, tier_cfg),
        "temperature_K": args.temperature, "friction_per_ps": friction,
        "dt_fs": args.dt_fs, "constraints": "HBonds (heavy atoms flexible)",
        "production_ns": args.ns, "equil_ns": args.equil_ns,
        "report_ps": args.report_ps, "n_atoms": n_atoms,
        "openmm_version": (mm.version.version if mm else None),
        "git_commit": git,
        "slow_dofs": SYSTEMS[args.residue]["slow_dofs"],
        "rare_processes": SYSTEMS[args.residue]["rare"],
    }


def run(args):
    if mm is None:
        sys.exit("OpenMM not importable. Install openmm (and openmmforcefields "
                 "if amber14/ff03 not bundled) on the run host.")
    if args.residue not in SYSTEMS:
        sys.exit(f"unknown residue {args.residue!r}; choose from {list(SYSTEMS)}")
    tier_cfg = TIERS[args.tier]
    friction = (args.friction if args.friction is not None
                else tier_cfg["friction_per_ps"])
    outdir = os.path.join("out", args.residue, args.tier, f"r{args.replica}")
    os.makedirs(outdir, exist_ok=True)

    topology, positions, system = _build_system(args.residue, tier_cfg, args)
    n_atoms = topology.getNumAtoms()
    # topology.pdb lets MDTraj/benchmark.curate load the DCD afterwards
    with open(os.path.join(outdir, "topology.pdb"), "w") as fh:
        app.PDBFile.writeFile(topology, positions, fh)
    print(f"[{args.residue}/{args.tier}/r{args.replica}] {n_atoms} atoms, "
          f"ff={_forcefield_files(args, tier_cfg)}, T={args.temperature}K, "
          f"friction={friction}/ps, seed={args.seed}")

    steps_per_ps = int(round(1000.0 / args.dt_fs))
    report_interval = max(1, int(round(args.report_ps * steps_per_ps)))
    equil_steps = int(round(args.equil_ns * 1000 * steps_per_ps))
    prod_steps = int(round(args.ns * 1000 * steps_per_ps))
    platform = (mm.Platform.getPlatformByName(args.platform)
                if args.platform else None)

    # ---- minimize + equilibrate ----------------------------------------
    # Explicit: barostat ON during equilibration (NPT) to set the box, then a
    # barostat-free System for NVT production. Implicit: no box; equilibrate at T.
    if tier_cfg["periodic"]:
        system.addForce(mm.MonteCarloBarostat(
            1.0 * unit.bar, args.temperature * unit.kelvin, 25))
    sim = (app.Simulation(topology, system, _make_integrator(args, friction),
                          platform) if platform else
           app.Simulation(topology, system, _make_integrator(args, friction)))
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
        _, _, system = _build_system(args.residue, tier_cfg, args)  # no barostat
    sim = (app.Simulation(topology, system, _make_integrator(args, friction),
                          platform) if platform else
           app.Simulation(topology, system, _make_integrator(args, friction)))
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
        json.dump(_provenance(args, tier_cfg, friction, n_atoms), f, indent=2)
    print(f"  wrote {outdir}/ (traj.dcd, scalars.csv, state.xml, meta.json)")


def build_parser():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--residue", required=True, help="3-letter id, e.g. ala")
    p.add_argument("--tier", required=True, choices=list(TIERS))
    p.add_argument("--replica", type=int, default=1)
    p.add_argument("--ns", type=float, default=100.0, help="production ns")
    p.add_argument("--equil-ns", type=float, default=0.5)
    p.add_argument("--report-ps", type=float, default=0.1,
                   help="trajectory/scalars stride (ps); 0.1 = 100 fs (MDGen)")
    p.add_argument("--dt-fs", type=float, default=2.0,
                   help="timestep (fs); 2 fs with H-bond constraints (field std)")
    p.add_argument("--temperature", type=float, default=310.0,
                   help="K; field uses 310-350 to accelerate sampling")
    p.add_argument("--friction", type=float, default=None,
                   help="Langevin friction (1/ps); default per tier "
                        "(0.3 explicit / 0.1 implicit). Use 10 for GROMACS check.")
    p.add_argument("--forcefield", nargs="+", default=None,
                   help="override ff xml files, e.g. amber03.xml tip3p.xml "
                        "(for the GROMACS ff03 cross-engine check)")
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

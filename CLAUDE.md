# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

A GROMACS molecular-dynamics pipeline for capped dipeptides ("mini-proteins", Ace-X-Nme). It runs unbiased Langevin (stochastic) MD and **extracts atomic forces** to build training datasets for ML potentials, while also serving as a validation ensemble against Ramachandran/J-coupling literature. The 9 curated input structures live in `data/<id>.pdb` (three-letter codes: `ala`, `gly`, `ile`, `leu`, `met`, `phe`, `pro`, `trp`, `val`).

## Running the pipeline

Everything is driven by `ID=<code>` (the three-letter residue) as an environment variable. The numbered scripts are meant to run in order; `run.sh` chains all five.

```bash
ID=ala sh scripts/run.sh           # full pipeline: preprocess -> equil -> production -> post -> analyze
GPU=1 ID=ala sh scripts/run.sh     # offload MD stages to a CUDA GPU
```

Individual stages (same `ID=` convention):
- `scripts/0_preprocess.sh` — pdb2gmx → editconf → solvate → genion. Writes to `out/<id>/raw/`.
- `scripts/1_equil.sh` — EM (steepest descent) → NVT (100 ps) → NPT (200 ps).
- `scripts/2_md_lang.sh` — `NREP` independent production replicas (default 5 × 100 ns). Override: `NREP=3 NSTEPS=50000000 ID=ala sh scripts/2_md_lang.sh` (`NSTEPS` = 1 fs steps, overrides the `.mdp`).
- `scripts/3_post.sh` — per replica: extract solute forces (`*_forces.xvg`, the ML labels) + PBC-corrected trajectory (`<id>-sim-rN.pdb`); runs `plot.py`; promotes products from `raw/` to `out/<id>/data/`.
- `scripts/4_analyze.sh` — `gmx rama` per replica then `analyze.py`, which emits a single validation + dataset-QC report (`out/<id>/data/validation.txt` and `validation.json`): provenance, convergence stats (statistical inefficiency / effective samples / basin-transition counts / half-half stationarity), FES + per-basin free energies, a side-by-side comparison to the published ff03/TIP3P populations, ³J(HN,Hα) under several Karplus sets, glycine symmetry, and force-label QC. Force-label and energy-drift sections read the per-replica `*_forces.xvg` / `*_etot.xvg` / `*_temp.xvg` from stage 3 (in `data/`); without them analyze.py still runs and prints the skip.
- `scripts/bench.sh` — times production throughput and prints estimated `$/protein`. Requires `out/<id>/raw/npt.gro` (run 0 + 1 first). Env: `GPU`, `PRICE` ($/hr), `TARGET_NS`, `BENCH_NSTEPS`.

The Python scripts are invoked by the shell scripts and also read `ID` from the environment: `ID=ala python scripts/plot.py`. Deps: `pip install -r requirements.txt` (numpy, matplotlib).

There are no tests, linters, or a build step — this is a simulation pipeline, not an application.

### Storage tiers (RunPod is ephemeral — nothing valuable lives only there)

- **git:** code, configs, input PDBs, docs, env lock (`requirements-ml.lock`), and **small certified reports** committed under `results/` (the `out/` tree is gitignored, so `scripts/pull_results.sh` copies `validation.*`/`curation.*`/figures into `results/` for committing — run it before stopping the pod).
- **RunPod (scratch):** raw/intermediate trajectories + the GPU env; regenerable from code + the seeds in each run's `meta.json`.
- **HuggingFace** (`hheiden/mini-proteins-bench`, dataset repo): `curated.npz` and any trajectories worth keeping — too big for git, the durable/shareable deliverable.

### Benchmark suite (OpenMM, newer work — see `docs/benchmark.md`)

A second, in-progress track is a **learning + interoperable framework** (`docs/benchmark.md`) for the craft of generating MD data and benchmarking ML methods (reduced-representation, dynamics propagators, generative/diffusion) on these dipeptides. Goal is learning + being comparable to the field, **not** novelty: four tracks — generate (OpenMM), certify (convergence/ground-truth rigor, our lens), interoperate (use MDGen's released data + TICA/JSD/MSM metrics so results are comparable), replicate (re-implement canonical methods as milestones). It is **single-engine OpenMM** with two solvent tiers — explicit TIP3P/PME (realism) and implicit **GBn2** (generative — tractable `u(x)` for Boltzmann generators). Force field is **amber14/ff14SB** across both, chosen for comparability with the literature this benchmark targets (MDGen, Timewarp, Transferable BG; all OpenMM/amber14). `ff03` is kept **only** for a one-off GROMACS cross-engine check (`--forcefield amber03.xml tip3p.xml --temperature 298 --friction 10` → match the GROMACS ff03 FES; a thermodynamics check). Field-aligned defaults: 310 K, 2 fs + HBond constraints, friction 0.1/0.3 ps⁻¹, 100 fs save stride. The differentiator vs existing released datasets is **rigor** (certified-converged ground truth, per-DOF convergence gating, shipped MSM/FES answer key, honest kinetic ceilings) not system size; tetrapeptides are Phase 5.

- `benchmark/systems.py` — per-residue metadata (slow DOFs incl. side-chain χ / Pro ω+pucker, difficulty, rare processes, transferability split). Single source of truth for what each system's ground truth must cover; φ/ψ convergence ≠ full convergence.
- `benchmark/simulate.py` — Phase-1 OpenMM generation skeleton (plain MD, both tiers): `python -m benchmark.simulate --residue ala --tier explicit --ns 100 --replica 1 --seed 1`. Explicit tier mirrors the GROMACS protocol (rigid water, flexible solute, Langevin friction = 1/tau_t) for the cross-validation. Outputs to `out/<id>/<tier>/r<n>/`.
- `benchmark/run_shakeout.sh` — drives the Ala/Gly/Pro trio × both tiers × replicas.
- **Phase 0 curation/analysis (engine-agnostic, the foundation):** `benchmark/features.py` (torsions → sin/cos; pure-numpy `dihedral` + MDTraj loader), `benchmark/certify.py` (convergence stats: statistical inefficiency, effective samples, transitions, stationarity — pure numpy, ported from `analyze.py`), `benchmark/observables.py` (FES + JSD in numpy; TICA + MSM via deeptime, guarded), `benchmark/curate.py` (CLI: load replicas → featurize → certify → observables → write `out/<id>/<tier>/curated.npz` + `curation.{json,txt}`). Works on OpenMM DCD **and** any MDTraj-readable trajectory (e.g. MDGen's released data) — point `--traj/--top` at it. The numpy cores are unit-tested; TICA/MSM (deeptime) and trajectory loading (mdtraj) need the run host. Deps: `requirements-ml.txt` (openmm, openmmforcefields, mdtraj, deeptime).
- Phases 2+ (REST2 for rare-event equilibrium, adaptive-sampling→MSM for kinetics, curation/splits/metrics) are specified in `docs/benchmark.md` but not yet built. `cis/trans ω` is scoped to thermodynamics only (seconds-timescale rate is beyond reach).

## Architecture & conventions

- **Output layout** (all gitignored under `out/`): `out/<id>/raw/` holds GROMACS working files (`.tpr`, `.gro`, `.trr`, `.cpt`, intermediate `.xvg`); `out/<id>/data/` holds final products; `out/<id>/figs/` holds PNG plots. `3_post.sh` is what moves files from `raw/` → `data/`.
- **Replica discovery is by glob**, not by re-reading `NREP`. `2_md_lang.sh` writes `md_lang_r<N>.trr`; `3_post.sh` and `4_analyze.sh` loop over `md_lang_r*.trr`, so they stay correct for any replica count without coordination. Preserve the `md_lang_r<N>` naming if you touch these.
- **Independent replicas, shared start.** Each replica is `grompp`'d separately from the same `npt.gro` but with `gen_vel=yes` and `gen_seed=ld_seed=-1` (random per run). Crucially, production does **not** pass `-t npt.cpt` — that is deliberate, so replicas diverge instead of being identical. Don't "fix" this by adding the checkpoint.

## Physics invariants (do not break these silently)

These are the scientific point of the repo; changing them changes what the data means.

- **One Hamiltonian throughout.** EM → NVT → NPT → production all use PME at 1.0 nm cutoff with rigid TIP3P water. Keep electrostatics/cutoffs consistent across the four `.mdp` files in `config/`.
- **Solute is unconstrained in production** (only water is rigid via SETTLE) so the written forces are complete `-∇U` labels. The `.trr` writes forces (`nstfout`) because force extraction is the primary purpose.
- **`sd` (Langevin) integrator** for canonical sampling, not velocity-rescaling. This has GPU consequences:
  - EM (steepest descent) does **not** support PME-on-GPU → EM gets `-nb gpu` only; NVT/NPT get `-nb gpu -pme gpu`.
  - Production cannot use `-update gpu` (GROMACS allows it only with the `md` integrator, not `sd`), so integration stays on CPU. This limits GPU speedup for these tiny (~2,300-atom) systems — a known limitation, see `docs/cloud-run.md`.
- Force field / water model (`amber03` / `tip3p`) are set at the top of `0_preprocess.sh`. README notes that modern ML work may want a newer FF (ff19SB / CHARMM36m); changing it is intentional and localized there.

## GROMACS interaction patterns

- `gmx energy`/`gmx traj`/`gmx trjconv`/`gmx rama` need an interactive group/term selection on stdin. Scripts pipe these with `printf 'Name\n'` (selecting **by name**, e.g. `Potential`, `Protein`) rather than by index number — this is intentional robustness against force-field/term-layout changes. Keep using names.
- `read_xvg` in both Python scripts filters header lines by prefix (`#`, `@`, `&`) rather than skipping a fixed count, because the header length varies. Reuse this when parsing new `.xvg` outputs.

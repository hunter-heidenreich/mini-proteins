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

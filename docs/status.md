# Status & next actions

The single "where are we / what's next" file — open this to pick up without
re-reading chat history. Keep it current. (As of 2026-06-21.)

For the *why* behind the design see `docs/benchmark.md`; for the box procedure
see `docs/cloud-run.md`; for commands see `CLAUDE.md`.

## Where we are

**GROMACS track — done & certified.**
- Pipeline (`scripts/0_*`–`4_*`, `run.sh`, `bench.sh`) validated.
- `scripts/analyze.py` produces a full validation + dataset-QC report
  (`validation.txt`/`.json`): convergence, FES, populations vs Vymětal ff03/TIP3P,
  ³J(HN,Hα) under multiple Karplus sets (Graf 2007 value pinned).
- **Alanine result certified** (the `validation.txt` you pasted): converged
  (~1,134 effective samples), populations within ~1σ of the ff03/TIP3P reference,
  forces clean. Findings: ff03's "helix" is **bridge/C7eq-dominated** (canonical
  αR≈0.21, bridge≈0.25), and ⟨³J⟩ is **+1.15 Hz** vs experiment under Graf's own
  Karplus set (a real ensemble difference, ff03's known over-helical bias).
- **Alanine fully persisted & raw reclaimed (2026-06-21).** Certified report
  (`validation.*` + figures) committed + pushed to git under `results/ala/`. The
  big data products — gzipped force labels + PBC trajectories (5 replicas, 1.3 GB
  → 249 MB) + a `README.md` manifest — uploaded to HF
  `hheiden/mini-proteins-bench/ala/gromacs-ff03/` (commit `dd39360`). The ~27 GB
  of raw `.trr`/`.edr`/`.log` and the local `out/ala/data/` are **deleted**;
  `out/ala` is now ~34 MB (figs + small restart/rama files). Regenerable from the
  seeds in `validation.json`. **`gh` and the modern `hf` CLI are installed on the
  box** (container disk).

**Benchmark track (OpenMM, single-engine, ff14SB, two tiers) — written, partially verified.**
- `benchmark/systems.py` — 9-residue metadata (slow DOFs incl. χ / Pro ω+pucker),
  shake-out trio, transferability split.
- `benchmark/simulate.py` — Phase-1 plain-MD generation skeleton (both tiers).
- Phase-0 curation: `features.py` (torsions→sin/cos), `certify.py` (convergence
  stats), `observables.py` (FES/JSD + TICA/MSM via deeptime), `curate.py` (CLI).
- Design + landscape + roadmap in `docs/benchmark.md`.

**Tooling & storage.**
- `uv` (`pyproject.toml` + `uv.lock`): `uv sync` (core) / `uv sync --extra ml`.
- git = code + small certified reports (`results/`, via `scripts/pull_results.sh`);
  RunPod = scratch; HF dataset `hheiden/mini-proteins-bench` = `curated.npz` + big data.

## Verified locally vs pending the box

| | status |
|---|---|
| numpy cores (dihedral, certify stats, FES, JSD) | ✅ unit-tested locally |
| uv sync / uv run / imports / `uv lock --check` | ✅ locally |
| OpenMM generation (`simulate.py`) | ⏳ never run — needs openmm + GPU |
| TICA/MSM (`observables.py`, deeptime) | ⏳ guarded; deeptime API unverified |
| MDTraj trajectory loading (`features.py`) | ⏳ unverified on real DCD/TRR |

## Next actions — on RunPod, in order

```bash
git pull
# 1. environment
uv sync --extra ml
# 2. confirm the openmm PyPI wheel sees the GPU
uv run --extra ml python -c "import openmm; \
  print(openmm.Platform.getPlatformByName('CUDA').getName())"
# 3. shake-out generation (Ala/Gly/Pro, both tiers, short)
NS=20 PLATFORM=CUDA sh benchmark/run_shakeout.sh
# 4. ff03 cross-engine check vs the certified GROMACS alanine FES
uv run --extra ml python -m benchmark.simulate --residue ala --tier explicit \
  --forcefield amber03.xml tip3p.xml --temperature 298 --friction 10 --ns 20
# 5. curate alanine -> first replication milestone (MSM self-check)
uv run --extra ml python -m benchmark.curate --residue ala --tier explicit
# 6. before stopping the pod: persist
sh scripts/pull_results.sh && git add results && git commit -m "results: shake-out" && git push
# HF uses the modern `hf` CLI (huggingface-cli is deprecated); hf is installed
# standalone at /root/.local/bin/hf, no `uv run` needed. Xet is on by default.
HF_XET_HIGH_PERFORMANCE=1 hf upload hheiden/mini-proteins-bench \
  out/ala/explicit/curated.npz ala/explicit/curated.npz --repo-type=dataset
```

**Self-check for step 5:** the alanine slow φ/ψ implied timescale should land in
the **hundreds-of-ps** range (GROMACS gave statistical inefficiency g≈441 ps). If
the MSM timescales are wildly off, the TICA/MSM lag or clustering needs tuning —
not a data problem.

## Open confirmations / risks (resolve as you hit them)

- **openmm CUDA from PyPI** — if the wheel's CUDA build mismatches the box, fall
  back to a conda openmm; everything else is pure-PyPI.
- **deeptime API drift** — `observables.tica_project` / `msm_timescales` use
  `fit_fetch` / `timescales`; verify against the installed deeptime version (calls
  are wrapped in try/except in `curate.py`).
- **topology for MDTraj** — OpenMM runs write `topology.pdb` (handled). For the
  GROMACS cross-check on the existing `.trr`, MDTraj may need a `.gro`/`.pdb`
  rather than `.tpr` as `--top`.

## Reclaiming disk (persist, *then* delete)

> **Done for alanine (2026-06-21)** following exactly this procedure — see "Where
> we are". The pattern below is the template for the remaining residues.
> Packaging that worked well: **gzip the verbose text before HF upload** (PDB
> ~7.9×, force `.xvg` ~2.2×; 1.3 GB → 249 MB) and ship a `README.md` manifest
> (provenance/units) alongside. Keep `validation.*` uncompressed so it renders on
> the HF web UI. Don't hand-roll an `.npz` here — that's `benchmark/curate.py`'s job.

Raw trajectories are scratch and regenerable — `validation.json` records the
seeds + GROMACS version + git commit, so a run can be reproduced. Safe to delete
to reclaim space, **but only after the small certified derivatives are in git**
(they otherwise live only on the ephemeral pod):

```bash
sh scripts/pull_results.sh && git add results && git commit -m "results: <id>" && git push
du -sh out/<id>/raw/*.trr                 # the space hog (production trajectories)
rm -f out/<id>/raw/md_lang_r*.trr out/<id>/raw/*.cpt \
      out/<id>/raw/em.* out/<id>/raw/nvt.* out/<id>/raw/npt.*
```

Nothing downstream needs the `.trr`: the ff03 cross-check (step 4) compares
*results* (FES/populations in `validation.json` + the FES figure), not raw
frames. Optional insurance: push the small `md_lang_r*_rama.xvg` (φ/ψ series, a
few MB) to HF if you want to recompute dihedral observables without re-simulating.

## After the shake-out checks out

- Phase 2 (kinetics: adaptive-sampling/MSM), Phase 3 (all 9 dipeptides, both
  tiers, per-DOF convergence gating), Phase 4 (replication milestones), Phase 5
  (tetrapeptides — ingest MDGen's released data). See `docs/benchmark.md`.

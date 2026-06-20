# Status & next actions

The single "where are we / what's next" file — open this to pick up without
re-reading chat history. Keep it current. (As of commit `eb3faf2`, 2026-06.)

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
uv run --extra ml huggingface-cli upload --repo-type dataset \
  hheiden/mini-proteins-bench out/ala/explicit/curated.npz ala/explicit/curated.npz
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

## After the shake-out checks out

- Phase 2 (kinetics: adaptive-sampling/MSM), Phase 3 (all 9 dipeptides, both
  tiers, per-DOF convergence gating), Phase 4 (replication milestones), Phase 5
  (tetrapeptides — ingest MDGen's released data). See `docs/benchmark.md`.

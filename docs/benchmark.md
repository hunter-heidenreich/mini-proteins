# Capped-dipeptide benchmark suite — design

A benchmark for **method development** on small biomolecular systems, not a
production dataset for a better force field. The systems are deliberately tiny
and their equilibrium/kinetic ground truth is knowable, so a method can be
*scored*, not just demonstrated. Target method classes:

1. **Reduced-representation modeling** — autoencoders / manifold learning that
   collapse conformations into a low-dimensional latent.
2. **Dynamics propagation & sampling** — transfer-operator / latent propagators
   (VAMPnets, ITO-style) and equilibrium samplers (Boltzmann generators).
3. **Generative / diffusion** — conformation diffusion and *video-diffusion*
   (trajectory-segment generation, transition-path inpainting).

This sits in a well-established line: dipeptide suites are the de-facto
transferability benchmark for learned dynamics (Timewarp), latent propagators
(Implicit Transfer Operators), and Boltzmann generators; torsional diffusion
works on exactly the torsion torus these systems live on.

## Systems — the difficulty curriculum

Nine capped dipeptides (Ace-X-Nme), chosen so the suite spans intrinsic
dimension and topology. The benchmark's value is partly this gradient: a method
that aces Ala should be *challenged* by Trp/Pro.

| id | residue | extra slow DOFs beyond φ/ψ | what it tests |
|----|---------|----------------------------|---------------|
| gly | Gly | none (symmetric) | symmetric ~2D manifold; (φ,ψ)→(−φ,−ψ) symmetry is a free correctness check |
| ala | Ala | none | canonical 2D Ramachandran; the reference everything is calibrated against |
| val | Val | χ1 (β-branched) | rotamer coupling to backbone |
| ile | Ile | χ1, χ2 (β-branched) | 2 coupled side-chain torsions |
| leu | Leu | χ1, χ2 | flexible aliphatic side chain |
| met | Met | χ1, χ2, χ3 | long, floppy side chain → highest side-chain dimension |
| phe | Phe | χ1, χ2 (ring flip) | aromatic rotamer + symmetry |
| trp | Trp | χ1, χ2 (indole) | largest side chain; slow χ relaxation |
| pro | Pro | ω cis/trans, ring pucker (endo/exo) | restricted φ; *different topology*; the hard case |

Implication for analysis: **φ/ψ convergence ≠ full convergence.** Each system's
ground truth must include its side-chain χ angles and (for Pro) ω and ring
pucker. Per-system run length is gated by the *slowest relevant DOF*, not a flat
wall-clock (see Convergence gating).

## Engine: OpenMM, single engine

One engine for both solvent tiers. Rationale:

- **Removes an engine confound from the headline comparison.** The benchmark's
  key axis is implicit-vs-explicit solvent. If the two tiers came from different
  engines, thermostat/PME/constraint/RNG differences would contaminate that
  comparison — you could not attribute a tier difference to *solvent* vs
  *engine*. Single engine eliminates this.
- **Energy-tractable in-loop.** OpenMM evaluates `u(x)` and forces for any
  configuration in Python — exactly the API Boltzmann generators and
  energy-guided diffusion need at train time.
- **Ecosystem fit.** Timewarp / bgflow / mdgen / ITO all live in OpenMM, so the
  benchmark is reproducible and comparable for its intended audience.
- At ~2,300 atoms (explicit) GROMACS's speed edge is irrelevant.

Force field is **ff03 across both tiers** (consistency + comparability with the
prior GROMACS validation). ff03 in OpenMM comes from `amber03.xml` if bundled,
otherwise via `openmmforcefields`.

### GROMACS is retained as a cross-engine validation anchor
The earlier GROMACS alanine runs (ff03/TIP3P, validated against Vymětal &
Vondrášek 2010) are **not discarded**. Re-running alanine-explicit in OpenMM and
confirming it reproduces the GROMACS FES / populations / timescales is an
*independent-engine agreement check* — evidence the OpenMM pipeline is correct
that a single-engine project normally cannot provide. The GROMACS analysis logic
(`scripts/analyze.py`: basin partitions, convergence stats, literature
comparison, J-couplings, force QC) is engine-agnostic and carries over with new
I/O adapters (MDTraj instead of `gmx rama`/`gmx energy`).

## Two tiers

| | Tier 1 — explicit (realism) | Tier 2 — implicit (generative) |
|---|---|---|
| solvent | TIP3P, PME, 1.0 nm cutoff | GB (OBC2 / GBn2), no box |
| `u(solute)` | not tractable (solvent-marginalized) | **tractable** → BG / energy-guided / reweighting |
| cost | the compute sink (×9 REST2) | cheap → aggressive rare-event sampling |
| role | physical reference; cross-validates Tier 2 | primary substrate for energy-based generative methods |

Expect the tiers to **differ** in populations and kinetics — that difference is
itself a benchmark axis ("how much does the solvent model distort the learned
object"), not an error to reconcile.

## Scope: kinetically faithful, with an honest ceiling

Ground truth must reproduce equilibrium *and* kinetics including rare states —
but rates and free energies are different targets across timescales, so the
benchmark tiers the **processes**:

| process | timescale | free energy | rate / kinetics |
|---|---|---|---|
| αR/bridge/PPII/β | 10²–10³ ps | plain MD | plain MD + MSM |
| αL | barrier-limited | enhanced sampling | seeded-MSM |
| side-chain χ flips | 1–100 ns | longer MD / REST2 | seeded-MSM |
| **cis/trans ω** (Pro esp.) | **seconds–min** | enhanced sampling | **out of scope as a rate** — reported as thermodynamics only, documented |

A benchmark that ships a rate it cannot defend is broken; the seconds-timescale
ω rate is explicitly declared beyond the kinetic ceiling.

### Method stack
- **Equilibrium (all DOFs, all rare states, both tiers): REST2** (replica
  exchange with solute tempering). CV-free, so it works uniformly across all 9
  residues — including Pro pucker and bulky χ — without hand-picking collective
  variables; the base replica is unbiased. (Metadynamics on φ,ψ is the classic
  alanine route but needs per-residue CVs.)
- **Kinetics (rates/timescales): adaptive sampling → MSM.** Biased trajectories
  distort time, so rates come from *unbiased* short trajectories seeded from the
  states REST2 discovered; build an MSM, read off implied timescales and rates.
  The MSM is what ships as the kinetic answer key.

## Curation: dataset → benchmark

The connective tissue that makes it scorable and reproducible:

1. **ML-ready export** — per (residue, tier, replica): arrays of
   `(coords, time, replica, [forces, energy])` in `.npz`/`.h5`, plus a fixed
   **internal-coordinate** convention (φ,ψ,χ…,ω,pucker on the torus) — diffusion
   and reduced-rep both prefer torsions; sidesteps SE(3) equivariance.
2. **Reference observables (the answer key)** — FES, per-DOF marginals, MSM
   implied timescales, metastable states (PCCA+) — shipped per system/tier.
3. **Standardized splits**
   - *within-system:* by replica or time block (never random — autocorrelation
     ~200+ ps means random splits leak near-identical frames → inflated metrics).
   - *transferability:* train-residues vs held-out-residues (Timewarp-style),
     e.g. train {ala,gly,val,leu,phe} → test {ile,met,trp,pro}.
4. **Per-system property sheet** — intrinsic dimension, slow DOFs, what each
   residue stresses.

## Convergence gating (replaces flat run length)

Each system runs until its **slowest relevant DOF** is converged — quantified,
not assumed:
- REST2 free energy of every basin stable to < ~0.2 kcal/mol across the second
  half of the run;
- MSM implied timescales flat vs lag time (Chapman–Kolmogorov pass);
- per-DOF (φ,ψ,ω,χ,pucker) statistical inefficiency → effective-sample and
  transition counts adequate (carried over from `analyze.py`).

Documented per residue, so the kinetic claims are auditable.

## Benchmark tasks & metrics (the leaderboard)

| task | method class | metric vs ground truth |
|---|---|---|
| latent embedding | reduced-rep / AE | reconstruction error; recovered intrinsic dim; FES preserved under embedding |
| propagator | transfer operator | implied-timescale error; Chapman–Kolmogorov; rate error |
| equilibrium sampler | Boltzmann generator | FES / marginal match; reweighted free-energy error (Tier 2, needs `u(x)`) |
| trajectory generation | (video-)diffusion | distributional+kinetic match of generated vs MD: FES, timescales, n-step transition densities (frame stride is a parameter; ~10–50 ps shows motion, 1 ps is jitter) |

Metrics are **distributional/kinetic**, not per-frame: a generated trajectory is
good iff its implied FES and timescales match the MD's.

## Roadmap

- **Phase 0 — curation layer (engine-agnostic).** Export + reference-observable
  computer + correlation-aware/transferability splits, runnable on the existing
  GROMACS alanine data. Unblocks ML immediately.
- **Phase 1 — OpenMM generation skeleton.** Plain MD, both tiers, Ala/Gly/Pro
  shake-out (canonical / symmetric / ring+cis-trans). Cross-validate
  Ala-explicit against GROMACS. *(This document ships with the Phase-1
  skeleton: `benchmark/simulate.py`, `benchmark/systems.py`.)*
- **Phase 2 — kinetics.** Adaptive-sampling/MSM; ship MSM timescales; ω
  thermodynamics-only.
- **Phase 3 — scale to all 9, both tiers**, with per-DOF convergence gating.
- **Phase 4 — define tasks/metrics/splits as a released benchmark.**

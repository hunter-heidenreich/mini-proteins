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

## Prior art & how this differs (verified against the papers)

| work | systems | engine / FF / solvent | what it does |
|---|---|---|---|
| **Timewarp** (Klein 2023) | AD, 2AA, **4AA** | OpenMM / amber14 / implicit | transferable time-coarsened MD acceleration (normalizing flow) |
| **MDGen** (Jing 2024) | **tetra/penta-peptides** | OpenMM / amber14 / gbn2 + tip3pfb | generative trajectories: forward sim, TPS, upsampling, inpainting, design |
| **ITO** (Schreiner 2023) | alanine dipeptide + CG | OpenMM / — | multi-time-resolution latent propagators (SE(3) diffusion) |
| **Transferable BG** (Klein 2024) | dipeptides | — / — | transferable zero-shot Boltzmann generators + reweighting |

Two consequences this benchmark is designed around:

1. **Dipeptides alone are largely solved** — Timewarp/MDGen have moved to
   tetra/pentapeptides. So the dipeptide suite is positioned as the
   **interpretable diagnostic floor**, not the frontier, with tetrapeptides
   planned (Phase 5). The frontier-relevance comes from Phase 5; the dipeptides'
   value is being *controlled and known*.
2. **The differentiation is rigor, not systems.** The released datasets above are
   "we ran it." None ship *certified-converged ground truth*: audited
   convergence (effective samples, transition counts), per-DOF convergence
   gating, a shipped MSM/FES answer key, honest kinetic ceilings (e.g. cis/trans
   ω declared out of scope), literature-anchored validation, and a *per-residue
   diagnostic curriculum with known slow DOFs*. That is the niche — the
   **validated/diagnostic** benchmark methods can be *scored and debugged*
   against, not just trained on.

To stay cross-comparable we adopt their conventions where sensible: ff14SB,
gbn2/TIP3P, OpenMM, the transferability split, torsion (sin/cos) featurization →
TICA, and Jensen–Shannon divergence as a distributional metric (see Metrics).

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

Force field is **amber14 / ff14SB across both tiers** — the de-facto standard in
the comparable literature (MDGen, Timewarp, Transferable Boltzmann Generators),
chosen for cross-comparability. ff03 (the earlier GROMACS choice) is older and
documented over-helical; carrying it into the benchmark would make the data
incomparable. ff14SB comes from `amber14-all.xml` (or `openmmforcefields`).

### GROMACS is retained as a cross-engine validation anchor
The earlier GROMACS alanine runs (ff03/TIP3P, validated against the literature)
are **not discarded** — but because the benchmark FF is now ff14SB, they do *not*
validate the benchmark FES directly (different FF, different surface). The check
instead is: **run ff03 in OpenMM** (`--forcefield amber03.xml tip3p.xml
--temperature 298 --friction 10`) and confirm it reproduces the GROMACS ff03
FES/populations. Equilibrium populations are unaffected by
thermostat/friction/constraint differences, so this *thermodynamics* check
cleanly validates the OpenMM pipeline against an independent engine. The GROMACS
analysis logic (`scripts/analyze.py`: basin partitions, convergence stats,
literature comparison, J-couplings) is engine-agnostic and carries over with new
I/O adapters (MDTraj instead of `gmx rama`/`gmx energy`).

## Two tiers

| | Tier 1 — explicit (realism) | Tier 2 — implicit (generative) |
|---|---|---|
| solvent | TIP3P, PME, 1.0 nm cutoff | **GBn2**, no box |
| `u(solute)` | not tractable (solvent-marginalized) | **tractable** → BG / energy-guided / reweighting |
| cost | the compute sink (×9 REST2) | cheap → aggressive rare-event sampling |
| role | physical reference | primary substrate for energy-based generative methods |

Shared settings (field-aligned with MDGen/Timewarp): ff14SB, Langevin **310 K**,
**2 fs** with H-bond constraints, low friction (**0.3 ps⁻¹** explicit /
**0.1 ps⁻¹** implicit), **100 fs** save stride (fine enough for the
trajectory-upsampling task; coarser strides are derived by subsampling).

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

Shared featurization (MDGen convention, for cross-comparability): sin/cos of all
torsions (φ,ψ,χ,ω) → **TICA** (deeptime/PyEMMA). Primary distributional metric:
**Jensen–Shannon divergence** on the TICA components and per-DOF marginals.

| task | method class | metric vs ground truth |
|---|---|---|
| latent embedding | reduced-rep / AE | reconstruction error; recovered intrinsic dim; FES/TICA preserved under embedding |
| propagator | transfer operator | implied-timescale error; Chapman–Kolmogorov; rate error |
| equilibrium sampler | Boltzmann generator | JSD on TICA/marginals; reweighted free-energy error (Tier 2, needs `u(x)`) |
| trajectory generation | (video-)diffusion | distributional+kinetic match of generated vs MD: JSD on TICA, implied timescales, n-step transition densities (frame stride is a parameter; 100 fs supports upsampling, ~10–50 ps shows macro-motion) |

Metrics are **distributional/kinetic**, not per-frame: a generated trajectory is
good iff its implied FES/TICA and timescales match the MD's.

## Roadmap

- **Phase 0 — curation layer (engine-agnostic).** Export + reference-observable
  computer + correlation-aware/transferability splits, runnable on the existing
  GROMACS alanine data. Unblocks ML immediately.
- **Phase 1 — OpenMM generation skeleton.** Plain MD, both tiers (ff14SB), the
  Ala/Gly/Pro shake-out. Plus the one-off ff03 cross-engine check against
  GROMACS. *(Ships with this doc: `benchmark/simulate.py`,
  `benchmark/systems.py`.)*
- **Phase 2 — kinetics.** Adaptive-sampling/MSM; ship MSM timescales; ω
  thermodynamics-only.
- **Phase 3 — scale to all 9 dipeptides, both tiers**, with per-DOF convergence
  gating.
- **Phase 4 — define tasks/metrics/splits as a released benchmark.**
- **Phase 5 — tetrapeptides.** Extend to 4-mers for frontier relevance
  (Timewarp/MDGen scale), using the dipeptides as the interpretable floor.
  Requires a sequence-sampling design (which 4-mers) on top of the proven
  pipeline.

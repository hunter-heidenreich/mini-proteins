# Mini-Proteins

A repo of simple molecular dynamics simulations of small proteins with GROMACS.

## Overview

This repo contains a pipeline to generate molecular dynamics (MD) trajectories and **extract atomic forces** for capped dipeptides (often called "mini-proteins"). It is designed to create training datasets for Machine Learning potentials or to study simple folding dynamics.

Key features include:

- **Langevin Dynamics:** Uses stochastic dynamics (SD) for proper canonical sampling (NVT/NPT), replacing the standard velocity-rescaling often found in tutorials. The forces written to the `.trr` are the conservative (force-field) forces, i.e. valid `-∇U` labels.
- **Force Extraction:** Configured to write uncompressed trajectory files (`.trr`), allowing for the extraction of atomic forces essential for force-matching/ML applications. The solute is left **unconstrained** in production so its force labels are complete.
- **Consistent electrostatics:** All stages (EM → NVT → NPT → production) use **PME** at a 1.0 nm cutoff with **rigid TIP3P** water, matching how Amber was parameterized — one Hamiltonian throughout.
- **Diverse Residues:** Extends the standard Alanine Dipeptide model to include bulky (Trp), sulfur-containing (Met), and flexible (Gly) residues.

**Note on Physics:** By default, this pipeline uses the **Amber03** force field and **TIP3P** water. Users intending to use this data for modern production-grade ML models should consider updating the `0_preprocess.sh` script to use newer force fields (e.g., CHARMM36m or Amber ff19SB). The pipeline extracts forces on the **solute only** (rigid water as environment); training an ML potential on explicit-water forces would require making water flexible (`-DFLEXIBLE`) and reducing the timestep to ~0.5 fs.

The scripts are written building off of the [GROMACS tutorial](https://cbp-unitn.gitlab.io/qcb22-23/QCB/tutorial2_gromacs) 
by Luca Tubiana at the University of Trento.

## Setup

The Python environment is managed with [uv](https://docs.astral.sh/uv/)
(`pyproject.toml` + committed `uv.lock`):

```bash
uv sync                 # core: numpy + matplotlib (GROMACS analysis/plots)
uv sync --extra ml      # + the OpenMM benchmark track (openmm, mdtraj, deeptime)
```

You also need a working **GROMACS** (`gmx`) install for the simulation pipeline
below. The pipeline scripts call `uv run` internally, so the commands shown
work directly once `uv sync` has run.

## Usage

### 0. Prepare the simulation structure

The first step is to prepare the simulation structure.
This is done by running the `0_preprocess.sh` script:
```
ID=ala sh scripts/0_preprocess.sh
```
where `ID` is the three-letter amino acid code of the protein to simulate.

This script will:
- Build the protein topology from the PDB file
- Build the box
- Solvate the protein in water 
- Add ions to neutralize the system

Additional parameters can be found at the top of the script.

### 1. Energy Minimization & Equilibration

The next step is to perform energy minimization and equilibration.
This is done by running the `1_equil.sh` script:
```bash
ID=ala sh scripts/1_equil.sh
```
where `ID` is the three-letter amino acid code of the protein to simulate.

This script will:
- Perform energy minimization (using steepest descent, see `config/minim.mdp` for all parameters)
- Equilibrate the system with constant volume (NVT, T=298K, see `config/nvt_langevin.mdp` for all parameters) for 100 ps 
- Equilibrate the system with constant pressure (NPT, T=298K, P=1bar, see `config/npt_langevin.mdp` for all parameters) for 200 ps

Additional parameters can be found at the top of the script.

### 2. Production Simulation

The next step is to run the production simulation.
This is done by running the `2_md_lang.sh` script:
```bash
ID=ala sh scripts/2_md_lang.sh
```
where `ID` is the three-letter amino acid code of the protein to simulate.

This script will:
- Run `NREP` **independent** production replicas (Langevin Dynamics, T=298K, see `config/md_langevin.mdp`), each 100 ns by default
    - Replicas start from the equilibrated NPT structure with fresh Maxwell velocities and independent random seeds, so they diverge immediately. This is plain unbiased MD — the trajectories serve as both ML force-label data and the validation ensemble.
    - Run length and replica count are env-overridable: `NREP=3 NSTEPS=50000000 ID=ala sh scripts/2_md_lang.sh` (default `NREP=5`, `NSTEPS=100000000` ≈ 100 ns).

Additional parameters can be found at the top of the script.

### 3. Post-Process Simulation

This is done by running the `3_post.sh` script:
```bash
ID=ala sh scripts/3_post.sh
```
where `ID` is the three-letter amino acid code of the protein to simulate.

This script will, **for every replica**:
- Extract the forces as an xvg file (the ML training labels)
- Extract the PBC-corrected solute trajectory as a PDB file

and generate sanity plots of the energy minimization, and of potential energy,
total energy, and temperature over the NVT→NPT→production timeline (replicas overlaid).

### 4. Validation Analysis

This is done by running the `4_analyze.sh` script:
```bash
ID=ala sh scripts/4_analyze.sh
```

This script computes backbone dihedrals (`gmx rama`) across all replicas and
writes a single validation + dataset-QC report (`out/${ID}/data/validation.txt`,
plus a machine-readable `validation.json`) designed to be credible to a
computational chemist. It contains:

- **Provenance** — engine version, force field/water, integrator, constraints,
  system size, sampling, seeds, and git commit.
- **Sampling & convergence** — statistical inefficiency `g` and effective
  independent sample count (φ/ψ are autocorrelated, so raw frame counts
  overstate the data), helical↔extended **basin-transition counts** (Hu 2003
  ergodicity check), and a first-half vs second-half stationarity test.
- **Backbone ensemble** — φ/ψ **free-energy surface** (`figs/rama_fes.png`,
  kcal/mol, contoured and basin-labeled), marginals, a cumulative-population
  **convergence plot** (`figs/convergence.png`), and basin populations + per-basin
  free energies under a refined partition.
- **Direct literature comparison** — populations under the published Wang–Duan
  partition printed side-by-side with Vymětal & Vondrášek 2010 (the same
  ff03/TIP3P force field + water model).
- **J-couplings** — ³J(HN,Hα) under several Karplus parameterizations (the choice
  swings ~1 Hz; Best 2008) vs the experimental range for short Ala peptides.
- **Force-label QC** — the actual ML product: |F| distribution, outlier/max-force
  and non-finite checks, and the lag-1 force autocorrelation that justifies the
  1 ps sampling stride (forces decorrelate fast; configurations do not).
- **Thermodynamic stability** — temperature/energy drift slopes (note: `sd` is
  thermostatted, so the check is stationarity, not energy conservation).
- For **glycine**, a (φ,ψ)→(−φ,−ψ) **symmetry index** (→ 0 at convergence).

Comparison references:
Vymětal & Vondrášek 2010 ([ff03 FES](https://doi.org/10.1021/jp100950w)),
Graf 2007 / Best 2008 (J-couplings: [JACS](https://doi.org/10.1021/ja0660406),
[Biophys. J.](https://doi.org/10.1529/biophysj.108.132696)),
Hu 2003 ([Ace-Ala/Gly-Nme maps](https://doi.org/10.1002/prot.10279)),
Duan 2003 ([ff03 definition](https://doi.org/10.1002/jcc.10349)).

**Note on Amber03's Ramachandran balance:** ff03 gives a prominent right-handed
helical basin for the alanine dipeptide (αR + α′ ≈ 42% under the reference's
partition, matching Vymětal & Vondrášek 2010), larger than the near-pure-PPII
picture from peptide NMR. This is the documented
[over-helical](https://doi.org/10.1529/biophysj.108.132696) bias of the Amber
line, not a simulation error — a healthy ff03 run is *expected* to reproduce it.

The refined partition in `analyze.py` adds a sharper detail the coarse
literature box hides: that "helical" density is **C7eq/bridge-dominated**
(ψ ≈ 0), not canonical α-helix (ψ ≈ −45). In our 5 × 100 ns alanine run, true
αR is only ~21% while the bridge basin is ~25% — they sum to the ~46% that a
single broad αR box (like the reference's) reports as one number. Likewise the
backbone ³J(HN,Hα) lands ~1.1 Hz above experiment under the *same* Karplus
parameterization Graf 2007 used (Hu & Bax), i.e. a real ensemble difference
consistent with that helical/bridge excess, not an analysis artifact.

### All-in-one

Alternatively, all of the above steps can be run at once by running the `run.sh` script:
```bash
ID=ala sh scripts/run.sh
```
where `ID` is the three-letter amino acid code of the protein to simulate.

### Running on a GPU / in the cloud

Every MD step accepts a `GPU=1` toggle to offload to a CUDA GPU (default `0` = CPU):
```bash
GPU=1 ID=ala sh scripts/run.sh
```
`scripts/bench.sh` reports production throughput and an estimated `$/protein` so
you can size a rental. Since these systems are small, a consumer GPU (e.g. RTX
4090 on Vast.ai/RunPod) is the cost-optimal choice — roughly a few dollars per
protein. See [`docs/cloud-run.md`](docs/cloud-run.md) for a full setup checklist.

## Included Proteins (And Provenance)

### Alanine Dipeptide

![Alanine Dipeptide](https://hunterheidenreich.com/img/alanine-dipeptide-molecular-dynamics.webp)

- `data/ala.pdb`: 
- Alanine Dipeptide (Ace-Ala-Nme) 
- PubChem CID: 5484387 ([URL](https://pubchem.ncbi.nlm.nih.gov/compound/5484387))
- ATB: [URL](https://atb.uq.edu.au/molecule.py?molid=757282)

### Glycine Dipeptide

![Glycine Dipeptide](https://hunterheidenreich.com/img/glycine-dipeptide-molecular-dynamics.webp)

- `data/gly.pdb`: 
- Glycine Dipeptide (Ace-Gly-Nme) 
- PubChem CID: 439506 ([URL](https://pubchem.ncbi.nlm.nih.gov/compound/2-acetamido-N-methylacetamide))
- ATB: [URL](https://atb.uq.edu.au/molecule.py?molid=22639)

### Isoleucine Dipeptide

![Isoleucine Dipeptide](https://hunterheidenreich.com/img/isoleucine-dipeptide-molecular-dynamics.webp)

- `data/ile.pdb`
- Isoleucine Dipeptide (Ace-Ile-Nme)
- PubChem CID: 7019852 ([URL](https://pubchem.ncbi.nlm.nih.gov/compound/7019852))
- ATB: [URL](https://atb.uq.edu.au/molecule.py?molid=40061)

### Leucine Dipeptide

![Leucine Dipeptide](https://hunterheidenreich.com/img/leucine-dipeptide-molecular-dynamics.webp)

- `data/leu.pdb`
- Leucine Dipeptide (Ace-Leu-Nme)
- PubChem CID: 6950977 ([URL](https://pubchem.ncbi.nlm.nih.gov/compound/6950977))
- ATB: [URL](https://atb.uq.edu.au/molecule.py?molid=40062)

### Methionine Dipeptide

![Methionine Dipeptide](https://hunterheidenreich.com/img/methionine-dipeptide-molecular-dynamics.webp)

- `data/met.pdb`
- Methionine Dipeptide (Ace-Met-Nme)
- PubChem CID: 13875186 ([URL](https://pubchem.ncbi.nlm.nih.gov/compound/13875186))
- ATB: [URL](https://atb.uq.edu.au/molecule.py?molid=618964)
- Contains a thioether side chain (sulfur)

### Phenylalanine Dipeptide

![Phenylalanine Dipeptide](https://hunterheidenreich.com/img/phenylalanine-dipeptide-molecular-dynamics.webp)

- `data/phe.pdb`
- Phenylalanine Dipeptide (Ace-Phe-Nme)
- PubChem CID: 7019860 ([URL](https://pubchem.ncbi.nlm.nih.gov/compound/7019860))
- ATB: [URL](https://atb.uq.edu.au/molecule.py?molid=459390)

### Proline Dipeptide

![Proline Dipeptide](https://hunterheidenreich.com/img/proline-dipeptide-molecular-dynamics.webp)

- `data/pro.pdb`
- Proline Dipeptide (Ace-Pro-Nme)
- PubChem CID: 5245806 ([URL](https://pubchem.ncbi.nlm.nih.gov/compound/5245806))
- ATB: [URL](https://atb.uq.edu.au/molecule.py?molid=1175407)

### Tryptophan Dipeptide

![Tryptophan Dipeptide](https://hunterheidenreich.com/img/tryptophan-dipeptide-molecular-dynamics.webp)

- `data/trp.pdb`
- Tryptophan Dipeptide (Ace-Trp-Nme)
- PubChem CID: 151412 ([URL](https://pubchem.ncbi.nlm.nih.gov/compound/151412))
- ATB: [URL](https://atb.uq.edu.au/molecule.py?molid=40063)

### Valine Dipeptide

![Valine Dipeptide](https://hunterheidenreich.com/img/valine-dipeptide-molecular-dynamics.webp)

- `data/val.pdb`
- Valine Dipeptide (Ace-Val-Nme)
- PubChem CID: 13875188 ([URL](https://pubchem.ncbi.nlm.nih.gov/compound/13875188))
- ATB: [URL](https://atb.uq.edu.au/molecule.py?molid=40060)

## Citation

If you use this repo in your research, please cite:

```
@misc{Heidenreich_Mini-proteins_2023,
author = {Heidenreich, Hunter},
month = sep,
title = {{Mini-proteins}},
url = {https://github.com/hunter-heidenreich/mini-proteins},
year = {2023}
}
```

## License & Attribution

The original work in this repository (the Python analysis script, the Bash pipeline as assembled here, and the curated dipeptide dataset under `data/`) is released under the [MIT License](LICENSE), copyright (c) 2023 Hunter Heidenreich.

Components with their own terms, not covered by that license:

* **Pipeline scripts**: written building off [Luca Tubiana's GROMACS tutorial](https://cbp-unitn.gitlab.io/qcb22-23/QCB/tutorial2_gromacs) (University of Trento); please respect the original author's terms.
* **Force field and water model**: Amber03 and TIP3P are third-party and distributed with GROMACS under their own licenses.
# Mini-Proteins

A repo of simple molecular dynamics simulations of small proteins with GROMACS.

## Overview

This repo contains a pipeline to generate molecular dynamics (MD) trajectories and **extract atomic forces** for capped dipeptides (often called "mini-proteins"). It is designed to create training datasets for Machine Learning potentials or to study simple folding dynamics.

Key features include:

- **Langevin Dynamics:** Uses stochastic dynamics (SD) for proper canonical sampling (NVT/NPT), replacing the standard velocity-rescaling often found in tutorials.
- **Force Extraction:** Configured to write uncompressed trajectory files (`.trr`), allowing for the extraction of atomic forces essential for force-matching/ML applications.
- **Diverse Residues:** Extends the standard Alanine Dipeptide model to include bulky (Trp), sulfur-containing (Met), and flexible (Gly) residues.

**Note on Physics:** By default, this pipeline uses the **Amber03** force field and **TIP3P** water. Users intending to use this data for modern production-grade ML models should consider updating the `0_preprocess.sh` script to use newer force fields (e.g., CHARMM36m or Amber ff19SB).

The scripts are written building off of the [GROMACS tutorial](https://cbp-unitn.gitlab.io/qcb22-23/QCB/tutorial2_gromacs) 
by Luca Tubiana at the University of Trento.

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
- Run the production simulation (Langevin Dynamics, T=298K, see `config/md_langevin.mdp` for all parameters) for 1 ns
    - A full simulation would be much longer, but this is sufficient for a demonstration

Additional parameters can be found at the top of the script.

### 3. Post-Process Simulation

The final step is to post-process the simulation.

This is done by running the `3_post.sh` script:
```bash
ID=ala sh scripts/3_post.sh
```
where `ID` is the three-letter amino acid code of the protein to simulate.

This script will:
- Generate a plot of the potential energy over time
- Generate a plot of the total energy over time
- Generate a plot of the temperature over time
- Extract the trajectory as a PDB file
- Extract the forces as a xvg file

Additional parameters can be found at the top of the script.

### All-in-one

Alternatively, all of the above steps can be run at once by running the `run.sh` script:
```bash
ID=ala sh scripts/run.sh
```
where `ID` is the three-letter amino acid code of the protein to simulate.

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
- Contains a disulfide bond

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
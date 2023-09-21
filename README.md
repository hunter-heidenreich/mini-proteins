# Mini-Proteins

A repo of simple molecular dynamics simulations of small proteins with GROMACS.

## Overview

This repo contains scripts to: 
- Perform energy minimization
- Solvate the protein
- Add ions to neutralize the system
- Equilibrate the system (NVT)
- Equilibrate the system (NPT)
- Run a production simulation
- Post-process the simulation
of a mini-protein using GROMACS.

In this repo, we consider a "mini-protein" to be 
a non-technical designation for a single amino acid residue (or a dipeptide), 
capped with an acetyl group and an N-methyl group.

Frequently, alanine dipeptide (Ace-Ala-Nme) is used as a model system for
protein folding studies.
It's especially enjoyed by machine learning researchers,
because it's small enough to be simulated quickly,
but large enough to exhibit interesting folding behavior.

This repo extends a typical data generation of alanine dipeptide
to include other amino acids.
While not all amino acids are included, these scripts could allow for easy generation
of multiple so-called dipeptide "mini-proteins" for machine learning studies
to add slight diversity to the models considered.

For example, the addition of a disulfide bond in methionine dipeptide
could be used to study the effects of disulfide bonds on protein folding.
Or the addition of a tryptophan residue could be used to study the effects
of aromatic residues on protein folding.
Furthermore, glycine dipeptide could be used to study the effects of
a residue with a small side chain on protein folding, inducing more flexibility. 

The scripts are written building off of the [GROMACS tutorial](https://cbp-unitn.gitlab.io/qcb22-23/QCB/tutorial2_gromacs) 
by Luca Tubiana at the University of Trento.
We make several key deviations: 
- langevin dynamics is used instead of velocity rescaling
- the production simulation is run for a longer time
- the production simulation writes uncompressed trajectory files, 
    which are much larger but allow for force extraction

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
- Equilibrate the system with constant volume (NVT, T=298K, see `config/nvt.mdp` for all parameters) for 100 ps 
- Equilibrate the system with constant pressure (NPT, T=298K, P=1bar, see `config/npt.mdp` for all parameters) for 200 ps

Additional parameters can be found at the top of the script.

### 2. Production Simulation

The next step is to run the production simulation.
This is done by running the `2_prod.sh` script:
```bash
ID=ala sh scripts/2_prod.sh
```
where `ID` is the three-letter amino acid code of the protein to simulate.

This script will:
- Run the production simulation (NVT, T=298K, see `config/prod.mdp` for all parameters) for 1 ns
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

## Included Proteins (And Providence)

### Alanine Dipeptide

![Alanine Dipetide](imgs/ala.gif)

- `data/ala.pdb`: 
- Alanine Dipeptide (Ace-Ala-Nme) 
- PubChem CID: 5484387 ([URL](https://pubchem.ncbi.nlm.nih.gov/compound/5484387))
- ATB: [URL](https://atb.uq.edu.au/molecule.py?molid=757282)

### Glycine Dipeptide

![Glycine Dipetide](imgs/gly.gif)

- `data/gly.pdb`: 
- Glycine Dipeptide (Ace-Gly-Nme) 
- PubChem CID: 439506 ([URL](https://pubchem.ncbi.nlm.nih.gov/compound/2-acetamido-N-methylacetamide))
- ATB: [URL](https://atb.uq.edu.au/molecule.py?molid=22639)

### Isoleucine Dipeptide

![Isoleucine Dipetide](imgs/ile.gif)

- `data/ile.pdb`
- Isoleucine Dipeptide (Ace-Ile-Nme)
- PubChem CID: 7019852 ([URL](https://pubchem.ncbi.nlm.nih.gov/compound/7019852))
- ATB: [URL](https://atb.uq.edu.au/molecule.py?molid=40061)

### Leucine Dipeptide

![Leucine Dipetide](imgs/leu.gif)

- `data/leu.pdb`
- Leucine Dipeptide (Ace-Leu-Nme)
- PubChem CID: 6950977 ([URL](https://pubchem.ncbi.nlm.nih.gov/compound/6950977))
- ATB: [URL](https://atb.uq.edu.au/molecule.py?molid=40062)

### Methionine Dipeptide

![Methionine Dipetide](imgs/met.gif)

- `data/met.pdb`
- Methionine Dipeptide (Ace-Met-Nme)
- PubChem CID: 13875186 ([URL](https://pubchem.ncbi.nlm.nih.gov/compound/13875186))
- ATB: [URL](https://atb.uq.edu.au/molecule.py?molid=618964)
- Contains a disulfide bond

### Phenylalanine Dipeptide

![Phenylalanine Dipetide](imgs/phe.gif)

- `data/phe.pdb`
- Phenylalanine Dipeptide (Ace-Phe-Nme)
- PubChem CID: 7019860 ([URL](https://pubchem.ncbi.nlm.nih.gov/compound/7019860))
- ATB: [URL](https://atb.uq.edu.au/molecule.py?molid=459390)

### Proline Dipeptide

![Proline Dipetide](imgs/pro.gif)

- `data/pro.pdb`
- Proline Dipeptide (Ace-Pro-Nme)
- PubChem CID: 5245806 ([URL](https://pubchem.ncbi.nlm.nih.gov/compound/5245806))
- ATB: [URL](https://atb.uq.edu.au/molecule.py?molid=1175407)

### Tryptophan Dipeptide

![Tryptophan Dipetide](imgs/trp.gif)

- `data/trp.pdb`
- Tryptophan Dipeptide (Ace-Trp-Nme)
- PubChem CID: 151412 ([URL](https://pubchem.ncbi.nlm.nih.gov/compound/151412))
- ATB: [URL](https://atb.uq.edu.au/molecule.py?molid=40063)

### Valine Dipeptide

![Valine Dipetide](imgs/val.gif)

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
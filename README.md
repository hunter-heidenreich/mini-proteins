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

## Included Proteins (And Providence)

### Alanine Dipeptide

- `data/ala.pdb`: 
- Alanine Dipeptide (Ace-Ala-Nme) 
- PubChem CID: 5484387 ([URL](https://pubchem.ncbi.nlm.nih.gov/compound/5484387))
- ATB: [URL](https://atb.uq.edu.au/molecule.py?molid=757282)

### Glycine Dipeptide

- `data/gly.pdb`: 
- Glycine Dipeptide (Ace-Gly-Nme) 
- PubChem CID: 439506 ([URL](https://pubchem.ncbi.nlm.nih.gov/compound/2-acetamido-N-methylacetamide))
- ATB: [URL](https://atb.uq.edu.au/molecule.py?molid=22639)

### Isoleucine Dipeptide

- `data/ile.pdb`
- Isoleucine Dipeptide (Ace-Ile-Nme)
- PubChem CID: 7019852 ([URL](https://pubchem.ncbi.nlm.nih.gov/compound/7019852))
- ATB: [URL](https://atb.uq.edu.au/molecule.py?molid=40061)

### Leucine Dipeptide

- `data/leu.pdb`
- Leucine Dipeptide (Ace-Leu-Nme)
- PubChem CID: 6950977 ([URL](https://pubchem.ncbi.nlm.nih.gov/compound/6950977))
- ATB: [URL](https://atb.uq.edu.au/molecule.py?molid=40062)

### Methionine Dipeptide

- `data/met.pdb`
- Methionine Dipeptide (Ace-Met-Nme)
- PubChem CID: 13875186 ([URL](https://pubchem.ncbi.nlm.nih.gov/compound/13875186))
- ATB: [URL](https://atb.uq.edu.au/molecule.py?molid=618964)
- Contains a disulfide bond

### Phenylalanine Dipeptide

- `data/phe.pdb`
- Phenylalanine Dipeptide (Ace-Phe-Nme)
- PubChem CID: 7019860 ([URL](https://pubchem.ncbi.nlm.nih.gov/compound/7019860))
- ATB: [URL](https://atb.uq.edu.au/molecule.py?molid=459390)

### Proline Dipeptide

- `data/pro.pdb`
- Proline Dipeptide (Ace-Pro-Nme)
- PubChem CID: 5245806 ([URL](https://pubchem.ncbi.nlm.nih.gov/compound/5245806))
- ATB: [URL](https://atb.uq.edu.au/molecule.py?molid=1175407)

### Tryptophan Dipeptide

- `data/trp.pdb`
- Tryptophan Dipeptide (Ace-Trp-Nme)
- PubChem CID: 151412 ([URL](https://pubchem.ncbi.nlm.nih.gov/compound/151412))
- ATB: [URL](https://atb.uq.edu.au/molecule.py?molid=40063)

### Valine Dipeptide

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
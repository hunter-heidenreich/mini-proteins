"""Per-system metadata for the capped-dipeptide benchmark suite.

The benchmark's ground truth must cover each system's *slow* degrees of freedom,
which for several residues extend beyond the backbone (phi, psi) into side-chain
chi angles, and for proline into omega (cis/trans) and ring pucker. This table
is the single source of truth for:
  - which dihedrals to track in convergence / reference-observable analysis,
  - the difficulty ordering of the suite,
  - per-residue notes that gate run length (slowest DOF) and scope (rare events).

Chi-angle atom names follow the IUPAC/Amber side-chain convention and match the
ACE-X-NME PDBs in data/. omega is the C(i-1)-N-CA-C... here the Ace-X peptide
bond, measured CH3(ACE)-C(ACE)-N(X)-CA(X).
"""

# difficulty: rough intrinsic-dimension / hardness rank (1 = easiest)
SYSTEMS = {
    "gly": {
        "resname": "GLY", "difficulty": 1,
        "slow_dofs": ["phi", "psi"],
        "chi": [],
        "rare": [],
        "notes": "symmetric; samples 4 Ramachandran corners; "
                 "(phi,psi)->(-phi,-psi) symmetry is a free correctness check",
    },
    "ala": {
        "resname": "ALA", "difficulty": 1,
        "slow_dofs": ["phi", "psi"],
        "chi": [],
        "rare": ["alphaL"],
        "notes": "canonical 2D Ramachandran; calibration reference; "
                 "alphaL is barrier-limited (needs enhanced sampling)",
    },
    "val": {
        "resname": "VAL", "difficulty": 2,
        "slow_dofs": ["phi", "psi", "chi1"],
        "chi": [["N", "CA", "CB", "CG1"]],
        "rare": ["alphaL", "chi1_flip"],
        "notes": "beta-branched; chi1 couples to backbone",
    },
    "leu": {
        "resname": "LEU", "difficulty": 3,
        "slow_dofs": ["phi", "psi", "chi1", "chi2"],
        "chi": [["N", "CA", "CB", "CG"], ["CA", "CB", "CG", "CD1"]],
        "rare": ["alphaL", "chi_flip"],
        "notes": "flexible aliphatic side chain",
    },
    "ile": {
        "resname": "ILE", "difficulty": 3,
        "slow_dofs": ["phi", "psi", "chi1", "chi2"],
        "chi": [["N", "CA", "CB", "CG1"], ["CA", "CB", "CG1", "CD1"]],
        "rare": ["alphaL", "chi_flip"],
        "notes": "beta-branched; two coupled side-chain torsions",
    },
    "phe": {
        "resname": "PHE", "difficulty": 3,
        "slow_dofs": ["phi", "psi", "chi1", "chi2"],
        "chi": [["N", "CA", "CB", "CG"], ["CA", "CB", "CG", "CD1"]],
        "rare": ["alphaL", "ring_flip"],
        "notes": "aromatic; chi2 ring flip is 2-fold symmetric",
    },
    "met": {
        "resname": "MET", "difficulty": 4,
        "slow_dofs": ["phi", "psi", "chi1", "chi2", "chi3"],
        "chi": [["N", "CA", "CB", "CG"], ["CA", "CB", "CG", "SD"],
                ["CB", "CG", "SD", "CE"]],
        "rare": ["alphaL", "chi_flip"],
        "notes": "long floppy thioether side chain; highest side-chain dimension",
    },
    "trp": {
        "resname": "TRP", "difficulty": 4,
        "slow_dofs": ["phi", "psi", "chi1", "chi2"],
        "chi": [["N", "CA", "CB", "CG"], ["CA", "CB", "CG", "CD1"]],
        "rare": ["alphaL", "chi_flip"],
        "notes": "largest side chain; slow chi relaxation -> longest runs",
    },
    "pro": {
        "resname": "PRO", "difficulty": 5,
        "slow_dofs": ["phi", "psi", "omega", "pucker"],
        "chi": [],
        "rare": ["cis_trans_omega", "pucker"],
        "notes": "restricted phi; ring pucker (endo/exo) + cis/trans omega; "
                 "different topology; omega rate is OUT OF SCOPE (seconds) -- "
                 "reported as thermodynamics only",
    },
}

# Backbone dihedral atom specs (Ace-X-Nme). phi/psi reference the central
# residue X; omega is the Ace-X peptide bond.
BACKBONE = {
    "phi": ["C", "N", "CA", "C"],      # C(ACE)-N(X)-CA(X)-C(X)
    "psi": ["N", "CA", "C", "N"],      # N(X)-CA(X)-C(X)-N(NME)
    "omega": ["CH3", "C", "N", "CA"],  # CH3(ACE)-C(ACE)-N(X)-CA(X)
}

# Shake-out trio for Phase 1: canonical / symmetric / ring+cis-trans.
SHAKEOUT = ["ala", "gly", "pro"]

# Suggested transferability split (Timewarp-style): train on a diverse subset,
# test generalization to held-out residues.
TRANSFER_SPLIT = {
    "train": ["ala", "gly", "val", "leu", "phe"],
    "test": ["ile", "met", "trp", "pro"],
}

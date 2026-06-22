"""Generate a per-residue README manifest for the OpenMM ff14SB benchmark data
uploaded to HuggingFace (hheiden/mini-proteins-bench/<id>/openmm-ff14sb/).
Writes out/<id>/openmm-ff14sb-README.md for each residue. Run in the omm env so
benchmark.systems is importable."""
import json
import os

from benchmark.systems import SYSTEMS

GIT = os.popen("git rev-parse --short HEAD").read().strip()
RESID = ["ala", "gly", "pro"]

HEADER = """# {resname} ({rid}) — OpenMM ff14SB Phase-1 shake-out

Capped dipeptide Ace-{rid_up}-Nme. Unbiased Langevin MD generated with
`benchmark/simulate.py`; curated with `benchmark/curate.py`. Part of
`mini-proteins-bench` (see the repo for the certified GROMACS ff03 track under
`{rid}/gromacs-ff03/`). Force field **amber14/ff14SB** for comparability with
MDGen / Timewarp / Transferable-BG.

## Layout
- `explicit/curated.npz` + `curation.json` — TIP3P + PME (1.0 nm cutoff), realism tier
- `implicit/curated.npz` + `curation.json` — GBn2, generative tier (tractable u(x))

## Generation (both tiers)
- 3 independent replicas (seeds 1–3, distinct `gen` per replica), 20 ns each + 0.5 ns equil
- 310 K, 2 fs timestep, H-bond constraints (heavy atoms flexible), rigid water (explicit)
- Langevin friction: 0.3 ps⁻¹ (explicit) / 0.1 ps⁻¹ (implicit)
- Save stride 100 fs (0.1 ps) → 200k frames/replica
- OpenMM 8.2.0, CUDA. Provenance (seeds, versions, git) in each `meta.json` / `curation.json`.

## curated.npz contents
`angles_deg` (torsions, deg), `features` (sin/cos of each torsion), `replica_id`,
`time_ps`, `torsion_names`, `fes` + `fes_edges` (φ/ψ free energy, kcal/mol),
`split='by_replica'` (correlation-aware; never random-split frames).

## Slow degrees of freedom (this residue)
{slow}

{note}
_Generated from git {git}._
"""

PRO_NOTE = (
    "## ⚠️ Kinetic ceiling (proline)\n"
    "Proline's φ is ring-locked (statistical inefficiency g≈1.5, near-iid) and its\n"
    "defining slow process — **cis/trans ω isomerization — is seconds-timescale** and\n"
    "is NOT sampled by 20 ns plain MD. This data is a **trans-locked, PPII-dominant\n"
    "thermodynamic** ensemble; it is intentionally *not* kinetically complete for ω.\n"
    "The helical/extended transition counts in `curation.json` use an alanine-centric\n"
    "state definition and are not meaningful for proline (they read 0).\n"
)

for rid in RESID:
    s = SYSTEMS[rid]
    slow = "\n".join(f"- {d}" for d in s["slow_dofs"])
    note = PRO_NOTE if rid == "pro" else ""
    txt = HEADER.format(resname=s["resname"], rid=rid, rid_up=s["resname"],
                        slow=slow, note=note, git=GIT)
    out = f"out/{rid}/openmm-ff14sb-README.md"
    with open(out, "w") as fh:
        fh.write(txt)
    print("wrote", out)

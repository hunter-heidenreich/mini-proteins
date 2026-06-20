"""Validation analysis for plain (unbiased) long MD of capped dipeptides.

Reads the per-replica gmx rama output and produces:
  - a phi/psi free-energy surface (FES = -kT ln P), figs/rama_fes.png
  - phi/psi marginal distributions, figs/rama_marginals.png
  - basin populations with per-replica mean +/- std (convergence/error bars)
  - for glycine: a (phi,psi)->(-phi,-psi) symmetry index (must -> 0 at convergence)
  - Karplus 3J(HN,Ha) predicted from phi, for comparison with NMR experiment

Comparison targets (plain MD at convergence is method-independent, so these
enhanced-sampling / experimental references still apply):
  - phi/psi FES & populations across force fields incl. Amber ff03:
      Vymetal & Vondrasek (2010) J. Phys. Chem. B, doi:10.1021/jp100950w
  - experimental 3J couplings for alanine peptides:
      Graf, Nguyen, Stock & Schwalbe (2007) JACS, doi:10.1021/ja0660406
      Best, Buchete & Hummer (2008) Biophys. J., doi:10.1529/biophysj.108.132696
  - Ace-Ala-Nme / Ace-Gly-Nme reference maps in water:
      Hu, Elstner & Hermans (2003) Proteins, doi:10.1002/prot.10279

Karplus coefficients for 3J(HN,Ha) are from Hu & Bax (1997) JACS 119:6360,
with theta = phi - 60 deg. Other parameter sets exist; swap KARPLUS if needed.
"""

import glob
import os

import matplotlib.pyplot as plt
import numpy as np

R_KJ = 8.314462618e-3      # kJ/(mol K)
TEMP = 298.0               # K (matches ref_t in the .mdp files)
KT = R_KJ * TEMP           # ~2.478 kJ/mol
KARPLUS = (7.09, -1.42, 1.55)   # (A, B, C) for 3J(HN,Ha), Hu & Bax 1997
CAPS = {"ACE", "NME", "NAC", "NH2", "NHE"}

# Approximate Ramachandran basins (degrees). Boundaries are conventional and
# coarse; the FES figure is the primary artifact, these just give populations.
# Checked in order; first match wins, otherwise 'other'.
#
# alphaR is kept narrow (psi ~ -120..-20) so it counts *canonical* right-handed
# helix only. The psi~0 density at helical phi -- the alphaR/C7eq "bridge" that
# lies between alphaR and PPII -- gets its own region instead of inflating the
# helix count (an earlier, broad alphaR box reached up to psi=+45 and absorbed
# it). alphaR and alphaL are exact mirrors so the glycine symmetry check (which
# compares their populations) stays fair; bridge is the L-residue region only.
REGIONS = [
    ("alphaR", lambda p, s: -160 <= p <= -20 and -120 <= s <= -20),
    ("alphaL", lambda p, s: 20 <= p <= 160 and 20 <= s <= 120),
    ("PPII", lambda p, s: -110 <= p <= -20 and (s > 90 or s < -160)),
    ("beta", lambda p, s: -180 <= p < -110 and (s > 90 or s < -160)),
    ("bridge", lambda p, s: -160 <= p <= -20 and -20 < s <= 90),
]


def read_rama(filename):
    """Return (phi, psi) arrays in degrees for the central (non-cap) residue."""
    by_res = {}
    with open(filename) as f:
        for line in f:
            line = line.strip()
            if not line or line[0] in ("#", "@", "&"):
                continue
            parts = line.split()
            if len(parts) < 3:
                continue
            phi, psi, label = float(parts[0]), float(parts[1]), parts[2]
            by_res.setdefault(label, []).append((phi, psi))
    # Prefer a residue that is not a cap; fall back to the most-populated label.
    candidates = {k: v for k, v in by_res.items()
                  if k.split("-")[0].upper() not in CAPS}
    if not candidates:
        candidates = by_res
    label = max(candidates, key=lambda k: len(candidates[k]))
    arr = np.array(candidates[label], dtype=float)
    return arr[:, 0], arr[:, 1]


def classify(phi, psi):
    for name, test in REGIONS:
        if test(phi, psi):
            return name
    return "other"


def populations(phi, psi):
    names = [n for n, _ in REGIONS] + ["other"]
    counts = dict.fromkeys(names, 0)
    for p, s in zip(phi, psi):
        counts[classify(p, s)] += 1
    n = len(phi)
    return {k: counts[k] / n for k in names}


_ID = os.environ.get("ID")
print(f"=== validation analysis: {_ID} ===")

RAW = f"out/{_ID}/raw"
FIGS = f"out/{_ID}/figs"
DATA = f"out/{_ID}/data"
os.makedirs(FIGS, exist_ok=True)
os.makedirs(DATA, exist_ok=True)

rama_files = sorted(glob.glob(f"{RAW}/md_lang_r*_rama.xvg"))
if not rama_files:
    raise SystemExit(f"no rama files found in {RAW} (run 4_analyze.sh after 2/3)")

per_rep_phi, per_rep_psi = [], []
for path in rama_files:
    phi, psi = read_rama(path)
    per_rep_phi.append(phi)
    per_rep_psi.append(psi)

all_phi = np.concatenate(per_rep_phi)
all_psi = np.concatenate(per_rep_psi)
print(f"replicas: {len(rama_files)}   frames total: {all_phi.size}")

# ---- Free-energy surface -------------------------------------------------
edges = np.linspace(-180, 180, 73)  # 5-degree bins
H, _, _ = np.histogram2d(all_phi, all_psi, bins=[edges, edges])
P = H / H.sum()
with np.errstate(divide="ignore"):
    F = -KT * np.log(P)
F -= np.nanmin(F[np.isfinite(F)])
F[~np.isfinite(F)] = np.nan

plt.figure(figsize=(5.5, 4.5))
im = plt.pcolormesh(edges, edges, F.T, cmap="viridis_r", shading="auto",
                    vmin=0, vmax=30)
plt.colorbar(im, label="free energy (kJ/mol)")
plt.xlabel(r"$\phi$ (deg)")
plt.ylabel(r"$\psi$ (deg)")
plt.title(f"{_ID}: Ramachandran FES")
plt.xlim(-180, 180)
plt.ylim(-180, 180)
plt.xticks(range(-180, 181, 90))
plt.yticks(range(-180, 181, 90))
plt.tight_layout()
plt.savefig(f"{FIGS}/rama_fes.png", dpi=150)
plt.close()

# ---- Marginals -----------------------------------------------------------
fig, ax = plt.subplots(1, 2, figsize=(8, 3.2))
ax[0].hist(all_phi, bins=edges, density=True)
ax[0].set_xlabel(r"$\phi$ (deg)")
ax[0].set_ylabel("density")
ax[1].hist(all_psi, bins=edges, density=True)
ax[1].set_xlabel(r"$\psi$ (deg)")
fig.suptitle(f"{_ID}: backbone dihedral marginals")
fig.tight_layout()
fig.savefig(f"{FIGS}/rama_marginals.png", dpi=150)
plt.close(fig)

# ---- Basin populations with per-replica error bars -----------------------
rep_pops = [populations(p, s) for p, s in zip(per_rep_phi, per_rep_psi)]
names = [n for n, _ in REGIONS] + ["other"]

lines = [f"Validation analysis for {_ID}",
         f"replicas={len(rama_files)}  frames={all_phi.size}  T={TEMP} K", "",
         "Basin populations (mean +/- std across replicas):"]
for name in names:
    vals = np.array([rp[name] for rp in rep_pops])
    lines.append(f"  {name:8s} {vals.mean():6.3f} +/- {vals.std():.3f}")

# ---- Glycine symmetry check ---------------------------------------------
if (_ID or "").lower() == "gly":
    Hm = H[::-1, ::-1]  # mirror (phi,psi)->(-phi,-psi)
    denom = (H + Hm).sum()
    asym = np.abs(H - Hm).sum() / denom if denom else float("nan")
    aR = np.array([rp["alphaR"] for rp in rep_pops])
    aL = np.array([rp["alphaL"] for rp in rep_pops])
    lines += ["",
              "Glycine symmetry (should -> 0 / equal at convergence):",
              f"  asymmetry index           {asym:.3f}",
              f"  alphaR pop {aR.mean():.3f}  vs  alphaL pop {aL.mean():.3f}"]

# ---- Karplus 3J(HN,Ha) ---------------------------------------------------
A, B, C = KARPLUS
theta = np.deg2rad(all_phi - 60.0)
J = A * np.cos(theta) ** 2 + B * np.cos(theta) + C
rep_J = []
for phi in per_rep_phi:
    th = np.deg2rad(phi - 60.0)
    rep_J.append((A * np.cos(th) ** 2 + B * np.cos(th) + C).mean())
rep_J = np.array(rep_J)
lines += ["",
          "Karplus 3J(HN,Ha) [Hz] (Hu & Bax 1997; compare to NMR, Graf 2007):",
          f"  <3J> {rep_J.mean():.2f} +/- {rep_J.std():.2f}"]

report = "\n".join(lines)
print(report)
with open(f"{DATA}/validation.txt", "w") as f:
    f.write(report + "\n")
print(f"\nwrote {DATA}/validation.txt, {FIGS}/rama_fes.png, {FIGS}/rama_marginals.png")

"""ff03 cross-engine check: compare OpenMM ff03 alanine basin populations to the
certified GROMACS ff03 reference (results/ala/validation.json). A thermodynamics
check -- thermostat/friction differences should not matter; the equilibrium
populations should agree within sampling error.

Reads phi/psi from a curated.npz (default the ff03 run's) and applies the SAME
refined basin partition as scripts/analyze.py (lines 108-112).

  micromamba run -n omm python scripts/_ff03_compare.py [curated.npz]
"""
import json
import sys

import numpy as np

# Refined partition, copied verbatim from scripts/analyze.py REFINED_PARTITION.
REFINED = [
    ("alphaR", lambda p, s: -160 <= p <= -20 and -120 <= s <= -20),
    ("alphaL", lambda p, s: 20 <= p <= 160 and 20 <= s <= 120),
    ("PPII",   lambda p, s: -110 <= p <= -20 and (s > 90 or s < -160)),
    ("beta",   lambda p, s: -180 <= p < -110 and (s > 90 or s < -160)),
    ("bridge", lambda p, s: -160 <= p <= -20 and -20 < s <= 90),
]
NAMES = [n for n, _ in REFINED] + ["other"]


def populations(phi, psi):
    counts = {n: 0 for n in NAMES}
    for p, s in zip(phi, psi):
        for n, fn in REFINED:
            if fn(p, s):
                counts[n] += 1
                break
        else:
            counts["other"] += 1
    tot = len(phi)
    return {n: counts[n] / tot for n in NAMES}


def main():
    npz = sys.argv[1] if len(sys.argv) > 1 else "out/ala/explicit/curated.npz"
    d = np.load(npz, allow_pickle=True)
    names = list(d["torsion_names"])
    ang = d["angles_deg"]
    phi = ang[:, names.index("phi")]
    psi = ang[:, names.index("psi")]
    pops = populations(phi, psi)

    ref = json.load(open("results/ala/validation.json"))["populations_refined"]
    print(f"ff03 cross-engine check  ({npz}, {len(phi)} frames)")
    print("=" * 64)
    print(f"{'basin':8s} {'OpenMM ff03':>12s} {'GROMACS ff03':>14s} "
          f"{'delta':>8s} {'GROMACS sigma':>10s}")
    maxdev = 0.0
    for n in NAMES:
        omm = pops[n]
        g = ref.get(n, {}).get("mean")
        gs = ref.get(n, {}).get("std")
        if g is None:
            print(f"{n:8s} {omm:12.4f} {'n/a':>14s}")
            continue
        delta = omm - g
        nsig = abs(delta) / gs if gs else float("nan")
        flag = "" if (gs and abs(delta) <= 3 * gs) else "  <-- >3sigma"
        maxdev = max(maxdev, nsig if gs else 0)
        print(f"{n:8s} {omm:12.4f} {g:14.4f} {delta:+8.4f} "
              f"{nsig:8.1f}sig{flag}")
    print("=" * 64)
    verdict = ("PASS: all basins within 3 GROMACS sigma"
               if maxdev <= 3 else
               f"CHECK: max deviation {maxdev:.1f} sigma (single replica; "
               "some sampling spread expected)")
    print(verdict)


if __name__ == "__main__":
    main()

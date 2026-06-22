# ff03 cross-engine check — OpenMM vs certified GROMACS (alanine)

Date 2026-06-22. A *thermodynamics* sanity check: does OpenMM reproduce the
GROMACS ff03/TIP3P alanine result (`results/ala/validation.json`)? Run:

```
benchmark.simulate --residue ala --tier explicit --forcefield amber03.xml tip3p.xml \
  --temperature 298 --friction 10 --ns 20 --replica 1   # then benchmark.curate
scripts/_ff03_compare.py   # applies analyze.py's refined basin partition
```

## Verdict: engines agree — kinetics cleanly, thermodynamics within error

**Kinetics — clean pass.** At matched friction (10 ps⁻¹) the OpenMM ff03 MSM slow
φ/ψ timescale is **300 ps**, same order as the GROMACS statistical inefficiency
g≈441 ps. The two engines agree on the slow process.

**ff03 signature reproduced.** αR+bridge ("helix region") dominant (0.57 vs 0.46),
αL dead (~0), bridge prominent — the known ff03 over-helical/bridge-heavy
character the GROMACS validation flagged. Same physics.

**Thermodynamics — consistent but underpowered.** Single 20 ns @ friction 10 gave
only **24 effective ψ samples** (ψ decorrelates in 474 ps at that friction). Every
basin agrees with GROMACS within **~1σ of the run's own sampling error** (≈0.09 per
population):

| basin | OpenMM ff03 | GROMACS ff03 | Δ | run's own 1σ | deviation |
|---|---|---|---|---|---|
| αR | 0.262 | 0.206 | +0.056 | 0.090 | 0.6σ |
| PPII | 0.260 | 0.337 | −0.077 | 0.090 | 0.9σ |
| β | 0.143 | 0.176 | −0.034 | 0.071 | 0.5σ |
| bridge | 0.305 | 0.252 | +0.053 | 0.094 | 0.6σ |
| αL | 0.000 | 0.001 | −0.000 | — | ~0σ |

(The `_ff03_compare.py` script prints 4–5σ because it divides by the GROMACS
*converged 5-replica* σ ≈0.01 — the wrong yardstick for a single short replica.
Against the run's own error, everything is <1σ.)

## Methodology fix for a conclusive thermo check

`--friction 10` is right for a *kinetics* match but **wrong for the populations/FES
check** it is labelled as: equilibrium populations are friction-independent, so the
thermo check should use *low* friction (e.g. 0.3–1 ps⁻¹) to sample 10–30× more per
wall-second, and/or multiple replicas. friction=10 crippled the ψ statistics
(24 eff samples). Re-run at low friction for a tight thermodynamic comparison.

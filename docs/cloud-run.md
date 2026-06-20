# Running on a rented GPU box

Goal: validate **alanine (`ala`)** end-to-end on a cheap cloud GPU, confirm the
simulation against the literature, *then* decide whether to run the rest.

This system is tiny (~2,300 atoms), so a consumer GPU (e.g. **RTX 4090** on
Vast.ai or RunPod) is the cost-optimal choice — far cheaper than datacenter
A100/H100 and not meaningfully slower at this size. CPU-only is *more*
expensive here because small systems scale poorly across cores.

## 1. Rent a box

- **Vast.ai or RunPod**, pick a single **RTX 4090** (or 3090/A4000 if cheaper).
- Use a **CUDA image** (RunPod "PyTorch/CUDA" template, or any `nvidia/cuda`
  base). Free egress on both providers — relevant because the dataset is
  ~25–30 GB/protein of `.trr`.

## 2. Get a GPU build of GROMACS

Easiest is the NGC container (already built with CUDA):

```bash
# inside the instance (Docker available on most templates)
docker run --gpus all -it -v $PWD:/work nvcr.io/hpc/gromacs:2024  # check NGC for the latest tag
```

Or build/install directly if you prefer (conda-forge `gromacs` CUDA variant, or
source build with `-DGMX_GPU=CUDA`). Confirm GPU is visible:

```bash
gmx --version | grep -i gpu   # should say "GPU support: CUDA"
nvidia-smi
```

## 3. Clone the repo + Python deps

```bash
git clone https://github.com/hunter-heidenreich/mini-proteins.git
cd mini-proteins
python3 -m venv .venv && . .venv/bin/activate
pip install -r requirements.txt        # numpy, matplotlib (for plots/analysis)
```

## 4. Benchmark first (≈ a few minutes, ~$0.05)

Equilibrate, then time the production stage to get a real ns/day and lock in the
per-protein cost:

```bash
GPU=1 ID=ala sh scripts/0_preprocess.sh
GPU=1 ID=ala sh scripts/1_equil.sh
GPU=1 PRICE=0.40 ID=ala sh scripts/bench.sh   # set PRICE to your $/hr
```

`bench.sh` prints throughput and `$X/protein`. If that's comfortably under
budget, proceed; if not, reduce `NREP` or per-replica length.

## 5. Full alanine run

`bench.sh` already did EM+NVT+NPT, so just run production → post → analysis.
Defaults are 5 replicas × 100 ns:

```bash
GPU=1 ID=ala sh scripts/2_md_lang.sh      # or: NREP=5 NSTEPS=100000000
GPU=1 ID=ala sh scripts/3_post.sh
ID=ala sh scripts/4_analyze.sh
```

## 6. Confirm the simulation is correct

Check `out/ala/figs/rama_fes.png` and `out/ala/data/validation.txt` against the
references in `scripts/analyze.py` (Vymětal & Vondrášek 2010 for the ff03 FES;
Graf 2007 / Best 2008 for ³J couplings). The report prints the populations
side-by-side with the Vymětal ff03/TIP3P numbers, so agreement is read directly.
Expected for **Amber03** alanine in water: a prominent right-handed helical basin
(their αR + α′ ≈ 42%) co-dominant with the PPII-like β basin (≈ 41%), an extended
C5 minority (≈ 17%), and αL essentially unsampled (≈ 0.1% — a documented ff03
trait plus a plain-MD limit, not a bug). Note this is *more* helical than the
near-pure-PPII picture from peptide NMR; that gap is Amber's known over-helical
bias, expected for a correct ff03 run. The convergence section (effective sample
count, basin-transition counts, half/half stationarity) and replica-to-replica
agreement are the checks that the run itself is sound.

## 7. Pull results

```bash
tar czf ala-data.tgz out/ala/data out/ala/figs
# then scp/rsync ala-data.tgz down (or use the provider's file browser)
```

Only after alanine checks out do we scale to gly + the rest (same commands,
different `ID`).

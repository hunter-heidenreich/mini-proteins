# results/

Versioned, durable home for the **small certified artifacts** — the project's
actual deliverable — so they survive the ephemeral RunPod box and stay in git
history.

Populated by `sh scripts/pull_results.sh` (run on the box before stopping the
pod), which copies from the gitignored `out/` tree:

```
results/<id>/validation.{txt,json}        # GROMACS validation (scripts/analyze.py)
results/<id>/<tier>/curation.{txt,json}   # benchmark curation (benchmark/curate.py)
results/<id>/figs/*.png                    # FES, convergence, marginals
```

## What does NOT go here

- **Raw / intermediate trajectories** (`.trr`, `.dcd`, `.tpr`, logs) — large and
  regenerable from code + the seeds in `meta.json`. They stay in `out/`
  (gitignored) on the pod; keep only if you push them to HuggingFace.
- **`curated.npz`** (the ML-ready arrays) — too big for git. → **HuggingFace.**

See `docs/cloud-run.md` ("Before you stop the pod") and `docs/benchmark.md`
(storage tiers).

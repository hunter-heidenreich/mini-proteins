#!/bin/sh
# Copy the small, durable artifacts out of the gitignored, ephemeral out/ tree
# into the versioned results/ dir, so the certified reports survive the (RunPod)
# pod and get committed to git. Run this BEFORE stopping the pod, then
# `git add results && git commit && git push`.
#
# Big data (curated.npz, raw trajectories) does NOT go here -- it goes to
# HuggingFace. See docs/cloud-run.md and docs/benchmark.md.

mkdir -p results

# GROMACS validation reports (scripts/analyze.py): out/<id>/data/validation.*
for f in out/*/data/validation.txt out/*/data/validation.json; do
    [ -f "$f" ] || continue
    id=$(echo "$f" | cut -d/ -f2)
    mkdir -p "results/${id}"
    cp "$f" "results/${id}/"
done

# Benchmark curation reports (benchmark/curate.py): out/<id>/<tier>/curation.*
for f in out/*/*/curation.txt out/*/*/curation.json; do
    [ -f "$f" ] || continue
    id=$(echo "$f" | cut -d/ -f2)
    tier=$(echo "$f" | cut -d/ -f3)
    mkdir -p "results/${id}/${tier}"
    cp "$f" "results/${id}/${tier}/"
done

# Figures (small PNGs): out/<id>/figs/*.png
for d in out/*/figs; do
    [ -d "$d" ] || continue
    id=$(echo "$d" | cut -d/ -f2)
    mkdir -p "results/${id}/figs"
    cp "$d"/*.png "results/${id}/figs/" 2>/dev/null || true
done

echo "Pulled small reports into results/ (commit these)."
echo "Big data (curated.npz, trajectories) -> HuggingFace, not git."

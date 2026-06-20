#!/bin/sh
# Phase-1 shake-out: plain MD for the Ala/Gly/Pro trio in both solvent tiers.
# Run from the repo root. Each replica gets a distinct seed for independence.
#
# Env overrides:
#   NREP   replicas per (residue,tier)        (default 3)
#   NS     production ns per replica           (default 20; bump to 100 for real)
#   TIERS  space-separated tiers               (default "implicit explicit")
#   PLATFORM  OpenMM platform                  (default CUDA)
#
# Cross-validation check: compare out/ala/explicit/* against the GROMACS Ala
# FES/populations (scripts/analyze.py logic) -- they should agree.

NREP=${NREP:-3}
NS=${NS:-20}
TIERS=${TIERS:-"implicit explicit"}
PLATFORM=${PLATFORM:-CUDA}

for RES in ala gly pro; do
    for TIER in ${TIERS}; do
        REP=1
        while [ ${REP} -le ${NREP} ]; do
            echo "=== ${RES} / ${TIER} / replica ${REP}/${NREP} (${NS} ns) ==="
            uv run --extra ml python -m benchmark.simulate \
                --residue ${RES} --tier ${TIER} --replica ${REP} \
                --ns ${NS} --seed ${REP} --platform ${PLATFORM}
            REP=$((REP + 1))
        done
    done
done

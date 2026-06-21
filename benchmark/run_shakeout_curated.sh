#!/bin/sh
# Disk-safe Phase-1 shake-out: generate -> curate -> prune, per (residue, tier).
#
# Unlike run_shakeout.sh (which only generates), this curates each (residue,tier)
# trio as soon as its replicas finish and then DELETES the raw traj.dcd files,
# keeping only the durable curated.npz + reports (+ small scalars/state/meta).
# Peak disk stays ~one trio of DCDs (~9 GB) instead of the full ~28 GB, so the
# whole shake-out fits a 50 GB volume with headroom. DCDs are scratch anyway --
# regenerable from the seeds recorded in each meta.json.
#
# Env overrides:
#   NS     production ns per replica   (default 20)
#   NREP   replicas per (res,tier)     (default 3)
#   RESID  residues                    (default "ala gly pro")
#   TIERS  tiers                       (default "implicit explicit")
#   KEEP_DCD=1  skip the prune step (keep raw trajectories)
#
# Runtime ~9.5 GPU-hr on CUDA for the default 18-run shake-out. Run from repo root:
#   MAMBA_ROOT_PREFIX=/root/micromamba sh benchmark/run_shakeout_curated.sh
export MAMBA_ROOT_PREFIX=${MAMBA_ROOT_PREFIX:-/root/micromamba}
PY=${PY:-"/root/bin/micromamba run -n omm python"}
NS=${NS:-20}
NREP=${NREP:-3}
RESID=${RESID:-"ala gly pro"}
TIERS=${TIERS:-"implicit explicit"}

for RES in ${RESID}; do
  for TIER in ${TIERS}; do
    echo "===== GENERATE ${RES}/${TIER} (${NREP}x${NS}ns) $(date -u +%H:%M:%S) ====="
    REP=1
    while [ ${REP} -le ${NREP} ]; do
      ${PY} -m benchmark.simulate --residue ${RES} --tier ${TIER} \
            --replica ${REP} --ns ${NS} --equil-ns 0.5 --seed ${REP} \
            --platform CUDA || echo "WARN: ${RES}/${TIER}/r${REP} generation failed"
      REP=$((REP + 1))
    done
    echo "===== CURATE ${RES}/${TIER} $(date -u +%H:%M:%S) ====="
    ${PY} -m benchmark.curate --residue ${RES} --tier ${TIER} \
      || echo "WARN: ${RES}/${TIER} curate failed"
    if [ -z "${KEEP_DCD}" ] && [ -f "out/${RES}/${TIER}/curated.npz" ]; then
      echo "===== PRUNE DCDs ${RES}/${TIER} ====="
      rm -f out/${RES}/${TIER}/r*/traj.dcd
    fi
    du -sh "out/${RES}/${TIER}" 2>/dev/null
  done
done
echo "===== SHAKE-OUT DONE $(date -u +%H:%M:%S) ====="
du -sh out/* 2>/dev/null

#!/bin/sh
# Quick throughput benchmark for the production stage, so cost/run can be
# computed from a real number instead of an estimate. Requires that
# equilibration has already produced out/${ID}/raw/npt.gro (run scripts 0 and 1
# first -- on a GPU box that is only a few minutes).
#
# Usage:
#   GPU=1 ID=ala sh scripts/bench.sh
#   GPU=1 PRICE=0.40 TARGET_NS=500 BENCH_NSTEPS=50000 ID=ala sh scripts/bench.sh
#
# Env:
#   GPU          1 to offload to GPU (default 0 = CPU)
#   BENCH_NSTEPS steps to time (default 50000 = 50 ps; long enough past warmup)
#   TARGET_NS    production ns per protein to extrapolate to (default 500)
#   PRICE        instance $/hour for the cost estimate (default 0.40)
#   NTO          OpenMP threads (default 8)

NTO=${NTO:-8}
BENCH_NSTEPS=${BENCH_NSTEPS:-50000}
TARGET_NS=${TARGET_NS:-500}
PRICE=${PRICE:-0.40}

if [ "${GPU:-0}" = "1" ]; then PROD_GPU="-nb gpu -pme gpu -bonded gpu -update gpu"; else PROD_GPU=""; fi

OUTDIR=out/${ID}/raw
if [ ! -f "${OUTDIR}/npt.gro" ]; then
    echo "ERROR: ${OUTDIR}/npt.gro not found. Run scripts 0 and 1 first." >&2
    exit 1
fi

echo "=== benchmark: ID=${ID} GPU=${GPU:-0} steps=${BENCH_NSTEPS} threads=${NTO} ==="
gmx grompp -f config/md_langevin.mdp -c ${OUTDIR}/npt.gro -p ${OUTDIR}/topol.top -o ${OUTDIR}/bench.tpr
gmx mdrun -ntmpi 1 -ntomp ${NTO} ${PROD_GPU} -nsteps ${BENCH_NSTEPS} -deffnm ${OUTDIR}/bench

# Parse "Performance:   <ns/day>  ..." from the run log.
NSDAY=$(awk '/Performance:/ {print $2}' ${OUTDIR}/bench.log | tail -1)
if [ -z "${NSDAY}" ]; then
    echo "ERROR: could not read Performance from ${OUTDIR}/bench.log" >&2
    exit 1
fi

# hours for one protein = TARGET_NS / ns_per_day * 24 ; cost = hours * PRICE
echo ""
echo "------------------------------------------------------------"
awk -v ns="${NSDAY}" -v tgt="${TARGET_NS}" -v price="${PRICE}" 'BEGIN {
    hrs_per_rep = (tgt/5.0)/ns*24.0;   # assuming 5 replicas of tgt/5 ns each
    hrs = tgt/ns*24.0;                 # total GPU-hours (replicas run sequentially)
    printf "Throughput : %.1f ns/day (single replica)\n", ns;
    printf "Target     : %.0f ns/protein\n", tgt;
    printf "Wall time  : %.1f h/protein (5 replicas, sequential)\n", hrs;
    printf "Cost       : $%.2f/protein at $%.2f/h\n", hrs*price, price;
}'
echo "------------------------------------------------------------"
echo "(rerun with PRICE=<your \$/hr> to match the instance you rented)"

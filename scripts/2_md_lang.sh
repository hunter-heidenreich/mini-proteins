NTM=1
NTO=8

# Independent production replicas. Each is grompp'd separately so it gets its
# own random gen_vel/ld_seed (config/md_langevin.mdp) -> independent trajectory
# from the shared NPT structure. Override at the command line, e.g.
#   NREP=3 NSTEPS=50000000 ID=ala sh scripts/2_md_lang.sh
NREP=${NREP:-5}            # number of independent replicas
NSTEPS=${NSTEPS:-100000000} # steps per replica (1 fs * 1e8 = 100 ns); overrides the .mdp

OUTDIR=out/${ID}/raw

REP=1
while [ ${REP} -le ${NREP} ]; do
    TAG=md_lang_r${REP}
    echo "=== production replica ${REP}/${NREP} (${TAG}) ==="

    # Configure production Langevin dynamics. We start from the equilibrated NPT
    # coordinates/box but regenerate velocities (gen_vel = yes), so we do NOT pass
    # -t npt.cpt here -- that keeps the replicas independent rather than identical.
    gmx grompp -f config/md_langevin.mdp -c ${OUTDIR}/npt.gro -p ${OUTDIR}/topol.top -o ${OUTDIR}/${TAG}.tpr

    # Run production. -nsteps overrides the value baked into the .tpr.
    gmx mdrun -ntmpi ${NTM} -ntomp ${NTO} -nsteps ${NSTEPS} -deffnm ${OUTDIR}/${TAG}

    # Thermodynamic quantities (selected by name; robust to term-layout changes)
    printf 'Potential\n'    | gmx energy -f ${OUTDIR}/${TAG} -o ${OUTDIR}/${TAG}_pot.xvg
    printf 'Total-Energy\n' | gmx energy -f ${OUTDIR}/${TAG} -o ${OUTDIR}/${TAG}_etot.xvg
    printf 'Temperature\n'  | gmx energy -f ${OUTDIR}/${TAG} -o ${OUTDIR}/${TAG}_temp.xvg

    REP=$((REP + 1))
done

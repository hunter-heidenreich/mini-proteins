NTM=1
NTO=8

OUTDIR=out/${ID}/raw

# Configure production Langevin dynamics simulation
# -f: input file
# -c: input structure
# -t: input checkpoint file
# -p: input topology
# -o: output file
gmx grompp -f config/md_langevin.mdp -c ${OUTDIR}/npt.gro -t ${OUTDIR}/npt.cpt -p ${OUTDIR}/topol.top -o ${OUTDIR}/md_lang.tpr

# Run production Langevin dynamics simulation
# -ntmpi: number of MPI threads
# -ntomp: number of OpenMP threads
# -deffnm: default file name
gmx mdrun -ntmpi $NTM -ntomp $NTO -deffnm ${OUTDIR}/md_lang

# Calculate thermodynamic quantity during production Langevin dynamics simulation
# -f: input file
# -o: output file
echo "10\n"| gmx energy -f ${OUTDIR}/md_lang -o ${OUTDIR}/md_lang_pot.xvg
echo "12\n"| gmx energy -f ${OUTDIR}/md_lang -o ${OUTDIR}/md_lang_etot.xvg
echo "13\n"| gmx energy -f ${OUTDIR}/md_lang -o ${OUTDIR}/md_lang_temp.xvg

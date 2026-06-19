NTM=1
NTO=8

# GPU offload toggle. GPU=1 offloads nonbonded + PME to the GPU; this subset is
# safe with the position restraints used during equilibration (we deliberately
# do NOT force -update gpu here, which can conflict with posres). Default 0 = CPU.
if [ "${GPU:-0}" = "1" ]; then EQ_GPU="-nb gpu -pme gpu"; else EQ_GPU=""; fi

OUTDIR=out/${ID}/raw

# Configure relaxation
# -f: input file
# -c: input structure
# -p: input topology
# -o: output file
gmx grompp -f config/minim.mdp -c ${OUTDIR}/solv_ions.gro -p ${OUTDIR}/topol.top -o ${OUTDIR}/em.tpr

# Relax
# -ntmpi: number of MPI threads
# -ntomp: number of OpenMP threads
# -v: verbose
# -deffnm: default file name
gmx mdrun -ntmpi $NTM -ntomp $NTO ${EQ_GPU} -v -deffnm ${OUTDIR}/em

# Calculate energy during relaxation
# -f: input file
# -o: output file
#
# Note: terms are selected by name (robust to FF/electrostatics changes).
# printf (not echo) so the newline is interpreted across shells.
printf 'Potential\n' | gmx energy -f ${OUTDIR}/em -o ${OUTDIR}/em_pot.xvg

# Configure NVT equilibration
# -f: input file
# -r: input structure
# -p: input topology
# -o: output file
gmx grompp -f config/nvt_langevin.mdp -c ${OUTDIR}/em.gro -r ${OUTDIR}/em.gro -p ${OUTDIR}/topol.top -o ${OUTDIR}/nvt.tpr

# Equilibrate NVT
# -ntmpi: number of MPI threads
# -ntomp: number of OpenMP threads
# -deffnm: default file name
gmx mdrun -ntmpi $NTM -ntomp $NTO ${EQ_GPU} -deffnm ${OUTDIR}/nvt

# Calculate thermodynamic quantity during NVT equilibration
# -f: input file
# -o: output file
#
# Note: the echo command is used to select the term
printf 'Potential\n' | gmx energy -f ${OUTDIR}/nvt -o ${OUTDIR}/nvt_pot.xvg
printf 'Total-Energy\n' | gmx energy -f ${OUTDIR}/nvt -o ${OUTDIR}/nvt_etot.xvg
printf 'Temperature\n' | gmx energy -f ${OUTDIR}/nvt -o ${OUTDIR}/nvt_temp.xvg

# Configure NPT equilibration
# -f: input file
# -c: input structure
# -r: input structure
# -t: input checkpoint file
# -p: input topology
# -o: output file
gmx grompp -f config/npt_langevin.mdp -c ${OUTDIR}/nvt.gro -r ${OUTDIR}/nvt.gro -t ${OUTDIR}/nvt.cpt -p ${OUTDIR}/topol.top -o ${OUTDIR}/npt.tpr

# Equilibrate NPT
# -ntmpi: number of MPI threads
# -ntomp: number of OpenMP threads
# -deffnm: default file name
gmx mdrun -ntmpi $NTM -ntomp $NTO ${EQ_GPU} -deffnm ${OUTDIR}/npt

# Calculate thermodynamic quantity during NPT equilibration
# -f: input file
# -o: output file
printf 'Potential\n' | gmx energy -f ${OUTDIR}/npt -o ${OUTDIR}/npt_pot.xvg
printf 'Total-Energy\n' | gmx energy -f ${OUTDIR}/npt -o ${OUTDIR}/npt_etot.xvg
printf 'Temperature\n' | gmx energy -f ${OUTDIR}/npt -o ${OUTDIR}/npt_temp.xvg
#printf 'Density\n' | gmx energy -f ${OUTDIR}/npt -o ${OUTDIR}/npt_dens.xvg

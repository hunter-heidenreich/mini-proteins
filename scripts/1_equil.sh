NTM=1
NTO=8

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
gmx mdrun -ntmpi $NTM -ntomp $NTO -v -deffnm ${OUTDIR}/em

# Calculate energy during relaxation
# -f: input file
# -o: output file
#
# Note: the echo command is used to select the potential energy term
echo "10\n" | gmx energy -f ${OUTDIR}/em -o ${OUTDIR}/em_pot.xvg

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
gmx mdrun -ntmpi $NTM -ntomp $NTO -deffnm ${OUTDIR}/nvt

# Calculate thermodynamic quantity during NVT equilibration
# -f: input file
# -o: output file
#
# Note: the echo command is used to select the term
echo "11\n" | gmx energy -f ${OUTDIR}/nvt -o ${OUTDIR}/nvt_pot.xvg
echo "13\n" | gmx energy -f ${OUTDIR}/nvt -o ${OUTDIR}/nvt_etot.xvg
echo "14\n" | gmx energy -f ${OUTDIR}/nvt -o ${OUTDIR}/nvt_temp.xvg

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
gmx mdrun -ntmpi $NTM -ntomp $NTO -deffnm ${OUTDIR}/npt

# Calculate thermodynamic quantity during NPT equilibration
# -f: input file
# -o: output file
echo "11\n"| gmx energy -f ${OUTDIR}/npt -o ${OUTDIR}/npt_pot.xvg
echo "13\n"| gmx energy -f ${OUTDIR}/npt -o ${OUTDIR}/npt_etot.xvg
echo "14\n"| gmx energy -f ${OUTDIR}/npt -o ${OUTDIR}/npt_temp.xvg
#echo "23\n"| gmx energy -f ${OUTDIR}/npt -o ${OUTDIR}/npt_dens.xvg

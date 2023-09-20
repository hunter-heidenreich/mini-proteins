FF=amber03  # force field
WATER=tip3p  # water model
BOX=1.0  # distance between molecule and box edge
BOX_TYPE=cubic  # box type
PNAME=NA  # positive ion
CNAME=CL  # negative ion
CONC=0.1  # ion concentration (mol/L)

OUTDIR=out/${ID}/raw

mkdir -p ${OUTDIR}

# Build the topology
# -f: input coordinate file
# -o: output coordinate file
# -p: output topology file
# -i: output position restraint file
# -ff: force field
# -water: water model
gmx pdb2gmx -f data/${ID}.pdb -o ${OUTDIR}/processed.gro -p ${OUTDIR}/topol.top -i ${PWD}/${OUTDIR}/posre.itp -ff ${FF} -water ${WATER}

# Center the molecule in the box
# -f: input coordinate file
# -o: output coordinate file
# -c: center the molecule in the box
# -d: distance between molecule and box edge
# -bt: box type
gmx editconf -f ${OUTDIR}/processed.gro -o ${OUTDIR}/newbox.gro -c -d ${BOX} -bt ${BOX_TYPE}

# Solvate
# -cp: input coordinate file
# -cs: input coordinate file
# -o: output coordinate file
# -p: output topology file
gmx solvate -cp ${OUTDIR}/newbox.gro -cs spc216.gro -o ${OUTDIR}/solv.gro -p ${OUTDIR}/topol.top

# Generate parameter files
# -f: input mdp file
# -c: input coordinate file
# -p: input topology file
# -o: output tpr file
gmx grompp -f config/minim.mdp -c ${OUTDIR}/solv.gro -p ${OUTDIR}/topol.top -o ${OUTDIR}/ions.tpr

# Add ions
# -s: input tpr file
# -o: output gro file
# -p: output top file
# -pname: name of positive ion
# -nname: name of negative ion
# - neutral: add ions to neutralize the system
# - conc: concentration of ions in mol/l
#
# Note: automatically selecting solution group 'SOL' to be replaced by ions
echo SOL | gmx genion -s ${OUTDIR}/ions.tpr -o ${OUTDIR}/solv_ions.gro -p ${OUTDIR}/topol.top -pname ${PNAME} -nname ${CNAME} -conc ${CONC} -neutral

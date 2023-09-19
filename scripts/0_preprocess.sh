FF=amber03  # force field
WATER=tip3p  # water model
BOX=1.0  # distance between molecule and box edge
BOX_TYPE=cubic  # box type
PNAME=NA  # positive ion
CNAME=CL  # negative ion
CONC=0.1  # ion concentration (mol/L)

mkdir -p out/${ID}

# Build the topology
# -f: input coordinate file
# -o: output coordinate file
# -p: output topology file
# -i: output position restraint file
# -ff: force field
# -water: water model
gmx pdb2gmx -f data/${ID}.pdb -o out/${ID}/processed.gro -p out/${ID}/topol.top -i out/${ID}/posre.itp -ff ${FF} -water ${WATER}

# Center the molecule in the box
# -f: input coordinate file
# -o: output coordinate file
# -c: center the molecule in the box
# -d: distance between molecule and box edge
# -bt: box type
gmx editconf -f out/${ID}/processed.gro -o out/${ID}/newbox.gro -c -d ${BOX} -bt ${BOX_TYPE}

# Solvate
# -cp: input coordinate file
# -cs: input coordinate file
# -o: output coordinate file
# -p: output topology file
gmx solvate -cp out/${ID}/newbox.gro -cs spc216.gro -o out/${ID}/solv.gro -p out/${ID}/topol.top

# Generate parameter files
# -f: input mdp file
# -c: input coordinate file
# -p: input topology file
# -o: output tpr file
gmx grompp -f config/minim.mdp -c out/${ID}/solv.gro -p out/${ID}/topol.top -o out/${ID}/ions.tpr

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
echo SOL | gmx genion -s out/${ID}/ions.tpr -o out/${ID}/solv_ions.gro -p out/${ID}/topol.top -pname ${PNAME} -nname ${CNAME} -conc ${CONC} -neutral

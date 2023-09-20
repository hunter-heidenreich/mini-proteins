OUTDIR=out/${ID}/raw

# Post-process production Langevin dynamics simulation to extract forces
# -f: input file
# -s: input structure
# -of: output file
#
# Note: the echo command is used to select the forces
echo "1\n" | gmx traj -f ${OUTDIR}/md_lang.trr -s ${OUTDIR}/md_lang.tpr -of ${OUTDIR}/md_lang_forces.xvg

# Post-process production Langevin dynamics simulation to extract trajectory
# -s: input structure
# -f: input file
# -o: output file
# -pbc: periodic boundary conditions
# -center: center the system
#
# Note: the echo command is used to select the protein
echo "1\n1\n" | gmx trjconv -s ${OUTDIR}/md_lang.tpr -f ${OUTDIR}/md_lang.trr -o ${OUTDIR}/md_lang_noPBC.pdb -pbc mol -center

ID=$ID python scripts/plot.py

mkdir -p out/${ID}/data
mv out/${ID}/raw/md_lang_*.xvg out/${ID}/data/
mv out/${ID}/raw/md_lang_noPBC.pdb out/${ID}/data/${ID}-sim.pdb

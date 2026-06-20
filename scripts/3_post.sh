OUTDIR=out/${ID}/raw

# Extract forces and a PBC-corrected trajectory for every production replica.
# Replicas are discovered by globbing, so this stays in sync with 2_md_lang.sh
# regardless of NREP.
for TRR in ${OUTDIR}/md_lang_r*.trr; do
    TAG=$(basename "${TRR}" .trr)
    echo "=== post-processing ${TAG} ==="

    # Forces on the solute (the ML training labels). Group selected by name.
    printf 'Protein\n' | gmx traj -f "${TRR}" -s ${OUTDIR}/${TAG}.tpr -of ${OUTDIR}/${TAG}_forces.xvg

    # PBC-corrected solute trajectory: center on Protein, output Protein.
    printf 'Protein\nProtein\n' | gmx trjconv -s ${OUTDIR}/${TAG}.tpr -f "${TRR}" \
        -o ${OUTDIR}/${TAG}_noPBC.pdb -pbc mol -center
done

# Energy/temperature sanity plots (reads the per-replica xvg files in raw/)
ID=$ID uv run python scripts/plot.py

# Move final products into data/
mkdir -p out/${ID}/data
mv out/${ID}/raw/md_lang_r*_*.xvg out/${ID}/data/
for PDB in out/${ID}/raw/md_lang_r*_noPBC.pdb; do
    TAG=$(basename "${PDB}" _noPBC.pdb)   # e.g. md_lang_r1
    REP=${TAG#md_lang_}                   # e.g. r1
    mv "${PDB}" out/${ID}/data/${ID}-sim-${REP}.pdb
done

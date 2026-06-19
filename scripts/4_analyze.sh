OUTDIR=out/${ID}/raw

# Backbone phi/psi for every production replica, via gmx rama (non-interactive;
# emits one (phi, psi, residue) record per frame for each residue that has both
# dihedrals defined -- for a capped dipeptide that is just the central residue).
for TRR in ${OUTDIR}/md_lang_r*.trr; do
    TAG=$(basename "${TRR}" .trr)
    echo "=== rama ${TAG} ==="
    gmx rama -s ${OUTDIR}/${TAG}.tpr -f "${TRR}" -o ${OUTDIR}/${TAG}_rama.xvg
done

# Ramachandran free-energy surface, basin populations (per-replica error bars),
# glycine symmetry check, and Karplus 3J(HN,Ha) for comparison with experiment.
ID=$ID python scripts/analyze.py

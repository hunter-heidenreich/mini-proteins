OUTDIR=out/${ID}/raw

# Backbone phi/psi for every production replica, via gmx rama (non-interactive;
# emits one (phi, psi, residue) record per frame for each residue that has both
# dihedrals defined -- for a capped dipeptide that is just the central residue).
for TRR in ${OUTDIR}/md_lang_r*.trr; do
    TAG=$(basename "${TRR}" .trr)
    echo "=== rama ${TAG} ==="
    gmx rama -s ${OUTDIR}/${TAG}.tpr -f "${TRR}" -o ${OUTDIR}/${TAG}_rama.xvg
done

# Full validation + dataset-QC report (writes data/validation.txt + .json):
#   convergence (statistical inefficiency, effective samples, basin-transition
#   counts, half/half stationarity), FES + per-basin free energies (kcal/mol),
#   populations vs the published ff03/TIP3P reference, 3J(HN,Ha) under several
#   Karplus sets, glycine symmetry, and force-label QC. The force-label and
#   energy-drift sections read the per-replica xvg files produced by 3_post.sh,
#   so run stage 3 first; without them analyze.py still runs and notes the skip.
ID=$ID python scripts/analyze.py

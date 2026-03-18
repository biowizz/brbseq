#!/bin/bash
# ===============================
# Step 6: Count matrix assembly
# Merges ReadsPerGene.out.tab
# from all 12 STAR outputs into
# a single TSV count matrix
#
# BRB-seq = forward stranded
# --> use column 2 of STAR output
#     (1-indexed col 2 = forward)
#
# STAR ReadsPerGene.out.tab cols:
#   col1 = gene_id
#   col2 = unstranded
#   col3 = forward strand  <-- USE
#   col4 = reverse strand
# ===============================

echo ""
echo "============================================"
echo " STEP 6: Count matrix assembly"
echo " $(date)"
echo "============================================"

OUTFILE="${COUNTS}/counts_raw.tsv"

# ---- Define sample order (matches metadata) ----
SAMPLES=(
    "1B_ctrl_1_F03"
    "1B_ctrl_2_G03"
    "1B_ctrl_3_H03"
    "1B_KD4h_1_A04"
    "1B_KD4h_2_B04"
    "1B_KD4h_3_C04"
    "D2_ctrl_1_G04"
    "D2_ctrl_2_H04"
    "D2_ctrl_3_A05"
    "D2_KD4h_1_B05"
    "D2_KD4h_2_C05"
    "D2_KD4h_3_D05"
)

# ---- Verify all count files exist ----
echo "Checking count files..."
for s in "${SAMPLES[@]}"; do
    f="${ALIGNED}/${s}/ReadsPerGene.out.tab"
    if [ ! -f "$f" ]; then
        echo "[ERROR] Missing: $f" >&2
        exit 1
    fi
    echo "  Found: $f"
done

# ---- Build header line ----
HEADER="gene_id"
for s in "${SAMPLES[@]}"; do
    HEADER="${HEADER}\t${s}"
done

# ---- Extract gene IDs from first sample (skip 4 STAR summary rows) ----
FIRST="${ALIGNED}/${SAMPLES[0]}/ReadsPerGene.out.tab"

# ---- Paste all forward-strand count columns (col 3) ----
# tail -n +5 skips the first 4 summary lines (N_unmapped, N_multimapping, etc.)
echo -e "$HEADER" > $OUTFILE

paste \
    <(tail -n +5 ${ALIGNED}/${SAMPLES[0]}/ReadsPerGene.out.tab  | cut -f1) \
    <(tail -n +5 ${ALIGNED}/${SAMPLES[0]}/ReadsPerGene.out.tab  | cut -f3) \
    <(tail -n +5 ${ALIGNED}/${SAMPLES[1]}/ReadsPerGene.out.tab  | cut -f3) \
    <(tail -n +5 ${ALIGNED}/${SAMPLES[2]}/ReadsPerGene.out.tab  | cut -f3) \
    <(tail -n +5 ${ALIGNED}/${SAMPLES[3]}/ReadsPerGene.out.tab  | cut -f3) \
    <(tail -n +5 ${ALIGNED}/${SAMPLES[4]}/ReadsPerGene.out.tab  | cut -f3) \
    <(tail -n +5 ${ALIGNED}/${SAMPLES[5]}/ReadsPerGene.out.tab  | cut -f3) \
    <(tail -n +5 ${ALIGNED}/${SAMPLES[6]}/ReadsPerGene.out.tab  | cut -f3) \
    <(tail -n +5 ${ALIGNED}/${SAMPLES[7]}/ReadsPerGene.out.tab  | cut -f3) \
    <(tail -n +5 ${ALIGNED}/${SAMPLES[8]}/ReadsPerGene.out.tab  | cut -f3) \
    <(tail -n +5 ${ALIGNED}/${SAMPLES[9]}/ReadsPerGene.out.tab  | cut -f3) \
    <(tail -n +5 ${ALIGNED}/${SAMPLES[10]}/ReadsPerGene.out.tab | cut -f3) \
    <(tail -n +5 ${ALIGNED}/${SAMPLES[11]}/ReadsPerGene.out.tab | cut -f3) \
    >> $OUTFILE

echo "[$(date)] STEP 6 COMPLETE"
echo " Count matrix: $OUTFILE"
echo " Dimensions  : $(tail -n +2 $OUTFILE | wc -l) genes x 12 samples"
echo ""

# ---- Library size per sample ----
echo " === Library sizes (forward-strand mapped counts) ==="
printf "%-35s %s\n" "Sample" "Total counts"
i=2  # column index starts at 2 (col1 = gene_id)
for s in "${SAMPLES[@]}"; do
    total=$(tail -n +2 $OUTFILE | awk -v col=$i '{sum+=$col} END {print sum}' OFS='\t')
    printf "%-35s %s\n" "$s" "$total"
    ((i++))
done

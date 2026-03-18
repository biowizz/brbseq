#!/bin/bash
# ===============================
# Step 3: Post-trim QC
# FastQC on trimmed reads
# MultiQC aggregation
# ===============================

echo ""
echo "============================================"
echo " STEP 3: Post-trim QC"
echo " $(date)"
echo "============================================"

fastqc \
    --threads $THREADS \
    --outdir $POSTQC \
    --format fastq \
    $POSTQC/*_trimmed.fq.gz

echo "[$(date)] FastQC done. Running MultiQC..."

multiqc \
    $POSTQC \
    --outdir $POSTQC \
    --filename multiqc_postqc \
    --title "Post-trim QC — Pawel BRB-seq" \
    --force

echo "[$(date)] STEP 3 COMPLETE"
echo " Report: $POSTQC/multiqc_postqc.html"
echo ""
echo " === Checklist before proceeding to alignment ==="
echo "   - Adapter content should be near 0%"
echo "   - Mean quality should be > Q28 across all samples"
echo "   - No sample should have fewer than 1M reads after trimming"

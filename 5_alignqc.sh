#!/bin/bash
# ===============================
# Step 5: Post-alignment QC
# Qualimap rnaseq per sample
# MultiQC over STAR logs +
# Qualimap reports
# ===============================

echo ""
echo "============================================"
echo " STEP 5: Post-alignment QC"
echo " $(date)"
echo "============================================"

for bam in $ALIGNED/*/Aligned.sortedByCoord.out.bam; do

    sample=$(basename $(dirname $bam))
    qm_out="${ALIGNED}/${sample}/qualimap"

    echo "[$(date)] Qualimap: $sample"

    qualimap rnaseq \
        -bam $bam \
        -gtf $GTF \
        -outdir $qm_out \
        -outformat HTML \
        -p strand-specific-forward \
        --java-mem-size=16G \
        -nt $THREADS

    echo "[$(date)] Done: $sample"

done

echo "[$(date)] All Qualimap jobs done. Running MultiQC..."

multiqc \
    $ALIGNED \
    --outdir $ALIGNED/multiqc_aligned \
    --filename multiqc_aligned \
    --title "Post-alignment QC — Pawel BRB-seq" \
    --force

echo "[$(date)] STEP 5 COMPLETE"
echo " Report: $ALIGNED/multiqc_aligned/multiqc_aligned.html"
echo ""
echo " === QC thresholds to verify in MultiQC report ==="
echo "   Uniquely mapped reads  : target > 80%"
echo "   Exonic rate            : target > 70%"
echo "   3-prime transcript bias: expected HIGH — normal for BRB-seq"
echo "   Duplication rate       : can be high — normal for 3-end protocol"

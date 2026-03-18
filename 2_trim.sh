#!/bin/bash
# ===============================
# Step 2: Adapter Trimming
# Trim Galore — Illumina TruSeq
# Single-end, 100 bp reads
# ===============================

echo ""
echo "============================================"
echo " STEP 2: Adapter Trimming"
echo " $(date)"
echo "============================================"

for fq in $RAW_DATA/*.fastq.gz; do

    sample=$(basename $fq .fastq.gz)
    echo "[$(date)] Trimming: $sample"

    trim_galore \
        --illumina \
        --cores 8 \
        --quality 20 \
        --length 30 \
        --output_dir $POSTQC \
        $fq

    echo "[$(date)] Done: $sample"
done

echo "[$(date)] STEP 2 COMPLETE"
echo " Trimmed files: $(ls $POSTQC/*_trimmed.fq.gz 2>/dev/null | wc -l) / 12 expected"

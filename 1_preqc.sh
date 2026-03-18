#!/bin/bash
# ===============================
# Step 1: Pre-alignment QC
# FastQC on all 12 raw FASTQs
# MultiQC aggregation
# ===============================

echo ""
echo "============================================"
echo " STEP 1: Pre-alignment QC"
echo " $(date)"
echo "============================================"

fastqc \
    --threads $THREADS \
    --outdir $PREQC \
    --format fastq \
    $RAW_DATA/*.fastq.gz

echo "[$(date)] FastQC done. Running MultiQC..."

multiqc \
    $PREQC \
    --outdir $PREQC \
    --filename multiqc_preqc \
    --title "Pre-alignment QC — Pawel BRB-seq" \
    --force

echo "[$(date)] STEP 1 COMPLETE"
echo " Report: $PREQC/multiqc_preqc.html"
echo " Files processed: $(ls $RAW_DATA/*.fastq.gz | wc -l)"

#!/bin/bash
#SBATCH -J pawel_rnaseq
#SBATCH -p standard
#SBATCH -c 40
#SBATCH --mem=248G
#SBATCH -t 70:00:00
#SBATCH -o pawel_rnaseq.%j.out
#SBATCH -e pawel_rnaseq.%j.err

# ===============================
# BRB-seq Pipeline — Pawel samples
# RPRD1B and RPRD2 knockdown (4h)
# 12 samples, single-end, 100 bp
# GRCh38 + GENCODE v48
# ===============================
# USAGE:
#   Comment/uncomment steps below
#   as you progress through the
#   pipeline. Always keep step 0
#   (paths) sourced.
#
#   sbatch run_pipeline.sh
# ===============================

### ---- Load project paths ----
### Always keep enabled
source /scratch/sandeepm/pawel/0_paths.sh

### ---- Step 1: Pre-alignment QC ----
#source /scratch/sandeepm/pawel/1_preqc.sh

### ---- Step 2: Adapter trimming ----
#source /scratch/sandeepm/pawel/2_trim.sh

### ---- Step 3: Post-trim QC ----
#source /scratch/sandeepm/pawel/3_postqc.sh

### ---- Step 4: STAR alignment ----
#source /scratch/sandeepm/pawel/4_align.sh

### ---- Step 5: Post-alignment QC ----
#source /scratch/sandeepm/pawel/5_alignqc.sh

### ---- Step 6: Count matrix assembly ----
#source /scratch/sandeepm/pawel/6_counts.sh

### ---- Step 7: DESeq2 (run interactively in R) ----
# Rscript /scratch/sandeepm/pawel/7_deseq2.R

### ---- Step 8: Enrichment (run interactively in R) ----
# Rscript /scratch/sandeepm/pawel/8_enrichment.R

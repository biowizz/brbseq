#!/bin/bash
# ===============================
# RNA-seq project paths
# Pawel BRB-seq — RPRD1B / RPRD2 KD
# ===============================

# ---- Project root ----
PROJECT_ROOT="/scratch/sandeepm/pawel"

# ---- Raw data ----
RAW_DATA="${PROJECT_ROOT}/raw"

# ---- Reference genome ----
GENOME_DIR="/scratch/mallya/mihirn/lncRNA/genome"

GENOME_FASTA="${GENOME_DIR}/GRCh38.primary_assembly.genome.fa"
GENOME_FAI="${GENOME_DIR}/GRCh38.primary_assembly.genome.fa.fai"

# ---- Annotation ----
GTF="${GENOME_DIR}/gencode.v48.annotation.gtf"

# ---- STAR index ----
STAR_INDEX="${GENOME_DIR}/star_index"

# ---- Output structure ----
RESULTS="${PROJECT_ROOT}/results"

PREQC="${RESULTS}/preqc"
POSTQC="${RESULTS}/postqc"
ALIGNED="${RESULTS}/aligned"
COUNTS="${RESULTS}/counts"

# ---- Threads ----
THREADS=40

# ===============================
# Create directories if missing
# ===============================
mkdir -p \
  "$PREQC"   \
  "$POSTQC"  \
  "$ALIGNED" \
  "$COUNTS"

echo "[paths] Project root : $PROJECT_ROOT"
echo "[paths] Raw data     : $RAW_DATA"
echo "[paths] Results      : $RESULTS"
echo "[paths] STAR index   : $STAR_INDEX"
echo "[paths] GTF          : $GTF"

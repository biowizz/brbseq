#!/bin/bash
# ===============================
# Step 4: STAR Alignment
# Single-end, 100 bp reads
# GRCh38 + GENCODE v48
# --quantMode GeneCounts for
# per-gene count output
# ===============================

echo ""
echo "============================================"
echo " STEP 4: STAR Alignment"
echo " $(date)"
echo "============================================"

for fq in $POSTQC/*_trimmed.fq.gz; do

    sample=$(basename $fq _trimmed.fq.gz)
    outdir="${ALIGNED}/${sample}"
    mkdir -p $outdir

    echo "[$(date)] Aligning: $sample"

    STAR \
        --runThreadN $THREADS \
        --genomeDir $STAR_INDEX \
        --sjdbGTFfile $GTF \
        --readFilesIn $fq \
        --readFilesCommand zcat \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMattributes NH HI AS NM MD \
        --outBAMsortingThreadN 4 \
        --quantMode GeneCounts \
        --outFileNamePrefix ${outdir}/ \
        --sjdbOverhang 99 \
        --outFilterMultimapNmax 10 \
        --outFilterMismatchNmax 999 \
        --outFilterMismatchNoverReadLmax 0.04 \
        --alignSJoverhangMin 8 \
        --alignSJDBoverhangMin 1 \
        --alignIntronMin 20 \
        --alignIntronMax 1000000 \
        --alignMatesGapMax 1000000 \
        --outFilterType BySJout \
        --outSAMstrandField intronMotif \
        --limitBAMsortRAM 30000000000

    # Index BAM for downstream tools
    samtools index -@ $THREADS ${outdir}/Aligned.sortedByCoord.out.bam

    echo "[$(date)] Done: $sample"

done

echo ""
echo "[$(date)] STEP 4 COMPLETE"
echo ""
echo " === Mapping rate summary ==="
printf "%-35s %-15s %-18s\n" "Sample" "Input reads" "Uniquely mapped%"
for logfile in $ALIGNED/*/Log.final.out; do
    sample=$(basename $(dirname $logfile))
    input=$(grep "Number of input reads"      $logfile | awk '{print $NF}')
    uniq=$( grep "Uniquely mapped reads %"    $logfile | awk '{print $NF}')
    printf "%-35s %-15s %-18s\n" "$sample" "$input" "$uniq"
done
echo " Target: >80% uniquely mapped per sample"

#!/bin/bash

# Directory for input and output
INPUT_DIR="/home/viji_msc_student_group01/Suryadip/GSE197829/GSE197829"
OUTPUT_DIR="/home/viji_msc_student_group01/Suryadip/GSE197829/star_aln_out"
INDEX_DIR="/home/viji_msc_student_group01/Suryadip/GSE197829/hg38_main/hg38_idx"

# Create output directory if it doesn't exist
mkdir -p $OUTPUT_DIR

# Loop over all paired-end FASTQ files in the input directory
for R1 in $INPUT_DIR/*_1.fastq.gz; do
    # Get the base filename (without _1.fastq.gz)
    BASE=$(basename $R1 _1.fastq.gz)
    
    # Define R2 file
    R2="${INPUT_DIR}/${BASE}_2.fastq.gz"
    
    # Define output file name prefix
    OUTPUT_PREFIX="$OUTPUT_DIR/${BASE}_"

    # Align with STAR
    echo "Aligning ${R1} and ${R2} to HG38 genome......."
    STAR --runThreadN 30 --genomeDir $INDEX_DIR --readFilesIn $R1 $R2 \
         --outFileNamePrefix $OUTPUT_PREFIX --outSAMtype BAM SortedByCoordinate \
         --readFilesCommand gunzip -c --limitBAMsortRAM 1756793772
    
    # Index the sorted BAM file
    echo "Indexing ${BASE}_Aligned.sortedByCoord.out.bam....."
    samtools index ${OUTPUT_PREFIX}Aligned.sortedByCoord.out.bam
done

echo "RNA-seq alignment complete!"


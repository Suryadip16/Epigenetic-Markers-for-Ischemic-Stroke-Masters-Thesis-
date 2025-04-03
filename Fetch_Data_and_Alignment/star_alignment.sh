#!/bin/bash

# Directory for input and output
INPUT_DIR="/home/suryadip/fastq_files"
OUTPUT_DIR="/home/suryadip/star_aln_out"
INDEX_DIR="/home/suryadip/hg38_index"

# Create output directory if it doesn't exist
mkdir -p $OUTPUT_DIR

# Loop over all single-end FASTQ files in the input directory
for R1 in $INPUT_DIR/*.fastq.gz; do
    # Get the base filename (without extension)
    BASE=$(basename $R1 .fastq.gz)
    
    # Define output file name prefix
    OUTPUT_PREFIX="$OUTPUT_DIR/${BASE}_"

    # Align with STAR
    echo "Aligning ${R1} to HG38 genome......."
    STAR --runThreadN 30 --genomeDir $INDEX_DIR --readFilesIn $R1 \
         --outFileNamePrefix $OUTPUT_PREFIX --outSAMtype BAM SortedByCoordinate \
         --readFilesCommand gunzip -c --limitBAMsortRAM 1756793772
    
    # Index the sorted BAM file
    echo "Indexing ${BASE}_Aligned.sortedByCoord.out.bam....."
    samtools index ${OUTPUT_PREFIX}Aligned.sortedByCoord.out.bam
done

echo "RNA-seq alignment complete!"


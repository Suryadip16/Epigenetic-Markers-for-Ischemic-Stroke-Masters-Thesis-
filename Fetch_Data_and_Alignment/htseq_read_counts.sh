#!/bin/bash

# Check if a directory is provided
if [ $# -ne 1 ]; then
    echo "Usage: $0 <directory_with_bam_files>"
    exit 1
fi

# Input directory
BAM_DIR="$1"

# Output directory
OUTPUT_DIR="star_htseq_counts"
mkdir -p "$OUTPUT_DIR"

# Path to the GTF file
GTF_FILE="/home/viji_msc_student_group01/Suryadip/GSE197829/hg38_main/hg38.ensGene.gtf"

# Loop through each BAM file in the directory
for BAM_FILE in "$BAM_DIR"/*.bam; do
    # Extract the filename without extension
    SAMPLE_NAME=$(basename "$BAM_FILE" .bam)
    
    # Run htseq-count
    htseq-count -f bam -r pos -a 27 -m intersection-strict --stranded=no -t exon -i gene_id "$BAM_FILE" "$GTF_FILE" > "$OUTPUT_DIR/counts_${SAMPLE_NAME}.txt"
    
    echo "Processed $BAM_FILE -> $OUTPUT_DIR/counts_${SAMPLE_NAME}.txt"
done

echo "All BAM files processed."


#!/bin/bash

# Check if the SRR list file is provided as an argument
if [ $# -ne 1 ]; then
    echo "Usage: $0 SRR_Acc_List.txt"
    exit 1
fi

# File containing the list of SRR numbers
SRR_LIST=$1

# Check if the file exists
if [ ! -f "$SRR_LIST" ]; then
    echo "File $SRR_LIST does not exist."
    exit 1
fi

# Loop through each SRR number in the list and download the data
while IFS= read -r run; do
    echo "Processing $run..."
    
    # Prefetch the SRA data
    prefetch "$run"
    
    if [ $? -ne 0 ]; then
        echo "Failed to prefetch $run. Skipping."
        continue
    fi
    
    # Convert the SRA file to FASTQ format
    fastq-dump "$run"
    
    if [ $? -ne 0 ]; then
        echo "Failed to dump FASTQ for $run."
    else
        echo "$run downloaded successfully."
        
        # Gzip the FASTQ file to save storage
        echo "Compressing $run.fastq..."
        gzip "$run.fastq"
        
        if [ $? -ne 0 ]; then
            echo "Failed to compress $run.fastq."
        else
            echo "$run.fastq.gz created successfully."
            
            # Remove the original FASTQ file
            rm -f "$run.fastq"
            
            # Remove the directory created by prefetch
            rm -rf "$run"
        fi
    fi
done < "$SRR_LIST"

echo "All SRR numbers processed."


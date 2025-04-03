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
    
    # Convert the SRA file to FASTQ format (paired-end)
    fastq-dump --split-files "$run"
    
    if [ $? -ne 0 ]; then
        echo "Failed to dump FASTQ for $run."
    else
        echo "$run downloaded successfully."
        
        # Compress the FASTQ files to save storage
        echo "Compressing $run\_1.fastq and $run\_2.fastq..."
        gzip "$run"_1.fastq "$run"_2.fastq
        
        if [ $? -ne 0 ]; then
            echo "Failed to compress FASTQ files for $run."
        else
            echo "$run\_1.fastq.gz and $run\_2.fastq.gz created successfully."
            
            # Remove the original FASTQ files
            rm -f "$run"_1.fastq "$run"_2.fastq
            
            # Remove the directory created by prefetch
            rm -rf "$run"
        fi
    fi
done < "$SRR_LIST"

echo "All SRR numbers processed."


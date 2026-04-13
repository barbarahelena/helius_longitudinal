#!/bin/bash

# Script to generate samplesheet.csv for HELIUS paired-end fastq files
# Usage: ./generate_samplesheet.sh

# Set the data directory (current directory where script is located)
DATA_DIR="/projects/prjs1784/heliuspaired/data/fastqsraw"
OUTPUT_FILE="/projects/prjs1784/heliuspaired/data/samplesheet.csv"

# Create the CSV header
echo "sample,group,short_reads_1,short_reads_2,long_reads,short_reads_platform" > "${OUTPUT_FILE}"

# First, collect all unique IDs and assign group numbers
declare -A id_to_group
group_counter=0

# Extract all unique IDs from the files
for file in "${DATA_DIR}"/HELIUS.*.1.fq.gz; do
    basename_file=$(basename "${file}")
    if [[ "${basename_file}" =~ ^HELIUS\.Metagenome\.([0-9]+)\.Fecal\.([A-Z]+)\.NA\.1\.fq\.gz$ ]]; then
        id="${BASH_REMATCH[1]}"
        
        # If this ID hasn't been seen before, assign it a group number
        if [[ -z "${id_to_group[$id]}" ]]; then
            id_to_group[$id]=$group_counter
            ((group_counter++))
        fi
    fi
done

# Now process all _1.fq.gz files and create sample entries
find "${DATA_DIR}" -name "HELIUS.*.1.fq.gz" | sort | while read -r file1; do
    # Extract the base filename without the _1.fq.gz suffix
    basename_file=$(basename "${file1}")
    
    # Check if this is a valid HELIUS file pattern
    if [[ "${basename_file}" =~ ^HELIUS\.Metagenome\.([0-9]+)\.Fecal\.([A-Z]+)\.NA\.1\.fq\.gz$ ]]; then
        # Extract ID and timepoint from the filename
        id="${BASH_REMATCH[1]}"
        timepoint="${BASH_REMATCH[2]}"
        
        # Create sample name: HELI + timepoint + underscore + ID
        sample_name="HELI${timepoint}_${id}"
        
        # Get the group number for this ID
        group_num="${id_to_group[$id]}"
        
        # Construct the corresponding _2.fq.gz filename
        file2="${file1/1.fq.gz/2.fq.gz}"
        
        # Check if the paired file exists
        if [[ -f "${file2}" ]]; then
            # Get relative paths from the data directory
            rel_path1="data/fastqsraw/$(basename "${file1}")"
            rel_path2="data/fastqsraw/$(basename "${file2}")"
            
            # Write the sample entry to CSV with the correct group number
            echo "${sample_name},${group_num},${rel_path1},${rel_path2},,ILLUMINA" >> "${OUTPUT_FILE}"
            
            echo "Added sample: ${sample_name} (group ${group_num})"
        else
            echo "Warning: Paired file not found for ${file1}"
        fi
    else
        echo "Warning: File ${basename_file} does not match expected HELIUS pattern"
    fi
done

echo "Samplesheet generated: ${OUTPUT_FILE}"
echo "Total samples processed: $(tail -n +2 "${OUTPUT_FILE}" | wc -l)"

#!/bin/bash

# might have to use dos2unix compute_coverage.sh to convert to Unix line endings first

# Define input files
bed_file="AS2_merged.bed"
bam_file="BAMS_output/aligned_sorted_samtools.bam"
average_coverage_output="average_coverage.txt"
region_coverage_output="region_coverage.txt"
genome_size_file="/mysrael_ONT/cheetah_ncbi_dataset/ncbi_dataset/data/GCF_027475565.1_chr_size.txt"  # File containing the size of each chromosome

# Ensure the BAM file is indexed
echo "Indexing BAM file..."
samtools index "$bam_file"
if [ $? -ne 0 ]; then
    echo "Error indexing BAM file. Exiting."
    exit 1
fi
echo "BAM file indexing complete."

# Initialize total coverage and regions count
total_coverage=0
total_regions=0

# Create or clear the region coverage output file
echo "Region,Average Coverage" > "$region_coverage_output"

# Function to clean up temporary files
cleanup() {
    rm -f temp_bed_file.bed temp_coverage.txt
}

# Trap to ensure cleanup happens on exit
trap cleanup EXIT

# Compute the coverage for each region in the BED file
echo "Computing coverage for each region in the BED file..."
while read -r line; do
    chrom=$(echo "$line" | awk '{print $1}')
    start=$(echo "$line" | awk '{print $2}')
    end=$(echo "$line" | awk '{print $3}')
    samtools depth -r "${chrom}:${start}-${end}" "$bam_file" > temp_coverage.txt
    if [ $? -ne 0 ]; then
        echo "Coverage computation failed for $line."
        exit 1
    fi

    if [ -s temp_coverage.txt ]; then
        coverage=$(awk '{sum+=$3} END {if (NR > 0) print sum/NR; else print 0}' temp_coverage.txt)
        total_coverage=$(echo "$total_coverage + $coverage" | bc)
        total_regions=$((total_regions + 1))
        echo "Processed region: $line"
        echo "$line,$coverage" >> "$region_coverage_output"
    else
        echo "No coverage data for region: $line"
        echo "$line,0" >> "$region_coverage_output"
    fi
    
    # Cleanup temporary coverage file for each iteration
    rm -f temp_coverage.txt
done < "$bed_file"

if [ "$total_regions" -gt 0 ]; then
    average_coverage=$(echo "scale=2; $total_coverage / $total_regions" | bc)
else
    average_coverage=0
fi
echo "Average coverage for regions in BED file calculated."

# Calculate the total size of the genome
echo "Calculating total genome size..."
total_genome_size=0
while read -r line; do
    size=$(echo "$line" | awk '{print $NF}')
    total_genome_size=$(echo "$total_genome_size + $size" | bc)
done < "$genome_size_file"
echo "Total genome size calculated."

# Calculate the size of the regions covered by the BED file
echo "Calculating size of regions covered by BED file..."
total_bed_size=0
while read -r line; do
    start=$(echo "$line" | awk '{print $2}')
    end=$(echo "$line" | awk '{print $3}')
    region_size=$((end - start))
    total_bed_size=$(echo "$total_bed_size + $region_size" | bc)
done < "$bed_file"
echo "Size of regions covered by BED file calculated."

# Calculate the coverage for the entire genome
echo "Calculating genome coverage..."
genome_coverage=$(samtools depth -a "$bam_file" | awk '{sum+=$3} END {print sum}')
genome_coverage=$(echo "scale=2; $genome_coverage / $total_genome_size" | bc)
echo "Genome coverage calculated."

# Calculate the average coverage for regions not in the BED file
echo "Calculating average coverage for regions not in BED file..."
uncovered_genome_size=$(echo "$total_genome_size - $total_bed_size" | bc)
if [ "$uncovered_genome_size" -gt 0 ]; then
    uncovered_genome_coverage=$(echo "($genome_coverage * $total_genome_size - $total_coverage)" | bc)
    average_uncovered_coverage=$(echo "scale=2; $uncovered_genome_coverage / $uncovered_genome_size" | bc)
else
    average_uncovered_coverage=0
fi
echo "Average coverage for regions not in BED file calculated."

# Write the average coverage to the output file
echo "Writing average coverage to output file..."
echo "Average Coverage (regions in BED file): $average_coverage" > "$average_coverage_output"
echo "Average Coverage (regions not in BED file): $average_uncovered_coverage" >> "$average_coverage_output"
echo "Average coverage values written to $average_coverage_output."

# Print completion message
echo "Coverage computation is complete. Results are saved in $average_coverage_output and $region_coverage_output"

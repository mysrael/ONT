#!/bin/bash

# Give permissions to use script before
# chmod +x post_sequencing_analysis.sh

# Run using
# ./post_sequencing_analysis.sh

# If you edit this script scp updated file to USS
# scp -i mysrael.pem C:\Users\mysrael\Documents\ONT\post_sequencing_analysis.sh ubuntu@132.249.242.229:/mysrael_ONT

# Always open tmux terminal
#tmux new -s session_name
#tmux attach-session -t session_name

# Activate conda environment for pycoQC
#conda activate pycoQC

# Variables: replace with path on USS to these folders (use pwd)
POD5s="/mysrael_ONT/tox443_AS_minion/POD5s_combined" # folder with pod5s
BAMs="/mysrael_ONT/tox443_AS_minion/BAMS_output" # folder for output bams
GENOMIC_FILES="/mysrael_ONT/cheetah_ncbi_dataset/ncbi_dataset/data/GCF_027475565.1" # folder with genomic files
#BED_FILE="/mysrael_ONT/tox443_AS_minion/AS2_merged.bed" # .bed file used for AS
SUMMARY_FILES="/mysrael_ONT/tox443_AS_minion/summary_files" # folder for output QC and summary files
METHYLATION_FILES="/mysrael_ONT/tox443_AS_minion/methylation_files" # folder for methylation analysis files
DORADO="/home/ubuntu/ONT_tools/dorado-0.5.3-linux-x64/bin" # path to dorado bin

# Dorado for basecalling with modified base calling
echo "Executing Dorado basecaller..."
"$DORADO/dorado" basecaller --min-qscore 8 "$DORADO/models/dna_r10.4.1_e8.2_400bps_sup@v4.2.0" "$POD5s" --modified-bases 5mCG_5hmCG > "$BAMs/calls_dorado.bam"
echo "Basecalling finished"

# Dorado for alignment to reference genome
echo "Executing Dorado alignment..."
"$DORADO/dorado" aligner "$GENOMIC_FILES/GCF_027475565.1_VMU_Ajub_asm_v1.0_genomic.fna" "$BAMs/calls_dorado.bam" > "$BAMs/aligned_dorado.bam"
echo "Alignment finished"

# Sort and index with samtools
echo "Executing Samtools sorting and indexing..."
samtools sort -o "$BAMs/aligned_sorted_samtools.bam" "$BAMs/aligned_dorado.bam"
samtools index "$BAMs/aligned_sorted_samtools.bam"
echo "Sorting and indexing finished"

# Create summary file
echo "Creating Dorado summary file..."
"$DORADO/dorado" summary "$BAMs/aligned_dorado.bam" > "$SUMMARY_FILES/summary_dorado.txt"
echo "Summary file created"

# Initialize conda (if not already initialized)
conda init bash

# Activate the pycoQC environment
source activate pycoQC

# Ensure the bam file exists before running PycoQC analysis
if [ ! -f "$BAMs/aligned_sorted_samtools.bam" ]; then
    echo "Error: Input bam file not found!"
    exit 1
fi

# Run PycoQC analysis
echo "Running PycoQC analysis..."
pycoQC --summary_file "$SUMMARY_FILES/summary_dorado.txt" --bam_file "$BAMs/aligned_sorted_samtools.bam" --html_outfile "$SUMMARY_FILES/summary_pycoQC.html"

# Deactivate the conda environment
# conda deactivate
echo "PycoQC analysis finished"

# Modkit for bedmethyl file to view in IGV
echo "Running modkit analysis..."
#modkit pileup "$BAMs/aligned_sorted_samtools.bam" "$METHYLATION_FILES/pileup_merged.bed" --log-filepath "$METHYLATION_FILES/pileup_merged.log"

modkit pileup "$BAMs/BAMS_output/aligned_sorted_samtools.bam" "$METHYLATION_FILES/pileup.bed" --ref "$GENOMIC_FILES//GCF_027475565.1_VMU_Ajub_asm_v1.0_genomic.fna" --cpg --combine-strands --combine-mods --log-filepath "$METHYLATION_FILES/pileup.log"


echo "modkit finished"

echo "use scripts pileup_to_csv3.sh for .csv conversion and compute_coverage.sh for coverage analysis" 
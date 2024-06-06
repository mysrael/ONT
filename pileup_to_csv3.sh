#!/bin/bash

#might need to make Unix
# dos2unix pileup_to_csv3.sh

# Give permissions to use script before
# chmod +x pileup_to_csv3.sh

# Run using
# ./pileup_to_csv3.sh

# If you edit this script scp updated file to USS
# scp -i mysrael.pem C:\Users\mysrael\Documents\ONT\cheetah_AS\AS_run2_tox00443\pileup_to_csv3.sh ubuntu@132.249.242.229:/mysrael_ONT/tox443_AS_minion

# Input and output files
input_csv="AS2_cpgs_postseq.csv"
bed_file="methylation_files/pileup2.bed"
output_csv="pileup3_output.csv"
log_file="no_matches_log3.txt"

# Create the output file and add headers
echo "ID,liftover_chr,liftover_start,liftover_end,note_for_bedfile,chr,start,end,base,depth,A,B,C,D,E,F,G,H,I,J,K" > "$output_csv"
# Create or clear the log file
> "$log_file"

# Read the input CSV file line by line
tail -n +2 "$input_csv" | while IFS=, read -r id liftover_chr liftover_start liftover_end note_for_bedfile
do
  # Debugging: Print the variables to ensure they are read correctly
  echo "Processing: ${id}, ${liftover_chr}, ${liftover_start}, ${liftover_end}, ${note_for_bedfile}"
  
  match_found=false
  
  # Loop through the range +/- 2bp for the start and end positions
  for (( i = liftover_start - 2; i <= liftover_start + 2; i++ )); do
    for (( j = liftover_end - 2; j <= liftover_end + 2; j++ )); do
      # Construct the grep pattern for each combination
      pattern="${liftover_chr}[[:space:]]+${i}[[:space:]]+${j}"

      # Perform the grep command
      bed_line=$(grep -P "^${pattern}" "$bed_file")

      # If a match is found, append the data to the output CSV and exit the loop
      if [[ ! -z "$bed_line" ]]; then
        match_found=true
        # Split the bed_line into columns
        IFS=$'\t' read -r chr start end base depth A B C D E F G H I J K <<< "$bed_line"
        echo "${id},${liftover_chr},${liftover_start},${liftover_end},${note_for_bedfile},${chr},${start},${end},${base},${depth},${A},${B},${C},${D},${E},${F},${G},${H},${I},${J},${K}" >> "$output_csv"
        break 2
      fi
    done
  done
  
  if [ "$match_found" = false ]; then
    # Print a message if no match is found and append to the log file
    echo "No match found for: ${id}, ${liftover_chr}, ${liftover_start}, ${liftover_end}" | tee -a "$log_file"
  fi
done

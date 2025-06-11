#!/bin/bash

# Define input and output directories
INPUT_DIR="../data/splitbychip"
OUTPUT_DIR="../data/SNPmiss005"

# List of chip names (match your filenames, without extensions)
chips=("unknown_chip" "HTS_iSelect_HD" "Illumina_GSAs" "OmniExpress" "OmniExpress_plus")

# Loop through each chip and run the PLINK command
for chip in "${chips[@]}"
do
    echo "Running PLINK on $chip..."
    plink --bfile "${INPUT_DIR}/${chip}" \
          --geno 0.05 \
          --make-bed \
          --out "${OUTPUT_DIR}/${chip}_geno005"
done

echo "All done!"

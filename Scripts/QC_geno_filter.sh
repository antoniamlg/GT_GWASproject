#!/bin/bash

# Define input and output directories
INPUT_DIR="../data/splitbychip"
OUTPUT_DIR="../data/SNPmiss1"

# List of chip names (match your filenames, without extensions)
chips=("unknown_chip_data" "HTS_iSelect_HD_data" "Illumina_GSAs_data" "OmniExpress_data" "OmniExpress_plus_data")

# Loop through each chip and run the PLINK command
for chip in "${chips[@]}"
do
    echo "Running PLINK on $chip..."
    plink --bfile "${INPUT_DIR}/${chip}" \
          --geno 1 \
          --make-bed \
          --out "${OUTPUT_DIR}/${chip}_geno1_filtered"
done

echo "All done!"

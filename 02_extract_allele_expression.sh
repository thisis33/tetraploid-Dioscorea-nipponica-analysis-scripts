#!/bin/bash

# Define input and output files; replace as needed
allele_file="two_allele_onetoone.txt" 
fpkm_file="hapABCD_average_fpkm_log2.txt"
output_file="two_allele_expression.txt"

# Create the output file and ensure any existing content is cleared
> "$output_file"

# Skip the header and process the allele_expression.txt file
tail -n +2 "$allele_file" | while IFS=$'\t' read -r AlleleA AlleleB AlleleC AlleleD; do
    # Create a new row
    new_row=""
    
    # Process each allele (A, B, C, D)
    for allele in "$AlleleA" "$AlleleB" "$AlleleC" "$AlleleD"; do
        # Use awk to find expression data for the gene ID in fpkm_file
        expression=$(awk -v gene="$allele" '$1 == gene {print $2, $3, $4, $5, $6, $7}' "$fpkm_file")
        
        # If expression data are found
        if [[ -n "$expression" ]]; then
            new_row="$new_row$allele $expression "
        else
            new_row="$new_row$allele NA NA NA NA NA NA "
        fi
    done
    
    # Remove trailing spaces and write the row to the output file
    echo -e "${new_row%"${new_row##*[![:space:]]}"}" >> "$output_file"
done

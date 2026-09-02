#!/usr/bin/bash

echo -n "Enter renamed FASTA output directory: "
read output_dir

echo
echo "========================================"
echo "Checking renamed FASTA files"
echo "========================================"

for gene_file in "$output_dir"/*.genes.fna
do
    [ -e "$gene_file" ] || continue

    sample=$(basename "$gene_file" .genes.fna)
    protein_file="$output_dir/${sample}.proteins.faa"

    echo
    echo "----------------------------------------"
    echo "Sample: $sample"
    echo "----------------------------------------"

    # Check gene FASTA
    gene_header=$(grep '^>' "$gene_file" | head -1)

    if [[ "$gene_header" == ">${sample}|"* ]]; then
        echo "GENE header : PASS"
        echo "  $gene_header"
    else
        echo "GENE header : FAIL"
        echo "  $gene_header"
    fi

    # Check protein FASTA
    if [ -f "$protein_file" ]; then
        protein_header=$(grep '^>' "$protein_file" | head -1)

        if [[ "$protein_header" == ">${sample}|"* ]]; then
            echo "PROTEIN header : PASS"
            echo "  $protein_header"
        else
            echo "PROTEIN header : FAIL"
            echo "  $protein_header"
        fi
    else
        echo "PROTEIN file : MISSING"
    fi

done

echo
echo "========================================"
echo "Checking completed."
echo "========================================"

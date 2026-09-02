#!/usr/bin/bash

echo -n "Enter Pyrodigal input directory: "
read input_dir

echo -n "Enter renamed FASTA output directory: "
read output_dir

echo
echo "========================================"
echo "Checking FASTA header renaming"
echo "========================================"

for orf in "$input_dir"/*_final.contigs_orf.fna
do
    [ -e "$orf" ] || continue

    sample=$(basename "$orf" _final.contigs_orf.fna)

    protein="$input_dir/${sample}_final.contigs.proteins.faa"

    gene_out="$output_dir/${sample}.genes.fna"
    protein_out="$output_dir/${sample}.proteins.faa"

    echo
    echo "----------------------------------------"
    echo "Sample: $sample"
    echo "----------------------------------------"

    # Check gene FASTA
    if [ -f "$gene_out" ]; then
        input_header=$(grep '^>' "$orf" | head -1)
        output_header=$(grep '^>' "$gene_out" | head -1)

        expected_header=$(echo "$input_header" | sed -E "s/^>([^[:space:]]+)/>${sample}|\1/")

        echo "GENE input    : $input_header"
        echo "GENE output   : $output_header"
        echo "GENE expected : $expected_header"

        if [ "$output_header" = "$expected_header" ]; then
            echo "GENE CHECK    : PASS"
        else
            echo "GENE CHECK    : FAIL"
        fi
    else
        echo "GENE output   : MISSING"
    fi

    # Check protein FASTA
    if [ -f "$protein" ] && [ -f "$protein_out" ]; then
        input_header=$(grep '^>' "$protein" | head -1)
        output_header=$(grep '^>' "$protein_out" | head -1)

        expected_header=$(echo "$input_header" | sed -E "s/^>([^[:space:]]+)/>${sample}|\1/")

        echo "PROTEIN input    : $input_header"
        echo "PROTEIN output   : $output_header"
        echo "PROTEIN expected : $expected_header"

        if [ "$output_header" = "$expected_header" ]; then
            echo "PROTEIN CHECK    : PASS"
        else
            echo "PROTEIN CHECK    : FAIL"
        fi
    else
        echo "PROTEIN output   : MISSING"
    fi

done

echo
echo "========================================"
echo "FASTA renaming check completed."
echo "========================================"

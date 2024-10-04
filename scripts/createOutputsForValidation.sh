#!/bin/bash

# Check if the user provided the number of iterations
if [ "$#" -ne 5 ]; then
    echo "Usage: $0 <number_of_iterations> <coverage> <path_to_genome> <genome_size> <path_to_read_length_distribution>"
    exit 1
fi

# Get the number of iterations from the command line
NUMITER=$1
COVERAGE=$2
GENOME=$3
GENOMESIZE=$4
DISTRIBUTION=$5

SEQREQUESTER=/path/to/seqrequester

# Loop n times
for (( i=1; i<=NUMITER; i++ )); do
    # 1. Run the seqrequester command
    $SEQREQUESTER simulate -truncate -genome $GENOME -genomesize $GENOMESIZE -coverage $COVERAGE -distribution $DISTRIBUTION > out1.fa
    
    # 2. Append the iteration header to output.txt
    echo "!BEGIN iteration=$i" >> outputHap1.txt
    
    # 3. Process out.fasta and append results to output.txt
    awk '/^>/ {
        # Extracting the read ID
        split($0, a, ",");
        for (i in a) {
            if (match(a[i], /read=([0-9]+)/, read_id)) {
                read = read_id[1];
            }
            if (match(a[i], /position=([0-9]+)-([0-9]+)/, positions)) {
                start = positions[1];
                end = positions[2];
            }
            if (match(a[i], /length=([0-9]+)/, len_data)) {
                len_value = len_data[1];
            }
        }
        
        # Print the extracted values
        printf "%s\t%s\t%s\t%s\n", read, start, end, len_value;
    }' out1.fa >> outputHap1.txt
    
    # 4. Append the iteration footer to output.txt
    echo "!END" >> outputHap1.txt
done

for (( i=1; i<=NUMITER; i++ )); do
    # 1. Run the seqrequester command
    $SEQREQUESTER simulate -truncate -genome $GENOME -genomesize $GENOMESIZE -coverage $COVERAGE -distribution $DISTRIBUTION > out2.fa

    # 2. Append the iteration header to output.txt
    echo "!BEGIN iteration=$i" >> outputHap2.txt

    # 3. Process out.fasta and append results to output.txt
    awk '/^>/ {
        # Extracting the read ID
        split($0, a, ",");
        for (i in a) {
            if (match(a[i], /read=([0-9]+)/, read_id)) {
                read = read_id[1];
            }
            if (match(a[i], /position=([0-9]+)-([0-9]+)/, positions)) {
                start = positions[1];
                end = positions[2];
            }
            if (match(a[i], /length=([0-9]+)/, len_data)) {
                len_value = len_data[1];
            }
        }

        # Print the extracted values
        printf "%s\t%s\t%s\t%s\n", read, start, end, len_value;
    }' out2.fa >> outputHap2.txt

    # 4. Append the iteration footer to output.txt
    echo "!END" >> outputHap2.txt
done

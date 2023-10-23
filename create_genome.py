import sys
import random

# Check for the correct number of command-line arguments
if len(sys.argv) != 2:
    print("Usage: python3 create_genome.py <sequence_length>")
    sys.exit(1)

# Get the sequence length from the command-line argument
sequence_length = int(sys.argv[1])

# Generate a random DNA sequence
sequence = ''.join(random.choice('ATCG') for _ in range(sequence_length))

# Create and write the FASTA file
with open('simGenome.fasta', 'w') as fasta_file:
    # Write the header
    fasta_file.write(">chr1\n")
    # Write the sequence
    fasta_file.write(sequence + '\n')

print(f"Simulated genome with length {sequence_length} written to simGenome.fasta.")
## <a name="intro"></a>Introduction

CGProb estimates the probability of observing a coverage gap after a heterozygous locus on the second haplotype of a circular diploid genome. CGProb takes the genome length, coverage on each haplotype, the sequencing read length distribution as input. It estimates the probability of observing a coverage gap after a heterozygous locus on the second haplotype by counting the number of read sequencing outputs which have a coverage gap and dividing it by the total number of read sequencing outputs.

For ease of computation CGProb allows users to scale down the longer read lengths by a constant factor. This factor can be chosen by the user. 

## <a name="started"></a>How to Use

```sh
# Install seqrequester
# Dependencies: GCC (tested with 11.3.0), git 2.25.1 or higher, gmp library for c++

cd /install/dir/for/seqrequester

git clone https://github.com/marbl/seqrequester.git
git checkout 31141c14d80fe0dde2766bcc3e622bf600e2bba4

cd src && make -j 12

# Install CGProb
# Dependencies: GCC 11.2.0 or higher, C++ version standard c++2b
cd /install/dir/for/CGProb

git clone https://github.com/ssk5e/covGapProb.git
cd covGapProb && make

# Create artifical genome for seqrequester. 
# Input desired length of chromosome after createGenome.py 
# Output in "simGenome.fasta"

python3 /scripts/create_genome.py 1000000

# Simulate and sumarize reads from read length distribution for two haplotypes

/install/dir/for/seqrequester/build/bin/seqrequester simulate \
-truncate \
-genome simGenome.fasta \
-genomesize 1000000 \
-coverage 50 \
-distribution readLengthDistribution.txt \
> readsHap1.fasta

/install/dir/for/seqrequester/build/bin/seqrequester summarize -simple readsHap1.fasta > originalReadLengthsHap1.txt

/install/dir/for/seqrequester/build/bin/seqrequester simulate \
-truncate \
-genome simGenome.fasta \
-genomesize 1000000 \
-coverage 50 \
-distribution readLengthDistribution.txt \
> readsHap2.fasta

/install/dir/for/seqrequester/build/bin/seqrequester summarize -simple readsHap2.fasta > originalReadLengthsHap2.txt

# Choose scaling factor to scale down read lengths
# For computation new genome size and read lengths will be ceiling value of initial_size / scaling_factor

python3 /scripts/scale_down.py originalReadLengthsHap1.txt shortLengthsHap1.txt 1000
python3 /scripts/scale_down.py originalReadLengthsHap2.txt shortLengthsHap2.txt 1000

# Run computation
tail -2 readsHap1.fasta | head -1 | awk -F'[=,]' '{print $2}' > readCountHap1.txt
tail -2 readsHap2.fasta | head -1 | awk -F'[=,]' '{print $2}' > readCountHap2.txt

/install/dir/for/CGProb/compute \
-g 1000 \
-R $(cat readCountHap1.txt) \
-r $(cat readCountHap2.txt) \
-h 200 \
-D shortLengthsHap1.txt \
-d shortLengthsHap2.txt \
-p 128 \
-t 256 \
> output.txt

tail output.txt | grep "probability = " | awk '{print $3}' > probability.txt
```

The output is stored in `probability.txt`. It is recommended that genome size is set to about 1,000 bp.

## <a name="use"></a>Usage Details

The required inputs are:
1. Genome Length [INT]
    - Artificial Genome (created using create_genome.py)
2. Read Length Distribution of sequencing instrument
    - Format: read_length number_of_reads
3. Coverage on first haplotype [INT]
4. Coverage on second haplotype [INT]
5. Reads on haplotypes
    - Simulated using seqrequester
    - readsHap1.fasta, originalReadLengthsHap1.txt, readCountHap1.txt
    - readsHap2.fasta, originalReadLengthsHap2.txt, readCountHap1.txt
6. Scaling Factor [INT]
    - Ideally, genome length is a multiple of this parameter, not necessary
    - shortLengthsHap1.txt, generated using scale_down.py
    - shortLengthsHap2.txt, generated using scale_down.py

Inputs to the executable `compute` are:

    -g INT: genome size
    -R INT: number of reads on first haplotype
    -r INT: number of reads on second haplotype
    -h INT: location of heterozygous locus; 1 < heterozygous locus < genome size;
    -D INT: read length distribution on first haplotype
    -d INT: read legnth distribution on second haplotype
    -p INT: number of bits of precision
    -t INT: number of threads

Output File:
Currently, the executable `compute` prints to stdout. This can be redirected as evinced in the How to Use section. The last line in `output.txt` contains the probability of observing a coverage gap on the second haplotype near a heterozygous locus due to contained read deletion.

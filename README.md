## <a name="intro"></a>Introduction

The [string graph formulation](https://doi.org/10.1093/bioinformatics/bti1114) for genome assembly involves deleting reads that are entirely contained in longer reads. This is known to be an unsafe operation because this heuristic [occasionally disconnects the walks](https://doi.org/10.1093/bioinformatics/btad124) corresponding to true chromosome sequences.  

CGProb computes the frequency of the occurrence of an assembly gap due to contained read deletion. CGProb takes the genome length, coverage on each haplotype, and read-length distribution as input. It computes the frequency of the occurrence of an assembly gap after a heterozygous locus due to contained read deletion by counting the number of read sequencing outputs where an assembly gap appears after contained read deletion and dividing it by the total number of read sequencing outputs where the heterozygous locus is sequenced on each haplotype.

The size of the sample space, i.e., the number of all possible read sequencing outputs grows combinatorially with the number of distinct read lengths, coverage and genome length. Simulating and testing each individual read sequencing output for an assembly gap due to contained read deletion would be a Sisyphean endeavour. Instead, CGProb uses efficient partitioning of the sample space and [ordinary generating functions](https://en.wikipedia.org/wiki/Generating_function) to calculate the frequency in polynomial time. See our [preprint](#papers) to know more about the method.

## <a name="install"></a>Installation

1. Install [seqrequester](https://github.com/marbl/seqrequester)

Dependencies: GCC (tested with 11.3.0), git 2.25.1 or higher, gmp library for c++

```
git clone https://github.com/marbl/seqrequester.git
cd seqrequester
git checkout 31141c14d80fe0dde2766bcc3e622bf600e2bba4
cd src && make -j 12
```

2. Install gmp (https://anaconda.org/conda-forge/gmp)

3. Activate an environment containing gmp

```
conda activate <myenv>
```

3. Install CGProb

Dependencies: gcc version 11.2.0 or higher; GMP version 6.2.1 or higher

```
git clone https://github.com/at-cg/CGProb.git
cd CGProb
python3 create_makefile.py /home/anaconda3/envs/myenv/include /home/anaconda3/envs/myenv/lib > makefile
make
```

## <a name="started"></a>How to Use

```sh
# Activate conda environment containing gmp (https://anaconda.org/conda-forge/gmp)
conda activate myenv

# Create an artificial genome for seqrequester. 
# Input desired length of chromosome after createGenome.py 
# Output in "simGenome.fasta"

python3 create_genome.py 1000000

# Simulate and summarise reads from read length distribution for two haplotypes

/install/dir/for/seqrequester/build/bin/seqrequester simulate \
-circular \
-genome simGenome.fasta \
-genomesize 1000000 \
-coverage 10 \
-distribution CGProb/distributions/fixed.txt \
> readsHap1.fasta

/install/dir/for/seqrequester/build/bin/seqrequester summarize -simple readsHap1.fasta > originalReadLengthsHap1.txt

/install/dir/for/seqrequester/build/bin/seqrequester simulate \
-circular \
-genome simGenome.fasta \
-genomesize 1000000 \
-coverage 10 \
-distribution CGProb/distributions/fixed.txt \
> readsHap2.fasta

/install/dir/for/seqrequester/build/bin/seqrequester summarize -simple readsHap2.fasta > originalReadLengthsHap2.txt

# Choose a scaling factor to scale down read lengths
# For computation, new genome size and read lengths will be the ceiling value of initial_size / scaling_factor

python3 scale_down.py originalReadLengthsHap1.txt shortLengthsHap1.txt 1000
python3 scale_down.py originalReadLengthsHap2.txt shortLengthsHap2.txt 1000

# Run computation
tail -2 readsHap1.fasta | head -1 | awk -F'[=,]' '{print $2}' > readCountHap1.txt
tail -2 readsHap2.fasta | head -1 | awk -F'[=,]' '{print $2}' > readCountHap2.txt

/install/dir/for/CGProb/compute \
-R $(cat readCountHap1.txt) -r $(cat readCountHap2.txt) -D shortLengthsHap1.txt -d shortLengthsHap2.txt > output.txt

tail output.txt | grep "probability = " | awk '{print $3}' > frequency.txt
```

The output is stored in `frequency.txt`. It is recommended that the initial genome size be set to about 1 Mbp. For ease of computation, CGProb allows users to scale down the longer read lengths by a constant factor. This factor can be chosen by the user. We recommend scaling read lengths down by a factor of 1000.

## <a name="use"></a>Usage Details

The following options can be used to customise the behaviour of CGProb:

1. Genome Length [INT]
    - Artificial Genome (created using create_genome.py)
2. Read Length Distribution of the sequencing instrument
    - Format: read_length number_of_reads
3. Coverage on the first haplotype [INT]
4. Coverage on the second haplotype [INT]
5. Reads on haplotypes
    - Simulated using seqrequester
    - readsHap1.fasta, originalReadLengthsHap1.txt, readCountHap1.txt
    - readsHap2.fasta, originalReadLengthsHap2.txt, readCountHap1.txt
6. Scaling Factor [INT]
    - Ideally, genome length is a multiple of this parameter, not necessary
    - shortLengthsHap1.txt, generated using scale_down.py
    - shortLengthsHap2.txt, generated using scale_down.py

Inputs to the executable `compute` are:

```
-g INT: genome size [1000]
-R INT: number of reads on the first haplotype (required input)
-r INT: number of reads on the second haplotype (required input)
-h INT: location of heterozygous locus [200]; 1 <= heterozygous locus <= genome size / 2
-D INT: read length distribution on the first haplotype (required input)
-d INT: read length distribution on the second haplotype (required input)
-p INT: number of bits of precision [128]
-t INT: number of threads [32]
```

Output File:
Currently, the executable `compute` prints to stdout. This can be redirected as evinced in the How to Use section. The last line in `output.txt` contains the frequency of observing an assembly gap after a heterozygous locus due to contained read deletion.

## <a name="examples"></a>Examples

1. Computing the frequency of observing an assembly gap after a heterozygous locus due to contained read deletion: Coverage 15x on each haplotype; PacBio HiFi read length distribution; Chromosome of size 100000; Reads scaled down by a factor of 100; Precision 256 bits; Number of threads 128

```
python3 create_genome.py 100000
seqrequester simulate -circular -genome simGenome.fasta -genomesize 100000 -coverage 15 -distribution CGProb/distributions/pacbio_hifi.txt \
> readsHap1.fasta
seqrequester summarize -simple readsHap1.fasta > originalReadLengthsHap1.txt

seqrequester simulate -circular -genome simGenome.fasta -genomesize 100000 -coverage 15 -distribution CGProb/distributions/pacbio_hifi.txt \
> readsHap2.fasta
seqrequester summarize -simple readsHap2.fasta > originalReadLengthsHap2.txt

# Scaling factor 100; compressed genome size 1000.
python3 scale_down.py originalReadLengthsHap1.txt shortLengthsHap1.txt 100
python3 scale_down.py originalReadLengthsHap2.txt shortLengthsHap2.txt 100
tail -2 readsHap1.fasta | head -1 | awk -F'[=,]' '{print $2}' > readCountHap1.txt
tail -2 readsHap2.fasta | head -1 | awk -F'[=,]' '{print $2}' > readCountHap2.txt

CGProb/compute -R $(cat readCountHap1.txt) -r $(cat readCountHap2.txt) -D shortLengthsHap1.txt -d shortLengthsHap2.txt -p 256 -t 128 > output.txt
tail output.txt | grep "probability = " | awk '{print $3}' > frequency.txt
```

2. Computing the frequency of observing an assembly gap after a somatic mutation locus: Coverage 20x on haplotype 1 and 5x on haplotype 2; ONT Simplex read length distribution; Chromosome of size 1000000; Reads scaled down by a factor of 10000; Somatic mutation locus 50; Precision 256 bits; Number of threads 64.

```
python3 create_genome.py 1000000
seqrequester simulate -circular -genome simGenome.fasta -genomesize 1000000 -coverage 20 \
-distribution CGProb/distributions/nanopore_simplex.txt > readsHap1.fasta
seqrequester summarize -simple readsHap1.fasta > originalReadLengthsHap1.txt

seqrequester simulate -circular -genome simGenome.fasta -genomesize 100000 -coverage 5 \
-distribution CGProb/distributions/nanopore_simplex.txt > readsHap2.fasta
seqrequester summarize -simple readsHap2.fasta > originalReadLengthsHap2.txt

# Scaling factor 10000; compressed genome size 100.
python3 scale_down.py originalReadLengthsHap1.txt shortLengthsHap1.txt 100
python3 scale_down.py originalReadLengthsHap2.txt shortLengthsHap2.txt 100
tail -2 readsHap1.fasta | head -1 | awk -F'[=,]' '{print $2}' > readCountHap1.txt
tail -2 readsHap2.fasta | head -1 | awk -F'[=,]' '{print $2}' > readCountHap2.txt

CGProb/compute -g 100 -h 50 \
-R $(cat readCountHap1.txt) -r $(cat readCountHap2.txt) \
-D shortLengthsHap1.txt -d shortLengthsHap2.txt \
-p 256 -t 128 > output.txt
tail output.txt | grep "probability = " | awk '{print $3}' > frequency.txt
```

## <a name="papers"></a>Preprint
Sudhanva Shyam Kamath, Mehak Bindra, Debnath Pal, Chirag Jain. [Telomere-to-telomere assembly by preserving contained reads](https://doi.org/10.1101/2023.11.07.565066). Biorxiv (November 2023).

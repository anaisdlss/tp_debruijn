# De Bruijn Graph-Based Assembler

## Introduction

This project implements a simplified genome assembler based on a De Bruijn graph, a method widely used in modern assembly pipelines (e.g., Velvet, SPAdes).

The goal is to reconstruct the genome of a virus (Enterovirus A71) from short Illumina reads.
The pipeline allows you to:
- generate k-mers and build the De Bruijn graph
- clean the graph (tips, bubbles)
- extract contigs
- compare these contigs to a reference genome using BLAST

This practical exercise illustrates how modern assemblers work internally: how they manage sequencing errors, how graph structures represent genome overlaps, and how complete sequences can be reconstructed from fragmented reads.


## Prerequisites

You must have **Python 3.10+** installed:

```bash
python3 --version
```

You also need the **uv** environment manager:
```
curl -LsSf https://astral.sh/uv/install.sh | sh
uv --version
```
And you must have **BLAST+** installed:
```
# Linux
sudo apt install ncbi-blast+

# macOS (Homebrew)
brew install blast              
```

## Usage


First, clone the repository and move into the created directory:
```
git clone https://github.com/anaisdlss/tp_debruijn.git
cd tp_debruijn
```

Synchronize the environment:
```
uv sync
```
Create an output directory:
```
mkdir resultats
```

Then run the assembler script:
```
uv run python debruijn/debruijn.py \
    -i data/eva71_plus_perfect.fq \
    -k 22 \
    -o resultats/contigs.fasta \
    -f resultats/graph.png
```
- -i : input FASTQ file containing sequencing reads. You may replace eva71_plus_perfect.fq with any FASTQ dataset you want to assemble.
- -k : k-mer size used to build the De Bruijn graph (default = 22).
Larger k-mers give more specific assemblies; smaller ones are more error-tolerant.
- -o : output FASTA file where the assembled contigs will be saved.
- -f : optional PNG image showing the cleaned De Bruijn graph.

The generated contigs and graph are saved respectively in:
- resultats/contigs.fasta
- resultats/graph.png



After that, create a BLAST database from the reference data/eva71.fna:
```
makeblastdb -in data/eva71.fna -dbtype nucl
```
- makeblastdb : creates a BLAST-searchable database from a reference sequence.
- -in data/eva71.fna : input reference genome (you can replace it with any FASTA file).
- -dbtype nucl : specifies that the database is nucleotide (use prot for protein databases).

Then compare your assembled contigs (resultats/contigs.fasta) to the reference genome:
```
blastn -query resultats/contigs.fasta -db data/eva71.fna -out resultats/blast_result.txt -outfmt 6
```
- blastn compares your assembled nucleotide contigs to the reference genome.
- -query resultats/contigs.fasta is the FASTA file containing the contigs produced by the assembler.
- -db data/eva71.fna is the BLAST database created from the reference genome.
- -out resultats/blast_result.txt saves the alignment results to a text file.
- -outfmt 6 outputs results in a simple tabular format (easy to parse or inspect).

This step verifies how well the reconstructed contigs align with the known reference, confirming that the assembly reflects the genome from which the reads originally came.

The BLAST results are written to:```resultats/blast_result.txt```.

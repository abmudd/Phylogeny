#!/bin/bash

# Make directory
mkdir -p prep/

# Download all amphibian genes
scripts/esearch -db nucleotide -query "txid8342[Organism] biomol_genomic[PROP]" | scripts/efetch -format gb > prep/NCBI_full.gb
scripts/esearch -db nucleotide -query "txid8296[Organism] biomol_genomic[PROP]" | scripts/efetch -format gb >> prep/NCBI_full.gb

# Initial query for gene names/counts and feature counts
scripts/allgenesingb.py prep/NCBI_full.gb

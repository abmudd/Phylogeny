#!/usr/bin/env python

import argparse, os, subprocess, sys
from Bio import SeqIO

# Parse arguments
parser = argparse.ArgumentParser(description='This script finds any gene names that are not in the key files, allowing the user to check for potential mistaken gene names in the input and to expand the key files with other gene names.')
parser.add_argument("--output", metavar='STR', help="output prefix [input file name]", type=str)
required = parser.add_argument_group('required arguments')
required.add_argument("gene_key", help="gene names key", type=str)
required.add_argument("input", help="input gb file", type=str)
args = parser.parse_args()

# Open input file and set output names
in_gb = open(args.input, 'r')
if not args.output:
    args.output = args.input
out_gene_name = args.output+'.gene.names'
if os.path.exists(out_gene_name):
    subprocess.call(['rm', out_gene_name])

# Set variables
list_features = ['gene','CDS','mRNA','rRNA','tRNA']
gene_names_seen = []
gene_names = []

# Set existing and ignore names
for line in open(args.gene_key, 'r'):
    if not line.startswith('GB_name') and ';' in line:
        line = line.upper().rstrip().split(';')
        gene_names.append(line[0])

# Check if gene name is in the key file or already seen; if not, then write out id and name
def writename(test_id, test_name, outfile_name, list_name, list_name_seen):
    if test_name not in list_name and test_name not in list_name_seen:
        outfile = open(outfile_name, 'a')
        outfile.write(test_id+'\t'+test_name+'\n')
        outfile.close()
        list_name_seen.append(test_name)
    return list_name_seen

# Parse input gb file
for record in SeqIO.parse(in_gb, "gb"):
    for feature in record.features:
        if feature.type in list_features:
            if 'gene' in feature.qualifiers:
                gene_names_seen = writename(record.id, feature.qualifiers['gene'][0].upper(), out_gene_name, gene_names, gene_names_seen)
            elif 'product' in feature.qualifiers:
                gene_names_seen = writename(record.id, feature.qualifiers['product'][0].upper(), out_gene_name, gene_names, gene_names_seen)

# Close input file
in_gb.close()

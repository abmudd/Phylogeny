#!/usr/bin/env python

# Copyright (c)2017. The Regents of the University of California (Regents).
# All Rights Reserved. Permission to use, copy, modify, and distribute this
# software and its documentation for educational, research, and
# not-for-profit purposes, without fee and without a signed licensing
# agreement, is hereby granted, provided that the above copyright notice,
# this paragraph and the following two paragraphs appear in all copies,
# modifications, and distributions. Contact the Office of Technology
# Licensing, UC Berkeley, 2150 Shattuck Avenue, Suite 510, Berkeley, CA
# 94720-1620, (510) 643-7201, for commercial licensing opportunities.

# IN NO EVENT SHALL REGENTS BE LIABLE TO ANY PARTY FOR DIRECT, INDIRECT,
# SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES, INCLUDING LOST PROFITS,
# ARISING OUT OF THE USE OF THIS SOFTWARE AND ITS DOCUMENTATION, EVEN IF
# REGENTS HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGE

# REGENTS SPECIFICALLY DISCLAIMS ANY WARRANTIES, INCLUDING, BUT NOT LIMITED
# TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
# PURPOSE. THE SOFTWARE AND ACCOMPANYING DOCUMENTATION, IF ANY, PROVIDED
# HEREUNDER IS PROVIDED "AS IS". REGENTS HAS NO OBLIGATION TO PROVIDE
# MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR MODIFICATIONS.


###############################################################################
# Setup
###############################################################################

# Import modules
# ============================================================

import argparse, os, sys
from Bio import SeqIO


# Parse arguments
# ============================================================

parser = argparse.ArgumentParser(description='This script reads through a list of GenBank entries, extracts the gene entry for each species matching the key with (1) the longest gene sequence and (2) the longest total sequence, and then writes out the fasta file for the gene. GenBank entries are skipped if no gene entries matching the gene sequence are found. All entries in the key are assumed to be from the same gene.')
parser.add_argument("-s", "--silent", help="run script in silent mode", action='store_true')
parser.add_argument("-l", "--log", metavar='STR', help="output file name for log of species extracted [stderr]", type=str, default='stderr')
parser.add_argument("-o", "--output", metavar='STR', help="output prefix for fasta files [sequence]", type=str, default='sequence')
parser.add_argument("-v", "--version", help="show version info and exit", action='version', version='%(prog)s 0.1')
required = parser.add_argument_group('required arguments')
required.add_argument("gene_key", help="gene names key", type=str)
required.add_argument("input", help="input gb file", type=str)
args = parser.parse_args()


# Set input and output files
# ============================================================

in_gb = open(args.input, 'r')
out_log = sys.stderr
if args.log != 'stderr':
    out_log = open(args.log, 'w')


# Set variables
# ============================================================

gene_names = []
gene_dict = {}
full_len_dict = {}
gene_len_dict = {}
ids_dict = {}
id_seen = []
species_seen = []
final_species_seen = []
list_features = ['gene','CDS','mRNA','rRNA','tRNA']


# Function to extract gene fasta
# ============================================================

def geneextract(temp_count, file_name, test_name, fasta_header, test_seq, start_pos, stop_pos, strand_value):
    temp_count += 1
    outfile = open(file_name+"."+test_name+".fa", 'w')
    if int(strand_value) < 0:
        gene_seq = test_seq[int(start_pos):int(stop_pos)].reverse_complement()
    else:
        gene_seq = test_seq[int(start_pos):int(stop_pos)]
    outfile.write('>'+fasta_header+'\n')
    outfile.write(str(gene_seq)+'\n')
    outfile.close()
    return temp_count


###############################################################################
# Run
###############################################################################

# Set existing and ignore names
# ============================================================

First_line = True
for line in open(args.gene_key, 'r'):
    if not line.startswith('GB_name'):
        line = line.upper().rstrip()
        if First_line:
            initial_name = line
            First_line = False
        gene_names.append(line)
        gene_dict[line] = initial_name


# Parse input gb file
# ============================================================

for record in SeqIO.parse(in_gb, "gb"):

    # Skip entries with duplicate accession numbers.
    # ============================================================

    if record.id in id_seen:
        if not args.silent:
            out_log.write("Skipping duplicate entry for accession "+record.id+".\n")

    # Skip entries without a species name.
    # ============================================================

    elif 'organism' not in record.annotations:
        if not args.silent:
            out_log.write("Skipping entry for accession "+record.id+" as species name is missing.\n")

    # Skip entries with cf., sp., aff., or environmental in the species name.
    # ============================================================

    elif 'cf.' in record.annotations['organism'] or 'sp.' in record.annotations['organism'] or 'aff.' in record.annotations['organism'] or 'environmental' in record.annotations['organism']:
        if not args.silent:
            out_log.write("Skipping entry for species "+record.annotations['organism']+" with accession "+record.id+".\n")

    # Skip entries that only contain N's.
    # ============================================================

    elif record.seq.count('N') == len(record.seq):
        if not args.silent: 
            out_log.write("Skipping entry for accession "+record.id+" as the sequence only contains N's.\n")

    else:
        id_seen.append(record.id)
        species_name = '_'.join(record.annotations['organism'].split(' ')[0:2])

        # Count number of annotations matching genes
        # ============================================================

        query_length = False
        for feature in record.features:
            if feature.type in list_features:
                if 'gene' in feature.qualifiers and feature.qualifiers['gene'][0].upper() in gene_names:
                    query_length = feature.location.nofuzzy_end - feature.location.nofuzzy_start
                    break
                elif 'product' in feature.qualifiers and feature.qualifiers['product'][0].upper() in gene_names:
                    query_length = feature.location.nofuzzy_end - feature.location.nofuzzy_start
                    break

        # Replace old entry if new has (1) a longer gene sequence and (2) a longer total sequence
        # ============================================================

        length_full_seq_query = len(record.seq)
        if species_name in species_seen and (query_length > gene_len_dict[species_name] or (query_length == gene_len_dict[species_name] and length_full_seq_query > full_len_dict[species_name])):
            gene_len_dict[species_name] = query_length
            full_len_dict[species_name] = len(record.seq)
            ids_dict[species_name] = record.id

        # Add new entry if species has not been previously seen
        # ============================================================

        elif species_name not in species_seen and query_length:
            species_seen.append(species_name)
            gene_len_dict[species_name] = query_length
            full_len_dict[species_name] = len(record.seq)
            ids_dict[species_name] = record.id


# Write list of final accession numbers for gene extraction
# ============================================================

final_ids = []
out_log.write("Extracting a total of "+str(len(species_seen))+" species.\n")
for species_name in species_seen:
    final_ids.append(ids_dict[species_name])


# Go to top of input gb file
# ============================================================

in_gb.seek(0)
count = 0


# Parse input gb file
# ============================================================

for record in SeqIO.parse(in_gb, "gb"):

    # If accession number is in the list of final accession numbers for gene extraction
    # ============================================================

    if record.id in final_ids:
        count += 1
        species_name = '_'.join(record.annotations['organism'].split(' ')[0:2])
        out_log.write(str(count)+'\t'+species_name+'\t'+ids_dict[species_name]+'\n')

        # Write out the fasta sequence for each gene
        # ============================================================

        temp_features = 0
        for feature in record.features:
            run_geneextract = False
            if feature.type in list_features:
                if 'gene' in feature.qualifiers and feature.qualifiers['gene'][0].upper() in gene_names:
                    temp_gene = feature.qualifiers['gene'][0].upper()
                    if temp_gene in gene_dict:
                        temp_gene = gene_dict[temp_gene]
                    run_geneextract = True
                elif 'product' in feature.qualifiers and feature.qualifiers['product'][0].upper() in gene_names:
                    temp_gene = feature.qualifiers['product'][0].upper()
                    if temp_gene in gene_dict:
                        temp_gene = gene_dict[temp_gene]
                    run_geneextract = True
            
            # Running the extract gene fasta function
            # ============================================================

            if run_geneextract:
                temp_features = geneextract(temp_features, args.output+str(count), temp_gene, species_name, record.seq, feature.location.start, feature.location.end, feature.strand)
                break


# Close input and output files
# ============================================================

in_gb.close()
if args.log != 'stderr':
    out_log.close()

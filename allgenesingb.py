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
from collections import defaultdict


# Parse arguments
# ============================================================

parser = argparse.ArgumentParser(description='This script finds all synonymous gene names (output.genes_names), counts all of the gene names (output.genes_count), and counts all of the feature types (output.feature_count) in a gb file.')
parser.add_argument("-o", "--output", metavar='STR', help="output prefix [input file name]", type=str)
parser.add_argument("-v", "--version", help="show version info and exit", action='version', version='%(prog)s 0.1')
required = parser.add_argument_group('required arguments')
required.add_argument("input", help="input gb file", type=str)
args = parser.parse_args()


# Set input and output names
# ============================================================

in_gb = open(args.input, 'r')
if not args.output:
    args.output = args.input
out_gene_name = open(args.output+'.genes_name', 'w')
out_gene_count = open(args.output+'.genes_count', 'w')
out_feature_count = open(args.output+'.feature_count', 'w')


# Set variables
# ============================================================

list_features = ['gene','CDS','mRNA','rRNA','tRNA']
seen_gene_name_dictlist = defaultdict(list)
primary_gene_name_list = []
alternate_gene_name_list = []
unable_to_join_list = []
genes_total_seen_count_dict = {}
species_seen_list = []
species_seen_count = 0
species_genes_seen_dictlist = defaultdict(list)
record_id_seen_list = []
feature_count_dict = {}


# Function to search dictionary for values and export key
# ============================================================

def search_dict(test_dict, test_value):
    for key, value in test_dict.iteritems():
        if test_value in value:
            return key


# Function to check if gene name is in the key file or already seen; if not, then write out id and name
# ============================================================

def writename(test_name, names_seen, count_dict, species, species_dict):
    if test_name not in names_seen and test_name not in species_dict[species]:
        if test_name in count_dict:
            count_dict[test_name] = count_dict[test_name] + 1
        else:
            count_dict[test_name] = 1
        names_seen.append(test_name)
        species_dict[species].append(test_name)
    return names_seen, count_dict, species_dict


###############################################################################
# Run
###############################################################################

# Parse input gb file
# ============================================================

for record in SeqIO.parse(in_gb, "gb"):

    # Skip records already seen, without a species name, or with cf., sp., or aff. in species name
    # ============================================================

    if record.id in record_id_seen_list or 'organism' not in record.annotations or 'cf.' in record.annotations['organism'] or 'sp.' in record.annotations['organism'] or 'aff.' in record.annotations['organism']:
        continue

    else:
        old_start = False
        old_end = False
        old_strand = False
        old_name = False
        name = False
        genes_record_seen_list = []
        record_id_seen_list.append(record.id)

        # Count number of species
        # ============================================================

        species_name = '_'.join(record.annotations['organism'].split(' ')[0:2])
        if not species_name in species_seen_list:
            species_seen_count += 1
            species_seen_list.append(species_name)

        # Parse through annotated features
        # ============================================================

        for feature in record.features:
            next = False
            if feature.type in feature_count_dict:
                feature_count_dict[feature.type] += 1
            else:
                feature_count_dict[feature.type] = 1
            if feature.type in list_features:
                if 'gene' in feature.qualifiers:
                    name = feature.qualifiers['gene'][0].upper()
                    genes_record_seen_list, genes_total_seen_count_dict, species_genes_seen_dictlist = writename(name, genes_record_seen_list, genes_total_seen_count_dict, species_name, species_genes_seen_dictlist)
                    if name in primary_gene_name_list: # Primary gene name; previously seen
                        continue
                    elif name in alternate_gene_name_list: # Alternate gene name; previously seen
                        name = search_dict(seen_gene_name_dictlist, name)
                    else: # Not previously seen
                        if feature.location.start == old_start and feature.location.end == old_end and feature.location.strand == old_strand and name != old_name: # Matches last parsed annotated feature
                            if old_name not in primary_gene_name_list:
                                primary_gene_name_list.append(old_name)
                            if name not in alternate_gene_name_list:
                                alternate_gene_name_list.append(name)
                            seen_gene_name_dictlist[old_name].append(name)
                            name = old_name
                        else: # Doesn't match last parsed annotated feature
                            unable_to_join_list.append(name)
                    old_start = feature.location.start
                    old_end = feature.location.end
                    old_strand = feature.location.strand
                    old_name = name
                elif 'product' in feature.qualifiers:
                    name = feature.qualifiers['product'][0].upper()
                    genes_record_seen_list, genes_total_seen_count_dict, species_genes_seen_dictlist = writename(name, genes_record_seen_list, genes_total_seen_count_dict, species_name, species_genes_seen_dictlist)
                    if name in primary_gene_name_list: # Primary gene name; previously seen
                        continue
                    elif name in alternate_gene_name_list: # Alternate gene name; previously seen
                        name = search_dict(seen_gene_name_dictlist, name)
                    else: # Not previously seen
                        if feature.location.start == old_start and feature.location.end == old_end and feature.location.strand == old_strand and name != old_name: # Matches last parsed annotated feature
                            if old_name not in primary_gene_name_list:
                                primary_gene_name_list.append(old_name)
                            if name not in alternate_gene_name_list:
                                alternate_gene_name_list.append(name)
                            seen_gene_name_dictlist[old_name].append(name)
                            name = old_name
                        else: # Doesn't match last parsed annotated feature
                            unable_to_join_list.append(name)
                    old_start = feature.location.start
                    old_end = feature.location.end
                    old_strand = feature.location.strand
                    old_name = name


# Write out seen_gene_name_dictlist, species count, gene names, and percent of species with each gene name
# ============================================================

for key, value in sorted(seen_gene_name_dictlist.iteritems(), key=lambda (k,v): (k,v)):
    out_gene_name.write(key+'\t'+'\t'.join(value)+'\n')
out_gene_count.write("A total of "+str(species_seen_count)+" species were counted.\n")
for key in genes_total_seen_count_dict:
    genes_total_seen_count_dict[key] = int(genes_total_seen_count_dict[key]*100/species_seen_count)
for key, value in sorted(genes_total_seen_count_dict.iteritems(), key=lambda (k,v): (-v,k)):
    out_gene_count.write(str(value)+'\t'+key+'\n')
for key, value in sorted(feature_count_dict.iteritems(), key=lambda (k,v): (-v,k)):
    out_feature_count.write(key+'\t'+str(value)+'\n')


# Close input and output files
# ============================================================

in_gb.close()
out_gene_name.close()
out_gene_count.close()
out_feature_count.close()

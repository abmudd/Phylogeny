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

import argparse, os, subprocess, sys
from Bio import SeqIO


# Parse arguments
# ============================================================

parser = argparse.ArgumentParser(description='This script finds any gene names that are not in the key files, allowing the user to check for potential mistaken gene names in the input and to expand the key files with other gene names.')
parser.add_argument("-o", "--output", metavar='STR', help="output prefix [input file name]", type=str)
parser.add_argument("-v", "--version", help="show version info and exit", action='version', version='%(prog)s 0.1')
required = parser.add_argument_group('required arguments')
required.add_argument("gene_key", help="gene names key", type=str)
required.add_argument("input", help="input gb file", type=str)
args = parser.parse_args()


# Set input and output files
# ============================================================

in_gb = open(args.input, 'r')
if not args.output:
    args.output = args.input
out_gene_name = args.output+'.gene.names'
if os.path.exists(out_gene_name):
    subprocess.call(['rm', out_gene_name])


# Set variables
# ============================================================

list_features = ['gene','CDS','mRNA','rRNA','tRNA']
gene_names_seen = []
gene_names = []


# Function to check if gene name is in the key file or already seen; if not, then write out id and name
# ============================================================

def writename(test_id, test_name, outfile_name, list_name, list_name_seen):
    if test_name not in list_name and test_name not in list_name_seen:
        outfile = open(outfile_name, 'a')
        outfile.write(test_id+'\t'+test_name+'\n')
        outfile.close()
        list_name_seen.append(test_name)
    return list_name_seen


###############################################################################
# Run
###############################################################################

# Set existing and ignore names
# ============================================================

for line in open(args.gene_key, 'r'):
    if not line.startswith('GB_name') and ';' in line:
        line = line.upper().rstrip().split(';')
        gene_names.append(line[0])


# Parse input gb file
# ============================================================

for record in SeqIO.parse(in_gb, "gb"):
    for feature in record.features:
        if feature.type in list_features:
            if 'gene' in feature.qualifiers:
                gene_names_seen = writename(record.id, feature.qualifiers['gene'][0].upper(), out_gene_name, gene_names, gene_names_seen)
            elif 'product' in feature.qualifiers:
                gene_names_seen = writename(record.id, feature.qualifiers['product'][0].upper(), out_gene_name, gene_names, gene_names_seen)


# Close input and output files
# ============================================================

in_gb.close()

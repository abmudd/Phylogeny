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

parser = argparse.ArgumentParser(description='This script reads through the gb file with all individuals and replaces the input tree name with the species name.')
parser.add_argument("-l", "--log", metavar='STR', help="output file name for log of species extracted [stderr]", type=str, default='stderr')
parser.add_argument("-n", "--names", action='store_true', default=False, help="write out final names to file names.replaced.txt")
parser.add_argument("-o", "--output", metavar='STR', help="output file name for new tree [stdout]", type=str, default='stdout')
required = parser.add_argument_group('required arguments')
required.add_argument("input_tree", help="input tree file", type=str)
required.add_argument("input_gb", help="input gb file", type=str)
args = parser.parse_args()


# Set input and output files
# ============================================================

in_tree = open(args.input_tree, 'r')
in_gb = open(args.input_gb, 'r')
out_file = sys.stdout
if args.output != 'stdout':
    out_file = open(args.output, 'w')
out_log = sys.stderr
if args.output != 'stderr':
    out_log = open(args.log, 'w')


# Set variables
# ============================================================

id_dict = {}
id_seen = []
count = 0


###############################################################################
# Run
###############################################################################

# Parse input gb file
# ============================================================

for record in SeqIO.parse(in_gb, "gb"):

    # Skip entries with duplicate accession numbers.
    # ============================================================

    if record.id in id_seen:
        out_log.write("Skipping duplicate entry for accession "+record.id+".\n")

    # Skip entries without a species name.
    # ============================================================

    elif 'organism' not in record.annotations:
        out_log.write("Skipping entry for accession "+record.id+" as species name is missing.\n")

    # Skip entries with cf., sp., or aff. in the species name.
    # ============================================================

    elif 'cf.' in record.annotations['organism'] or 'sp.' in record.annotations['organism'] or 'aff.' in record.annotations['organism']:
        out_log.write("Skipping entry for species "+record.annotations['organism']+" with accession "+record.id+".\n")

    # Skip entries that only contain N's.
    # ============================================================

    elif record.seq.count('N') == len(record.seq):
        out_log.write("Skipping entry for accession "+record.id+" as the sequence only contains N's.\n")

    else:
        id_seen.append(record.id)
        id_dict[record.id] = '_'.join(record.annotations['organism'].split(' ')[0:2])
        count += 1

out_log.write("Identified "+str(count)+" possible sequence names to replace.\n")
count = 0


# Parse input tree file
# ============================================================

for item in in_tree:

    # Replace the record id with the species name using the dictionary
    # ============================================================

    for i in range(0,len(id_seen)):
        if id_seen[i] in item:
            item = item.replace(ids_seen[i], id_dict[ids_seen[i]])
            if args.names:
                names_out = open("names.replaced.txt", 'a')
                names_out.write(id_dict[ids_seen[i]]+'\n')
                names_out.close()
            count += 1

    # Write out the new tree
    # ============================================================

    out_file.write(item)

out_log.write("Replaced "+str(count)+" sequence names.\n")


# Close input and output files
# ============================================================

in_tree.close()
in_gb.close()
if args.log != 'stderr':
    out_log.close()
if args.output != 'stdout':
    out_file.close()

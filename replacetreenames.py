#!/usr/bin/env python

import argparse, os, sys
from Bio import SeqIO

# Parse arguments
parser = argparse.ArgumentParser(description='This script reads through the gb file with all individuals and replaces the input tree name with the species name.')
parser.add_argument("--output", metavar='STR', help="output file name for new tree [stdout]", type=str, default='stdout')
parser.add_argument("--log", metavar='STR', help="output file name for log of species extracted [stderr]", type=str, default='stderr')
parser.add_argument("--names", action='store_true', default=False, help="write out final names to file names.replaced.txt")
required = parser.add_argument_group('required arguments')
required.add_argument("input_tree", help="input tree file", type=str)
required.add_argument("input_gb", help="input gb file", type=str)
args = parser.parse_args()

# Open input and output files
in_tree = open(args.input_tree, 'r')
in_gb = open(args.input_gb, 'r')
out_file = sys.stdout
if args.output != 'stdout':
    out_file = open(args.output, 'w')
out_log = sys.stderr
if args.output != 'stderr':
    out_log = open(args.log, 'w')

# Set variables
id_dict = {}
id_seen = []
count = 0

# Parse input gb file
for record in SeqIO.parse(in_gb, "gb"):

    # Skip entries with duplicate accession numbers.
    if record.id in id_seen:
        out_log.write("Skipping duplicate entry for accession "+record.id+".\n")

    # Skip entries without a species name.
    elif 'organism' not in record.annotations:
        out_log.write("Skipping entry for accession "+record.id+" as species name is missing.\n")

    # Skip entries with cf., sp., or aff. in the species name.
    elif 'cf.' in record.annotations['organism'] or 'sp.' in record.annotations['organism'] or 'aff.' in record.annotations['organism']:
        out_log.write("Skipping entry for species "+record.annotations['organism']+" with accession "+record.id+".\n")

    else:
        id_seen.append(record.id)
        id_dict[record.id] = '_'.join(record.annotations['organism'].split(' ')[0:2])
        count += 1

out_log.write("Identified "+str(count)+" possible sequence names to replace.\n")
count = 0

# Parse input tree file
for item in in_tree:

    # Replace the record id with the species name using the dictionary
    for i in range(0,len(id_seen)):
        if id_seen[i] in item:
            item = item.replace(ids_seen[i], id_dict[ids_seen[i]])
            if args.names:
                names_out = open("names.replaced.txt", 'a')
                names_out.write(id_dict[ids_seen[i]]+'\n')
                names_out.close()
            count += 1

    # Write out the new tree
    out_file.write(item)

out_log.write("Replaced "+str(count)+" sequence names.\n")

#Close input and output files
in_tree.close()
in_gb.close()
if args.log != 'stderr':
    out_log.close()
if args.output != 'stdout':
    out_file.close()

#!/usr/bin/env python

import argparse, os, sys
from Bio import SeqIO

# Parse arguments
parser = argparse.ArgumentParser(description='This script removes individuals without any data and duplicate individuals.')
parser.add_argument("-l", "--log", metavar='STR', help="output file name for log of species extracted [stderr]", type=str, default='stderr')
parser.add_argument("-o", "--output", metavar='STR', help="output prefix [stdout]", type=str, default='stdout')
parser.add_argument("-s", "--silent", help="run script in silent mode", action='store_true')
required = parser.add_argument_group('required arguments')
required.add_argument("input", help="input phylip file", type=str)
args = parser.parse_args()

# Open input file and set output names
in_phy = open(args.input, 'r')
out_phy = sys.stdout
if args.output != 'stdout':
    out_phy = open(args.output, 'w')
out_log = sys.stderr
if args.log != 'stderr':
    out_log = open(args.log, 'w')

# Set variables
keep_indv = []
seen_seq = []
indv_seq = {}
seq_indv = {}

# Parse input phy file
for record in SeqIO.parse(in_phy, "phylip"):

    # Skip records with only missing data
    len_record_seq = len(record.seq)
    if record.seq.count('-') == len_record_seq:
        if not args.silent:
            out_log.write("Skipping record "+record.id+" due to it only containing -'s.\n")

    # Skip duplicate records
    elif record.seq in seen_seq:
        rm_old_name = seq_indv[record.seq]
        if rm_old_name in keep_indv:
            keep_indv.remove(rm_old_name)
        if not args.silent:
            out_log.write("Skipping duplicate records "+record.id+" and "+rm_old_name+".\n")

    else:
        keep_indv.append(record.id)
        seen_seq.append(record.seq)
        indv_seq[record.id] = record.seq
        seq_indv[record.seq] = record.id

# Write out keep_indv
out_phy.write(' '+str(len(keep_indv))+' '+str(len_record_seq)+'\n')
for item in keep_indv:
    rec_item = item
    for i in range(len(item),13):
        item = item+' '
    out_phy.write(str(item)+str(indv_seq[rec_item])+'\n')

# Close input file
in_phy.close()
if args.output != 'stdout':
    out_phy.close()
if args.log != 'stderr':
    out_log.close()

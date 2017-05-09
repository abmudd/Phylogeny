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

parser = argparse.ArgumentParser(description='This script removes individuals without any data and duplicate individuals.')
parser.add_argument("-l", "--log", metavar='STR', help="output file name for log of species extracted [stderr]", type=str, default='stderr')
parser.add_argument("-o", "--output", metavar='STR', help="output prefix [stdout]", type=str, default='stdout')
parser.add_argument("-s", "--silent", help="run script in silent mode", action='store_true')
parser.add_argument("-v", "--version", help="show version info and exit", action='version', version='%(prog)s 0.1')
required = parser.add_argument_group('required arguments')
required.add_argument("input", help="input phylip file", type=str)
args = parser.parse_args()


# Set input and output files
# ============================================================

in_phy = open(args.input, 'r')
out_phy = sys.stdout
if args.output != 'stdout':
    out_phy = open(args.output, 'w')
out_log = sys.stderr
if args.log != 'stderr':
    out_log = open(args.log, 'w')


# Set variables
# ============================================================

keep_indv = []
seen_seq = []
indv_seq = {}
seq_indv = {}


###############################################################################
# Run
###############################################################################

# Parse input phy file
# ============================================================

for record in SeqIO.parse(in_phy, "phylip"):

    # Skip records with only missing data
    # ============================================================

    len_record_seq = len(record.seq)
    if record.seq.count('-') == len_record_seq:
        if not args.silent:
            out_log.write("Skipping record "+record.id+" due to it only containing -'s.\n")

    # Skip duplicate records
    # ============================================================

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
# ============================================================

out_phy.write(' '+str(len(keep_indv))+' '+str(len_record_seq)+'\n')
for item in keep_indv:
    rec_item = item
    for i in range(len(item),13):
        item = item+' '
    out_phy.write(str(item)+str(indv_seq[rec_item])+'\n')


# Close input and output files
# ============================================================

in_phy.close()
if args.output != 'stdout':
    out_phy.close()
if args.log != 'stderr':
    out_log.close()

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


# Parse arguments
# ============================================================

parser = argparse.ArgumentParser(description='This script creates shell scripts for downloading and analyzing a GenBank record, tblastx confirmation, MAFFT or MUSCLE alignment, gBlocks filtering, and RAxML tree creation and bootstrap.')
parser.add_argument("-e", "--edirect-dir", metavar='STR', help="path for edirect directory [edirect]", type=str, default='edirect')
parser.add_argument("-g", "--gblocks", metavar='STR', help="path for gBlocks [Gblocks]", type=str, default='Gblocks')
parser.add_argument("-ma", "--mafft", metavar='STR', help="aligns sequences with MAFFT and designates path for MAFFT", type=str)
parser.add_argument("-mu", "--muscle", metavar='STR', help="aligns sequences with MUSCLE and designates path for MUSCLE", type=str)
parser.add_argument("-o", "--output", metavar='STR', help="output script name (excluding .sh) [gene]", type=str, default='gene')
parser.add_argument("-p", "--prep-dir", metavar='STR', help="directory created by prep.sh [prep/]", type=str, default='prep/')
parser.add_argument("-r", "--raxml", metavar='STR', help="path for RAxML [raxmlHPC-HYBRID-SSE3]", type=str, default='raxmlHPC-HYBRID-SSE3')
parser.add_argument("-s", "--scripts-dir", metavar='STR', help="path for directory with GitHub Phylogeny scripts [Phylogeny/]", type=str, default='Phylogeny/')
parser.add_argument("-t", "--tblastx", metavar='STR', help="path for tblastx [tblastx]", type=str, default='tblastx')
parser.add_argument("-v", "--version", help="show version info and exit", action='version', version='%(prog)s 0.1')
parser.add_argument("-w", "--workdir", metavar='STR', help="path for working directory [./]", type=str, default='./')
required = parser.add_argument_group('required arguments')
required.add_argument("gene-name", help="gene name", type=str)
required.add_argument("txid", help="taxonomy ID from NCBI for classification", type=str)
args = parser.parse_args()


# Error function
def error(output, number):
    sys.stderr.write('['+os.path.basename(sys.argv[0])+']: '+output+'\n')
    sys.exit(number)


# Check executable function
def check_exe(function):
    try:
        subprocess.call(['which', function])
    except:
        error(function+' not in PATH env / specified location or not executable', 127)

# Add slash after directories function
def add_slash(directory):
    if not directory.endswith('/'):
        directory = directory+'/'
    return directory

args.edirect-dir = add_slash(args.edirect-dir)
args.scripts-dir = add_slash(args.scripts-dir)
args.workdir = add_slash(args.workdir)


# Check that external tools are accessible
check_exe(args.scripts-dir+'/esearch')
check_exe(args.scripts-dir+'/efetch')
if args.mafft and not args.muscle:
    check_exe(args.mafft)
elif args.muscle and not args.mafft:
    check_exe(args.muscle)
elif args.muscle and args.mafft:
    error('Please select and identify path for either MAFFT or MUSCLE, not both', 127)
else:
    error('Must select and identify path for either MAFFT or MUSCLE', 127)
check_exe(args.raxml)
check_exe(args.tblastx)
check_exe(args.scripts-dir+'genenamesfromgb.py')

# Make gene directory
subprocess.call(['mkdir', '-p', args.workdir+"/"+args.gene-name+"/"])

# Open output files
out_1_sh = open(args.workdir+"/"+args.gene-name+"/"+args.gene-name+'.part1.sh', 'w')
out_2_sh = open(args.workdir+"/"+args.gene-name+"/"+args.gene-name+'.part2.sh', 'w')
out_3_sh = open(args.workdir+"/"+args.gene-name+"/"+args.gene-name+'.part3.sh', 'w')
out_4_sh = open(args.workdir+"/"+args.gene-name+"/"+args.gene-name+'.part4.sh', 'w')
out_key = open(args.workdir+"/"+args.gene-name+"/"+args.gene-name+'.key', 'w')

# Write output 1: download GenBank records; check key file
out_1_sh.write('#!/bin/bash\n\n'
               +'# Download 'args.gene-name' GenBank records from NCBI\n'
               +args.scripts-dir+'/esearch -db nucleotide -query \'txid'+args.txid+'[Organism] biomol_'
               +'genomic[PROP]'+args.gene-name+'[All Fields]\' | 'args.scripts-dir+'efetch -format gb > '
               +'NCBI_query.gb\n\n'
               +'# Create list of genes not in key file\n'
               +args.scripts-dir+'genenamesfromgb.py '+args.gene-name+'.key NCBI_query.gb\n'
               +'cut -f2 NCBI_query.gb.gene.names | sort | uniq -c >NCBI_query.gb.uniq.names\n')

# Write output 2: recheck key file
out_2_sh.write('#!/bin/bash\n\n'
               +'# Check list of genes not in key file\n'
               +args.scripts-dir+'genenamesfromgb.py '+args.gene-name+'.key NCBI_query.gb\n'
               +'cut -f2 NCBI_query.gb.gene.names | sort | uniq -c >NCBI_query.gb.uniq.names\n')

# Write output 3: extract GenBank records and make alignment
out_3_sh.write('#!/bin/bash\n\n'
               +'# Extract gene annotations from GenBank records\n'
               +'mkdir -p extract/\n'
               +args.scripts-dir+'extractgb.py --silent --log log.extractgb --output extract/sequence '
               +args.gene-name+'.key '+args.prep-dir+'NCBI_full.gb\n\n'
               +'# Merge sequences for all species\n'
               +'cat extract/sequence*.'+args.gene-name+'.fa > extract/'+args.gene-name+'.cat.fa\n\n'
               +'# Run mafft and muscle on merged gene file\n'
               +'module load meraculous\n'
               +'module load muscle\n'
               +'memtime -f memtime.mafft '+args.mafft+' --maxiterate 1000000 --genafpair --thread 20 '
               +'extract/'+args.gene-name+'.cat.fa >../mafft/'+args.gene-name+'.cat.mafft.fa 2>log.mafft &\n'
               +args.muscle+' -in extract/'+args.gene-name+'.cat.fa -out ../muscle/'+args.gene-name
               +'.cat.muscle.fa 2>log.muscle &\nwait\n')

# Make mafft and muscle directories and write analysis files
if not os.path.isdir(args.workdir+'/mafft/') and not os.path.isdir(args.workdir+'/muscle/'):
    def analysis_sh(name):
        subprocess.call(['mkdir', '-p', name])
        out_5_sh = open(args.workdir+'/'+name+'/analysis.sh', 'w')
        out_5_sh.write('#!/bin/bash\n\n'
                       +'cat '+args.workdir+'/*/log.extractgb | cut -f2 | sort | uniq | awk \'{print $0 '
                       +'"\\t" "taxon" NR}\' >translate.dict\n'
                       +args.scripts-dir+'BeforePhylo.pl -type=dna -conc=raxml -sort -Gblocks='+args.gblocks
                       +' -trim -translate=translate.dict -output=phylip *.fa >log.BeforePhylo\n'
                       +args.scripts-dir+'phyfilter.py -l output.filter.log -o output.filter.phy output.phy'
                       +'\n'
                       +'awk \'{print gsub(/-/,"")}\' output.filter.phy | tail -n +2 | hist - 1 100 > '
                       +'output.filter.missing.hist\n'
                       +'module load RAxML/8.2.2\n'
                       +'module load meraculous\n'
                       +'memtime -f memtime.RAxML '+args.raxml+' -f a -T 15 -m GTRGAMMA -n RAxML -# '
                       +'autoMRE_IGN -x $RANDOM -p $RANDOM -q output_partitions.txt -s output.filter.phy '
                       +'-o taxon131\n')
        out_5_sh.close()
        subprocess.call(['chmod', 'u+x', './'+name+'/analysis.sh'])
    analysis_sh('mafft')
    analysis_sh('muscle')

# Write output key file
out_key.write('GB_name\n'+args.gene-name+'\n')

# Close output files
out_1_sh.close()
out_2_sh.close()
out_3_sh.close()
out_4_sh.close()
out_key.close()

# Chmod u+x for sh files
subprocess.call(['chmod', 'u+x', args.workdir+"/"+args.gene-name+"/"+args.gene-name+'.part1.sh'])
subprocess.call(['chmod', 'u+x', args.workdir+"/"+args.gene-name+"/"+args.gene-name+'.part2.sh'])
subprocess.call(['chmod', 'u+x', args.workdir+"/"+args.gene-name+"/"+args.gene-name+'.part3.sh'])
subprocess.call(['chmod', 'u+x', args.workdir+"/"+args.gene-name+"/"+args.gene-name+'.part4.sh'])

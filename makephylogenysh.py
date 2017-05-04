#!/usr/bin/env python

import argparse, os, subprocess, sys

# Parse arguments
parser = argparse.ArgumentParser(description='This script creates a shell script for running a GenBank record extraction, MAFFT and MUSCLE alignment, and RAxML tree creation and bootstrap.')
parser.add_argument("-d", "--output_dir", metavar='STR', help="directory for script output [./gene/]", type=str, default='./gene/')
parser.add_argument("-g", "--gblocks", metavar='STR', help="path for gBlocks [/usr/common/jgi/phylogenetics/Gblocks/0.91b/Gblocks]", type=str, default='/usr/common/jgi/phylogenetics/Gblocks/0.91b/Gblocks')
parser.add_argument("-ma", "--mafft", metavar='STR', help="path for MAFFT [~/tools/mafft-7.307-without-extensions/bin/mafft]", type=str, default='~/tools/mafft-7.307-without-extensions/bin/mafft')
parser.add_argument("-mu", "--muscle", metavar='STR', help="path for MUSCLE [muscle]", type=str, default='muscle')
parser.add_argument("-o", "--output", metavar='STR', help="output script name (excluding .sh) [gene]", type=str, default='gene')
parser.add_argument("-p", "--prep_dir", metavar='STR', help="directory for prep file [./prep/]", type=str, default='../prep/')
parser.add_argument("-r", "--raxml", metavar='STR', help="path for RAxML [/usr/common/jgi/phylogenetics/RAxML/8.2.2/bin/raxmlHPC-HYBRID-SSE3]", type=str, default='/usr/common/jgi/phylogenetics/RAxML/8.2.2/bin/raxmlHPC-HYBRID-SSE3')
parser.add_argument("-s", "--scripts_dir", metavar='STR', help="directory for phylogeny scripts [./scripts/]", type=str, default='../scripts/')
required = parser.add_argument_group('required arguments')
required.add_argument("NCBI_query", help="NCBI query", type=str)
args = parser.parse_args()

# Make gene directory
subprocess.call(['mkdir', '-p', args.output_dir])

# Open output files
if not args.output_dir.endswith('/'):
    args.output_dir = args.output_dir+'/'
out_1_sh = open(args.output_dir+args.output+'.part1.sh', 'w')
out_2_sh = open(args.output_dir+args.output+'.part2.sh', 'w')
out_3_sh = open(args.output_dir+args.output+'.part3.sh', 'w')
out_key = open(args.output_dir+args.output+'.key', 'w')

# Write output 1: download GenBank records; check key file
out_1_sh.write('#!/bin/bash\n\n'
               +'# Download GenBank records from NCBI\n'
               +args.scripts_dir+'esearch -db nucleotide -query \''+args.NCBI_query+'\' | '
               +args.scripts_dir+'efetch -format gb > NCBI_query.gb\n\n'
               +'# Create list of genes not in key file\n'
               +args.scripts_dir+'genenamesfromgb.py '+args.output+'.key NCBI_query.gb\n'
               +'cut -f2 NCBI_query.gb.gene.names | sort | uniq -c >NCBI_query.gb.uniq.names\n')

# Write output 2: recheck key file
out_2_sh.write('#!/bin/bash\n\n'
               +'# Check list of genes not in key file\n'
               +args.scripts_dir+'genenamesfromgb.py '+args.output+'.key NCBI_query.gb\n'
               +'cut -f2 NCBI_query.gb.gene.names | sort | uniq -c >NCBI_query.gb.uniq.names\n')

# Write output 3: extract GenBank records and make alignment
out_3_sh.write('#!/bin/bash\n\n'
               +'# Extract gene annotations from GenBank records\n'
               +'mkdir -p extract/\n'
               +args.scripts_dir+'extractgb.py --silent --log log.extractgb --output extract/sequence '
               +args.output+'.key '+args.prep_dir+'NCBI_full.gb\n\n'
               +'# Merge sequences for all species\n'
               +'cat extract/sequence*.'+args.output+'.fa > extract/'+args.output+'.cat.fa\n\n'
               +'# Run mafft and muscle on merged gene file\n'
               +'module load meraculous\n'
               +'module load muscle\n'
               +'memtime -f memtime.mafft '+args.mafft+' --maxiterate 1000000 --genafpair --thread 20 '
               +'extract/'+args.output+'.cat.fa >../mafft/'+args.output+'.cat.mafft.fa 2>log.mafft &\n'
               +args.muscle+' -in extract/'+args.output+'.cat.fa -out ../muscle/'+args.output
               +'.cat.muscle.fa 2>log.muscle &\nwait\n')

# Make mafft and muscle directories and write analysis files
if not os.path.isdir('./mafft/') and not os.path.isdir('./muscle/'):
    def analysis_sh(name):
        subprocess.call(['mkdir', '-p', name])
        out_4_sh = open('./'+name+'/analysis.sh', 'w')
        out_4_sh.write('#!/bin/bash\n\n'
                       +'cat ../*/log.extractgb | cut -f2 | sort | uniq | awk \'{print $0 "\\t" "taxon" '
                       +'NR}\' >translate.dict\n'
                       +args.scripts_dir+'BeforePhylo.pl -type=dna -conc=raxml -sort -Gblocks='+args.gblocks
                       +' -trim -translate=translate.dict -output=phylip *.fa >log.BeforePhylo\n'
                       +args.scripts_dir+'phyfilter.py -l output.filter.log -o output.filter.phy output.phy'
                       +'\n'
                       +'awk \'{print gsub(/-/,"")}\' output.filter.phy | tail -n +2 | hist - 1 100 > '
                       +'output.filter.missing.hist\n'
                       +'module load RAxML/8.2.2\n'
                       +'module load meraculous\n'
                       +'memtime -f memtime.RAxML '+args.raxml+' -f a -T 15 -m GTRGAMMA -n RAxML -# '
                       +'autoMRE_IGN -x $RANDOM -p $RANDOM -q output_partitions.txt -s output.filter.phy '
                       +'-o taxon131\n')
        out_4_sh.close()
        subprocess.call(['chmod', 'u+x', './'+name+'/analysis.sh'])
    analysis_sh('mafft')
    analysis_sh('muscle')

# Write output key file
out_key.write('GB_name;Consensus_name\n'+args.output+';'+args.output+'\n')

# Close output files
out_1_sh.close()
out_2_sh.close()
out_3_sh.close()
out_key.close()

# Chmod u+x for sh files
subprocess.call(['chmod', 'u+x', args.output_dir+args.output+'.part1.sh'])
subprocess.call(['chmod', 'u+x', args.output_dir+args.output+'.part2.sh'])
subprocess.call(['chmod', 'u+x', args.output_dir+args.output+'.part3.sh'])

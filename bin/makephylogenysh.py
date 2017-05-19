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

parser = argparse.ArgumentParser(description='This script creates shell scripts for downloading and analyzing a GenBank record, tblastx confirmation, MAFFT alignment, gBlocks filtering, and RAxML tree creation and bootstrap.', epilog='It is highly recommended that absolute paths be passed via the path flags.')
parser.add_argument("-e", "--esearch", metavar='STR', help="path for esearch [esearch]", type=str, default='esearch')
parser.add_argument("-g", "--gblocks", metavar='STR', help="path for Gblocks [Gblocks]", type=str, default='Gblocks')
parser.add_argument("-m", "--mafft", metavar='STR', help="path for MAFFT [mafft]", type=str, default='mafft')
parser.add_argument("-r", "--raxml", metavar='STR', help="path for RAxML [raxmlHPC-HYBRID-SSE3]", type=str, default='raxmlHPC-HYBRID-SSE3')
parser.add_argument("-s", "--samtools", metavar='STR', help="path for samtools [samtools]", type=str, default='samtools')
parser.add_argument("-t", "--tblastx", metavar='STR', help="path for tblastx [tblastx]", type=str, default='tblastx')
parser.add_argument("-u", "--threads", metavar='INT', help="threads for MAFFT and RAxML [10]", type=str, default='10')
parser.add_argument("-v", "--version", help="show version info and exit", action='version', version='%(prog)s 0.1')
parser.add_argument("-w", "--workdir", metavar='STR', help="path for working directory [pwd]", type=str, default=os.getcwd())
required = parser.add_argument_group('required arguments')
required.add_argument("genename", help="gene name", type=str)
required.add_argument("txid", help="taxonomy ID from NCBI for classification", type=str)
args = parser.parse_args()


# Error function
# ============================================================

def error(output, number):
    sys.stderr.write('['+os.path.basename(sys.argv[0])+']: '+output+'\n')
    sys.exit(number)


# Function to check directories
# ============================================================

def check_dir(directory):
    try:
        os.path.isdir(directory)
    except:
        error(directory+' not found', 1)


# Function to check executables
# ============================================================

if sys.version_info[:2] < (3,3):
    def which(pgm):
        path=os.getenv('PATH')
        for p in path.split(os.path.pathsep):
            p=os.path.join(p,pgm)
            if os.path.exists(p) and os.access(p,os.X_OK):
                return p
else:
    from shutil import which
def check_exe_return(function):
    if not which(function) or not os.access(which(function), os.X_OK):
        error(function+' not found in specified location or not executable', 127)
    else:
        return which(function)
def check_exe(function):
    if not which(function) or not os.access(which(function), os.X_OK):
        error(function+' not found in specified location or not executable', 127)


# Set variables
# ============================================================

scriptsdir = os.path.dirname(sys.argv[0])
efetch = os.path.dirname('esearch')+'efetch'
makeblastdb = os.path.dirname('tblastx')+'makeblastdb'
prepdir = args.workdir+'/prep/'
genedir = args.workdir+'/'+args.genename+'/'
mafftdir = args.workdir+'/mafft/'


# Check that directories exist
# ============================================================

check_dir(prepdir)
check_dir(scriptsdir)
check_dir(args.workdir)


# Check that external tools are executable
# ============================================================

esearch = check_exe_return(args.esearch)
efetch = check_exe_return(efetch)
gblocks = check_exe_return(args.gblocks)
mafft = check_exe_return(args.mafft)
raxml = check_exe_return(args.raxml)
check_exe(scriptsdir+'/allgenesingb.py')
check_exe(scriptsdir+'/AfterPhylo.pl')
check_exe(scriptsdir+'/BeforePhylo.pl')
check_exe(scriptsdir+'/genenamesfromgb.py')
check_exe(scriptsdir+'/extractgb.py')
check_exe(scriptsdir+'/phyfilter.py')
makeblastdb = check_exe_return(makeblastdb)
tblastx = check_exe_return(args.tblastx)
samtools = check_exe_return(args.samtools)


# Make gene directory
# ============================================================

subprocess.call(['mkdir', '-p', genedir])


# Set output names
# ============================================================

out_1_sh = open(genedir+args.genename+'.part1.sh', 'w')
out_2_sh = open(genedir+args.genename+'.part2.sh', 'w')
out_3_sh = open(genedir+args.genename+'.part3.sh', 'w')
out_4_sh = open(genedir+args.genename+'.part4.sh', 'w')
out_key = open(genedir+args.genename+'.key', 'w')


###############################################################################
# Run
###############################################################################

# Write output 1: download GenBank records; check key file
# ============================================================

out_1_sh.write('#!/bin/bash\n\n'
               +'# Download '+args.genename+' GenBank records for taxID '+args.txid+' from NCBI\n'
               +esearch+' -db nucleotide -query \'txid'+args.txid+'[Organism] biomol_genomic[PROP] '
               +args.genename+'[All Fields]\' | '+efetch+' -format gb > '+genedir+'NCBI_query.gb\n\n'
               +'# Create list of genes not in key file\n'
               +scriptsdir+'/genenamesfromgb.py '+genedir+args.genename+'.key '+genedir+'NCBI_query.gb\n'
               +'cut -f2 '+genedir+'NCBI_query.gb.gene.names | sort | uniq >'+genedir
               +'NCBI_query.gb.uniq.names\n')


# Write output 2: recheck key file
# ============================================================

out_2_sh.write('#!/bin/bash\n\n'
               +'# Recheck list of genes not in key file\n'
               +scriptsdir+'/genenamesfromgb.py '+genedir+args.genename+'.key '+genedir+'NCBI_query.gb\n'
               +'cut -f2 '+genedir+'NCBI_query.gb.gene.names | sort | uniq >'+genedir
               +'NCBI_query.gb.uniq.names\n')


# Write output 3: extract GenBank records and confirm with BLAST
# ============================================================

out_3_sh.write('#!/bin/bash\n\n'
               +'# Extract gene annotations from GenBank records\n'
               +'mkdir -p '+genedir+'extract/\n'
               +scriptsdir+'/extractgb.py --silent --log '+genedir+'log.extractgb --output '+genedir
               +'extract/sequence '+genedir+args.genename+'.key '+prepdir+'/NCBI_full.gb\n\n'
               +'# Merge sequences for all species\n'
               +'cat '+genedir+'extract/sequence*.'+args.genename+'.fa > '+genedir+'extract/'+args.genename
               +'.cat.fa\n\n'
               +'# Identify possible incorrect sequences with blast\n'
               +makeblastdb+' -dbtype nucl -in '+genedir+'extract/'+args.genename+'.cat.fa &>'+genedir
               +'blast.log\n'
               +tblastx+' -query extract/'+args.genename+'.cat.fa -db extract/'+args.genename+'.cat.fa '
               +'-evalue '+'0.00001 -outfmt 6 -max_target_seqs 30 1>'+genedir+'blast.out 2>>'+genedir
               +'blast.log\n'
               +samtools+' faidx '+genedir+'extract/'+args.genename+'.cat.fa\n'
               +'cat /dev/null >'+genedir+'blast.incorrect \n'
               +'awk \'{if($1 != $2) {print $1}}\' '+genedir+'blast.out | sort | uniq >'+genedir+'temp\n'
               +'cut -f1 '+genedir+'extract/'+args.genename+'.cat.fa.fai >>'+genedir+'temp\n'
               +'sort '+genedir+'temp | uniq -u | while read z; do grep ${z} '+genedir+'/log.extractgb >>'
               +genedir+'blast.incorrect; done\n'
               +'rm '+genedir+'temp\n')


# Write output 4: make alignment
# ============================================================

out_4_sh.write('#!/bin/bash\n\n'
               +'# Run mafft on merged gene file\n'
               +mafft+' --maxiterate 1000000 --genafpair --thread '+args.threads+' '+genedir+'extract/'
               +args.genename+'.cat.fa >'+mafftdir+args.genename+'.cat.mafft.fa 2>'+genedir+'log.mafft\n')


# Make mafft directory and write analysis script
# ============================================================

if not os.path.isdir(mafftdir):
    subprocess.call(['mkdir', '-p', mafftdir])
    out_5_sh = open(mafftdir+'/analysis.sh', 'w')
    out_5_sh.write('#!/bin/bash\n\n'
                   +'# Concatenate mafft alignments, run Gblocks, and change names\n'
                   +'cat '+args.workdir+'/*/log.extractgb | grep -v \'Extracting a total of\' | cut -f2 | '
                   +'sort | uniq | awk \'{print $0 "\\t" "taxon" NR}\' >'+mafftdir+'translate.dict\n'
                   +scriptsdir+'/BeforePhylo.pl -type=dna -conc=raxml -sort -Gblocks='+gblocks
                   +' -trim -translate='+mafftdir+'translate.dict -output=phylip '+mafftdir+'*.fa >'
                   +mafftdir+'log.BeforePhylo\n\n'
                   +'# Remove individuals without any data and duplicate individuals\n'
                   +scriptsdir+'/phyfilter.py -l '+mafftdir+'output.filter.log -o '+mafftdir
                   +'output.filter.phy '+mafftdir+'output.phy\n\n'
                   +'# Run RAxML\n'
                   +raxml+' -f a -T '+args.threads+' -m GTRGAMMA -n RAxML -# autoMRE_IGN -x $RANDOM -p '
                   +'$RANDOM -q '+mafftdir+'output_partitions.txt -s '+mafftdir+'output.filter.phy\n\n'
                   +'# Revert names and analyze tree\n'
                   +'awk \'{print $2 "\\t" $1}\' '+mafftdir+'translate.dict >'+mafftdir+'replace.dict\n'
                   +scriptsdir+'/AfterPhylo.pl -format=newick -replace -annotate='+mafftdir+'replace.dict '
                   +mafftdir+'RAxML_bipartitions.RAxML >'+mafftdir+'log.AfterPhylo\n'
                   +scriptsdir+'/AfterPhylo.pl -format=newick -average '+mafftdir+'RAxML_bipartitions.RAxML '
                   +'>>'+mafftdir+'log.AfterPhylo\n')
    out_5_sh.close()
    subprocess.call(['chmod', 'u+x', mafftdir+'/analysis.sh'])


# Write output key file
# ============================================================

out_key.write('GB_name\n'+args.genename+'\n')


# Close output files
# ============================================================

out_1_sh.close()
out_2_sh.close()
out_3_sh.close()
out_4_sh.close()
out_key.close()


# Make output scripts executable and run part 1
# ============================================================

subprocess.call(['chmod', 'u+x', genedir+args.genename+'.part1.sh'])
subprocess.call(['chmod', 'u+x', genedir+args.genename+'.part2.sh'])
subprocess.call(['chmod', 'u+x', genedir+args.genename+'.part3.sh'])
subprocess.call(['chmod', 'u+x', genedir+args.genename+'.part4.sh'])
subprocess.call([genedir+args.genename+'.part1.sh'])

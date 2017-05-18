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
parser.add_argument("-a", "--afterphylo", metavar='STR', help="path for AfterPhylo.pl from https://github.com/qiyunzhu/AfterPhylo.git [AfterPhylo.pl]", type=str, default='AfterPhylo.pl')
parser.add_argument("-b", "--beforephylo", metavar='STR', help="path for BeforePhylo.pl from https://github.com/qiyunzhu/BeforePhylo.git [BeforePhylo.pl]", type=str, default='BeforePhylo.pl')
parser.add_argument("-e", "--esearch", metavar='STR', help="path for esearch [esearch]", type=str, default='esearch')
parser.add_argument("-g", "--gblocks", metavar='STR', help="path for gBlocks [Gblocks]", type=str, default='Gblocks')
parser.add_argument("-m", "--mafft", metavar='STR', help="path for MAFFT [mafft]", type=str, default='mafft')
parser.add_argument("-o", "--samtools", metavar='STR', help="path for samtools [samtools]", type=str, default='samtools')
parser.add_argument("-r", "--raxml", metavar='STR', help="path for RAxML [raxmlHPC-HYBRID-SSE3]", type=str, default='raxmlHPC-HYBRID-SSE3')
parser.add_argument("-s", "--scriptsdir", metavar='STR', help="path for directory with GitHub Phylogeny scripts [Phylogeny/]", type=str, default='Phylogeny/')
parser.add_argument("-t", "--tblastx", metavar='STR', help="path for tblastx [tblastx]", type=str, default='tblastx')
parser.add_argument("-v", "--version", help="show version info and exit", action='version', version='%(prog)s 0.1')
parser.add_argument("-w", "--workdir", metavar='STR', help="path for working directory [./]", type=str, default='./')
required = parser.add_argument_group('required arguments')
required.add_argument("genename", help="gene name", type=str)
required.add_argument("txid", help="taxonomy ID from NCBI for classification", type=str)
args = parser.parse_args()


# Error function
def error(output, number):
    sys.stderr.write('['+os.path.basename(sys.argv[0])+']: '+output+'\n')
    sys.exit(number)


# Check directories function
def check_dir(directory):
    try:
        os.path.isdir(directory)
    except:
        error(directory+' not found', 1)


# Check executable function
if sys.version_info[:2] < (3,3):
    def which(pgm):
        path=os.getenv('PATH')
        for p in path.split(os.path.pathsep):
            p=os.path.join(p,pgm)
            if os.path.exists(p) and os.access(p,os.X_OK):
                return p
else:
    from shutil import which
def check_exe(function):
    if not which(function) or not os.access(which(function), os.X_OK):
        error(function+' not found in specified location or not executable', 127)


# Define variables
efetch_loc = args.esearch.replace('esearch','efetch')
mkblastdb_loc = args.tblastx.replace('tblastx','makeblastdb')
prepdir = args.workdir+'/prep/'


# Check that directories exist
check_dir(prepdir)
check_dir(args.scriptsdir)
check_dir(args.workdir)


# Check that external tools are accessible
check_exe(args.afterphylo)
check_exe(args.beforephylo)
check_exe(args.esearch)
check_exe(efetch_loc)
check_exe(args.gblocks)
check_exe(args.mafft)
check_exe(args.raxml)
check_exe(args.scriptsdir+'/allgenesingb.py')
check_exe(args.scriptsdir+'/genenamesfromgb.py')
check_exe(args.scriptsdir+'/extractgb.py')
check_exe(mkblastdb_loc)
check_exe(args.tblastx)
check_exe(args.samtools)


# Make gene directory
subprocess.call(['mkdir', '-p', args.workdir+"/"+args.genename+"/"])


# Open output files
out_1_sh = open(args.workdir+"/"+args.genename+"/"+args.genename+'.part1.sh', 'w')
out_2_sh = open(args.workdir+"/"+args.genename+"/"+args.genename+'.part2.sh', 'w')
out_3_sh = open(args.workdir+"/"+args.genename+"/"+args.genename+'.part3.sh', 'w')
out_4_sh = open(args.workdir+"/"+args.genename+"/"+args.genename+'.part4.sh', 'w')
out_key = open(args.workdir+"/"+args.genename+"/"+args.genename+'.key', 'w')


# Write output 1: download GenBank records; check key file
out_1_sh.write('#!/bin/bash\n\n'
               +'# Download '+args.genename+' GenBank records for taxID '+args.txid+' from NCBI\n'
               +args.esearch+' -db nucleotide -query \'txid'+args.txid+'[Organism] biomol_genomic[PROP]'
               +args.genename+'[All Fields]\' | '+efetch_loc+' -format gb > '+args.workdir+"/"+args.genename+"/"+'NCBI_query.gb\n\n'
               +'# Create list of genes not in key file\n'
               +args.scriptsdir+'/genenamesfromgb.py '+args.workdir+"/"+args.genename+"/"+args.genename
               +'.key '+args.workdir+"/"+args.genename+"/"+'NCBI_query.gb\n'
               +'cut -f2 '+args.workdir+"/"+args.genename+"/"+'NCBI_query.gb.gene.names | sort | uniq >'
               +args.workdir+"/"+args.genename+"/"+'NCBI_query.gb.uniq.names\n')

# Write output 2: recheck key file
out_2_sh.write('#!/bin/bash\n\n'
               +'# Recheck list of genes not in key file\n'
               +args.scriptsdir+'/genenamesfromgb.py '+args.genename+'.key NCBI_query.gb\n'
               +'cut -f2 NCBI_query.gb.gene.names | sort | uniq >NCBI_query.gb.uniq.names\n')

# Write output 3: extract GenBank records and confirm with BLAST
out_3_sh.write('#!/bin/bash\n\n'
               +'# Extract gene annotations from GenBank records\n'
               +'mkdir -p extract/\n'
               +args.scriptsdir+'/extractgb.py --silent --log log.extractgb --output extract/sequence '
               +args.genename+'.key '+prepdir+'/NCBI_full.gb\n\n'
               +'# Merge sequences for all species\n'
               +'cat extract/sequence*.'+args.genename+'.fa > extract/'+args.genename+'.cat.fa\n\n'
               +'# Identify possible incorrect sequences with blast\n'
               +mkblastdb_loc+' -dbtype nucl -in extract/'+args.genename+'.cat.fa\n'
               +args.tblastx+' -query extract/'+args.genename+'.cat.fa -db extract/'
               +args.genename+'.cat.fa -evalue 0.00001 -outfmt 6 -max_target_seqs 30 >blast.out\n'
               +args.samtools+' faidx extract/'+args.genename+'.cat.fa\n'
               +'echo "Potential incorrect NCBI accessions:" >blast.incorrect \n'
               +'awk \'{if($1 != $2) {print $1}}\' extract/'+args.genename+'.cat.fa | sort | uniq >temp\n'
               +'cut -f1 extract/'+args.genename+'.cat.fa.fai >>temp\n'
               +'sort temp | uniq -u | while read z; do grep ${z} '+args.workdir+'/'+args.genename
               +'/log.extractgb >>blast.incorrect\n'
               +'rm temp\n')

# Write output 4: make alignment
out_4_sh.write('#!/bin/bash\n\n'
               +'# Run mafft on merged gene file\n'
               +args.mafft+' --maxiterate 1000000 --genafpair --thread 20 extract/'+args.genename
               +'.cat.fa >../mafft/'+args.genename+'.cat.mafft.fa 2>log.mafft &\n')

# Make mafft directory and write analysis files
if not os.path.isdir(args.workdir+'/mafft/'):
    def analysis_sh(name):
        subprocess.call(['mkdir', '-p', args.workdir+'/'+name])
        out_5_sh = open(args.workdir+'/'+name+'/analysis.sh', 'w')
        out_5_sh.write('#!/bin/bash\n\n'
                       +'cat '+args.workdir+'/*/log.extractgb | cut -f2 | sort | uniq | awk \'{print $0 '
                       +'"\\t" "taxon" NR}\' >translate.dict\n'
                       +args.scriptsdir+'BeforePhylo.pl -type=dna -conc=raxml -sort -Gblocks='+args.gblocks
                       +' -trim -translate=translate.dict -output=phylip *.fa >log.BeforePhylo\n'
                       +args.scriptsdir+'phyfilter.py -l output.filter.log -o output.filter.phy output.phy'
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

# Write output key file
out_key.write('GB_name\n'+args.genename+'\n')

# Close output files
out_1_sh.close()
out_2_sh.close()
out_3_sh.close()
out_4_sh.close()
out_key.close()

# Chmod u+x for sh files
subprocess.call(['chmod', 'u+x', args.workdir+"/"+args.genename+"/"+args.genename+'.part1.sh'])
subprocess.call(['chmod', 'u+x', args.workdir+"/"+args.genename+"/"+args.genename+'.part2.sh'])
subprocess.call(['chmod', 'u+x', args.workdir+"/"+args.genename+"/"+args.genename+'.part3.sh'])
subprocess.call(['chmod', 'u+x', args.workdir+"/"+args.genename+"/"+args.genename+'.part4.sh'])
subprocess.call([args.workdir+"/"+args.genename+"/"+args.genename+'.part1.sh'])

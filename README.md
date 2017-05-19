README


***** VERSION *****

0.1


***** DESCRIPTION *****

The phylogeny pipeline allows the user to rapidly mine GenBank for a particular taxa, curate and align the resulting hits, and produce a high-density phylogenetic tree. The steps are outlined in more detail in the prepared publication.


***** REQUIREMENTS *****

Linux system (tested on Debian GNU/Linux 6.0.6)
Python (tested on v2.7.4)
BLAST (tested on BLAST+ v2.3.0)
Entrez Direct (tested on v4.70)
Gblocks (tested on v0.91b)
MAFFT (tested on v7.221)
RAxML (tested on v7.5.4)
Samtools (tested on v0.1.19)


***** INSTALLATION *****

1. git clone git@github.com:abmudd/Phylogeny.git
2. cd Phylogeny
3. make


***** SYNOPSIS *****

1. cd workdir
2. Phylogeny/bin/prep.sh -txid 0000
3. Phylogeny/bin/makephylogenysh.py gene 0000
4. cd gene
5. examine NCBI_query.gb.gene.names and edit gene.key to include synonymous names on new lines
6. ./gene.part2.sh
7. examine NCBI_query.gb.gene.names and confirm synonymous names are removed
8. ./gene.part3.sh
9. examine blast.incorrect for any potential incorrect sequences and delete the fa files from extract subdirectory using log.extractgb to determine numbering
10. ./gene.part4.sh
11. repeat steps 3-10 for each desired gene
12. cd workdir/mafft
13. ./analysis.sh


***** prep.sh v0.1 *****

Usage: prep.sh [--workdir STR] [--scripts-dir STR] [--txid INT] [--edirect-dir STR]
       [--help] [-h]

Prep script to download GenBanks records for taxonomic classification of phylogeny.

Required arguments:
       -txid INT           taxonomy ID from NCBI for classification

Optional arguments:
       -efetch STR         path for efetch if not in PATH [efetch]
       -esearch STR        path for esearch if not in PATH [esearch]
       -workdir STR        path for working directory [pwd]


***** makephylogenysh.py v0.1 *****

usage: makephylogenysh.py [-h] [-e STR] [-g STR] [-m STR] [-r STR] [-s STR]
                          [-t STR] [-u INT] [-v] [-w STR]
                          genename txid

This script creates shell scripts for downloading and analyzing a GenBank
record, tblastx confirmation, MAFFT alignment, gBlocks filtering, and RAxML
tree creation and bootstrap.

optional arguments:
  -h, --help            show this help message and exit
  -e STR, --esearch STR
                        path for esearch [esearch]
  -g STR, --gblocks STR
                        path for Gblocks [Gblocks]
  -m STR, --mafft STR   path for MAFFT [mafft]
  -r STR, --raxml STR   path for RAxML [raxmlHPC-HYBRID-SSE3]
  -s STR, --samtools STR
                        path for samtools [samtools]
  -t STR, --tblastx STR
                        path for tblastx [tblastx]
  -u INT, --threads INT
                        threads for MAFFT and RAxML [10]
  -v, --version         show version info and exit
  -w STR, --workdir STR
                        path for working directory [pwd]

required arguments:
  genename              gene name
  txid                  taxonomy ID from NCBI for classification

It is highly recommended that absolute paths be passed via the path flags.

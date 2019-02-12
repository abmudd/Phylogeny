# Phylogeny Pipeline

The phylogeny pipeline allows the user to rapidly mine GenBank for a particular taxa, curate and align the resulting hits, and produce a high-density phylogenetic tree. The steps are outlined in more detail in the prepared manuscript.

Version 0.1


## Prerequisites

Download and install the following tools prior to running this pipeline:

```
Linux system (tested on Debian GNU/Linux 6.0.6)
Python (tested on v2.7.4)
BLAST (tested on BLAST+ v2.3.0)
Entrez Direct (tested on v4.70)
Gblocks (tested on v0.91b)
MAFFT (tested on v7.221)
RAxML (tested on v7.5.4 and v8.2.2)
```

Required Python packages:

```
argparse, Bio, collections, os, subprocess, sys
```

## Installing

Run the following to install this pipeline:

```
git clone git@github.com:abmudd/Phylogeny.git
cd Phylogeny
make
```

## Method

The pipeline uses the following approach:

1.	The user provides a taxonomy ID from NCBI. The pipeline (```prep.sh```) downloads and analyzes all GenBank sequences for that clade, returning the total number of unique species with available data and a list of all identified gene names with the percent of species containing that name.
2.	From this list, the user selects an initial gene name that is represented by a large fraction of species and creates the script files for each gene (```makephylogenysh.py```). The pipeline (```gene/gene.part1.sh```) queries the gene name, downloads all GenBank results, and analyzes the dataset to determine potential synonymous names. This step can be easily repeated if multiple synonymous names are known for a particular gene.
3.	The user selects the correct synonymous names from the provided list and outputs them into a synonym key file. The pipeline (```gene/gene.part2.sh```) can rerun the analysis of GenBank results to confirm that synonymous names are correctly written in the key. Once confirmed, the pipeline (```gene/gene.part3.sh```) extracts sequences matching the synonymous names and uses BLAST+ to flag potentially incorrect sequences.
4.	The user validates any incorrect sequences, and the pipeline (```gene/gene.part4.sh```) aligns sequences with MAFFT.
5.	The user inspects and curates the alignment and then repeats steps 2-4 for all desired genes. After all genes are aligned, the pipeline (```mafft/analysis.sh```) concatenates the alignments, filters with Gblocks, and generates a phylogenetic tree with RAxML.

## Test Data

Once prerequisites are installed and in the PATH environment, run the following to test this pipeline using the mitochondrial genes ATP6 and ATP8 with the frog genus Bufo:

```
mkdir workdir && cd workdir
/path/to/Phylogeny/bin/prep.sh -txid 8383
/path/to/Phylogeny/bin/makephylogenysh.py ATP6 8383 && cd ATP6 && cp /path/to/Phylogeny/publication/ATP6.key . && sh ATP6.part2.sh && sh ATP6.part3.sh && sh ATP6.part4.sh && cd ..
/path/to/Phylogeny/bin/makephylogenysh.py ATP8 8383 && cd ATP8 && cp /path/to/Phylogeny/publication/ATP8.key . && sh ATP8.part2.sh && sh ATP8.part3.sh && sh ATP8.part4.sh && cd ..
cd mafft && sh analysis.sh && cd ..
```

## Copyright

Copyright (c)2017. The Regents of the University of California (Regents). All Rights Reserved. Permission to use, copy, modify, and distribute this software and its documentation for educational, research, and not-for-profit purposes, without fee and without a signed licensing agreement, is hereby granted, provided that the above copyright notice, this paragraph and the following two paragraphs appear in all copies, modifications, and distributions. Contact the Office of Technology Licensing, UC Berkeley, 2150 Shattuck Avenue, Suite 510, Berkeley, CA 94720-1620, (510) 643-7201, for commercial licensing opportunities.

IN NO EVENT SHALL REGENTS BE LIABLE TO ANY PARTY FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES, INCLUDING LOST PROFITS, ARISING OUT OF THE USE OF THIS SOFTWARE AND ITS DOCUMENTATION, EVEN IF REGENTS HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGE

REGENTS SPECIFICALLY DISCLAIMS ANY WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. THE SOFTWARE AND ACCOMPANYING DOCUMENTATION, IF ANY, PROVIDED HEREUNDER IS PROVIDED "AS IS". REGENTS HAS NO OBLIGATION TO PROVIDE MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR MODIFICATIONS.

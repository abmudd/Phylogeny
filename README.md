# Phylogeny Pipeline

The phylogeny pipeline allows the user to rapidly mine GenBank for a particular taxa, curate and align the resulting hits, and produce a high-density phylogenetic tree. The steps are outlined in more detail in the prepared publication.

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
RAxML (tested on v7.5.4)
Samtools (tested on v0.1.19)
```

## Installing

Run the following to install this pipeline:

```
git clone git@github.com:abmudd/Phylogeny.git
cd Phylogeny
make
```

## Test Data

Once prerequisites are installed and in the PATH environment, run the following to test this pipeline using the mitochondrial genes ATP6 and ATP8 with the frog genus Bufo:

```
mkdir workdir && cd workdir
~/Phylogeny/bin/prep.sh -txid 8383
~/Phylogeny/bin/makephylogenysh.py ATP6 8383 && cd ATP6 && cp ~/Phylogeny/publication/ATP6.key . && sh ATP6.part2.sh && sh ATP6.part3.sh && sh ATP6.part4.sh && cd ..
~/Phylogeny/bin/makephylogenysh.py ATP8 8383 && cd ATP8 && cp ~/Phylogeny/publication/ATP8.key . && sh ATP8.part2.sh && sh ATP8.part3.sh && sh ATP8.part4.sh && cd ..
cd mafft && sh analysis.sh && cd ..
```

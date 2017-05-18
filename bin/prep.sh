#!/bin/bash

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

# Set hard and optional variables
# ============================================================

VERSION='0.1';
COMMAND="$0 $*";
DATE=`date +%F`;
ESEARCH=`which esearch`;
EFETCH=`which efetch`;
SCRIPTSDIR=`dirname $0`
WORKDIR=`pwd`


# Error function
# ============================================================

function error () {
    printf "[%s] ERROR: " `basename $0` >&2;
    echo "$2." >&2;
    exit $1;
}


# Help usage function
# ============================================================

function usage () {
    printf "\n" >&2;
    printf "%s v%s \n" `basename $0` $VERSION >&2;
    printf "\n" >&2;
    printf "Usage: %s [--workdir STR] [--scripts-dir STR] [--txid INT] [--edirect-dir STR]\n" `basename $0` >&2;
    printf "       [--help] [-h]\n" >&2;
    printf "\n" >&2;
    printf "Prep script to download GenBanks records for taxonomic classification of phylogeny.\n" >&2;
    printf "\n" >&2;
    printf "Required arguments:\n" >&2;
    printf "       -txid INT           taxonomy ID from NCBI for classification\n" >&2;
    printf "\n" >&2;
    printf "Optional arguments:\n" >&2;
    printf "       -efetch STR         path for efetch if not in PATH [efetch]\n" >&2;
    printf "       -esearch STR        path for esearch if not in PATH [esearch]\n" >&2;
    printf "       -workdir STR        path for working directory [pwd]\n" >&2;
    printf "\n" >&2;
}


## Retrieve options on the command line and check for errors
# ============================================================

HELP_MESSAGE=;

while [[ -n $@ ]]; do
    case "$1" in
        '-workdir') shift; WORKDIR=$1;;
        '-txid') shift; TXID=$1;;
        '-efetch') shift; EFETCH=$1;;
        '-esearch') shift; ESEARCH=$1;;
        '-help') HELP_MESSAGE=1;;
        '-h') HELP_MESSAGE=1;;
        -*) usage; error 2 "Invalid option: ${1}";;
        *) break;;
    esac;
    shift;
done

if [[ -n "${HELP_MESSAGE}" ]]; then
    usage;
    exit 1;

elif [[ -z "${WORKDIR}" || ! -d "${WORKDIR}" ]]; then
    usage; error 1 "WORKDIR not defined or does not exist";

elif [[ -z "${TXID}" ]]; then
    usage; error 1 "TXID not defined";

else
    printf "[%s] Starting %s %s %s %s %s %s\n" `basename $0` `date`     >&2;
    printf "[%s] Command-line: $COMMAND\n" `basename $0` >&2;
    printf "[%s] Version: $VERSION\n" `basename $0` >&2;
    printf "[%s] PARAM: %s = %s\n" `basename $0` "WORKDIR"     $WORKDIR   >&2;
    printf "[%s] PARAM: %s = %s\n" `basename $0` "TXID"        $TXID    >&2;
fi


# Check that external tools are accessible
# ============================================================

PYTHON=`which python`;
ALLGENESINGB="${SCRIPTSDIR}/allgenesingb.py";

if [[ -z "${ESEARCH}" || ! -x "${ESEARCH}" ]]; then
    error 127 "esearch not in PATH env variable or not executable";

elif [[ -z "${EFETCH}" || ! -x "${EFETCH}" ]]; then
    error 127 "efetch not in PATH env variable or not executable";

elif [[ -z "${PYTHON}" || ! -x "${PYTHON}" ]]; then
    error 127 "python not in PATH env variable or not executable";

elif [[ -z "${ALLGENESINGB}" || ! -x "${ALLGENESINGB}" ]]; then
    error 127 "allgenesingb.py not found or not executable";
fi

    
###############################################################################
# Run
###############################################################################

# Make prep directory
# ============================================================

mkdir -p $WORKDIR/prep;


# Download GenBank records
# ============================================================

printf "[%s] Downloading GenBank records \n" `basename $0` >&2;
$ESEARCH -db nucleotide -query "txid${TXID}[Organism] biomol_genomic[PROP]" | $EFETCH -format gb > $WORKDIR/prep/NCBI_full.gb;

if [[ "$(find ${WORKDIR}/prep/* -maxdepth 0 -type f -name 'NCBI_full.gb' -size +1c | wc -l)" -eq 0 ]]; then
    error 1 "Downloading GenBank records failed; please identify error and restart";
fi


# Initial query for gene names/counts and feature counts
# ============================================================

printf "[%s] Analyzing GenBank records \n" `basename $0` >&2;
$PYTHON $ALLGENESINGB $WORKDIR/prep/NCBI_full.gb;

if [[ "$(find ${WORKDIR}/prep/* -maxdepth 0 -type f -name 'NCBI_full.gb.*' -size +1c | wc -l)" -eq 0 ]]; then
    error 1 "Analyzing GenBank records failed; please identify error and restart";
else
    printf "[%s] Finished %s %s %s %s %s %s\n" `basename $0` `date` >&2;
fi

exit 0;

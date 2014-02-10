#!/usr/bin/env bash

# Usage example (just use one SRA ID as argument):
#
#   add_sra.sh SRX110652
#

# Terminate script if a command fails.
set -e

# Load biolite resources.
sh biolite_env.sh

# Get SRA ID from argument.
ID=$1

# Parse email from my python file email.py.
EMAIL=`grep email_bruno email.py | sed 's:^email_bruno = "\(.*\)"$:\1:g'`

echo "Processing $ID with $EMAIL..."

## Go to data directory.
cd /sysdev/s9/bruno/acoel_phylo/data/

## Download SRA file, convert to gunzip paired FASTQ, insert to catalog and
## clean downloaded files to save disk space.
sra import --clean --email $EMAIL --gzip $ID

# Needs to verify the quality offset of FASTQ.
echo "Finished importing. Check ASCII offset for FASTQ quality scores before continuing."
echo -n "Done? Proceed to assembly? [y/n]: "
read assembly

if [ $assembly == "y" ]; then
    # Run assembly script.
    sh assemble_sra.sh $ID
else
    echo "Take your time. Bye."
    exit 1
fi

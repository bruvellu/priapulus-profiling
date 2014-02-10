#!/usr/bin/env bash

# Usage example (just use one catalog ID as argument):
#
#   assemble_sra.sh SRX110652
#

# Terminate script if a command fails.
set -e

# Load biolite resources.
sh biolite_env.sh

# Get SRA ID from argument.
ID=$1

echo "Assembling transcriptome for $ID..."

# Go to scratch directory.
cd /sysdev/s9/bruno/acoel_phylo/scratch/

# Initiate assembly.
agalma transcriptome --id $ID > $ID.out 2>&1

# Load sample for phylogeny.
agalma load --id $ID --previous $ID >> $ID.out 2>&1

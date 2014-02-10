# Load biolite resources for the acoel phylogeny.

DATABASE=/sysdev/s9/bruno/acoel_phylo/biolite.sqlite
OUTDIR=/sysdev/s9/bruno/acoel_phylo/analyses
THREADS=24
MEMORY=800G

export BIOLITE_RESOURCES="database=$DATABASE,threads=$THREADS,memory=$MEMORY,outdir=$OUTDIR"

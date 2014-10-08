#!/bin/bash

# Usage: sh profiler.sh path/transcripts.r

set -e

#outdir=`dirname $1`

# Use EdgeR to normalize read counts.
R --no-save < $1
#R --no-save --args $1 $outdir < profiler.r
wait

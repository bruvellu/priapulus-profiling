#!/bin/bash

set -e

# Run R script that creates objects for analysis.
R --save < load_read_counts.r
wait

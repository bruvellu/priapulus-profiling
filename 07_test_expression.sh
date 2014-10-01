#!/bin/bash

set -e

# Test for differential gene expression.
R --save < test_dge.r

# Count differentially expressed genes per each interval.
R --save < count_dge.r

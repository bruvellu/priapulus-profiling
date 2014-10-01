#!/bin/bash

set -e

# Use EdgeR to normalize read counts.
R --save < edger_normalize.r

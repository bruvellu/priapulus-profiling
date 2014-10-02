#!/bin/bash

set -e

sh 00_prepare_stem.sh
sh 01_format_to_stem.sh
sh 02_format_from_stem.sh
sh 03_process_stem.sh

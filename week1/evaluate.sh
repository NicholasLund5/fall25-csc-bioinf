#!/bin/bash

# Enable strict error handling
set -euo pipefail

# Unzip all datasets
unzip -o week1/data/data1.zip >/dev/null
unzip -o week1/data/data2.zip >/dev/null
unzip -o week1/data/data3.zip >/dev/null
unzip -o week1/data/data4.zip >/dev/null

# Set stack size limit
ulimit -s 8192000 >/dev/null 2>&1

# Print header
echo "Dataset Language Runtime N50"
echo "------------------------------"

# Process each dataset
for dataset in data1 data2 data3 data4; do
    python3 week1/code/main.py "$dataset"
done
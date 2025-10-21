#!/bin/bash
export PATH="${HOME}/.codon/bin:$PATH"

# Temporarily disable strict error handling for debugging
# set -euo pipefail

# Detect and set CODON_PYTHON if not already set
if [ -z "${CODON_PYTHON:-}" ]; then
  CODON_PYTHON=$(
    curl -sL https://raw.githubusercontent.com/exaloop/codon/refs/heads/develop/test/python/find-python-library.py | python3
  )
  export CODON_PYTHON
  echo "Using CODON_PYTHON=$CODON_PYTHON"
  echo ""
fi

# Set locale to avoid formatting errors
export LC_ALL=C
export LANG=C

# Print header
echo "Method                  Language        Runtime   "
echo "--------------------------------------------------"

cd week4/code 
codon build -release main.codon -o main_exe
./main_exe
python3 main.py
#!/bin/bash
export PATH="${HOME}/.codon/bin:$PATH"

# Enable strict error handling
set -euo pipefail

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
echo "Language    Runtime"
echo "-------------------"

cd week3/code
python3 test_phylo.py
codon build -release test_phylo.codon -o test_phylo_codon_exe
./test_phylo_codon_exe
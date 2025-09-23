#!/bin/bash
export PATH="${HOME}/.codon/bin:$PATH"

# Enable strict error handling
set -euo pipefail

echo "Testing BioPython Motifs Port to Codon"
echo "======================================="

echo "Compiling motifs tests..."
codon build -release week2/code/test_motifs.codon -o test_motifs_exe

echo "Running motifs tests..."

# Use Codon's helper to detect the right libpython
if [ -z "${CODON_PYTHON:-}" ]; then
  CODON_PYTHON=$(
    curl -sL https://raw.githubusercontent.com/exaloop/codon/refs/heads/develop/test/python/find-python-library.py | python3
  )
  export CODON_PYTHON
  echo "Using CODON_PYTHON=$CODON_PYTHON"
fi

./test_motifs_exe

echo "All tests completed successfully!"

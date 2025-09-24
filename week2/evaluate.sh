#!/bin/bash
export PATH="${HOME}/.codon/bin:$PATH"

# Enable strict error handling
set -euo pipefail

echo "Testing BioPython Motifs Port to Codon"
echo "======================================="

echo "Compiling motifs tests..."
codon build -release week2/code/matrix_test.codon -o test_matrix_exe
codon build -release week2/code/minimal_test.codon -o test_minimal_exe
codon build -release week2/code/thresholds_test.codon -o test_thresholds_exe
codon build -release week2/code/__init___test.codon -o test_init_exe

echo "Running motifs tests..."

# Use Codon's helper to detect the right libpython
if [ -z "${CODON_PYTHON:-}" ]; then
  CODON_PYTHON=$(
    curl -sL https://raw.githubusercontent.com/exaloop/codon/refs/heads/develop/test/python/find-python-library.py | python3
  )
  export CODON_PYTHON
  echo "Using CODON_PYTHON=$CODON_PYTHON"
fi

./test_matrix_exe
./test_minimal_exe
./test_thresholds_exe
./test_init_exe

echo "All tests completed successfully!"

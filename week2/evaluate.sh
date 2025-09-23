#!/bin/bash
export PATH="${HOME}/.codon/bin:$PATH"

# Enable strict error handling
set -euo pipefail

# Print header
echo "Testing BioPython Motifs Port to Codon"
echo "======================================="

# Compile and run tests
echo "Compiling motifs tests..."
codon build -release week2/tests/test_motifs.codon -o test_motifs_exe

echo "Running motifs tests..."
export CODON_PYTHON="${CODON_PYTHON}"
./test_motifs_exe

echo "All tests completed successfully!"
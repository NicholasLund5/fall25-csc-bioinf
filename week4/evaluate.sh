#!/bin/bash
export PATH="${HOME}/.codon/bin:$PATH"

if [ -z "${CODON_PYTHON:-}" ]; then
  CODON_PYTHON=$(
    curl -sL https://raw.githubusercontent.com/exaloop/codon/refs/heads/develop/test/python/find-python-library.py | python3
  )
  export CODON_PYTHON
  echo "Using CODON_PYTHON=$CODON_PYTHON"
  echo ""
fi

export LC_ALL=C
export LANG=C

echo "Method                  Language        Runtime   "
echo "--------------------------------------------------"

cd week4/code 
python3 main.py
codon build -release main.codon -o main_exe
./main_exe
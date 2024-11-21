#!/bin/bash

# usage: ./tester.sh <python program>

# Get the full path of the program
fullpath=$(dirname "$(realpath "$1")")  # Only the directory path
testsdir="${fullpath}/tests"

# Colors for output
RED='\033[0;31m'      # Red
GREEN='\033[0;32m'    # Green
YELLOW='\033[0;33m'   # Yellow
NC='\033[0m'          # No Color

# Iterate over test (.test) and answer (.ans) files
#cd "$testsdir" || { echo "Tests directory not found: $testsdir"; exit 1; }

for t in tests/*.test; do
    a="${t%.test}.ans"  # Match corresponding .ans file
    if [[ ! -f "$a" ]]; then
        echo -e "${YELLOW}Warning: No matching answer file for test file: $t${NC}"
        continue
    fi

    echo -e "${GREEN}Running Test: $t${NC}"
    echo -e "${RED}Program output:${NC}"
    python3 "$1" < "$t"
    echo -e "${GREEN}Expected output:${NC}"
    cat "$a"
    echo  # Print an empty line for clarity
done

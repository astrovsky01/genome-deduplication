#!/usr/bin/env bash

. "$(dirname "$0")/common.sh"

# Test: repository Bash scripts parse successfully.
# Expected: bash -n exits successfully for synopsis.sh and every test script.
bash -n "$ROOTDIR/code/synopsis.sh"
for script in "$ROOTDIR"/tests/synopsis/*.sh; do
    bash -n "$script"
done

printf '%s\n' 'PASS: Bash syntax'

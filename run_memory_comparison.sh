#!/bin/bash

# Memory comparison script for dedup2.py vs dedup2.cpp
# This script runs both versions with memory tracking and compares results

echo "==================================================================="
echo "Memory Comparison: dedup2.py vs dedup2.cpp"
echo "==================================================================="
echo ""

# Set variables
INPUT_FILES="tests/test-data/GCA_*.fna.gz"
SEED=123
PY_OUT="testout/py_memory"
CPP_OUT="testout/cpp_memory"

echo ""
echo "==================================================================="
echo "Running Python version with memory tracking..."
echo "==================================================================="
yes | /usr/bin/time -v python code/dedup2.py -o "$PY_OUT" --seed "$SEED" $INPUT_FILES 2>&1 | tee testout/py_memory_report.txt

echo ""
echo "==================================================================="
echo "Running C++ version with memory tracking..."
echo "==================================================================="
yes | /usr/bin/time -v ./code/dedup2 -o "$CPP_OUT" --seed "$SEED" $INPUT_FILES 2>&1 | tee testout/cpp_memory_report.txt

echo ""
echo "==================================================================="
echo "Extracting memory statistics..."
echo "==================================================================="

# Extract key metrics from Python run
PY_MAX_RSS=$(grep "Maximum resident set size" testout/py_memory_report.txt | awk '{print $6}')
PY_TIME=$(grep "Elapsed (wall clock)" testout/py_memory_report.txt | awk '{print $8}')
PY_CPU=$(grep "Percent of CPU" testout/py_memory_report.txt | awk '{print $7}')

# Extract key metrics from C++ run
CPP_MAX_RSS=$(grep "Maximum resident set size" testout/cpp_memory_report.txt | awk '{print $6}')
CPP_TIME=$(grep "Elapsed (wall clock)" testout/cpp_memory_report.txt | awk '{print $8}')
CPP_CPU=$(grep "Percent of CPU" testout/cpp_memory_report.txt | awk '{print $7}')

echo ""
echo "==================================================================="
echo "COMPARISON SUMMARY"
echo "==================================================================="
echo ""
echo "Python Version:"
echo "  Max Memory:     ${PY_MAX_RSS} KB"
echo "  Elapsed Time:   ${PY_TIME}"
echo "  CPU Usage:      ${PY_CPU}"
echo ""
echo "C++ Version:"
echo "  Max Memory:     ${CPP_MAX_RSS} KB"
echo "  Elapsed Time:   ${CPP_TIME}"
echo "  CPU Usage:      ${CPP_CPU}"
echo ""

# Calculate memory ratio if both values exist
if [ -n "$PY_MAX_RSS" ] && [ -n "$CPP_MAX_RSS" ]; then
    RATIO=$(echo "scale=2; $PY_MAX_RSS / $CPP_MAX_RSS" | bc)
    echo "Memory Ratio (Python/C++): ${RATIO}x"
    echo ""
fi

echo "==================================================================="
echo "Verifying output correctness..."
echo "==================================================================="

# Compare outputs
echo ""
echo "Line counts:"
wc -l "$PY_OUT"/*.bed "$CPP_OUT"/*.bed

echo ""
echo "Checking for differences:"
for base in GCA_001775245.1_ASM177524v1_genomic GCA_001786795.1_ASM178679v1_genomic GCA_001820145.1_ASM182014v1_genomic; do
    echo ""
    echo "  === $base ==="
    for ext in samples masks ignored ambiguous; do
        DIFF_COUNT=$(diff "$PY_OUT"/${base}.$ext.bed "$CPP_OUT"/${base}.$ext.bed 2>/dev/null | wc -l)
        if [ "$DIFF_COUNT" -eq 0 ]; then
            echo "    $ext.bed: ✓ IDENTICAL"
        else
            echo "    $ext.bed: ✗ DIFFERENT ($DIFF_COUNT lines differ)"
        fi
    done
done

echo ""
echo "==================================================================="
echo "Full reports saved to:"
echo "  Python:  testout/py_memory_report.txt"
echo "  C++:     testout/cpp_memory_report.txt"
echo "==================================================================="

#!/bin/bash
# Dummy GNINA script for testing on Linux

echo "[DUMMY GNINA] Called with args: $*"

# Find output file argument (after --out)
out_file=""
prev=0
for arg in "$@"; do
    if [ "$prev" = "1" ]; then
        out_file="$arg"
        break
    fi
    if [ "$arg" = "--out" ]; then
        prev=1
    else
        prev=0
    fi
done

# Create dummy output file if specified
if [ -n "$out_file" ]; then
    out_dir=$(dirname "$out_file")
    if [ ! -d "$out_dir" ]; then
        mkdir -p "$out_dir"
    fi
    echo "REMARK DUMMY GNINA OUTPUT" > "$out_file"
fi

exit 0 
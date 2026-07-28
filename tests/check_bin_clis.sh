#!/usr/bin/env bash
set -euo pipefail

echo "Testing CLI --help for all scripts in bin/..."

for f in bin/*; do
    if [ -f "$f" ]; then
        echo "Testing $f --help..."
        python3 "$f" --help >/dev/null 2>&1 || "$f" --help >/dev/null 2>&1 || { echo "FAIL $f"; exit 1; }
    fi
done

echo "CLI-OK"

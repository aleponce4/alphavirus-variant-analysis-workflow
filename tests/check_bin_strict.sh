#!/usr/bin/env bash
set -euo pipefail

echo "Testing strict CLI behavior for all scripts in bin/..."

for f in bin/*; do
    if [ -f "$f" ]; then
        "$f" --help >/dev/null 2>&1 || { echo "FAIL --help $f"; exit 1; }
        "$f" --version >/dev/null 2>&1 || { echo "FAIL --version $f"; exit 1; }
        if "$f" --badflag >/dev/null 2>&1; then
            echo "NO-STRICT $f (accepted invalid flag without exit code 2)"
            exit 1
        fi
    fi
done

echo "BIN-OK"

#!/usr/bin/env bash
set -euo pipefail

echo "Running container lint checks..."

ERRORS=0

# 1. Check that no process file contains inline hardcoded container strings (must use params.container_*)
for nf_file in $(find modules subworkflows main.nf -name "*.nf" 2>/dev/null); do
    if grep -E "container\s+['\"][^$]" "$nf_file" >/dev/null 2>&1; then
        echo "ERROR: Hardcoded container directive found in $nf_file. Use params.container_*."
        ERRORS=$((ERRORS + 1))
    fi
done

# 2. Check conf/containers.config for forbidden 'latest' tag or unpinned tags or invalid registries
if [ -f "conf/containers.config" ]; then
    if grep -E ":latest|:latest'" conf/containers.config >/dev/null 2>&1; then
        echo "ERROR: 'latest' tag found in conf/containers.config."
        ERRORS=$((ERRORS + 1))
    fi

    if grep -E "container_" conf/containers.config | grep -vE "quay\.io/biocontainers/|ghcr\.io/" >/dev/null 2>&1; then
        echo "ERROR: Container registry outside quay.io/biocontainers or ghcr.io found in conf/containers.config."
        ERRORS=$((ERRORS + 1))
    fi
fi

if [ "$ERRORS" -gt 0 ]; then
    echo "CONTAINERS LINT FAILED ($ERRORS errors)"
    exit 1
fi

echo "CONTAINERS OK"

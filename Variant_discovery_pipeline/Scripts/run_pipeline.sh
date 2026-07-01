#!/bin/bash

# Backward-compatible Stage 1 entrypoint.
# Delegates to the manifest-driven Stage 1 orchestrator.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"

bash "${REPO_DIR}/run_stage1.sh" "$@"

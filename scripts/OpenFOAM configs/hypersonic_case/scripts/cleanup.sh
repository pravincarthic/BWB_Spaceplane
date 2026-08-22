#!/bin/bash
# cleanup.sh - archive old run results before a fresh run
#
# Usage: bash scripts/cleanup.sh [archive_name]

set -e

CASE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../case" && pwd)"
ARCHIVE_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)/archive"
NAME="${1:-run_$(date +%Y%m%d_%H%M%S)}"

mkdir -p "$ARCHIVE_ROOT/$NAME"

cd "$CASE_DIR"

echo "Archiving time directories (excluding 0/) to: $ARCHIVE_ROOT/$NAME"
for d in [1-9]*; do
    if [ -d "$d" ]; then
        mv "$d" "$ARCHIVE_ROOT/$NAME/"
    fi
done

if [ -d processor0 ]; then
    echo "Archiving decomposed processor directories"
    mv processor* "$ARCHIVE_ROOT/$NAME/"
fi

if [ -d VTK ]; then
    echo "Archiving VTK output"
    mv VTK "$ARCHIVE_ROOT/$NAME/"
fi

echo "Cleanup complete. Case is ready for a fresh run."

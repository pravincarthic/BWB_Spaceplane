#!/bin/bash
# ============================================================================
# cleanup.sh - archive a finished run before starting a fresh one
#
# Usage: bash scripts/cleanup.sh [archive_name]
#
# With collated output the decomposed data lives in a single processors<N>
# directory rather than 384 processor<n> directories, so both layouts are
# handled here (the plain processor* form is still produced if someone runs
# with -fileHandler uncollated).
# ============================================================================

set -euo pipefail

PKG_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
CASE_DIR="$PKG_DIR/case"
ARCHIVE_ROOT="$PKG_DIR/archive"
NAME="${1:-run_$(date +%Y%m%d_%H%M%S)}"
DEST="$ARCHIVE_ROOT/$NAME"

mkdir -p "$DEST"
cd "$CASE_DIR"

echo "Archiving to: $DEST"

echo "  reconstructed time directories"
for d in [1-9]* 0.*; do
    [ -d "$d" ] && [ "$d" != "0" ] && mv "$d" "$DEST/"
done

for p in processors[0-9]* processor[0-9]*; do
    if [ -d "$p" ]; then
        echo "  decomposed data: $p"
        mv "$p" "$DEST/"
    fi
done

if [ -d postProcessing ]; then
    echo "  postProcessing (probes, surfaces, forces)"
    mv postProcessing "$DEST/"
fi

if [ -d VTK ]; then
    echo "  VTK"
    mv VTK "$DEST/"
fi

if [ -d "$PKG_DIR/logs" ]; then
    echo "  logs"
    mkdir -p "$DEST/logs"
    find "$PKG_DIR/logs" -maxdepth 1 -type f -name '*.log' -exec mv {} "$DEST/logs/" \;
fi

echo ""
echo "Kept in place: 0/, constant/ (mesh and thermo), system/."
echo "The case is ready for a fresh scripts/run.sh."

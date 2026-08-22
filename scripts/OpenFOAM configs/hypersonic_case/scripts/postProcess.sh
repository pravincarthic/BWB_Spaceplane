#!/bin/bash
# postProcess.sh - field extraction and basic monitoring
#
# Run after scripts/run.sh has produced time-step results.

set -e

CASE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../case" && pwd)"
cd "$CASE_DIR"

echo "Available time directories:"
foamListTimes

echo ""
echo "Species mass fraction sum check (should be approx 1.0 everywhere):"
echo "Manually verify in ParaView or with a custom function object -"
echo "this script does not compute it automatically."

echo ""
echo "Converting latest fields to VTK for visualization"
foamToVTK -latestTime

echo ""
echo "Done. VTK output in $CASE_DIR/VTK/"
echo "Load into ParaView and check per docs/VALIDATION.md:"
echo "  1. Tt shock-jump structure"
echo "  2. Tv lag behind Tt (non-equilibrium)"
echo "  3. Species dissociation front (N2, O2 decreasing; N, O increasing)"
echo "  4. sum(Y_i) conservation"

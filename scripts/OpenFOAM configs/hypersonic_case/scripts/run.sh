#!/bin/bash
# ============================================================================
# run.sh - decompose and run rhoCentralFoam
#
# Usage: bash scripts/run.sh [num_ranks]        (default 384)
#
# Collated output is used throughout (fileHandler collated, set in
# case/system/controlDict). Every parallel step must agree on that setting,
# so -fileHandler collated is passed explicitly here as well.
#
# There is deliberately NO reconstructPar at the end. With collated output the
# time directories are already single files that ParaView reads directly
# through case/case.foam - reconstructing a 31M cell case 20 times would cost
# hours and produce nothing new.
# ============================================================================

set -euo pipefail

NP="${1:-384}"
PKG_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
CASE_DIR="$PKG_DIR/case"
LOG_DIR="$PKG_DIR/logs"
mkdir -p "$LOG_DIR"

cd "$CASE_DIR"

if [ ! -f constant/polyMesh/owner ]; then
    echo "ERROR: no converted mesh in constant/polyMesh/. Run scripts/setup.sh first."
    exit 1
fi

DECOMP_N="$(foamDictionary -entry numberOfSubdomains -value system/decomposeParDict)"
if [ "$DECOMP_N" != "$NP" ]; then
    echo "ERROR: requested $NP ranks but system/decomposeParDict says $DECOMP_N."
    echo "Edit numberOfSubdomains in that file, or run with: bash scripts/run.sh $DECOMP_N"
    echo "The decomposition is not regenerated automatically - ptscotch on a"
    echo "31M cell mesh is a deliberate, reproducible choice and silently"
    echo "swapping it for a geometric split would change the load balance."
    exit 1
fi

STAMP="$(date +%Y%m%d_%H%M%S)"

if [ "$NP" -gt 1 ]; then
    if [ ! -d processor0 ] && [ ! -d processors"$NP" ]; then
        echo "Decomposing for $NP ranks (ptscotch, collated)"
        decomposePar -force -fileHandler collated 2>&1 | tee "$LOG_DIR/decomposePar_$STAMP.log"
    else
        echo "Existing decomposition found, reusing it."
        echo "Delete processors$NP/ (or processor*/) to force a fresh decomposePar."
    fi

    echo ""
    echo "Running rhoCentralFoam on $NP ranks"
    echo "  675,000 fixed steps of 45 ns, endTime 0.030375 s"
    echo "  200 volume writes, 675 surface frames, probes every 10 steps"
    LOG_FILE="$LOG_DIR/rhoCentralFoam_$STAMP.log"
    srun -np "$NP" rhoCentralFoam -parallel -fileHandler collated \
        2>&1 | tee "$LOG_FILE"
else
    echo "Running rhoCentralFoam in serial (this is a 31M cell case - expect"
    echo "this to be useful only for a smoke test on a cut-down mesh)"
    LOG_FILE="$LOG_DIR/rhoCentralFoam_$STAMP.log"
    rhoCentralFoam 2>&1 | tee "$LOG_FILE"
fi

echo ""
echo "Run finished. Log: $LOG_FILE"
echo ""
echo "Check the Courant monitor before trusting the result - the timestep is"
echo "fixed, so nothing adapted if it drifted up:"
echo "    grep -A2 'fieldMinMax courantMonitor' \"$LOG_FILE\" | grep max | tail"
echo "Design value is about 0.31. Anything approaching 1 means the run is"
echo "unstable and deltaT must be reduced."
echo ""
echo "Next: bash scripts/postProcess.sh"

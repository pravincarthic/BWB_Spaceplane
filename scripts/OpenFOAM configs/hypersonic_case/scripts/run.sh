#!/bin/bash
# run.sh - MPI execution wrapper for hy2Foam
#
# Usage: bash scripts/run.sh [num_ranks]
# Default num_ranks = 8 if not given.

set -e

NP="${1:-8}"
CASE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../case" && pwd)"
LOG_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../logs" && pwd)"

# --------------------------------------------------------------------------
# Extents used to drive the hierarchical decomposition below.
#
# These are the SPACEPLANE's own bounding box, not the outer farfield
# domain box. The farfield box (200 m cube) is just the artificial
# empty far-field cage the simulation runs inside - the mesh is almost
# entirely coarse and near-uniform out there, so its shape tells you
# nothing about how the actual mesh (and therefore load) is distributed.
# The cells that matter are clustered tightly around the vehicle.
#
# Measured from Geometry_precut_AoA_neg5_20260822.step (the vehicle
# solid, kept in the file as a reference body alongside the pre-cut
# domain):
#   X:  0.057 m to  49.643 m  (Lx = 49.585 m)  - nose to tail
#   Y: -22.459 m to 14.213 m  (Ly = 36.672 m)  - spanwise
#   Z: -0.022 m to  2.600 m   (Lz =  2.622 m)  - thickness
#
# Note: this case uses a y+<1 wall-resolved mesh (prism layers stacked
# in the thickness/wall-normal direction - see SETTINGS.md). Because Lz
# is so much smaller than Lx/Ly, the balanced-edge-length search below
# naturally comes out to nz=1 for realistic rank counts - i.e. it never
# places a decomposition boundary through the thin near-wall boundary
# layer stack. That is the physically correct behavior here, not a
# degenerate result: splitting through a wall-resolved viscous layer
# adds interprocessor communication exactly where the mesh is stiffest
# and finest, so leaving it whole and instead parallelizing across the
# streamwise (X, wake/shock direction) and spanwise (Y) extents - where
# the geometry and flow actually have room to spread across ranks - is
# the standard approach.
#
# If you re-derive the geometry, recompute LX/LY/LZ with a CAD kernel
# (e.g. cadquery: shape.BoundingBox() on the vehicle solid, not the
# domain solid) and update below. Once you have the actual generated
# mesh (after scripts/setup.sh), it's worth confirming real per-rank
# cell counts with `checkMesh -allTopology` or by inspecting
# processor*/constant/polyMesh cell counts after decomposePar, since
# this bounding-box weighting is still a geometric proxy, not a measure
# of actual cell density - the real mesh may cluster cells (bow shock,
# wake refinement) in ways this can't see.
# --------------------------------------------------------------------------
LX=49.585
LY=36.672
LZ=2.622

cd "$CASE_DIR"

if [ ! -d "constant/polyMesh" ] || [ ! -f "constant/polyMesh/owner" ]; then
    echo "ERROR: no converted mesh found in constant/polyMesh/"
    echo "Run bash scripts/setup.sh first."
    exit 1
fi

if [ "$NP" -gt 1 ]; then
    if [ ! -f "system/decomposeParDict" ]; then
        echo "NOTE: no system/decomposeParDict found."
        echo "Generating a hierarchical decomposition for $NP ranks."

        # Factor NP into three integers (nx ny nz) with nx*ny*nz = NP,
        # choosing the triple that keeps each subdomain's edge length
        # (Lx/nx, Ly/ny, Lz/nz) as close to equal as possible - i.e. the
        # split is weighted by the vehicle's own proportions (LX/LY/LZ
        # above), not the outer farfield box. For this vehicle (thin,
        # elongated, wall-resolved in Z) that means nz stays at 1 for
        # realistic rank counts and splitting happens in X/Y instead -
        # see the note above on why that's the physically correct call.
        read NX NY NZ <<< "$(awk -v np="$NP" -v lx="$LX" -v ly="$LY" -v lz="$LZ" '
            BEGIN {
                best_cost = -1
                for (a = 1; a <= np; a++) {
                    if (np % a != 0) continue
                    rem = np / a
                    for (b = 1; b <= rem; b++) {
                        if (rem % b != 0) continue
                        c = rem / b
                        sx = lx / a; sy = ly / b; sz = lz / c
                        mean = (sx + sy + sz) / 3
                        cost = (sx - mean)^2 + (sy - mean)^2 + (sz - mean)^2
                        if (best_cost < 0 || cost < best_cost) {
                            best_cost = cost; bnx = a; bny = b; bnz = c
                        }
                    }
                }
                print bnx, bny, bnz
            }
        ')"
        echo "hierarchical split (from vehicle bbox ${LX}x${LY}x${LZ} m): nx=$NX ny=$NY nz=$NZ (nx*ny*nz=$NP)"

        cat > system/decomposeParDict << EOF
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      decomposeParDict;
}
numberOfSubdomains  $NP;
method              hierarchical;

hierarchicalCoeffs
{
    n           ($NX $NY $NZ);
    delta       0.001;
    order       xyz;
}
EOF
    fi
    echo "Decomposing case for $NP ranks"
    decomposePar -force

    echo "Running hy2Foam in parallel on $NP ranks"
    LOG_FILE="$LOG_DIR/hy2Foam_run_$(date +%Y%m%d_%H%M%S).log"
    mpirun -np "$NP" hy2Foam -parallel | tee "$LOG_FILE"

    echo "Reconstructing decomposed case"
    reconstructPar
else
    echo "Running hy2Foam in serial"
    LOG_FILE="$LOG_DIR/hy2Foam_run_$(date +%Y%m%d_%H%M%S).log"
    hy2Foam | tee "$LOG_FILE"
fi

echo "Run complete. Log written to: $LOG_FILE"
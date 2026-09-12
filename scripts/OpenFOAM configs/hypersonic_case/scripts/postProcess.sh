#!/bin/bash
# ============================================================================
# postProcess.sh - collect the run output for ParaView and for the PSE solver
#
# Usage: bash scripts/postProcess.sh [num_ranks]     (default 384)
# ============================================================================

set -euo pipefail

NP="${1:-384}"
PKG_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
CASE_DIR="$PKG_DIR/case"
cd "$CASE_DIR"

echo "Time directories written:"
foamListTimes -withZero || true

echo ""
echo "=========================================================="
echo "1. ParaView"
echo "=========================================================="
echo "Nothing needs converting. Open case/case.foam directly - the OpenFOAM"
echo "reader understands the collated decomposition, so pick 'Decomposed Case'"
echo "in the reader panel if it does not select it automatically."
echo ""
echo "The shock-wave sampling is already in ParaView's native format. Open"
echo "these as file series (ParaView groups the numbered time directories"
echo "automatically):"
echo "    postProcessing/shockSurfaces/*/shockSchlierenStrong.vtp   shock surface"
echo "    postProcessing/shockSurfaces/*/shockSchlierenWeak.vtp     shock, downstream"
echo "    postProcessing/shockSurfaces/*/shockIsoP_outer.vtp        stand-off"
echo "    postProcessing/shockSurfaces/*/shockIsoP_mid.vtp          shock mid-jump"
echo "    postProcessing/shockSurfaces/*/sonicSurface.vtp           Ma = 1"
echo "    postProcessing/shockSurfaces/*/planeSymmetryY0.vtp        best single view"
echo "    postProcessing/shockSurfaces/*/planeX*.vtp                cross-flow"
echo "    postProcessing/wallSurface/*/Solid_Walls.vtp              surface p and heat flux"
echo "675 frames at 45 us spacing."
echo ""
echo "If the |grad rho| isosurfaces come out empty or noisy, retune isoValue"
echo "in case/system/surfacesShock - 5.0 and 1.0 kg/m^4 are estimates from the"
echo "normal-shock density jump, not measurements."

echo ""
echo "=========================================================="
echo "2. Volume fields to VTK (optional - only if you need them outside ParaView)"
echo "=========================================================="
echo "Skipping by default: 20 writes of a 31M cell mesh is a large conversion."
echo "To convert the last one:"
echo "    foamToVTK -latestTime -fileHandler collated"
echo "To convert one field only:"
echo "    foamToVTK -latestTime -fields '(p T U Ma magGradRho)' -fileHandler collated"

echo ""
echo "=========================================================="
echo "3. PSE input"
echo "=========================================================="
if [ -d postProcessing/pseWallWindward ]; then
    NROWS=$(find postProcessing/pseWallWindward -name 'p' | head -1 | xargs wc -l 2>/dev/null | awk '{print $1}')
    echo "Wall pressure time series found, about ${NROWS:-?} rows."
else
    echo "WARNING: postProcessing/pseWallWindward not found - the probes did not run."
fi
echo ""
echo "Sampling of the wall probes:"
echo "    dt_s = 10 * 45 ns = 450 ns      fs = 2.2222 MHz"
echo "    NYQUIST           = 1.1111 MHz"
echo "    record            = 0.030375 s, df = 32.92 Hz, 67,500 samples"
echo "The second Mack mode band across this vehicle is roughly 7.9 kHz at the"
echo "trailing edge to 170 kHz near the nose, so the Nyquist margin is at"
echo "least 6x and the series can be FFT'd as-is."
echo ""
echo "Mean baseflow for PSE (time-averaged over the second flow-through):"
echo "    UMean, pMean, TMean, rhoMean in the last time directory"
echo "    UPrime2Mean and pPrime2Mean show where the disturbance is growing -"
echo "    use them to confirm the probe stations are in the right place."
echo ""
echo "To point the PSE solver at this case, set in the pse_solver config:"
echo "    baseflow.source     = \"foam\""
echo "    baseflow.foam_case  = \"$CASE_DIR\""
echo "    baseflow.foam_time  = <latest time directory>"
echo "    baseflow.wall_patches = [\"Solid_Walls\"]"
echo "    baseflow.wall_temperature = 1500.0"
echo "and widen modes.f_min_hz / f_max_hz to cover the station you are"
echo "analysing - the defaults in config/mach10_130kft.json (40-320 kHz) are"
echo "sized for a 0.2-2.0 m cone, not for the 48 m vehicle."

echo ""
echo "=========================================================="
echo "4. Sanity checks"
echo "=========================================================="
echo "    postProcessing/courantMonitor/  max Co should sit near 0.31"
echo "    postProcessing/fieldMinMax/     T max must stay below 6000 K"
echo "                                    (janaf Thigh) and p must stay positive"
echo "    postProcessing/yPlus/           y+ below 1 over Solid_Walls"
echo "    postProcessing/forceCoeffs/     Cd, Cl settled after one flow-through"
echo "    postProcessing/wallHeatFlux/    peak heat flux at the nose"

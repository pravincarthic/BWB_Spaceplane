#!/bin/bash
# ============================================================================
# setup.sh - mesh conversion and case initialisation
#
# Stock OpenFOAM v2412. Nothing from the hyStrath / hy2Foam suite is needed
# any more: the chemistry, two-temperature and Blottner-Eucken transport
# dictionaries that used to live in case/constant/ have been removed, and the
# thermophysical model is now a single janaf + polynomial air pseudo-species
# built entirely from core OpenFOAM models.
#
# Run from the top-level extracted folder:  bash scripts/setup.sh
# ============================================================================

set -euo pipefail

# --------------------------------------------------------------------------
# PLACEHOLDER: point this at your actual Gmsh mesh file
# --------------------------------------------------------------------------
MESH_FILE="${MESH_FILE:-$HOME/plain_BWB/Gmsh_without_cavities_AoA_neg5_20260822.msh}"

CASE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../case" && pwd)"

echo "=========================================="
echo "rhoCentralFoam case setup (OpenFOAM v2412)"
echo "=========================================="
echo "Case directory : $CASE_DIR"
echo "Mesh file      : $MESH_FILE"
echo ""

if ! command -v rhoCentralFoam >/dev/null 2>&1; then
    echo "ERROR: rhoCentralFoam not on PATH. Source the OpenFOAM v2412"
    echo "environment first, for example:"
    echo "  source /opt/openfoam2412/etc/bashrc"
    exit 1
fi

FOAM_VER="$(foamVersion 2>/dev/null || echo unknown)"
echo "OpenFOAM version reported: $FOAM_VER"
case "$FOAM_VER" in
    *2412*) ;;
    *) echo "WARNING: this case is written against v2412. Other versions may"
       echo "         differ in function object keywords and sampling types." ;;
esac
echo ""

if [ ! -f "$MESH_FILE" ]; then
    echo "ERROR: MESH_FILE not found at: $MESH_FILE"
    echo "Edit scripts/setup.sh, or export MESH_FILE=/path/to/mesh.msh"
    exit 1
fi

cd "$CASE_DIR"

echo "Step 1: converting mesh with gmshToFoam"
gmshToFoam "$MESH_FILE"

echo ""
echo "Step 2: checking mesh integrity"
checkMesh -allTopology -allGeometry | tee ../logs/checkMesh.log

echo ""
echo "Step 3: boundary patch names in the converted mesh"
echo "-----------------------------------------------------------"
echo "These MUST match the patch names used in case/0/* and in the"
echo "function objects: Inlet, Outlet, Atmosphere, Solid_Walls."
echo "Rename them in Gmsh or edit the dictionaries if they differ."
echo "-----------------------------------------------------------"
foamDictionary -entry boundary -expand constant/polyMesh/boundary || true

echo ""
echo "Step 4: minimum cell size check"
echo "-----------------------------------------------------------"
echo "The run uses a FIXED 45 ns timestep. For the target CFL of 0.37 the"
echo "smallest cell in the domain must be at least 0.424 mm:"
echo "    dx_min = (u_inf + a_inf) * deltaT / Co"
echo "           = 3487 m/s * 45e-9 s / 0.37"
echo "           = 4.241e-4 m"
echo "The mesh generator floor (Mesh.MeshSizeMin in scripts/gmsh_config.json)"
echo "is 0.5 mm, which gives Co = 0.314. If checkMesh above reports a smaller"
echo "minimum edge length than 0.424 mm, REDUCE deltaT in"
echo "case/system/controlDict before running, or the explicit solver will"
echo "diverge - there is no adaptive timestep to save it."
echo "-----------------------------------------------------------"

echo ""
echo "Step 5: validating the case dictionaries"
foamDictionary -expand system/controlDict            > /dev/null && echo "  controlDict            ok"
foamDictionary -expand system/fvSchemes              > /dev/null && echo "  fvSchemes              ok"
foamDictionary -expand system/fvSolution             > /dev/null && echo "  fvSolution             ok"
foamDictionary -expand system/decomposeParDict       > /dev/null && echo "  decomposeParDict       ok"
foamDictionary -expand constant/thermophysicalProperties > /dev/null && echo "  thermophysicalProperties ok"
foamDictionary -expand constant/turbulenceProperties > /dev/null && echo "  turbulenceProperties   ok"

echo ""
echo "Setup complete."
echo ""
echo "Recommended next step: regenerate the PSE probe locations from the real"
echo "wall surface, since the ones shipped in case/system/probesPSE are"
echo "bounding-box estimates:"
echo "    python3 scripts/make_pse_probes.py"
echo ""
echo "Then run:  bash scripts/run.sh 384"

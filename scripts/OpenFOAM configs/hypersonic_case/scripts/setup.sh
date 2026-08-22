#!/bin/bash
# ============================================================================
# setup.sh - mesh conversion and case initialization
#
# THIS IS WHERE THE MESH FILE PLACEHOLDER LIVES. Edit MESH_FILE below to
# point to your actual .msh file path, then run this script from inside
# the top-level extracted folder (i.e. run: bash scripts/setup.sh).
#
# NOTE: the chemistry/thermo database is NOT a placeholder anymore - the
# actual files (Park 1993 11-species reactions, species thermo data, V-T
# relaxation model) were copied verbatim from the real hy2Foam repo
# (https://github.com/ivanZanardi/hypersonicfoam) into
# case/constant/chemDicts/, case/constant/thermoDEM, and
# case/constant/thermo2TModel, and are referenced via $FOAM_CASE-relative
# paths in case/constant/thermophysicalProperties. Nothing to edit there.
# ============================================================================

set -e

# --------------------------------------------------------------------------
# PLACEHOLDER: point this at your actual Gmsh mesh file
# --------------------------------------------------------------------------
MESH_FILE="~/plain_BWB/Gmsh_without_cavities_AoA_neg5_20260822.msh"

CASE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../case" && pwd)"

echo "=========================================="
echo "hy2Foam case setup"
echo "=========================================="
echo "Case directory : $CASE_DIR"
echo "Mesh file      : $MESH_FILE"
echo ""

if [ ! -f "$MESH_FILE" ]; then
    echo "ERROR: MESH_FILE not found at: $MESH_FILE"
    echo "Edit scripts/setup.sh and set MESH_FILE to your actual .msh path."
    exit 1
fi

echo "Step 1: converting mesh with gmshToFoam"
cd "$CASE_DIR"
gmshToFoam "$MESH_FILE"

echo "Step 2: checking mesh integrity"
checkMesh -allTopology -allGeometry

echo "Step 3: listing boundary patch names from the converted mesh"
echo "-----------------------------------------------------------"
echo "IMPORTANT: compare these against the placeholder patch names"
echo "(farfield, wall, outlet) used in case/0/* boundary conditions."
echo "Rename patches in Gmsh or edit the 0/ files if they do not match."
echo "The repo's own genericCase example uses different patch names"
echo "(inlet, cylinder) for its own geometry - patch names are always"
echo "specific to your mesh, never assume ours match yours."
echo "-----------------------------------------------------------"
foamDictionary -entry boundary -expand constant/polyMesh/boundary || true

echo ""
echo "Setup complete. Review the patch name comparison above before running."
echo "Next: bash scripts/run.sh"

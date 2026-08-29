#!/bin/bash
#SBATCH --job-name=OF_PreProc_CavityBWB
#SBATCH --partition=small
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=180G
#SBATCH --time=04:00:00          # Bumped to 4 hours due to slower uncollated file writing
#SBATCH --output=preproc-%j.out

# Load OpenFOAM v1706 
source $SCRATCH/$USER/OpenFOAM/OpenFOAM-v1706/etc/bashrc

cp -r ~/cavity_BWB $SCRATCH/$USER/cavity_BWB
cd $SCRATCH/$USER/cavity_BWB/hypersonic_case/case

echo "Starting gmshToFoam..."
gmshToFoam Gmsh_with_cavities_AoA_neg5_20260830.msh -scale 0.001 > ../../log.gmshToFoam.log 2>&1

echo "Starting uncollated decomposePar..."
# Stripped out the -fileHandler flag entirely to match v1706 capability
decomposePar > ../../log.decomposePar.log 2>&1

echo "Pre-processing complete!"

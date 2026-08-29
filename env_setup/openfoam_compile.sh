#!/bin/bash
#SBATCH --job-name=compile_foam
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --time=04:00:00
#SBATCH --partition=small
#SBATCH --output=compile_%j.out
#SBATCH --error=compile_%j.err

# 1. Load system compiler and MPI modules ???
module load gcc/9.3.0 openmpi/4.0.5 cmake/3.20.0
# (Adjust module names/versions based on available modules: 'module avail')

# 2. Source the OpenFOAM environment
source $SCRATCH/$USER/OpenFOAM/OpenFOAM-v1706/etc/bashrc

# 3. Compile ThirdParty libraries first
cd $WM_THIRD_PARTY_DIR
./Allwmake -j 16

# 4. Compile OpenFOAM core and solvers
cd $WM_PROJECT_DIR
./Allwmake -j 16
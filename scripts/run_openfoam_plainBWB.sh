#!/bin/bash
#SBATCH --job-name=of_plainBWB
#SBATCH --nodes=8
#SBATCH --ntasks-per-node=48
#SBATCH --cpus-per-task=1
#SBATCH --partition=medium
#SBATCH --exclusive
#SBATCH --time=48:00:00
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err

# 1. Environment Diagnostics
echo "Job ID: $SLURM_JOBID"
echo "Allocated Nodes: $SLURM_JOB_NODELIST"
echo "Total Tasks: $SLURM_NTASKS"

# Unlimit stack and core files
ulimit -s unlimited
ulimit -c unlimited

source $SCRATCH/$USER/OpenFOAM/OpenFOAM-v1706/etc/bashrc
cd $SCRATCH/$USER/plain_BWB/hypersonic_case/case

echo "Checking Mesh Quality..."
srun -n 384 checkMesh -allTopology -allGeometry -parallel 2>&1 | tee ../../log.checkmesh.log

#------------
#module load spack
#export SPACK_ROOT=/home/apps/SPACK
#. $SPACK_ROOT/share/spack/setup-env.sh
#spack load gcc@13.3.0
#spack load intel-oneapi-compilers@2024.2.1
#spack load intel-mpi@2021.11.0
#spack load openfoam
#------------

# Source system-wide OpenFOAM environment (modify path based on your module load setup)
# Example: source /home/apps/OpenFOAM/v2312/etc/bashrc
# If loaded via module: spack load openfoam

# Execution Workflow with Proper Logging and Pipe Redirection
# Launch OpenFOAM solver with srun on 384 MPI ranks
SOLVER="hy2Foam"  # Replace with your specific solver (e.g., scalarTransportFoam, hyStrath, etc.)

echo "Executing $SOLVER on $SLURM_NTASKS ranks..."
srun -n 384 --cpu-bind=cores $SOLVER -parallel 2>&1 | tee ../../log.$SOLVER.log
EXIT_CODE=$?
echo "$SOLVER exited with code $EXIT_CODE"
exit "$EXIT_CODE"

#!/bin/bash
# ============================================================================
# submit_384.slurm - batch submission for the 384-core production run
#
#   sbatch scripts/submit_384.slurm
#
# 675,000 fixed 45 ns steps on a 31M cell mesh across 384 ranks is a long job.
# Size --time from a short scaling run before committing: launch with
# endTime set to a few thousand steps, read the "ExecutionTime" per step from
# the log, and multiply. As a starting estimate, an explicit density-based
# solver on this mesh typically runs 1 to 3 seconds of wall time per 1000
# steps per 80k cells per rank, which puts the full run in the region of
# 24 to 72 hours. Set the walltime from your own measurement, not from that
# range, and use the restart note at the bottom if it does not fit in one
# allocation.
# ============================================================================
#SBATCH --job-name=bwb-rcf
#SBATCH --nodes=8
#SBATCH --ntasks-per-node=48          # 8 x 48 = 384 MPI ranks
#SBATCH --cpus-per-task=1
#SBATCH --time=48:00:00
#SBATCH --exclusive
#SBATCH --output=logs/bwb-rcf-%j.out
#SBATCH --error=logs/bwb-rcf-%j.err

set -euo pipefail

module purge
module load gcc openmpi                       # site-specific
source /opt/openfoam2412/etc/bashrc            # site-specific

# rhoCentralFoam is pure MPI. One thread per rank; any BLAS or OpenMP
# threading underneath would oversubscribe the cores.
export OMP_NUM_THREADS=1

PKG_DIR="$SLURM_SUBMIT_DIR"
CASE_DIR="$PKG_DIR/case"
mkdir -p "$PKG_DIR/logs"
cd "$CASE_DIR"

if [ ! -f constant/polyMesh/owner ]; then
    echo "ERROR: no mesh. Run scripts/setup.sh on a login node first."
    exit 1
fi

if [ ! -d processors384 ] && [ ! -d processor0 ]; then
    decomposePar -force -fileHandler collated
fi

srun --cpu-bind=cores rhoCentralFoam -parallel -fileHandler collated

# ---------------------------------------------------------------------------
# Restarting
# ---------------------------------------------------------------------------
# system/controlDict has startFrom latestTime, so resubmitting this script
# picks up from the newest of the 20 volume writes with no edit. Two things to
# know about a restart:
#   - the probe and surface output files are appended to, not overwritten, so
#     the time series stays continuous across the restart
#   - fieldAverage has restartOnRestart false, so the running mean survives;
#     it does NOT restart the averaging window
# The cost of a restart is up to 33,749 lost steps, since only the volume
# writes are restartable. If the queue limit forces frequent restarts, raise
# the volume write frequency (lower writeInterval in system/controlDict) and
# accept the extra disk.

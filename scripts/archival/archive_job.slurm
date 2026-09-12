#!/bin/bash
#SBATCH --job-name=foamarch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=01:00:00
#SBATCH --output=%x-%j.out

# Single-core housekeeping job. Runs one archive pass, then resubmits itself
# one hour later. Stops when the stop file appears or the solver job is gone.

export CASE=/scratch/$USER/bwb_mach10
export ARCH=$CASE/archive
export OFFLOAD=$HOME/foam_archive
export KEEP=2
export AGE=20
export THREADS=$SLURM_CPUS_PER_TASK
export LEVEL=3
export SCRATCH_PCT=80

INTERVAL=20          # minutes between passes
SOLVER_JOB=${SOLVER_JOB:-}   # optional job id of the running solver
STOP=$CASE/archive.stop

bash "$HOME/bin/foam_archive.sh"

if [ -f "$STOP" ]; then
    echo "stop file present, not resubmitting"
    exit 0
fi

if [ -n "$SOLVER_JOB" ]; then
    if ! squeue -h -j "$SOLVER_JOB" >/dev/null 2>&1 || \
       [ -z "$(squeue -h -j "$SOLVER_JOB" 2>/dev/null)" ]; then
        echo "solver job finished, final pass done, not resubmitting"
        exit 0
    fi
fi

sbatch --begin=now+${INTERVAL}minutes \
       --export=ALL,SOLVER_JOB="$SOLVER_JOB" \
       "$0"

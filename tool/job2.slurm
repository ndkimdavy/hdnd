#!/bin/bash
#SBATCH --job-name=genclosure1d
#SBATCH --partition=amd
#SBATCH --time=1-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --output=tool/%x.out
#SBATCH --error=tool/%x.err

set -euo pipefail

cd "$HOME/hdnd"

export PATH="$HOME/.juliaup/bin:$PATH"
julia --version

export JULIA_NUM_THREADS="$SLURM_CPUS_PER_TASK"
export OPENBLAS_NUM_THREADS=1

cd tool

REJECT=none LOG=genclosure1d julia -t "$JULIA_NUM_THREADS" -O3 -C native --project=.. genclosure1d.jl

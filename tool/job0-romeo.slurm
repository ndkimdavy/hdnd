#!/bin/bash
#SBATCH --job-name=hdnd
#SBATCH --account=r250142
#SBATCH --partition=long
#SBATCH --constraint=x64cpu
#SBATCH --time=2-00:00:00
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=8
#SBATCH --ntasks=16
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --output=tool/%x-%j.out
#SBATCH --error=tool/%x-%j.err

set -euo pipefail

XINPUT="${XINPUT:-input/config1d.txt}"
XWORK="${XWORK:-$HOME/hdnd}"
XOUT="${XOUT:-/gpfs/scratch/$USER}"

XOUT="$XOUT/$SLURM_JOB_NAME/$SLURM_JOB_ID"
XCONFIG="$XOUT/config1d.txt"

cd "$XWORK"

romeo_load_x64cpu_env
spack load openmpi@5.0.5 %aocc

unset OMPI_MCA_mtl
unset FI_PROVIDER

export PATH="$HOME/.juliaup/bin:$HOME/.julia/bin:$PATH"

export JULIA_NUM_THREADS="$SLURM_CPUS_PER_TASK"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

mkdir -p "$XOUT"

cp "$XINPUT" "$XCONFIG"

sed -i -E "s|^[[:space:]]*nprocx[[:space:]]*=.*|nprocx = $SLURM_NTASKS|" "$XCONFIG"
sed -i -E "s|^[[:space:]]*outpath[[:space:]]*=.*|outpath = $XOUT/out|" "$XCONFIG"

echo "XINPUT   = $XINPUT"
echo "XWORK    = $XWORK"
echo "XOUT     = $XOUT"
echo "XCONFIG  = $XCONFIG"
echo "NNODE    = $SLURM_NNODES"
echo "NPROC    = $SLURM_NTASKS"
echo "NTHREAD  = $JULIA_NUM_THREADS"
echo ""

srun julia -t "$JULIA_NUM_THREADS" -O3 -C native --project=. \
    launch/app.jl -f "$XCONFIG"

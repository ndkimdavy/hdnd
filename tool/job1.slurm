#!/bin/bash
#SBATCH --job-name=train-ai
#SBATCH --partition=visu
#SBATCH --time=1-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --output=tool/%x.out
#SBATCH --error=tool/%x.err

set -euo pipefail

cd "$HOME/hdnd"

source "$HOME/.conda/miniconda3/etc/profile.d/conda.sh"
conda activate dev

export OMP_NUM_THREADS="$SLURM_CPUS_PER_TASK"
export MKL_NUM_THREADS="$SLURM_CPUS_PER_TASK"

export CUDA_VISIBLE_DEVICES=0

python --version
python -c 'import os, torch; print(torch.__version__); print(torch.cuda.is_available()); print(os.environ["CUDA_VISIBLE_DEVICES"]); print(torch.cuda.get_device_name(0))'

cd tool

papermill train.ipynb train_executed.ipynb \
  --log-output \
  --progress-bar \
  --request-save-on-cell-execute \
  --autosave-cell-every 60

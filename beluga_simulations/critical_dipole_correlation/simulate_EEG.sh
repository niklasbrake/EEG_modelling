#!/bin/bash
#SBATCH --time=30:00:00
#SBATCH --account=def-akhadra
#SBATCH --mem=4G
#SBATCH --mail-user=niklas.brake@mail.mcgill.ca
#SBATCH --mail-type=FAIL,END
#SBATCH --output=simulation.log

folder=$1
S=$2

module load python/3.8.10
module load mpi4py
module load scipy-stack
module load matlab/2020a
module load neuron

virtualenv --no-download $SLURM_TMPDIR/env
source $SLURM_TMPDIR/env/bin/activate
pip install --no-index --upgrade pip --quiet
pip install --no-index LFPy --quiet

matlab -nodisplay -r "simulate_EEG('$folder',${SLURM_ARRAY_TASK_ID},$S)"


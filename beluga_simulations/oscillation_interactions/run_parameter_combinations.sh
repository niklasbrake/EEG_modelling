#!/bin/bash
#SBATCH --time=12:00:00
#SBATCH --account=def-akhadra
#SBATCH --mem=4G
#SBATCH --mail-user=niklas.brake@mail.mcgill.ca
#SBATCH --mail-type=FAIL,END
#SBATCH --output=simulation.log

module load python/3.8.10
module load mpi4py
module load scipy-stack
module load matlab/2020a
module load neuron

virtualenv --no-download $SLURM_TMPDIR/env
source $SLURM_TMPDIR/env/bin/activate
pip install --no-index --upgrade pip --quiet
# pip install --no-index -r requirements.txt --quiet
pip install --no-index LFPy --quiet
pip install --no-index umap-learn

matlab -nodisplay -r "run_parameter_combinations(${SLURM_ARRAY_TASK_ID})"


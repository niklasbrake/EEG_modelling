#!/bin/bash
#SBATCH --time=002:00:00
#SBATCH --account=def-akhadra
#SBATCH --mem=40G
#SBATCH --mail-user=niklas.brake@mail.mcgill.ca
#SBATCH --mail-type=FAIL
#SBATCH --cpus-per-task=10
#SBATCH --output=embedding.log
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

module load python/3.8.10
module load mpi4py
module load scipy-stack

virtualenv --no-download $SLURM_TMPDIR/env
source $SLURM_TMPDIR/env/bin/activate
pip install --no-index --upgrade pip --quiet
pip install --no-index umap-learn

folder=$1
subfolders=($(ls -d $folder/*))

for i in {0..9}; do
    b=$((10*(${SLURM_ARRAY_TASK_ID}-1)+i))
    tempFolder=${subfolders[b]}

    correlationFile=$tempFolder'/presynaptic_network/correlations.csv'
    embededFile=$tempFolder'/presynaptic_network/UMAP_embedding.csv'

    python /lustre04/scratch/nbrake/code/simulation_code/embed_data.py $correlationFile $embededFile
done


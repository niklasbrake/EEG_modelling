#!/bin/bash
#SBATCH --time=000:20:00
#SBATCH --account=def-akhadra
#SBATCH --mem=4G
#SBATCH --mail-user=niklas.brake@mail.mcgill.ca
#SBATCH --mail-type=FAIL,END

masterPath='/lustre04/scratch/nbrake/data/simulations/raw/synapse_embedding'
mValues='/home/nbrake/aperiodic_EEG_modelling/simulations/mValues'
m=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $mValues)
# DIR=$1

module load python/3.8.10
module load mpi4py
module load scipy-stack
module load matlab/2020a
module load neuron

virtualenv --no-download $SLURM_TMPDIR/env
source $SLURM_TMPDIR/env/bin/activate
pip install --no-index --upgrade pip --quiet
pip install --no-index -r requirements.txt --quiet
pip install --no-index LFPy --quiet
pip install --no-index umap-learn

for (( j=1; j<=4; j++ ))
do
    DIR="${masterPath}/m=${m}/rep${j}"
    for (( k=1; k<=10; k++ ))
    do
        matlab -nodisplay -r "runSimulationRep('$DIR',$k)"
    done
done
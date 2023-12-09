# Run simulations for m=0.98

cd /lustre04/scratch/nbrake/code/critical_dipole_correlation
folder='/lustre04/scratch/nbrake/data/simulations/embed_branching_network/m=0.98'

# Initialize network models in matlab and simulate crit dynamics
ID1=$(sbatch --parsable --array=1-10 initialize_network.sh $folder 0.98)

# Compute pairwise correlations in presynaptic network
ID2=$(sbatch --depend=afterok:$ID1 --array=1-10 --parsable compute_correlations.sh $folder)

# Embed network onto sphere using the UMAP algorithm
ID3=$(sbatch --depend=afterok:$ID2 --array=1-10 --parsable embed_data.sh $folder)

# Run simulations for varying embedding optimalities
sValues=("0.0" "0.2" "0.4" "0.6" "0.8" "1.0")
for i in {0..5}; do
    s=${sValues[i]}
    # Simulate single-neuron EEG signal
    ID1=$(sbatch --array=1-10 --parsable simulate_EEG.sh $folder $s)
    # Analyze results
    sbatch --depend=afterok:$ID1 analyze_simulations.sh $folder $s
done
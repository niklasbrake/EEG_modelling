# gcc compute_tiling_correlation.c -o compute_tiling_correlation.exe -lm -fopenmp

cd /lustre04/scratch/nbrake/code/critical_dipole_correlation

folder='/lustre04/scratch/nbrake/data/simulations/embed_branching_network/m=0.98'

# Initialize network models in matlab and simulate crit dynamics
ID1=$(sbatch --parsable --array=1-10 initialize_network.sh $folder 0.98)

# Compute pairwise correlations in presynaptic network
ID2=$(sbatch --depend=afterok:$ID1 --array=1-10 --parsable compute_correlations.sh $folder)

# Embed network onto sphere using the UMAP algorithm
ID3=$(sbatch --depend=afterok:$ID2 --array=1-10 --parsable embed_data.sh $folder)

# Simulate single-neuron EEG signal
ID4=$(sbatch --depend=afterok:$ID3 --array=1-10 --parsable simulate_EEG.sh $folder)


# Simulate single-neuron EEG signal
sbatch --depend=afterok:$ID4 analyze_simulations.sh $folder)


folder='/lustre04/scratch/nbrake/data/simulations/embed_branching_network/m=0.98'
sValues=("0.0" "0.2" "0.4" "0.6" "0.8" "1.0")

for i in {0..5}; do
    s=${sValues[i]}
    ID1=$(sbatch --array=1-10 --parsable simulate_EEG.sh $folder $s)
    sbatch --depend=afterok:$ID1 analyze_simulations.sh $folder $s
done




folder='/lustre04/scratch/nbrake/data/simulations/embed_branching_network/m=0.98'
sValues=("0.2" "0.4" "0.6" "0.8")

for i in {0..3}; do
    s=${sValues[i]}
    # sbatch --array=1-5 simulate_EEG.sh $folder $s
    sbatch analyze_simulations.sh $folder $s
done

folder='/lustre04/scratch/nbrake/data/simulations/embed_branching_network'
subfolders=($(ls -d $folder/*))

for i in {0,1,2,4,5,6}; do
    tempFolder=${subfolders[i]}
    echo $tempFolder
    sbatch compute_PSD.sh $tempFolder
done



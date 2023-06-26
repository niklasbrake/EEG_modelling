cd /lustre04/scratch/nbrake/code/parameter_sensitivity_analysis
ID1=$(sbatch --array=1-20 --parsable run_parameter_combinations.sh)
ID2=$(sbatch --array=1-20 --depend=afterok:$ID1 --parsable analyze_simulations.sh)
sbatch --depend=afterok:$ID2 combine_unitary_spectra.sh
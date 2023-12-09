ID1=$(sbatch --array=1-25 --parsable run_parameter_combinations.sh)
ID2=$(sbatch --depend=afterok:$ID1 --parsable analyze_simulations.sh)
sbatch --depend=afterok:$ID2 combine_unitary_spectra.sh
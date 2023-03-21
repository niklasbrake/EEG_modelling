set folder1=%1\presynaptic_network\spikeTimesLong.csv
set folder1=%folder1:\=/%
set folder2=%1\presynaptic_network\correlations.csv
set folder2=%folder2:\=/%
C:\Users\brake\Documents\GitHub\aperiodic_EEG_modelling\simulations\functions\compute_tiling_correlation.exe %folder1% %folder2%
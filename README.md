[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10359818.svg)](https://doi.org/10.5281/zenodo.10359818)

## Executing code in this repository

To reproduce the figures, data needs to be downloaded from https://doi.org/10.6084/m9.figshare.24777990. Once this data is downloaded, line 16 of EEG_modelling/model/network_simulation_beluga.m needs to be updated, such that the variable resourceFolder points to the data folder.

If you wish to simulate dipoles, the script EEG_modelling/simulation_examples/example_embedding.m outlines how the model can be simulated. To run this script, the file EEG_modelling/model/compute_tiling_correlation.c needs to be compiled, which can be accomplished with the following command:
````
gcc compute_tiling_correlation.c -o compute_tiling_correlation.exe -lm -fopenmp
````
The python packages in requirements.txt need to be installed, which can be accomplished with the following command:
````
pip install -r requirements.txt
````
Finally, the .MOD files need to be compiled for the neuron simulator. To do so, navigate to the folder EEG_modelling/model/mod_files and run the command
```
mknrndll
```

## Detrending EEG spectra
The file EEG_modelling/data_analysis/detrending.py contains the code used to fit Eq. 1, 5, and 6 from the manuscript to EEG spectra. The upper and lower bounds for the various parameters were optimized for our data and may need to be adjusted. EEG_modelling/data_analysis/synDetrend.m is a wrapper to run this function from MATLAB. 

## Acknowledgements
I completed this work as part of my PhD under the supervision of [Dr. Anmar Khadra](http://www.medicine.mcgill.ca/physio/khadralab/) and in collaboration with Dr. Gilles Plourde at McGill Univeristy. The EEG data used in this project was collected by Dr. Plourde and his team.

## License
<a rel="license" href="http://creativecommons.org/licenses/by/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by/4.0/88x31.png" /></a><br />This repository is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by/4.0/">Creative Commons Attribution 4.0 International License</a>.

For completeness, this repository includes fmriView, which was taken from the GitHub repository [ecog_fmri_visualization_matlab](https://github.com/edden-gerber/ecog_fmri_visualization_matlab) under a BSD 2-Clause "Simplified" License. This license is retained in the subdirectory EEG_modelling/auxiliary/fmriView.

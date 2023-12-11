## What is this?
In this repository is the code associated with the following paper:

**Brake N, Duc F, Rokos A, Arseneau F, Shahiri S, Khadra A, and Plourde G (2021). "A neurophysiological basis for aperiodic EEG and the background spectral trend."** *****Submitted.*****

This repository includes functions to generate all figures from the manuscript (under manuscript_figures), as well as the original code used to simulate the model (simulations).

To reproduce the figures, data needs to be downloaded from https://doi.org/10.6084/m9.figshare.24777990.v2. Once this data is downloaded, line 16 of model/network_simulation_beluga.m needs to be updated, such that the variable resourceFolder points to the data folder.

If you wish to simulate dipoles, the script simulation_examples/example_embedding.m outlines how the model is simulated. To run this script, the file model/compute_tiling_correlation.c needs to be compiled, which can be accomplished with the following command:

gcc compute_tiling_correlation.c -o compute_tiling_correlation.exe -lm -fopenmp

Furthermore, the python packages in requirements.txt need to be installed, which can be accomplished with the following command:

pip install -r requirements.txt

## Acknowledgements
I completed this work during 2020-2023 as part of my PhD under the supervision of [Dr. Anmar Khadra](http://www.medicine.mcgill.ca/physio/khadralab/) and in collaboration with Dr. Gilles Plourde at McGill Univeristy. The EEG data used in this project was collected by Dr. Plourde.

gcc compute_tiling_correlation.c -o compute_tiling_correlation.exe -lm -fopenmp

## License
<a rel="license" href="http://creativecommons.org/licenses/by/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by/4.0/88x31.png" /></a><br />This repository is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by/4.0/">Creative Commons Attribution 4.0 International License</a>.
c
To reference this code, please cite the article mentioned above.

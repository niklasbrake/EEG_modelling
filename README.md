## What is this?
In this repository is the code associated with the following paper:

**Brake N, Duc F, Rokos A, Arseneau F, Shahiri S, Khadra A, and Plourde G (2021). "Aperiodic EEG activity masks the dynamics of neural oscillations during loss of consciousness from propofol."** *****Submitted.*****

Preprint: https://doi.org/10.1101/2021.10.12.464109

This repository includes the functions used to generate the figures from the paper. Certain data is included in the data folder. Not all the data required to reporduce the figures are included in this repository. For raw EEG timeseries and spectrograms, please contact gilles.plourde@mcgill.ca.

## Dependencies
All the MATLAB code was developed in MATLAB 2017b. Spectral exponent estimation was done using the Python implementation of the FOOOF package (Donoghue et al. 2020) (FOOOF 1.0.0 executed from Python 3.8.10).

## Acknowledgements
I completed the work for this project in 2020-2021 during my PhD under the supervision of [Dr. Anmar Khadra](http://www.medicine.mcgill.ca/physio/khadralab/) and in collaboration with Dr. Gilles Plourde at McGill Univeristy. The EEG data used in this project was collected by Dr. Plourde. 

gcc compute_tiling_correlation.c -o compute_tiling_correlation.exe -lm -fopenmp

## License
<a rel="license" href="http://creativecommons.org/licenses/by/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by/4.0/88x31.png" /></a><br />This repository is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by/4.0/">Creative Commons Attribution 4.0 International License</a>.

To reference this code, please cite the article mentioned above.

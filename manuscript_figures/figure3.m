function figure2(dataFolder)

if(nargin<1)
    error('Path to data required as input argument. Data can be downloaded from link in README file.');
end

path = mfilename('fullpath');
basePath = fileparts(fileparts(path));
simulationPath2 = fullfile(basePath,'beluga_simulations','synthetic_spikes');


run(fullfile(simulationPath2,'plot_example_eeg.m'));
run(fullfile(simulationPath2,'plot_results.m'));


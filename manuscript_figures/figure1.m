function figure1(dataFolder)

if(nargin<1)
    error('Path to data required as input argument. Data can be downloaded from link in README file.');
end

path = mfilename('fullpath');
basePath = fileparts(fileparts(path));
simulationPath1 = fullfile(basePath,'beluga_simulations','parameter_sensitivity_analysis');

% run(fullfile(simulationPath1,'plot_example_morphology.m'));
run(fullfile(simulationPath1,'plot_example_dipoles.m'));
% run(fullfile(simulationPath1,'plot_example_source_locations.m'));
% run(fullfile(simulationPath1,'plot_fitted_tau.m'));
% run(fullfile(simulationPath1,'plot_parameter_distribution.m'));
% run(fullfile(simulationPath1,'plot_asynchronous_amplitude_estimate.m'));
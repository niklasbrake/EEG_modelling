function figure2(dataFolder)

if(nargin<1)
    error('Path to data required as input argument. Data can be downloaded from link in README file.');
end

path = mfilename('fullpath');
basePath = fileparts(fileparts(path));
simulationPath1 = fullfile(basePath,'beluga_simulations','oscillation_interactions','mixed_input');

run(fullfile(simulationPath1,'plot_results.m'));

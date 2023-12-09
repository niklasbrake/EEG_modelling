function figure2(dataFolder)

if(nargin<1)
    error('Path to data required as input argument. Data can be downloaded from link in README file.');
end

% Add paths for function execution
myPath = mfilename('fullpath');
basePath = fileparts(fileparts(myPath));
addpath(fullfile(basePath,'auxiliary'));
addpath(fullfile(basePath,'auxiliary','fmriView'));
addpath(fullfile(basePath,'model'));
addpath(fullfile(basePath,'data_analysis'));

% Run subscripts for plotting results
filePath = fullfile(basePath,'simulations','parameter_sensitivity_analysis','_plotting');
run(fullfile(filePath,'plot_parameter_distribution.m'));
run(fullfile(filePath,'plot_sensitivity_analysis.m'));
run(fullfile(filePath,'plot_slope_vs_EI.m'));


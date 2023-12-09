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
filePath = fullfile(basePath,'simulations','oscillation_interactions','single_input_unitary_spectra','_plotting');
run(fullfile(filePath,'plot_results.m'));
run(fullfile(filePath,'plot_supplemental.m'));

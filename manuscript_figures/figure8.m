function figure8(dataFolder)

if(nargin<1)
    error('Path to data required as input argument. Data can be downloaded from link in README file.');
end

% Add paths for function execution
myPath = mfilename('fullpath');
basePath = fileparts(fileparts(myPath));
addpath(fullfile(basePath,'auxiliary'));
addpath(fullfile(basePath,'model'));
addpath(fullfile(basePath,'data_analysis'));

% Run subscripts for plotting results
filePath = fullfile(basePath,'data_analysis','tau1_extraction');
% run(fullfile(filePath,'plot_results.m'));
run(fullfile(filePath,'plot_timeseries.m'));
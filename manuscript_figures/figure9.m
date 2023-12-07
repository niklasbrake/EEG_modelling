function figure9(dataFolder)

if(nargin<1)
    error('Path to data required as input argument. Data can be downloaded from link in README file.');
end

path = mfilename('fullpath');
basePath = fileparts(fileparts(path));
dataPath = fullfile(basePath,'data_analysis','detrending');

run(fullfile(dataPath,'plot_results.m'));
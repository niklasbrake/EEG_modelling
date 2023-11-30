function figure9(dataFolder)

if(nargin<1)
    error('Path to data required as input argument. Data can be downloaded from link in README file.');
end

run('../data_analysis/detrending/plot_results.m')
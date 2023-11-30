function figure8(dataFolder)

if(nargin<1)
    error('Path to data required as input argument. Data can be downloaded from link in README file.');
end

run('../data_analysis/tau1_extraction/plot_results.m')
run('../data_analysis/tau1_extraction/marsh_model_analysis.m')
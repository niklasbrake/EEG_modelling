function figure8(dataFolder)

if(nargin<1)
    error('Path to data required as input argument. Data can be downloaded from link in README file.');
end

% error('parameters need fixing');
run('../data_analysis/tau1_extraction/other_electrode_plots.m')

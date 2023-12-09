function combine_unitary_spectra

baseFolder = '/lustre04/scratch/nbrake/data/simulations/parameter_sensitivity_analysis2';
addpath('/lustre04/scratch/nbrake/code/simulation_code');

pars = zeros(9,2e4);
m = 1e3;

for iChunk = 1:20
    filename = fullfile(baseFolder,['analyzed_results_' int2str(iChunk) '.mat']);
    load(filename,'f','P','parameters');
    if(iChunk==1)
        psd = zeros(size(P,1),2e4);
    end
    psd(:,(iChunk-1)*m+1:iChunk*m) = P;
    pars(:,(iChunk-1)*m+1:iChunk*m) = parameters;
end

saveFile = '/lustre04/scratch/nbrake/data/simulation_analyzed/parameter_sensitivity_analysis2.mat';
save(saveFile,'f','psd','pars');

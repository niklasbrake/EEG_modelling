function analyze_simulations

% baseFolder = 'E:\Research_Projects\004_Propofol\data\simulations\raw\aperiodic_interactions';
baseFolder = '/lustre04/scratch/nbrake/data/simulations/parameter_sensitivity_analysis';
addpath('/lustre04/scratch/nbrake/code/simulation_code');


fs = 2e3;
dp = zeros(4001,3,10000);
pars = zeros(7,10000);
for iChunk = 1:20
    load(fullfile(baseFolder,['analyzed_results_' int2str(iChunk) '.mat']),'dipoles','parameters');

    idcs = [(iChunk-1)*500+1:iChunk*500];
    pars(:,idcs) = parameters;
    for k = 1:size(dipoles,3)
        dp(:,:,idcs(k)) = detrend(resample(dipoles(:,:,k),fs,16e3),'constant');
    end
end

save(fullfile(baseFolder,'analyzed_results_all.mat'),'dp','pars','fs');
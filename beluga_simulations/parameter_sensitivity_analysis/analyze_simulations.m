function analyze_simulations(iChunk)

baseFolder = '/lustre04/scratch/nbrake/data/simulations/parameter_sensitivity_analysis2';
addpath('/lustre04/scratch/nbrake/code/simulation_code');

[sa,X] = network_simulation_beluga.getHeadModel;
locations = sa.cortex2K.in_from_cortex75K;

filename = fullfile(baseFolder,['analyzed_results_' int2str(iChunk) '.mat']);
load(filename,'dipoles','parameters');

k = 100;
n = 10;

dp = [];
for k = 1:size(dipoles,3)
    dp(:,:,idcs(k)) = detrend(resample(dipoles(:,:,k),1e3,16e3),'constant');
end
dipoles = dp;

m = 100;
for i = 1:10
    dp = dipoles(:,:,m*(i-1)+1:m*i);
    eeg = zeros(size(dp,1),m*length(locations));
    G = zeros(1,m*length(locations));
    for k = 1:length(locations)
        eeg(:,m*(k-1)+1:m*k) = network_simulation_beluga.getEEG(dp,sa,locations(k));
        G(m*(k-1)+1:m*k) = 1:m;
    end
    psd = mypmtm(eeg,1e3,10);
    Ptemp = splitapply(@(x) mean(x,2),psd,G);
    if(i==1)
        P = nan(size(psd,1),size(dipoles,3));
    end
    P(:,m*(i-1)+1:m*i) = Ptemp;
end

f = [0.1:0.1:500]
save(filename,'f','P',"-append");


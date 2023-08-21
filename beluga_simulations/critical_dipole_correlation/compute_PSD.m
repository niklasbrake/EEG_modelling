function initialize_network(masterFolder)

addpath('/lustre04/scratch/nbrake/code/simulation_code');

[sa,X] = network_simulation_beluga.getHeadModel;
locations = sa.cortex2K.in_from_cortex75K;

% dp = zeros(2001,3,100);
% for arrayID = 1:10
%     for j = 1:10
%         i = 10*(arrayID-1)+j;
%         folder = fullfile(masterFolder,['run_' num2str(i,'%.3d')]);
%         for rep = 1:100
%             % saveFolder = fullfile(folder,['simulation' int2str(rep)]);
%             saveFolder = fullfile(folder,['simulation' int2str(rep)]);
%             load(fullfile(saveFolder,'simulation_data.mat'));
%             dp(:,:,rep) = dp(:,:,rep) + resample(sum(dipoles,3),1e3,16e3);
%         end
%     end
% end
% time = time(1:16:end);
% dipoles = dp;
load(fullfile(masterFolder,'dipoles.mat'),'dipoles','time');

m=10;
for i = 1:10
    dp = dipoles(:,:,m*(i-1)+1:m*i);
    eeg = zeros(size(dp,1),m*length(locations));
    G = zeros(1,m*length(locations));
    for k = 1:length(locations)
        eeg(:,m*(k-1)+1:m*k) = network_simulation_beluga.getEEG(dp,sa,locations(k));
        G(m*(k-1)+1:m*k) = 1:m;
    end
    psd = mypmtm(detrend(eeg,'constant'),1e3,2);
    Ptemp = splitapply(@(x) mean(x,2),psd,G);
    if(i==1)
        P = nan(size(psd,1),size(dp,3));
    end
    P(:,m*(i-1)+1:m*i) = Ptemp;
end

f = [0.5:0.5:500];
save(fullfile(masterFolder,'PSD.mat'),'f','P');


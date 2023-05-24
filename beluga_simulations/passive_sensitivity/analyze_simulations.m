[sa,X] = network_simulation_beluga.getHeadModel;

folder = '/lustre04/scratch/nbrake/data/simulations/passive_sensitivity';
saveFolder = '/lustre04/scratch/nbrake/data/simulation_analyzed/passive_sensitivity';

F = dir(folder); F = F(3:end);

% locations = randi(size(sa.cortex75K.vc,1),20,1);
locations = sa.cortex2K.in_from_cortex75K;
for i = 1:length(F)
    folder2 = fullfile(folder,F(i).name);
    F2 = dir(folder2); F2 = F2(3:end);
    for j = 1:length(F)
        load(fullfile(folder2,F2(j).name,'simulation','simulation_data.mat'));
        m = size(dipoles,3);
        dp = [];
        for k = 1:m
            dp(:,:,k) = detrend(resample(dipoles(:,:,k),2e3,16e3),'constant');
        end
        eeg = zeros(size(dp,1),m*length(locations));
        G = zeros(1,m*length(locations));
        for k = 1:length(locations)
            eeg(:,m*(k-1)+1:m*k) = network_simulation_beluga.getEEG(dp,sa,locations(k));
            G(m*(k-1)+1:m*k) = 1:m;
        end
        psd = mypmtm(eeg,2e3,2);
        P = splitapply(@(x) mean(x,2),psd,G);
        f = [0.5:0.5:1e3]
        fid = fopen(fullfile(folder2,F2(j).name,'simulation','_parameters.json'));
        str=char(fread(fid,inf)')
        fclose(fid); 
        pars = jsondecode(str);

        load(fullfile(folder2,F2(j).name,'model.mat'))
        mTypes = {network.neurons(:).mType};
        save(fullfile(saveFolder,['analyzed_' F(i).name '_' F2(j).name '.mat']),'f','P','pars','mType','i','j');
    end
end
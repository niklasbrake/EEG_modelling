load('/home/nbrake/data/resources/head_models/sa_nyhead.mat');
idcs = sa.cortex2K.in_from_cortex75K;

masterPath = '/lustre04/scratch/nbrake/data/simulations/raw/synapse_embedding_2';
F = dir(masterPath);
F = F(3:end);

m = sort([0.999,0.997,0.99,0.98,0.95,0.86,0.63,0]);

C = zeros(length(F),10);
m0 = C;
P = zeros(2e3,10,length(F));
for i = 1:length(F)
    folder = fullfile(masterPath,F(i).name);
    F2 = dir(folder);
    F2 = F2(3:end);
    for j = 1:length(F2)
        for rep = 1:10
            folder2 = fullfile(folder,F2(j).name,['LFPy_' int2str(rep)]);
            load(fullfile(folder2,'simulation_data.mat'));
            for k = 1:length(idcs)
                eeg = getEEG(dipoles,sa,idcs(k));
                C(i,j) = C(i,j) + corr(eeg(:,1),eeg(:,2));
            end
        end
        C(i,j) = C(i,j)/10/length(idcs);
        m0(i,j) = csvread(fullfile(folder,F2(j).name,'presynaptic_network','emperical_m.csv'));
    end
end

save(fullfile(masterPath,'correlations.mat'),'C');

function W = getEEG(Q,sa,idx)
    % Lead field for Cz recording
    L0 = squeeze(sa.cortex75K.V_fem(49,idx,:))'; % czIDX = 49

    % Orient dipole "up" direction to be normal to cortical surface
    vz = sa.cortex75K.normals(idx,:);
    [vx,vy] = getOrthBasis(vz);
    q = (Q(:,1,:).*vx+Q(:,2,:).*vy+Q(:,3,:).*vz)*1e-12; % nAum -> mAm

    W = squeeze(sum(L0.*q,2)*1e6); % V -> uV
end

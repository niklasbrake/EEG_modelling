load('/home/nbrake/data/resources/head_models/sa_nyhead.mat');
load('/home/nbrake/data/passive/dipoles.mat');
Q = gpuArray(dp);
idcs = 1:size(sa.cortex75K.normals,1);
L1 = squeeze(sa.cortex75K.V_fem(49,:,:)); % czIDX = 49
V = gpuArray(zeros(length(idcs),600));
for j = 1:length(idcs)
    L0 = L1(idcs(j),:);
    vz = sa.cortex75K.normals(idcs(j),:);
    [vx,vy] = getOrthBasis(vz);

    q = (Q(:,1,:).*vx+Q(:,2,:).*vy+Q(:,3,:).*vz)*1e-12; % nAum -> mAm

    W = squeeze(sum(L0.*q,2)*1e6); % V -> uV
    V(j,:) = var(W);
end
V = gather(V);
save('/home/nbrake/data/passive/totalVariance.mat','V');
function initialize_network(masterPath,branchNo)
% masterPath = '/lustre04/scratch/nbrake/data/simulations/raw/correlation_test';
addpath('/home/nbrake/aperiodic_EEG_modelling/simulations/functions');

tmax0 = 20e3;
tmax = 2e3;
nrnCount = 1;

network = network_simulation_beluga(masterPath);
network = network.initialize_postsynaptic_network(nrnCount,ones(nrnCount,1));
network.tmax = tmax;
% network = network.setsynapsecount(1e3);
N = network.getsynapsecount;

% Plane
elevation = rand(N,1);
azimuth = rand(N,1);
y1 = asin(2*elevation-1);
x1 = 2*pi*azimuth;
dist_metric = @(x,y)vecnorm(x-y,2,2);
dt = 4e-3;
nNeigh = 4;

idcs1 = find(y1>0);
N1 = floor(length(idcs1)*0.8);
idcs1E = idcs1(1:N1);
C = zeros(N1*nNeigh,3);
for i = 1:N1
    D0 = dist_metric([elevation(idcs1E(i)),azimuth(idcs1E(i))],[elevation(idcs1E),azimuth(idcs1E)]);
    [~,I] = sort(D0,'ascend'); I = I(2:nNeigh+1);
    C(nNeigh*(i-1)+1:nNeigh*i,1) = i*ones(nNeigh,1);
    C(nNeigh*(i-1)+1:nNeigh*i,2) = I;
end
C(:,3) = branchNo/nNeigh;
[ids0,ts0,ei0] = simulatespikes_osc(N1,branchNo,tmax0*1e-3,C,0);


idcs2 = find(y1<=0);
N2 = floor(length(idcs2)*0.8);
idcs2E = idcs2(1:N2);
C = zeros(N2*nNeigh,3);
for i = 1:N2
    D0 = dist_metric([elevation(idcs2E(i)),azimuth(idcs2E(i))],[elevation(idcs2E),azimuth(idcs2E)]);
    D0 = dist_metric([elevation(idcs2E(i)),azimuth(idcs2E(i))],[elevation(idcs2E),azimuth(idcs2E)]);
    [~,I] = sort(D0,'ascend'); I = I(2:nNeigh+1);
    C(nNeigh*(i-1)+1:nNeigh*i,1) = i*ones(nNeigh,1);
    C(nNeigh*(i-1)+1:nNeigh*i,2) = I;
end
C(:,3) = branchNo/nNeigh;
[ids1,ts1,ei1] = simulatespikes_osc(N2,branchNo,tmax0*1e-3,C,1);

% [ids,ts,ei,~,m0] = simulatespikes_det(N,branchNo,tmax0*1e-3,C,nNeigh);
ei = [ei0(:);ei1(:)];
ids = [idcs1(ids0);idcs2(ids1)];
ts = [ts0(:);ts1(:)];

file = fullfile(network.preNetwork,'spikeTimesLong.csv');
network_simulation_beluga.save_presynaptic_network(ids,ts,ei,N,file)

correlationFile = fullfile(network.preNetwork,'correlations.csv');
network = network.setCorrelationFile(correlationFile);

csvwrite(fullfile(network.preNetwork,'connections.csv'),C);
csvwrite(fullfile(network.preNetwork,'locations.csv'),[elevation(:),azimuth(:)]);
csvwrite(correlationFile,'');
% csvwrite(fullfile(network.preNetwork,'emperical_m.csv'),m0);

csvwrite(fullfile(network.preNetwork,'UMAP_embedding.csv'),[(1:N)',x1(:),y1(:)]);

network.save();
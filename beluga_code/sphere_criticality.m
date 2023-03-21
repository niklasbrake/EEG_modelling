function sphere_criticality(masterPath)
addpath('/home/nbrake/aperiodic_EEG_modelling/simulations/functions');

% Model parameters
branchNo = 0.98;
tmax0 = 2e3;
tmax = 2e3;
nrnCount = 2;

% Initalize network structure
network = network_simulation_beluga(masterPath);
network = network.initialize_postsynaptic_network(nrnCount,ones(nrnCount,1));
network.tmax = tmax;
N = network.getsynapsecount;

% Generate presyanptic network topology
elevation = asin(2*rand(N,1)-1);
azimuth = 2*pi*rand(N,1);
[x,y,z] = sph2cart(azimuth,elevation,1);
dist_metric = @(x,y) network_simulation.haversine_distance2(x,y);
dt = 4e-3;
nNeigh = 4;
pmax = branchNo/nNeigh;
C = zeros(N*nNeigh,4);
for i = 1:N
    D0 = dist_metric([elevation(i),azimuth(i)],[elevation,azimuth]);
    D = exp(-D0/0.02);
    D(i) = min(D);
    idcs = find(D>1e-9);
    I = [];
    for j = 1:nNeigh
        j0 = interp1(cumsum(D(idcs))/sum(D(idcs)),idcs,rand,'next','extrap');
        I(j) = j0;
        idcs = setdiff(idcs,j0);
    end
    C(nNeigh*(i-1)+1:nNeigh*i,1) = i*ones(nNeigh,1);
    C(nNeigh*(i-1)+1:nNeigh*i,2) = I;
    C(nNeigh*(i-1)+1:nNeigh*i,3) = rand(nNeigh,1)*pmax*2;
    C(nNeigh*(i-1)+1:nNeigh*i,4) = dt*2*rand(nNeigh,1); % Random transmission speed
end

% Simulate presynaptic spikes
[ids,ts,ei,C] = simulatespikes_det(N,branchNo,tmax0*1e-3,[],nNeigh);
file = fullfile(network.preNetwork,'spikeTimesLong.csv');
network_simulation_beluga.save_presynaptic_network(ids,ts,ei,N,file);

% Save subsample of presyanptic spikes for simulation
t0 = rand*max(tmax0-tmax,0);
idcs = find(and(ts>=t0,ts<t0+tmax));
N = network.getsynapsecount;
network_simulation_beluga.save_presynaptic_network(ids(idcs),ts(idcs)-t0,ei,N,network.spikingFile);

% Touch correlation file
correlationFile = fullfile(network.preNetwork,'correlations.csv');
network = network.setCorrelationFile(correlationFile);
csvwrite(correlationFile,'');

% Save sphereical coordinates and bypass umap embedding
umapFile = fullfile(network.preNetwork,'UMAP_embedding.csv');
csvwrite(umapFile, [(1:N)',azimuth,elevation-pi/2]);

% Form connections and simulate postsyanptic network
network.form_connections(1);
network = network.simulate();
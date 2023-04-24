fldr = '/lustre04/scratch/nbrake/data/simulations/active_crit/run1';
addpath('/lustre04/scratch/nbrake/code/simulation_code');

% eFiringRate = 0.5;
% iFiringRate = 2.5;

% Initialize network
network = network_simulation_beluga(fldr);
% network = network.setFiringRate(eFiringRate,iFiringRate);

% Initialize post network
nPostNeurons = 2;
network = network.initialize_postsynaptic_network(nPostNeurons,ones(nPostNeurons,1));

% Presyanptic network parameters
nPreNeurons = 30e3;
network.tmax = 2e3; % 2 seconds
network.branchNo = 0.98;
network.simulatespikes_critplane(nPreNeurons,network.tmax);

% Explicit mapping of plane onto sphere to save time (bypass UMAP; ignore seam)
copyfile(fullfile(network.preNetwork,'spikeTimes.csv'),fullfile(network.preNetwork,'spikeTimesLong.csv'));
x = csvread(fullfile(network.preNetwork,'locations.csv'));
x1 = 2*pi*x(:,1);
y1 = asin(2*x(:,2)-1);
data = [(1:nPreNeurons)'-1,x1(:),y1(:)];
dlmwrite(fullfile(network.preNetwork,'UMAP_embedding.csv'),data,'precision','%.4f');

% Form random connections
network = network.form_connections(1);

% Simulate dipoles
network.savePath = fullfile(fldr,'LFPy_passive');
network = network.simulate();

network = network.addActiveChannels();
network.savePath = fullfile(fldr,'LFPy_active');
network = network.simulate();

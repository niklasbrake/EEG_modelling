folder = fullfile(pwd,'example_simulation');

% Initialize network
network = network_simulation_beluga(folder);

% Initialize post network
nPostNeurons = 1;
mTypes = 1;
network = network.initialize_postsynaptic_network(nPostNeurons,mTypes);

% Presyanptic network parameters
N = network.getsynapsecount;
network.tmax = 2e3; % 2 seconds
network.branchNo = 0.98;
network.simulatespikes_critplane(m,network.tmax);

%% Compute pairwise correlations and UMAP embed onto sphere
% network = network.compute_presynaptic_correlations(network.spikingFile);
% network.embed_presyanptic_neurons;

%% Explicit mapping of plane onto sphere to save time (bypass UMAP; ignore seam)
% x = csvread(fullfile(folder,'presynaptic_network\locations.csv'));
% copyfile(fullfile(network.preNetwork,'spikeTimes.csv'),fullfile(network.preNetwork,'spikeTimesLong.csv'));
% x1 = 2*pi*x(:,1);
% y1 = asin(2*x(:,2)-1)-pi/2;
% data = [(1:N)'-1,x1(:),y1(:)];
% dlmwrite(fullfile(network.preNetwork,'UMAP_embedding.csv'),data,'precision','%.4f');

% Place syanpse optimally (input between 0 and 1, with 1 optimal and 0 random)
network = network.form_connections(0);

% Place syanpse optimally (input between 0 and 1, with 1 optimal and 0 random)
network.simulate();
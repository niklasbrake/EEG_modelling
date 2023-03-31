% Initialize network
network = network_simulation_beluga(fullfile(pwd,'avaETosc_example'));

% Initialize post network
nPostNeurons = 1;
network = network.initialize_postsynaptic_network(nPostNeurons,1);
N = 30000;
network.tmax = 2e3; % 2 seconds

[ids,ts,ei,x] = simulatespikes_avaETosc(preNetwork);
network_simulation_beluga.save_presynaptic_network(ids,ts,ei,N,network.spikingFile)

% Explicit mapping of plane onto sphere to save time (bypass UMAP; ignore seam)
x1 = 2*pi*x(:,1);
y1 = asin(2*x(:,2)-1);
data = [(1:N)',x1(:),y1(:)];
dlmwrite(fullfile(network.preNetwork,'UMAP_embedding.csv'),data,'precision','%.4f');
network.form_connections(1);

% Simulate dipoles
network = network.simulate();
[time,V,dipoles] = network.importSimulationResults;
folder = 'E:\Research_Projects\004_Propofol\data\simulations\raw';
% Initialize network
network = network_simulation_beluga(fullfile(folder,'EI_oscillations','WT'));

% Initialize post network
nPostNeurons = 1;
network = network.initialize_postsynaptic_network(nPostNeurons,1);

% Presyanptic network parameters
nPreNeurons = 30000;
network.tmax = 2e3; % 2 seconds
network.branchNo = 0;
oscillation_frequency = 10;
simulatespikes_osc(nPreNeurons,network.tmax,network,oscillation_frequency);

% Explicit mapping of plane onto sphere to save time (bypass UMAP; ignore seam)
copyfile(fullfile(network.preNetwork,'spikeTimes.csv'),fullfile(network.preNetwork,'spikeTimesLong.csv'));
x = csvread(fullfile(network.preNetwork,'locations.csv'));
x1 = 2*pi*x(:,1);
y1 = asin(2*x(:,2)-1);
data = [(1:nPreNeurons)'-1,x1(:),y1(:)];
dlmwrite(fullfile(network.preNetwork,'UMAP_embedding.csv'),data,'precision','%.4f');

% Place syanpse optimally (input between 0 and 1, with 1 optimal and 0 random)
network.form_connections(1);

% Simulate dipoles
network = network.simulate();
% [time,V,dipoles] = network.importSimulationResults;

% Compute EEG signal
% [sa,X] = network_simulation_beluga.getHeadModel;
% location = randi(size(sa.cortex75K.vc,1)); % Random location
% eeg = network_simulation_beluga.getEEG(dipoles,sa,location);
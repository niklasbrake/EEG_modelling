folder = 'E:\Research_Projects\004_Propofol\data\simulations\raw\test\critical_embedding';
% Initialize network
network = network_simulation_beluga(folder);

% Initialize post network
nPostNeurons = 2;
mTypes = 1;
network = network.initialize_postsynaptic_network(nPostNeurons,[1,7]);

% Presyanptic network parameters
nPreNeurons = 30e3;
network.tmax = 5e3; % 2 seconds
network.branchNo = 0.98;
network.simulatespikes_critplane(nPreNeurons,network.tmax);

network = network.compute_presynaptic_correlations(network.spikingFile)

network.embed_presyanptic_neurons;

% Place syanpse optimally (input between 0 and 1, with 1 optimal and 0 random)
network = network.form_connections(1);

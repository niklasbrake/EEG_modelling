folder = 'E:\Research_Projects\004_Propofol\data\simulations\raw\test\EI_network(m=0)';
% Initialize network
network = network_simulation_beluga(folder);

% Initialize post network
nPostNeurons = 2;
mTypes = 1;
network = network.initialize_postsynaptic_network(nPostNeurons,[1,1]);

% Presyanptic network parameters
N = 30e3;
network.tmax = 10e3; % 2 seconds

[ids,ts,ei,P,trueL] = run_EInetwork(N,0);
network_simulation_beluga.save_presynaptic_network(ids,ts,ei,N,network.spikingFile)
% network = network.compute_presynaptic_correlations(network.spikingFile)
% network.embed_presyanptic_neurons;

x1 = P(:,2);
y1 = P(:,1)+pi/2;
data = [(1:N)'-1,x1(:),y1(:)];
dlmwrite(fullfile(network.preNetwork,'UMAP_embedding.csv'),data,'precision','%.4f');

% Place syanpse optimally (input between 0 and 1, with 1 optimal and 0 random)
network = network.form_connections(1);

network.savePath = fullfile(network.outputPath,'simulation_gamE');
network.parameters.eSynParams.weight = 1.4e-3;
network = network.simulate();
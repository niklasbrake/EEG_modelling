
folder = 'E:\Research_Projects\004_Propofol\data\simulations\raw\test\EI_network';
% Initialize network
network = network_simulation_beluga(folder);

% Initialize post network
nPostNeurons = 1;
mTypes = 1;
network = network.initialize_postsynaptic_network(nPostNeurons,mTypes);

% Presyanptic network parameters
network.tmax = 2e3; % 2 seconds
load('E:\Research_Projects\004_Propofol\data\simulations\raw\EI_network\EI_network.mat')
preIdcs = find(postNeurons(1,:));
idsPre = [];
tsPre = [];
for i = 1:length(preIdcs)
    idcs = find(ids == preIdcs(i));
    t0 = ts(idcs);
    idsPre = [idsPre;i*ones(length(t0),1)];
    tsPre = [tsPre;t0(:)];
end
eiPre = ei(preIdcs);
network = network.setsynapsecount(length(preIdcs))

network_simulation_beluga.save_presynaptic_network(idsPre,tsPre,eiPre,network.getsynapsecount,network.spikingFile);

% Place syanpse optimally (input between 0 and 1, with 1 optimal and 0 random)
network = network.form_connections(0);

% Simulate dipoles
network = network.simulate();

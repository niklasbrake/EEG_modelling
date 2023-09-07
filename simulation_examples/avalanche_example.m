folder = 'E:\Research_Projects\004_Propofol\data\simulations\raw\example_embedding';
% Initialize network
network = network_simulation_beluga(folder);

% Initialize post network
nPostNeurons = 2;
mTypes = 1;
network = network.initialize_postsynaptic_network(nPostNeurons,[1,2]);

% Presyanptic network parameters
N = 30e3;
network.tmax = 6e3; % 2 seconds
network.branchNo = 0.98;

CC = 1;
m = floor(N/CC);
M = 0:m:N;
ids = {};
ts = {};
ei = {};
for i = 1:CC
    m = M(i+1)-M(i)+1;
    [ids{i},ts{i},ei{i}] = network.simulatespikes_critplane(m,network.tmax);
    ids{i} = ids{i} + sum(M(1:i));
end
ids = cat(1,ids{:});
ts = cat(1,ts{:});
ei = cat(1,ei{:});

csvwrite(fullfile(network.preNetwork,'multisynapse_IDs.csv'),-ones(N,1));
network_simulation_beluga.save_presynaptic_network(ids,ts,ei,N,network.spikingFile)


network = network.compute_presynaptic_correlations(network.spikingFile);
network.embed_presyanptic_neurons;

figureNB;
C = csvread('E:\Research_Projects\004_Propofol\data\simulations\raw\test\connected_components\CC1\presynaptic_network\UMAP_embedding.csv');
[x,y,z] = sph2cart(C(:,2),C(:,3)-pi/2,1+0*C(:,2));
scatter3(x,y,z,5,C(:,1),'filled')
set(gca,'DataAspectRatio',[1,1,1])

% Explicit mapping of plane onto sphere to save time (bypass UMAP; ignore seam)
% x = csvread(fullfile(folder,'presynaptic_network\locations.csv'));
% copyfile(fullfile(network.preNetwork,'spikeTimes.csv'),fullfile(network.preNetwork,'spikeTimesLong.csv'));
% x1 = 2*pi*x(:,1);
% y1 = asin(2*x(:,2)-1)-pi/2;
% data = [(1:N)'-1,x1(:),y1(:)];
% dlmwrite(fullfile(network.preNetwork,'UMAP_embedding.csv'),data,'precision','%.4f');

% Place syanpse optimally (input between 0 and 1, with 1 optimal and 0 random)
network = network.form_connections(1);


% Place syanpse optimally (input between 0 and 1, with 1 optimal and 0 random)
network.simulate();
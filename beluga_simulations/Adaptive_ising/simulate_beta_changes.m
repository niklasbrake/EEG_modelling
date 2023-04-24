addpath('/lustre04/scratch/nbrake/code/simulation_code');
folder = '/lustre04/scratch/nbrake/data/simulations/adaptive_ising_beta';
beta = [0.7,0.9,0.99,1.01,2];
for i = 1:length(beta)
    runSimulations(folder,beta(i));
end


function runSimulations(folder,beta)
    % Initialize network
    network = network_simulation_beluga(fullfile(folder,num2str(beta)));

    % Initialize post network
    nPostNeurons = 2;
    network = network.initialize_postsynaptic_network(nPostNeurons,ones(nPostNeurons,1));
    N = 30000;
    network.tmax = 2e3; % 2 seconds

    [ids,ts,ei,x] = simulatespikes_avaETosc(network,beta,N);
    network_simulation_beluga.save_presynaptic_network(ids,ts,ei,N,network.spikingFile)
    csvwrite(fullfile(network.preNetwork,'multisynapse_IDs.csv'),-ones(N,1));

    % Explicit mapping of plane onto sphere to save time (bypass UMAP; ignore seam)
    x1 = 2*pi*x(:,1);
    y1 = asin(2*x(:,2)-1);
    data = [(0:N-1)',x1(:),y1(:)];
    dlmwrite(fullfile(network.preNetwork,'UMAP_embedding.csv'),data,'precision','%.4f');
    network.form_connections(1);

    % Simulate dipoles
    network = network.simulate();
end
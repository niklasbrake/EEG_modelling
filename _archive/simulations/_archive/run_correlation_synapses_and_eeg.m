function main
    sim_dir = 'E:\Research_Projects\004_Propofol\Modelling\neuron_simulations\data\simulations\pairwise_correlation';
    tmax = 2e3;
    N = 10;
    network0 = network_simulation;
    network0.postNetwork = sim_dir;
    network0 = network0.initialize_postsyanptic_network(2*N,0.98);
    network0.preNetwork = fullfile(sim_dir,'spiking.csv');
    network0.outputPath = fullfile(sim_dir);
    network0 = network0.initialize_presyanptic_network(0,tmax);


    pair_corr = linspace(0,1,11);
    for i = 1:length(pair_corr)
        pc = pair_corr(i)
        network = network0;
        network.postNetwork = fullfile(sim_dir,num2str(pc,2));
        network.outputPath = fullfile(sim_dir,num2str(pc,2));
        network = network.initialize_postsyanptic_network(N,pc);
        network = network.simulate();
        network.save();
    end
end
[sa,X] = network_simulation.getHeadModel;

N = 30;
tmax = 2e3;

masterPath = 'E:\Research_Projects\004_Propofol\data\simulations\raw\Propofol\baseline';
network = network_simulation(masterPath);
% network = network.addPropofol;
network = network.initialize_postsynaptic_network(N);
network = network.initialize_presynaptic_network(0,tmax);
network.form_connections;
network = network.simulate;
network = network.importResults;

masterPath = 'E:\Research_Projects\004_Propofol\data\simulations\raw\Propofol\propofol_triple';
network = network_simulation(masterPath);
network = network.addPropofol;
network = network.initialize_postsynaptic_network(N);
network = network.initialize_presynaptic_network(0,tmax);
network.form_connections;
network = network.simulate;
network = network.importResults;


function main(N,tmax,corr_idx,sim_dir)
    cons = csvread('E:\Research_Projects\004_Propofol\data\simulations\raw\dyads\s=0\presynaptic_network\preConnections.csv');

    load('E:\Research_Projects\004_Propofol\data\simulations\raw\dyads\s=0\data.mat');
    i = 0;
    for m = [0,0.25,0.5,0.75,0.9,0.95,0.98,0.99,0.999];
        for k = 0:9
            newRunPath = ['E:\Research_Projects\004_Propofol\data\simulations\raw\dyad_network_criticality\m=' num2str(m,3) '\run' int2str(k)];
            fdr2 = ['s=' num2str(i,2)];
            masterPath3 = fullfile(newRunPath,fdr2);
            network2 = network.copy_network(masterPath3);
            N = network.getsynapsecount;
            [ids,ts,ei] = simulatespikes_det(N,m,network.tmax*1e-3,cons);
            network_simulation.save_presynaptic_network(ids,ts,ei,N,fullfile(network2.preNetwork,'spikeTimes.csv'))
            network2.form_connections(i);
            network2.simulate();
        end
    end
end


function run_correlation(m)
    sim_dir = ['E:\Research_Projects\004_Propofol\Modelling\neuron_simulations\data\simulations\pairwise_correlation\m' num2str(m,2)];

    tmax = 2e3;
    N = 10;
    network0 = network_simulation;
    network0.postNetwork = sim_dir;
    network0 = network0.initialize_postsyanptic_network(2*N,0.98);
    network0.preNetwork = fullfile(sim_dir,'spiking.csv');
    network0.outputPath = fullfile(sim_dir);
    network0 = network0.initialize_presyanptic_network(m,tmax);


    pair_corr = linspace(0,1,11);
    for i = 1:length(pair_corr)
        pc = pair_corr(i)
        network = network0;
        network.postNetwork = fullfile(sim_dir,['pc' num2str(pc,2)]);
        network.outputPath = fullfile(sim_dir,['pc' num2str(pc,2)]);
        network = network.initialize_postsyanptic_network(N,pc);
        network = network.simulate();
        network.save();
    end
end

function run_conductances()
    sim_dir = 'E:\Research_Projects\004_Propofol\Modelling\neuron_simulations\data\simulations\compare_conductances\active';

    tmax = 6e3;
    network = network_simulation;
    network.postNetwork = sim_dir;
    network = network.initialize_postsyanptic_network(100,0);
    network.preNetwork = fullfile(sim_dir,'spiking.csv');
    network.outputPath = fullfile(sim_dir);
    network = network.initialize_presyanptic_network(0,tmax);
    network = network.addActiveChannels;
    network = network.simulate();
    network.save();
end
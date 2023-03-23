function main(isWT)
    folder = '/lustre04/scratch/nbrake/data/simulations/passive/';
    eFiringRate = 0.5;
    iFiringRate = 2.5;
    if(isWT)
        fldr = fullfile(folder,'WT');
        passive_example(fldr,eFiringRate,iFiringRate)
    else
        fldr = fullfile(folder,'PMS');
        passive_example(fldr,1.5*eFiringRate,0.75*iFiringRate)
    end
end
function passive_example(fldr,eFiringRate,iFiringRate)

    addpath('/lustre04/scratch/nbrake/code/simulation_code');

    % Initialize network
    network = network_simulation_beluga(fldr);
    network = network.setFiringRate(eFiringRate,iFiringRate);

    % Initialize post network
    nPostNeurons = 100;
    network = network.initialize_postsynaptic_network(nPostNeurons);

    % Presyanptic network parameters
    nPreNeurons = network.getsynapsecount;
    network.tmax = 2e3; % 2 seconds
    network.branchNo = 0;

    % Set up Poissonian presynatpic network
    tmax = network.tmax+100;
    ME = floor(network.eiFraction*nPreNeurons);
    MI = nPreNeurons-ME;
    ei = [zeros(ME,1);ones(MI,1)];
    nE = poissrnd(network.eFiringRate*ME*tmax/1e3);
    nI = poissrnd(network.iFiringRate*MI*tmax/1e3);
    ts = tmax*rand(nE+nI,1);
    ids = [randi(ME,nE,1);ME+randi(MI,nI,1)];
    network_simulation_beluga.save_presynaptic_network(ids,ts,ei,nPreNeurons,network.spikingFile);

    % Form random connections
    network = network.form_connections(0);

    % Simulate dipoles
    network = network.simulate();
end
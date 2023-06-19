function main(isWT)
    folder = '/lustre04/scratch/nbrake/data/simulations/EI_crit/';
    eFiringRate = 0.5;
    iFiringRate = 2.5;
    for i = 1:5
        if(isWT)
            fldr = fullfile(folder,'WT',['run' int2str(i)]);
            passive_example(fldr,eFiringRate,iFiringRate)
        else
            fldr = fullfile(folder,'PMS',['run' int2str(i)]);
            passive_example(fldr,1.5*eFiringRate,0.75*iFiringRate)
        end
    end
end
function passive_example(fldr,eFiringRate,iFiringRate)

    addpath('/lustre04/scratch/nbrake/code/simulation_code');

    % Initialize network
    network = network_simulation_beluga(fldr);
    network = network.setFiringRate(eFiringRate,iFiringRate);

    % Initialize post network
    nPostNeurons = 1;
    network = network.initialize_postsynaptic_network(nPostNeurons,ones(nPostNeurons,1));

    % Presyanptic network parameters
    nPreNeurons = 30e3;
    network.tmax = 2e3; % 2 seconds
    network.branchNo = 0;
    network.eFiringRate = eFiringRate;
    network.iFiringRate = iFiringRate;
    network.simulatespikes_critplane(nPreNeurons,network.tmax);

    % Form random connections
    network = network.form_connections(0);

    % Simulate dipoles
    network = network.simulate();
end
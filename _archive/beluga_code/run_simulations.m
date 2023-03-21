function run_simulations(branchNo,prop,masterPath)
    addpath('/home/nbrake/aperiodic_EEG_modelling/simulations/functions');
    rng(1);
    N = 30e3;
    tmax0 = 4e3;
    tmax = 2e3;
    nrnCount = 2;

    network = network_simulation_beluga(masterPath);
    network = network.initialize_postsynaptic_network(nrnCount,ones(nrnCount,1));
    network.tmax = tmax;
    network.branchNo = branchNo;
    % N = network.getsynapsecount;
    % network = network.setsynapsecount(1e5);
    % network.simulatespikes_critplane(N,tmax0);
    [ids,ts,ei] = simulatespikes_osc(N,tmax0,network);

    % Explicit mapping of presynaptic connections onto sphere
    copyfile(fullfile(network.preNetwork,'spikeTimes.csv'),fullfile(network.preNetwork,'spikeTimesLong.csv'));
    x = csvread(fullfile(network.preNetwork,'locations.csv'));
    x1 = 2*pi*x(:,1);
    y1 = asin(2*x(:,2)-1);
    data = [(1:length(x1))',x1(:),y1(:)];
    dlmwrite(fullfile(network.preNetwork,'UMAP_embedding.csv'),data,'precision','%.4f');
    network.save();

    if(prop)
        propofol = [10,3,3]*[100;10;1];
    else
        propofol = [10,1,1]*[100;10;1];
    end
    network = network.addPropofol(propofol);
    network.form_connections(1);

    for repNo = 1:1
        [ids,ts,ei] = network.getprenetwork(fullfile(network.preNetwork,'spikeTimesLong.csv'));
        tmax = network.tmax;
        tmax0 = max(ts)+100;
        t0 = rand*max(tmax0-tmax,0);
        idcs = find(and(ts>=t0,ts<t0+tmax));

        network_simulation_beluga.save_presynaptic_network(ids(idcs),ts(idcs)-t0,ei,N,network.spikingFile)
        network.savePath = fullfile(masterPath,['LFPy_' int2str(repNo)]);
        network = network.simulate();
    end
end
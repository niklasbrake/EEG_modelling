function aperiodic_sensitivity(arrayID)

% baseFolder = '/lustre04/scratch/nbrake/data/simulations/aperiodic_sensitivity';
% funPath = '/lustre04/scratch/nbrake/code/simulation_code';
funPath = 'C:\Users\brake\Documents\GitHub\EEG_modelling\simulation_code';
baseFolder = 'E:\Research_Projects\004_Propofol\data\simulations\raw\test2';
addpath(funPath);

fid = fopen(fullfile(funPath,'mod_files','default_parameters.json'));
str = char(fread(fid,inf)');
fclose(fid);
pars = jsondecode(str);

simName = ['rhythm' int2str(arrayID)];
X = csvread([simName '.csv']);
lamFun = @(t0) interp1(X(:,1),X(:,2),t0);

for comboID = 1:1
    pars2 = changeParameters(pars,comboID);
    folder = fullfile(baseFolder,simName,['combo' int2str(comboID)]);
    runSimulation(folder,pars2,lamFun)
end


end
function pars = changeParameters(pars,comboID)
    switch comboID
    case 1
        % Baseline parameters (no change)
    case 2
        % Change tauI
        pars.iSynParams.tau2 = 20;
    case 3
        % Change EL
        pars.biophys_pars.pas_mem_pars.erev_leak = -45;
    case 4
        % Change firing rate without changing EI
        pars.eCellParams.firingRate = 2*pars.eCellParams.firingRate;
        pars.iCellParams.firingRate = 2*pars.iCellParams.firingRate;
    case 5
        % Change EI without changing firing rate
        meanFiringRate = pars.iCellParams.firingRate*(1-pars.eiFraction) + pars.eiFraction*pars.eCellParams.firingRate;
        pars.eCellParams.firingRate = meanFiringRate;
        pars.iCellParams.firingRate = meanFiringRate;
    end
end
function runSimulation(folder,pars,lamFun)
    % Initialize network
    network = network_simulation_beluga(folder,pars);

    % Initialize post network
    network = network.initialize_postsynaptic_network(10);

    % Presyanptic network parameters
    network.tmax = 1e3; % 2 seconds
    N = network.getsynapsecount;

    % Place syanpse randomly
    network = network.form_connections(0);
    cons = csvread(fullfile(network.postNetwork,'connections.csv'));
    synIDs = cons(:,3);
    mData = network.getmData;
    P = arrayfun(@(i,j) mData{i}.pos(j,:),cons(:,1),cons(:,2),'UniformOutput',false);
    P = cat(1,P{:});
    [thet,phi] = cart2sph(P(:,1),P(:,2),P(:,3));
    idcsUpper = synIDs(phi>0);
    idcsLower = synIDs(phi<=0);

    % Simulate coutnerphase rhythmic input into upper and lower hemispheres
    M1 = length(idcsLower);
    [ids1,ts1,ei1] = network.sample_spike_rate(lamFun,M1);
    M2 = length(idcsUpper);
    [ids2,ts2,ei2] = network.sample_spike_rate(@(x) 2-lamFun(x),M2);

    % Save network spike times
    ids = [idcsLower(ids1);idcsUpper(ids2)];
    ts = [ts1;ts2];
    ei = [ei1;ei2];
    network_simulation_beluga.save_presynaptic_network(ids,ts,ei,N,network.spikingFile)

    % Simulate dipoles
    network = network.simulate();
end

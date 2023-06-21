function aperiodic_sensitivity(arrayID)

baseFolder = '/lustre04/scratch/nbrake/data/simulations/trend_peak_interaction';
funPath = '/lustre04/scratch/nbrake/code/simulation_code';
% funPath = 'C:\Users\brake\Documents\GitHub\EEG_modelling\simulation_code';
% baseFolder = 'E:\Research_Projects\004_Propofol\data\simulations\raw\test2';
addpath(funPath);

fid = fopen(fullfile(funPath,'mod_files','default_parameters.json'));
str = char(fread(fid,inf)');
fclose(fid);
pars = jsondecode(str);

comboIDs = [1:5,1:5,1:5,1:5,1:5];
rhythmIDs = [zeros(1,5),ones(1,5),2*ones(1,5),3*ones(1,5),4*ones(1,5)];

rhythmID = rhythmIDs(arrayID);
comboID = comboIDs(arrayID);

simName = ['rhythm' int2str(rhythmID)];
if(rhythmID>0)
    X = csvread([simName '.csv']);
    lamFun = @(t0) interp1(X(:,1),X(:,2),t0,'linear','extrap');
else
    t = 0:1e4;
    lamFun = @(t0) interp1(t,ones(size(t)),t0,'linear','extrap');
end

pars2 = changeParameters(pars,comboID);
folder = fullfile(baseFolder,simName,['combo' int2str(comboID)]);
runSimulation(folder,pars2,lamFun)


end
function pars = changeParameters(pars,comboID)
    pars.eiFraction = 0.8;
    switch comboID
    case 1
        % Baseline parameters (no change)
    case 2
        % Change tauI
        pars.iSynParams.tau2 = 30;
    case 3
        % Change EL
        pars.biophys_pars.pas_mem_pars.erev_leak = -45;
    case 4
        % Change firing rate without changing EI
        pars.eCellParams.firingRate = 2*pars.eCellParams.firingRate;
        pars.iCellParams.firingRate = 2*pars.iCellParams.firingRate;
    case 5
        % Change EI without changing firing rate
        pars.iSynParams.weight = 14e-4;
    end
end
function runSimulation(folder,pars,lamFun)
    % Initialize network
    network = network_simulation_beluga(folder,pars);

    % Initialize post network
    network = network.initialize_postsynaptic_network(50);

    % Presyanptic network parameters
    network.tmax = 10e3; % 2 seconds
    N = network.getsynapsecount;

    % Place syanpse randomly
    network = network.form_connections(0);
    cons = csvread(fullfile(network.postNetwork,'connections.csv'));
    synIDs = cons(:,3);
    mData = network.getmData;
    P = arrayfun(@(i,j) mData{i}.pos(j,:),cons(:,1),cons(:,2),'UniformOutput',false);
    P = cat(1,P{:});
    [thet,phi] = cart2sph(P(:,1),P(:,2),P(:,3));
    [~,I] = sort(phi);
    [~,J] = sort(I);
    cons(:,3) = J;
    dlmwrite(fullfile(network.postNetwork,'connections.csv'),cons, 'precision', '%i');

    % ALTERNATIVE
    M1 = floor(N/2);
    M2 = N-M1;
    % Simulate coutnerphase rhythmic input into upper and lower hemispheres
    [ids1,ts1,ei1] = network.sample_spike_rate(lamFun,M1);
    [ids2,ts2,ei2] = network.sample_spike_rate(@(x) 2-lamFun(x),M2);

    % Save network spike times
    ids = [ids1;ids2+M1];
    ts = [ts1;ts2];
    ei = [ei1;ei2];
    network_simulation_beluga.save_presynaptic_network(ids,ts,ei,network.getsynapsecount,network.spikingFile)

    % Simulate dipoles
    network.parameters.eiFraction = 0.8;
    network = network.simulate();
end

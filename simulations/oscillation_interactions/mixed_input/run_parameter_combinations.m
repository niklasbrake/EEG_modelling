function aperiodic_sensitivity(model_case)

baseFolder = '/lustre04/scratch/nbrake/data/simulations/mixed_input';
funPath = '/lustre04/scratch/nbrake/code/simulation_code';
addpath(funPath);

fid = fopen(fullfile(funPath,'mod_files','default_parameters.json'));
str = char(fread(fid,inf)');
fclose(fid);
pars = jsondecode(str);

X1 = csvread('../rhythm1.csv');
X2 = csvread('../rhythm4.csv');


switch model_case
case 0
    % Baseline condition
    lamFun1 = @(t0) interp1(X1(:,1),max(1+1*(1-X1(:,2)),0),t0,'linear','extrap');
    lamFun2 = @(t0) interp1(X2(:,1),max(1+1*(1-X2(:,2)),0),t0,'linear','extrap');
    folder = fullfile(baseFolder,'baseline');
case 1
    % Increase oscillation
    lamFun1 = @(t0) interp1(X1(:,1),max(1+5*(1-X1(:,2)),0),t0,'linear','extrap');
    lamFun2 = @(t0) interp1(X2(:,1),max(1+1*(1-X2(:,2)),0),t0,'linear','extrap');
    folder = fullfile(baseFolder,'high_oscillation');
case 2
    % Increase aperiodic
    lamFun1 = @(t0) interp1(X1(:,1),max(1+1*(1-X1(:,2)),0),t0,'linear','extrap');
    lamFun2 = @(t0) interp1(X2(:,1),max(1+5*(1-X2(:,2)),0),t0,'linear','extrap');
    folder = fullfile(baseFolder,'high_aperiodic');
case 3
    % Increase both
    lamFun1 = @(t0) interp1(X1(:,1),max(1+5*(1-X1(:,2)),0),t0,'linear','extrap');
    lamFun2 = @(t0) interp1(X2(:,1),max(1+5*(1-X2(:,2)),0),t0,'linear','extrap');
    folder = fullfile(baseFolder,'high_both');
case 4
    % Oscillation only case
    lamFun1 = @(t0) interp1(X1(:,1),max(1+1*(1-X1(:,2)),0),t0,'linear','extrap');
    lamFun2 = @(t0) interp1(X2(:,1),max(1+1*(1-X1(:,2)),0),t0,'linear','extrap');
    folder = fullfile(baseFolder,'oscillation_only');
case 5
    % Aperiodic only case
    lamFun1 = @(t0) interp1(X1(:,1),max(1+1*(1-X2(:,2)),0),t0,'linear','extrap');
    lamFun2 = @(t0) interp1(X2(:,1),max(1+1*(1-X2(:,2)),0),t0,'linear','extrap');
    folder = fullfile(baseFolder,'aperiodic_only');
end

runSimulation(folder,pars,lamFun1,lamFun2);
end
function runSimulation(folder,pars,lamFun1,lamFun2)
    % Initialize network
    network = network_simulation_beluga(folder,pars);

    % Initialize post network
    network = network.initialize_postsynaptic_network(1);

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

    M11 = floor(M1/2);
    M12 = M1-M11;
    M21 = floor(M2/2);
    M22 = M2-M21;

    % Simulate coutnerphase rhythmic input into upper and lower hemispheres
    [ids11,ts11,ei11] = network.sample_spike_rate(lamFun1,M11);
    [ids12,ts12,ei12] = network.sample_spike_rate(lamFun2,M12);

    [ids21,ts21,ei21] = network.sample_spike_rate(@(x) 2-lamFun1(x),M21);
    [ids22,ts22,ei22] = network.sample_spike_rate(@(x) 2-lamFun2(x),M22);


    apicalIdcs = randperm(M1)';
    ids_apical = [apicalIdcs(ids11);apicalIdcs(ids12+M11)];

    basalIdcs = randperm(M2)';
    ids_basal = [basalIdcs(ids21);basalIdcs(ids22+M21)];
    % Save network spike times
    ids = [ids_apical;ids_basal+M1];
    ts = [ts11;ts12;ts21;ts22];
    ei = [ei11;ei12;ei21;ei22];
    network_simulation_beluga.save_presynaptic_network(ids,ts,ei,network.getsynapsecount,network.spikingFile)

    % Simulate dipoles
    network = network.simulate();
end

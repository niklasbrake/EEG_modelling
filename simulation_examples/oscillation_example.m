folder = 'E:\Research_Projects\004_Propofol\data\simulations\raw\oscillation_test3';
% Initialize network1
network1 = network_simulation_beluga(folder);

% Initialize post network1
nPostNeurons = 1;
network1 = network1.initialize_postsynaptic_network(nPostNeurons,1);
network1.tmax = 4e3; % 2 seconds
N = network1.getsynapsecount;

% Presyanptic network1 parameters
%{
nPreNeurons = 30000;
network1.branchNo = 0;
oscillation_frequency = 2;
% simulatespikes_osc(nPreNeurons,network1.tmax+100,network1,oscillation_frequency);
simulatespikes_osc(nPreNeurons,network1.tmax,network1,oscillation_frequency);

% Explicit mapping of plane onto sphere to save time (bypass UMAP; ignore seam)
copyfile(fullfile(network1.preNetwork,'spikeTimes.csv'),fullfile(network1.preNetwork,'spikeTimesLong.csv'));
x = csvread(fullfile(network1.preNetwork,'locations.csv'));
x1 = 2*pi*x(:,1);
y1 = asin(2*x(:,2)-1)-pi/2;
data = [(1:nPreNeurons)'-1,x1(:),y1(:)];
dlmwrite(fullfile(network1.preNetwork,'UMAP_embedding.csv'),data,'precision','%.4f');

% Place syanpse optimally (input between 0 and 1, with 1 optimal and 0 random)
network1.form_connections(1);
%}

X(:,1) = 0:1e4; X(:,2) = 1+0.7*sin(2*pi*X(:,1)*2*1e-3);
lamFun = @(t0) interp1(X(:,1),X(:,2),t0);

% Place syanpse randomly
network1 = network1.form_connections(0);
cons = csvread(fullfile(network1.postNetwork,'connections.csv'));
synIDs = cons(:,3);
mData = network1.getmData;
P = arrayfun(@(i,j) mData{i}.pos(j,:),cons(:,1),cons(:,2),'UniformOutput',false);
P = cat(1,P{:});
[thet,phi] = cart2sph(P(:,1),P(:,2),P(:,3));
[~,I] = sort(phi);
[~,J] = sort(I);
cons(:,3) = J;
csvwrite(fullfile(network1.postNetwork,'connections.csv'),cons);


% ALTERNATIVE
M1 = floor(N/2);
M2 = N-M1;
% Simulate coutnerphase rhythmic input into upper and lower hemispheres
[ids1,ts1,ei1] = network1.sample_spike_rate(lamFun,M1);
[ids2,ts2,ei2] = network1.sample_spike_rate(@(x) 2-lamFun(x),M2);

% Save network1 spike times
ids = [ids1;ids2+M1];
ts = [ts1;ts2];
ei = [ei1;ei2];
network_simulation_beluga.save_presynaptic_network(ids,ts,ei,network1.getsynapsecount,network1.spikingFile)

% Simulate dipoles
network1.parameters.eiFraction = 0.8;
network1.parameters.iSynParams.tau2 = 10;
network1.savePath = fullfile(folder,'simulations_tau10');
network1 = network1.simulate();


network1.parameters.iSynParams.tau2 = 20;
network1.savePath = fullfile(folder,'simulations_tau20');
network1 = network1.simulate();

% Set path and folder
% folder = '/lustre04/scratch/nbrake/data/simulations/synthetic_spikes';
folder = 'E:\Research_Projects\004_Propofol\data\simulations\raw\test\synth_spikes';
addpath('/lustre04/scratch/nbrake/code/simulation_code');

% Initialize network
nPostNeurons = 2;
mTypes = 5;
network = network_simulation_beluga(folder);
network = network.initialize_postsynaptic_network(nPostNeurons,5*ones(nPostNeurons,1));
network = network.setsynapsecount(500);
N = network.getsynapsecount;
network.tmax = 2e3;

% Distribute synapses across dendrites and get locations
network = network.form_connections(0);
cons = csvread(fullfile(network.postNetwork,'connections.csv'));
synIDs = cons(:,3);
mData = network.getmData;
P = arrayfun(@(i,j) mData{i}.pos(j,:),cons(:,1),cons(:,2),'UniformOutput',false);
P = cat(1,P{:});
[thet,phi] = cart2sph(P(:,1),P(:,2),P(:,3));
X = [thet(:),phi(:)];

% Define excitatory/inhibitory cells and avg. firing rates
dt = 1;
ei = rand(N,1)>network.parameters.eiFraction;
r = (network.parameters.eCellParams.firingRate*(1-ei) + network.parameters.iCellParams.firingRate*ei)*dt/1e3;
rVar = r.*(1-r);

% Max correlation
R = 0.2;

% Use haversine distance among synapses to construct covariance matrix
KxxCell = cell(N,3); % Store sparse cov matrix (i,j,v) for every synapse
hav = @(x) (1-cos(x))/2;
hav_d = @(p1,x) hav(p1(1)-x(:,1))+(1-hav(x(:,1)-p1(1))-hav(p1(1)+x(:,1))).*hav(p1(2)-x(:,2));
parfor i = 1:N
    d = hav_d(X(i,:),X(i+1:end,:));
    idcs = find(d<0.25);
    KxxCell(i,:) = {idcs+i,i+0*zeros(size(idcs)),R*exp(-10*d(idcs)).*sqrt(rVar(i)*rVar(i+idcs))};
end
i0 = cat(1,KxxCell{:,1});
j0 = cat(1,KxxCell{:,2});
S = cat(1,KxxCell{:,3});

% Find mu and sig of n-d Gaussian with which to sample spikes
gam = icdf('normal',r,0,1);
fun = @(x,kk,rr,gg) kk+rr-mvncdf(gg,0,[1,x;x,1]);
A = zeros(length(i0),1);
ibnds = [-1+1e-4,1-1e-4];
parfor k = 1:length(i0)
    i = i0(k);
    j = j0(k);
    A(k) = fzero(@(x) S(k)+r(i)*r(j)-mvncdf(gam([i,j]),0,[1,x;x,1]),ibnds);
end
L = sparse(i0,j0,A,N,N);
L = L+L'+eye(N);
L0 = nearcorr(L);
save(fullfile(network.preNetwork,'covariance.mat'),'L0','ei');
clearvars L S i0 j0 KxxCell A

% Sample spikes using a discretized Gaussian
M = ceil(network.tmax/dt);
Lc = chol(nearcorr(L0),'lower');
xCell = cell(M,2);
parfor i = 1:M
    ids0 = find((Lc*randn(N,1))<gam);
    ts0 = (i*ones(size(ids0))+rand(size(ids0)))*dt;
    xCell(i,:) = {ids0,ts0};
end
ids = cat(1,xCell{:,1});
ts = cat(1,xCell{:,2});

% Map ids onto synapse ids
ids = synIDs(ids);
[~,I] = sort(synIDs);
ei = ei(I);

network_simulation_beluga.save_presynaptic_network(ids,ts,ei,N,network.spikingFile);

network.simulate();
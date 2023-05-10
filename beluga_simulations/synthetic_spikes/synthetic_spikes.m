% function synthetic_spikes
% folder = '/lustre04/scratch/nbrake/data/simulations/synthetic_spikes';
folder = 'E:\Research_Projects\004_Propofol\data\simulations\raw\test\synth_spikes';
addpath('/lustre04/scratch/nbrake/code/simulation_code');

% Initialize network
nPostNeurons = 2;
mTypes = 5;
network = network_simulation_beluga(folder);
network = network.initialize_postsynaptic_network(nPostNeurons,[1,2]);
network = network.setsynapsecount(100);
network.getsynapsecount
network = network.form_connections(0);

%%%%% GET SYNAPSE LOCATIONS ON SPHERE %%%%%
mData = network.getmData;
cons = csvread(fullfile(network.postNetwork,'connections.csv'));
P = arrayfun(@(i,j) mData{i}.pos(j,:),cons(:,1),cons(:,2),'UniformOutput',false);
P = cat(1,P{:});
[thet,phi] = cart2sph(P(:,1),P(:,2),P(:,3));
synIDs = cons(:,3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate covariance matrix from haversine distances
R = 0.2;
N = network.getsynapsecount;
dt = 1;
ei = rand(N,1)>network.parameters.eiFraction;
r = (network.parameters.eCellParams.firingRate*(1-ei) + network.parameters.iCellParams.firingRate*ei)*dt/1e3;
rVar = r.*(1-r);
X = [thet(:),phi(:)];
hav = @(x) (1-cos(x))/2;
hav_d = @(p1,x) hav(p1(1)-x(:,1))+(1-hav(x(:,1)-p1(1))-hav(p1(1)+x(:,1))).*hav(p1(2)-x(:,2));

KxxCell = cell(N,3);
parfor i = 1:N
    d = hav_d(X(i,:),X(i+1:end,:));
    idcs = find(d<0.25);
    KxxCell(i,:) = {idcs+i,i+0*zeros(size(idcs)),R*exp(-10*d(idcs)).*sqrt(rVar(i)*rVar(i+idcs))};
end
i0 = cat(1,KxxCell{:,1});
j0 = cat(1,KxxCell{:,2});
S = cat(1,KxxCell{:,3});

% Set up covariance matrix for DG
gam = icdf('normal',r,0,1);
fun = @(x,kk,rr,gg) kk+rr-mvncdf(gg,0,[1,x;x,1]);
A = zeros(length(i0),1);
parfor k = 1:length(i0)
    i = i0(k);
    j = j0(k);
    A(k) = fzero(@(x) S(k)+r(i)*r(j)-mvncdf(gam([i,j]),0,[1,x;x,1]),[-1+1e-4,1-1e-4]);
end
L = sparse(i0,j0,A,N,N);
L = L+L'+eye(N);
L0 = sparse(nearcorr(L)); % (DG has unit variance)
save(fullfile(network.preNetwork,'covariance.mat'),'L0','ei');


network.tmax = 10e3;
M = ceil(network.tmax/dt);
% Simulate spikes
X = zeros(N,M);
parfor i = 1:M
    waitbar(i/M)
    X(:,i) = mvnrnd(gam(:),L0,1)>0;
end
[ids,ts] = find(X);

% Map ids onto synapse ids
ids2 = synIDs(ids);
[~,I] = sort(synIDs);
ei2 = ei(I);
ts2 = (ts+rand(size(ts)))*dt;

network_simulation_beluga.save_presynaptic_network(ids2,ts2,ei2,N,network.spikingFile);

% network.simulate();
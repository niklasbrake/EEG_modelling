folder = 'E:\Research_Projects\004_Propofol\data\simulations\raw\test\synth_spikes';

% Initialize network
nPostNeurons = 2;
mTypes = 5;
network = network_simulation_beluga(folder);
network = network.initialize_postsynaptic_network(2,5*ones(nPostNeurons,1));
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
dt = 50;
ei = rand(N,1)<network.parameters.eiFraction;
r = (network.parameters.eCellParams.firingRate + network.parameters.iCellParams.firingRate*(1-ei))*dt/1e3;
rVar = r.*(1-r);
X = [thet(:),phi(:)];
hav = @(x) (1-cos(x))/2;
hav_d = @(p1,x) hav(p1(1)-x(:,1))+(1-hav(x(:,1)-p1(1))-hav(p1(1)+x(:,1))).*hav(p1(2)-x(:,2));
D = zeros(N,N);
for i = 1:N
    waitbar(i/N);
    D(:,i) = hav_d(X(i,:),X);
end
S = R*exp(-10*D).*(1-eye(N));
S = S.*(S>0.05*R);
aux = eye(N).*sqrt(rVar);
Kxx = aux*S*aux;

% Set up covariance matrix for DG
idcs = find(triu(Kxx)~=0);
gam = icdf('normal',r,0,1);
fun = @(x,kk,rr,gg) kk+rr-mvncdf(gg,0,[1,x;x,1]);
L = sparse(zeros(N,N));
for j = idcs(:)'
    [i0,j0] = ind2sub([N,N],i);
    L(i) = fzero(@(x) fun(x,Kxx(i),r(i0)*r(j0),gam([i0,j0])),[-1+1e-4,1-1e-4]);
end
L = L+L'+eye(N);
L0 = sparse(nearcorr(L)); % (DG has unit variance)
save(fullfile(network.preNetwork,'covariance.mat'),'L0');


network.tmax = 10e3;
M = ceil(network.tmax/dt);
% Simulate spikes
X = zeros(N,M);
for i = 1:M
    waitbar(i/M)
    X(:,i) = mvnrnd(gam(:),L0,1)>0;
end
[ids,ts] = find(X);
ids = synIDs(ids);
ts = (ts+rand(size(ts)))*dt;

network_simulation_beluga.save_presynaptic_network(ids,ts,ei,N,network.spikingFile);
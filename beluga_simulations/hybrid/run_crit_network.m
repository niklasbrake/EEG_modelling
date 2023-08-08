folder = '/lustre04/scratch/nbrake/data/simulations/EI_network';

N = 30e3;
tmax = 10e3;

% Set up connectivity matrix
pars.g = 1;
pars.p = 0.17;
pars.alpha = 0.2;

A = 5;
L = 0;
pars.w = L*2/(N*pars.p*(1-pars.alpha) - pars.g*N*pars.p*pars.alpha/A);

ei = zeros(N,1);
ei(floor(N*(1-pars.alpha)):end) = 1;

iE = find(ei==0); nE = length(iE);
iI = find(ei==1); nI = length(iI);
k = round(normrnd(N*pars.p,sqrt(N*pars.p*(1-pars.p)),N,1));
J = zeros(N,N);


waitbar(0);
for i = 1:N
    waitbar(i/N);
    iTarget = randsample(N,k(i));

    if(ei(i)) % Inhibitory
        J(iTarget,i) = -pars.g*pars.w*rand(k(i),1)/A^2;
    else % excitatory
        eiTarget = ei(iTarget);
        k2 = sum(eiTarget);
        J(iTarget(find(eiTarget)),i) = A*pars.w*rand(k2,1);
        J(iTarget(find(~eiTarget)),i) = pars.w*rand(k(i)-k2,1);
    end
end

J = sparse(J);

% Set up simulations
dt = 4;
pars.pext_I = 0.002;
% pars.pext_I = 2.5*(1-L)*dt*1e-3 + 0.06*L*dt*1e-3;
pars.pext_I = 2.5*(1-L)*dt*1e-3 + 0.01*L*dt*1e-3;
pars.pext_E = 0.002;
% pars.pext_E = 0.5*(1-L)*dt*1e-3 + 0.06*L*dt*1e-3;
pars.pext_E = 0.5*(1-L)*dt*1e-3 + 0.01*L*dt*1e-3;
T = 0:dt:tmax;
M = length(T);

X = sparse(zeros(N,1));

nI_Ex = poissrnd(pars.pext_I*nI,1);
exIdcs = randsample(nI,nI_Ex,1)+nE;
X(exIdcs) = 1;

nE_Ex = poissrnd(pars.pext_E*nE,1);
exIdcs = randsample(nE,nE_Ex,1);
X(exIdcs) = 1;

ids = zeros(length(T)*N,1);
ts = zeros(length(T)*N,1);
count = 0;
spikeCount0 = sum(X);

trueL = zeros(length(T),1);
for i = 1:length(T)
    waitbar(i/length(T));
    pIn = max(min(J*X,1),0);
    X = (rand(N,1)<=pIn);
    spikeCount1 = sum(X);

    trueL(i) = spikeCount1/spikeCount0;

    nI_Ex = poissrnd(pars.pext_I*nI,1);
    exIdcs = randsample(nI,nI_Ex,1)+nE;
    X(exIdcs) = 1;

    nE_Ex = poissrnd(pars.pext_E*nE,1);
    exIdcs = randsample(nE,nE_Ex,1);
    X(exIdcs) = 1;

    K = sum(X);
    ids(count+1:count+K) = find(X);
    ts(count+1:count+K)=T(i)+dt*rand(K,1);
    count = count+K;

    spikeCount0 = sum(X);
end
ids(count+1:end) = [];
ts(count+1:end) = [];

sum(ei(ids)==1)/10/sum(ei==1)
sum(ei(ids)==0)/10/sum(ei==0)
nanmean(trueL)


folder = 'E:\Research_Projects\004_Propofol\data\simulations\raw\hybrid';
% Initialize network
network = network_simulation_beluga(folder);

% Initialize post network
nPostNeurons = 1;
mType = 4;
network = network.initialize_postsynaptic_network(nPostNeurons,mType);

network.tmax = 10e3; % 2 seconds

network_simulation_beluga.save_presynaptic_network(ids,ts,ei,N,network.spikingFile)

i0 = randi(N); % Choose random neuron for postsyanptic 
i1 = find(J(i0,:)); % Get presynaptic connections

X = csvread('E:\Research_Projects\004_Propofol\data\simulations\raw\hybrid\presynaptic_network\spikeTimes.csv',0,1);

% extract and format presyanptic network to the chosen postsynaptic neuron
ids = X(i1,:)*0+i1(:);
ts = X(i1,:);
ids = ids(:);
ts = ts(:);
[i,j,ts] = find(ts);
ids = ids(i);
[ids,j0] = findgroups(ids);
EI = ei(j0);



network.spikingFile = fullfile(network.preNetwork,'preSpikes.csv');
network_simulation_beluga.save_presynaptic_network(ids,ts,EI,length(i1),network.spikingFile)

network = network.setsynapsecount(length(EI));
network.form_connections(0);

network.simulate();

network.savePath = fullfile(network.outputPath,'simulation_tau');
network.parameters.iSynParams.tau2 = 20;
network.simulate();
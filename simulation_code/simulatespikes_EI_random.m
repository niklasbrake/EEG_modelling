
folder = '/lustre04/scratch/nbrake/data/simulations/EI_network';

N = 50e3;


% Set up connectivity matrix
pars.g = 1;
pars.p = 0.2;
pars.alpha = 0.2;

lTarget = 0.98;
pars.w = lTarget*2/(N*pars.p*(1-pars.alpha) - pars.g*N*pars.p*pars.alpha);

lam = pars.w/2*(N*pars.p*(1-pars.alpha) - pars.g*N*pars.p*pars.alpha);
sE = (pars.p/3-pars.p^2/4)*pars.w^2;
sI = (pars.p/3-pars.p^2/4)*(pars.g*pars.w)^2;
R = sqrt(N*(1-pars.alpha)*sE^2+pars.alpha*sI^2);
lMax = max(lam,R)

ei = zeros(N,1);
ei(floor(N*(1-pars.alpha)):end) = 1;

iE = find(ei==0);
iI = find(ei==1);
J = zeros(N,N);
k = binornd(N,pars.p,N,1);
waitbar(0);
for i = 1:N
    waitbar(i/N);
    iTarget = randsample(N,k(i));
    if(ei(i)) % Inhibitory
        J(iTarget,i) = -pars.g*pars.w*rand(k(i),1);
    else % excitatory
        J(iTarget,i) = pars.w*rand(k(i),1);
    end
end

J(iI,iE) = 5*J(iI,iE);
J(:,iI) = J(:,iI)/5;

J = sparse(J);

% Set up simulations
pars.pext = 0.005*10;
dt = 1;
tmax = 5e3;
T = 0:dt:tmax;
M = length(T);


X = sparse(zeros(N,1));
nEx = poissrnd(pars.pext,1);
exIdcs = randsample(N,nEx,1);
X(exIdcs) = 1;

ids = zeros(length(T)*N,1);
ts = zeros(length(T)*N,1);
count = 0;
for i = 1:length(T)
    waitbar(i/length(T));
    pIn = max(min(J*X,1),0);
    X = (rand(N,1)<=pIn);

    nEx = poissrnd(pars.pext,1);
    exIdcs = randsample(N,nEx,1);
    X(exIdcs) = 1;
    
    K = sum(X);
    ids(count+1:count+K) = find(X);
    ts(count+1:count+K)=i;
    count = count+K;
end
ids(count+1:end) = [];
ts(count+1:end) = [];


save(fullfile(folder,'EI_network.mat','ids','ts','ei','pars'))

return;
raster(ids,ts);


lamI = sum(ei(ids)==1)/sum(ei==1)/tmax*1e3
lamE = sum(ei(ids)==0)/sum(ei==0)/tmax*1e3


K = 100;
M = zeros(K,1);
X0 = sparse(zeros(N,1));
waitbar(0);
for i = 1:K
    waitbar(i/K);
    K0 = randi(N);
    for k = 1
        X1 = X0;
        X1(randsample(N,K0)) = 1;
        pIn = max(min(J*X1,1),0);
        M(i) = M(i) + sum((rand(N,1)<=pIn))/K0;
    end
end
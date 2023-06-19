% function [lamI,lamE,ei] = run_EInetwork(lTarget,A)

folder = '/lustre04/scratch/nbrake/data/simulations/EI_network';

N = 10e3;

% Set up connectivity matrix
pars.g = 1;
pars.p = 0.2;
pars.alpha = 0.2;

A = 5;
% pars.w = lTarget*2/(N*pars.p*(1-pars.alpha) - pars.g*N*pars.p*pars.alpha);
lTarget = 0.98;
pars.w = lTarget*2/(N*pars.p*(1-pars.alpha) - pars.g*N*pars.p*pars.alpha);
% pars.w = lTarget*2/(N*pars.p*(1-pars.alpha) - pars.g*N*pars.p*pars.alpha/A);

lam = pars.w/2*(N*pars.p*(1-pars.alpha) - pars.g*N*pars.p*pars.alpha);
% lam = pars.w/2*(N*pars.p*(1-pars.alpha)*(1-(1-A)*pars.alpha) - pars.g*N*pars.p*pars.alpha/A);

sE = (pars.p/3-pars.p^2/4)*pars.w^2;
sI = (pars.p/3-pars.p^2/4)*(pars.g*pars.w)^2;
R = sqrt(N*(1-pars.alpha)*sE^2+pars.alpha*sI^2);
lMax = max(lam,R);

ei = zeros(N,1);
ei(floor(N*(1-pars.alpha)):end) = 1;

iE = find(ei==0); nE = length(iE);
iI = find(ei==1); nI = length(iI);
k = round(normrnd(N*pars.p,sqrt(N*pars.p*(1-pars.p)),N,1));
J = zeros(N,N);

for i = 1:N
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
pars.pext_I = 0.002;%2.5*(1-lTarget)*dt*1e-3 + 0.06*lTarget*dt*1e-3;
pars.pext_E = 0.002;%0.5*(1-lTarget)*dt*1e-3 + 0.06*lTarget*dt*1e-3;
tmax = 5e3;
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
for i = 1:length(T)
    waitbar(i/length(T));
    pIn = max(min(J*X,1),0);
    X = (rand(N,1)<=pIn);

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
end
ids(count+1:end) = [];
ts(count+1:end) = [];


% save(fullfile(folder,'EI_network.mat','ids','ts','ei','pars'))

lamI = sum(ei(ids)==1)/sum(ei==1)/tmax*1e3
lamE = sum(ei(ids)==0)/sum(ei==0)/tmax*1e3
ei = lamI/lamE
raster(ids,ts);
return;




% K = 100;
% M = zeros(K,1);
% X0 = sparse(zeros(N,1));
% waitbar(0);
% for i = 1:K
%     waitbar(i/K);
%     K0 = randi(N);
%     for k = 1
%         X1 = X0;
%         X1(randsample(N,K0)) = 1;
%         pIn = max(min(J*X1,1),0);
%         M(i) = M(i) + sum((rand(N,1)<=pIn))/K0;
%     end
% end
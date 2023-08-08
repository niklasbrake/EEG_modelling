function [ids,ts,ei,P,trueL] = run_EInetwork(N,lTarget)
tmax = 1e4;

folder = '/lustre04/scratch/nbrake/data/simulations/EI_network';

% N = 30e3;
P = rand(N,2);
P(:,1) = acos(1-2*P(:,1))-pi/2;
P(:,2) = 2*P(:,2)*pi;

% Set up connectivity matrix
pars.g = 1;
pars.p = 0.2*0.1;
pars.alpha = 0.2;

A = 5;
% pars.w = lTarget*2/(N*pars.p*(1-pars.alpha) - pars.g*N*pars.p*pars.alpha);
% lTarget = 0.99;
% pars.w = lTarget*2/(N*pars.p*(1-pars.alpha) - pars.g*N*pars.p*pars.alpha);
pars.w = 5*lTarget*2/(N*pars.p*(1-pars.alpha) - pars.g*N*pars.p*pars.alpha/A);

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


hav = @(x) (1-cos(x))/2;
h_dist = @(p1,x) hav(p1(1)-x(:,1))+(1-hav(x(:,1)-p1(1))-hav(p1(1)+x(:,1))).*hav(p1(2)-x(:,2));

waitbar(0);
for i = 1:N
    waitbar(i/N);
    p0 = P(i,:);
    % d = max(exp(-sqrt(sum((p0-P).^2,2))/0.005),1e-9);
    % d = max(exp(-sqrt(sum((p0-P).^2,2))/0.006),1e-9);
    d = max(exp(-h_dist(p0,P)/8e-4),1e-9);
    idcs = [1:i-1,i+1:N];
    C = cumsum(d(idcs))/sum(d(idcs));
    iTarget = interp1(C,idcs,rand(k(i),1),'next','extrap');

    % iTarget = randsample(N,k(i));
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
pars.pext_I = 2.5*(1-lTarget)*dt*1e-3 + 0.06*lTarget*dt*1e-3;
pars.pext_E = 0.002;
pars.pext_E = 0.5*(1-lTarget)*dt*1e-3 + 0.06*lTarget*dt*1e-3;
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


nanmean(trueL)

% save(fullfile(folder,'EI_network.mat','ids','ts','ei','pars'))

% lamI = sum(ei(ids)==1)/sum(ei==1)/tmax*1e3
% lamE = sum(ei(ids)==0)/sum(ei==0)/tmax*1e3
% ei = lamI/lamE

% if(lamI<20)
    % [d,I] = sort(sqrt(sum((P-[0.5,0.5]).^2,2)));
    [d,I] = sort(P(:,2));
    [~,J] = sort(I);
    ids = J(ids);
    P = P(I,:);
%     raster(J(ids),ts);
% end
% raster(ids,ts);

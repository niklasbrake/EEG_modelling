function [ids,ts,ei,C,m0] = simulatespikes_det(N,branchNo,tmax,C,nNeigh)

dt = 16e-3;
if(nargin<5)
    nNeigh = 4;
end
if(nargin<4 || isempty(C))
    pmax = branchNo/nNeigh;
    idcs = randi(N,N*nNeigh,1); % Random start nodes
    C = zeros(N*nNeigh,4);
    for i = 1:length(idcs)
        C(i,1) = idcs(i);
        j0 = randi(N);
        while(j0==idcs(i))
            j0 = randi(N);
        end
        C(i,2) = j0; % Random end nodes
        C(i,3) = rand*pmax*2; % Random probability of propogation
        C(i,4) = dt*2*rand; % Random transmission speed
    end
    % C = randi(N,N,nNeigh);
    % C = [repmat((1:N)',[1,nNeigh]),C];
    % C = reshape(C,[],2);
end

% Look up table for node indices
for i = 1:N
    CLUT{i} = find(C(:,1)==i);
end

% tmax = 10;
t = 0:dt:tmax;
tN = length(t);
X = zeros(N,tN);

% Initialize
ids = -ones(ceil(2*N*tmax),1);
ts = -ones(ceil(2*N*tmax),1);
lamE = 1;
lamI = 1;

% Random external input
nTrans = poissrnd(lamE*dt*N);
ids(1:nTrans) = randperm(N,nTrans);
post = ids(1:nTrans);
k = nTrans;
ts(1:nTrans) = t(1)*ones(nTrans,1);
count = nTrans;
m0 = zeros(tN-1,1);
m0(1) = branchNo;

exN = poissrnd(lamE*dt*N*(1-branchNo),tN,1);
lamE = 1+0*0.6*sin(2*pi*t*10);

for i = 2:tN-1
    % Get spiking cells at previous time point
    % preI = unique(ids(count-nTrans+1:count));
    preI = ids(count-nTrans+1:count);
    % length(preI)
    jj = cat(1,CLUT{preI});

    % For each spiking cell, find neighbours...
    postI = C(jj,2);
    % ... and start flipping coins
    iTrans = find(rand(length(jj),1)<C(jj,3));
    postI = postI(iTrans);
    nTrans = length(postI);

    % Correct for redundent spike propogations
    % postI = unique(postI(iTrans));
    if(length(preI)>0)
        m0 = m0 + (min(nTrans/length(preI),branchNo)-m0)/i;
        m0(i) = nTrans/length(preI);
    end
    % exN(i) = poissrnd(dt*N*(1-min(m0(i),branchNo)));
    exN(i) = poissrnd(lamE(i)*dt*N*(1-branchNo));


    % Add propogated spikes
    ids(count+1:count+nTrans) = postI;
    ts(count+1:count+nTrans) = t(i)+ dt*rand(nTrans,1);;
    count = count+nTrans;

    % Add external noise
    ids(count+1:count+exN(i)) = randperm(N,exN(i));
    ts(count+1:count+exN(i)) = t(i)+ dt*rand(exN(i),1);
    count = count+exN(i);
    nTrans = nTrans + exN(i);
end
ts(count:end) = [];
ids(count:end) = [];


B = (lamI/lamE(1)-1)*histcounts(ts,t)/N/dt;
B = B*0+2;

ts = repmat(ts,[5,1]);
ids = repmat(ids,[5,1]);
NI = ceil(0.2*N);
for i = 1:tN-1
    nTrans = poissrnd(B(i)*dt*NI,1,1);
    ids(count+1:count+nTrans) = randperm(NI,nTrans)+N-NI;
    ts(count+1:count+nTrans) = t(i)+ dt*rand(nTrans,1);
    count = count+nTrans;
end
ts(count+1:end) = [];
ids(count+1:end) = [];
ei = [zeros(1,N-NI),ones(1,NI)];


ts = ts*1e3;
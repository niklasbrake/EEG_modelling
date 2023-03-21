function [ids,ts,ei,C] = simulatespikes_det(N,branchNo,tmax,C)

% branchNo = 0.98; 
% N = 20000; 
nNeigh = 4;
if(nargin<4)
    C = randi(N,N,nNeigh);
    C = [repmat((1:N)',[1,4]),C];
    C = reshape(C,[],2);
else
    if(size(C,1)~=nNeigh*N)
        error('Connection matrix is not the right size');
    end
end

% tmax = 10;
dt = 4e-3;
t = 0:dt:tmax;
tN = length(t);
X = zeros(N,tN);

% Initialize
ids = -ones(2*N*tmax,1);
ts = -ones(2*N*tmax,1);

% Random external input
h = dt*N*(1-branchNo); 
exN = poissrnd(h,tN,1);
nTrans = poissrnd(dt*N);
ids(1:nTrans) = randperm(N,nTrans);
ts(1:nTrans) = t(1)*ones(nTrans,1);
count = nTrans;

for i = 2:tN-1
    preI = ids(ts==t(i-1)); % Get spikes at time i
    % For each spike, stimulate neighbours
    for j = 1:length(preI)
        jj = find(C(:,1)==preI(j));
        postI = C(jj,2);
        nTrans = binornd(nNeigh,branchNo/nNeigh);
        ids(count+1:count+nTrans) = postI(randperm(nNeigh,nTrans));
        ts(count+1:count+nTrans) = t(i)*ones(nTrans,1);
        count = count+nTrans;
        % X(postI(randperm(nNeigh,nTrans)),i+1) = 1;
    end
    % Add external noise
    ids(count+1:count+exN(i)) = randperm(N,exN(i));
    ts(count+1:count+exN(i)) = t(i)*ones(exN(i),1);
    count = count+exN(i);
end
ts(count+1:end) = [];
ids(count+1:end) = [];


B = 4*histcounts(ts,t)/N/dt;

ts = repmat(ts,[5,1]);
ids = repmat(ids,[5,1]);
NI = ceil(0.2*N);
for i = 1:tN-1
    nTrans = poissrnd(B(i)*dt*NI,1,1);
    ids(count+1:count+nTrans) = randperm(NI,nTrans)+N-NI;
    ts(count+1:count+nTrans) = t(i)*ones(nTrans,1);
    count = count+nTrans;
end
ts(count+1:end) = [];
ids(count+1:end) = [];
ei = [zeros(1,N-NI),ones(1,NI)];


ts = (ts+dt*rand(size(ts)))*1e3;
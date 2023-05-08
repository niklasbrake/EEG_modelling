N = 100;
dt = 50; % ms

% 15% inhibitory
ei = rand(N,1)>0.2;

% e: 1Hz and i: 5Hz
r = (ei + 5*(1-ei))*dt/1e3;
% variance
rVar = r.*(1-r);

p = unique(r);
for i = 1:length(p)
    for j = i:length(p)
        Rlo(i,j) = max(-p(i)*p(j),-(1-p(i))*(1-p(j)));
        Rhi(i,j) = min((1-p(i))*p(j),(1-p(j))*p(i));
    end
end

R = 0.2;

% Get correlations based on haversine metrix
hav = @(x) (1-cos(x))/2;
hav_d = @(p1,x) hav(p1(1)-x(:,1))+(1-hav(x(:,1)-p1(1))-hav(p1(1)+x(:,1))).*hav(p1(2)-x(:,2));
thet = rand(N,1)*2*pi;
phi = asin(2*rand(N,1)-1);
X = [thet(:),phi(:)];
[x,y,z] = sph2cart(thet,phi,phi*0+1);
D = zeros(N,N);
for i = 1:N
    waitbar(i/N);
    D(:,i) = hav_d(X(i,:),X);
end
S = R*exp(-10*D).*(1-eye(N));
S = S.*(S>0.05*R);

% Convert correlation matrix to covariance matrix
aux = eye(N).*sqrt(rVar);
Kxx = aux*S*aux;

% Get nonzero indices for speed
idcs = find(triu(S)~=0);

% Mean and variance of nd Gaussian
gam = icdf('normal',r,0,1);
% fun = @(x,K,rr,gg) Kxx(i,j)+r(i)*r(j)-mvncdf(gam([i,j]),0,[1,x;x,1]);
fun = @(x,kk,rr,gg) kk+rr-mvncdf(gg,0,[1,x;x,1]);

L = sparse(zeros(N,N));
count = 1; K = 1000;
h = waitbar(0);
tic;
for j = 1:length(idcs)
    if(j > floor(count*length(idcs)/K))
        count = count+1;
        t = toc/j*(length(idcs)-j);
        waitbar(count/K,h,datestr(seconds(t),'HH:MM:SS'));
    end
    i = idcs(j);
    [i0,j0] = ind2sub([N,N],i);
    L(i)= fzero(@(x) fun(x,Kxx(i),r(i0)*r(j0),gam([i0,j0])),[-1+1e-4,1-1e-4]);
end
delete(h)
L = L+L'+eye(N);

L0 = nearcorr(L); % (DG has unit variance)


tmax = 20e3;
M = ceil(tmax/dt);
% Simulate spikes
X = zeros(N,M);
for i = 1:M
    waitbar(i/M)
    X(:,i) = mvnrnd(gam(:),L0,1)>0;
end
[ids,ts] = find(X);
ts = (ts+rand(size(ts)))*dt;
raster(ids,ts);

figureNB;
    h = scatter3(x,y,z,5*ones(N,1),zeros(N,1),'filled');
for i = 1:10:tmax
    idcs0 = find(and(ts>=i,ts<i+10));
    idcs = ids(idcs0);
    h.CData(idcs) = 1;
    h.SizeData(idcs) = 20;
    set(gca,'CLim',[0,1]);
    drawnow;
    pause(0.1);
    h.CData(idcs) = 0;
    h.SizeData(idcs) = 5;
end


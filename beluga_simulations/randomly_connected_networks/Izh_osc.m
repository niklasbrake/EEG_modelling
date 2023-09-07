
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Integration Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dt = 0.05;
tmax = 3e3;
N = tmax/dt;
nRS = 1e4;
nFS = ceil(nRS*0.15);
V = zeros(nRS+nFS,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Model parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CFS = 20;
kFS = 1;
vrFS = -55;
vtFS = -40;
vpFS = 25;
aFS = max(0.02*(1+3*rand(nFS,1)),0);
bFS = 8;
cFS = -55;
dFS = 200;
vFS = vrFS + (vtFS-vrFS)*rand(nFS,1); 
uFS = -vrFS/bFS.*ones(nFS,1); 

CRS = 100;
kRS = 3;
vrRS = -60;
vtRS = -50;
vpRS = 50;
aRS = max(0.002*(1+3*rand(nRS,1)),0);
bRS = 5;
cRS = -60;
dRS = 400;
vRS = vrRS + (vtRS-vrRS)*rand(nRS,1); 
uRS = -vrRS/bRS.*ones(nRS,1); 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Synaptic parameters 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gaFS = max(0,randn(nFS,1)/3);    ggFS = max(0,randn(nFS,1));
gaRS = zeros(nRS,1);    ggRS = max(0,randn(nRS,1));

IRS2 = max(100*(1+5*randn(nRS,1)),0);
IFS2 = 100*(1+0.15*randn(nFS,1));

tauG = 7;
tauE = 2;


RSxFS = randi(nRS,nFS,ceil(nRS*0.2));
RSxRS = randi(nRS,nRS,ceil(nRS*0.1));
FSxRS = randi(nFS,nRS,ceil(nFS*0.1));
FSxFS = randi(nFS,nFS,ceil(nFS*0.2));

ts = nan(N*nRS,1);
ids = nan(N*nRS,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Integrate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h= waitbar(0);
tic;
count = 0;
for i = 1:N
    if(~mod(i,100))
        waitbar(i/N,h,['Estimate time remaining: ' char(duration(0,0,toc/i*(N-i)))]);
    end

    % Update RS membrane voltage
    IRS = gaRS.*(vRS-0) + ggRS.*(vRS + 70) - IRS2; % input current
    % Update DEs
    vRS = vRS + 0.5*dt/CRS * (kRS*(vRS-vrRS).*(vRS-vtRS) - uRS - IRS);
    uRS = uRS + 0.5*dt*aRS.*(bRS*(vRS-vrRS)-uRS);
    % Detect spikes and perform resetting
    spRS = vRS>vpRS;
    vRS(spRS) = cRS;
    uRS(spRS) = uRS(spRS)+dRS;
    % Decay of synaptic input currents
    gaRS = gaRS - dt*gaRS/tauE; ggRS = ggRS - dt*ggRS/tauG;

    % Update FS membrane voltage
    IFS = gaFS.*(vFS-0) - IFS2; % input current
    % Update DEs
    vFS = vFS + dt/CFS * (kFS*(vFS-vrFS).*(vFS-vtFS) - uFS - IFS);
    uFS = uFS + dt*aFS.*(bFS*(vFS-vrFS)-uFS);
    % Detect spikes and perform resetting
    spFS = vFS>vpFS;
    vFS(spFS) = cFS;
    uFS(spFS) = uFS(spFS)+dFS;
    % Decay of synaptic input currents
    gaFS = gaFS - dt*gaFS/tauE; ggFS = ggFS - dt*ggFS/tauG;

    gaRS = gaRS + 75*mean(spRS(RSxRS),2);
    ggRS = ggRS + 25*mean(spFS(FSxRS),2);
    gaFS = gaFS + 50*mean(spRS(RSxFS),2);
    ggFS = ggFS + 10*mean(spFS(FSxFS),2);

    spikeIDs = [find(spRS); ...
            nRS+find(spFS)];
    k = length(spikeIDs);
    ids(count+1:count+k) = spikeIDs;
    ts(count+1:count+k) = dt*i*ones(size(spikeIDs));
    count = count+k;
end
ts(count:end) = [];
ids(count:end) = [];
raster(ids,ts);
%{



% Y = Y(X>2e3/dt);
% X = X(X>2e3/dt);
% T = tmax-2e3;

N = histcounts(X,'BinEdges',(0:10:T)/dt);

for i = 1:max(Y)
    n = histcounts(X(Y==i),'BinEdges',(0:10:T)/dt);
    C(i) = corr(n(:),N(:));
end

figure('color','w');
subplot(2,2,1)
    plot(X*dt/1e3,Y,'.k','MarkerSize',1);
    xlabel('Time (s)');
    ylim([1,size(V,1)]);
    gcaformat;
subplot(2,2,2);
    for i=  1:max(Y)
        F(i) = sum(Y==i)/T*1e3;
        line([0,F(i)],[i,i],'color','k','LineWidth',0.1);
    end
    ylim([1,size(V,1)]);
    xlabel('Firing rate (Hz)');
    gcaformat;
subplot(2,2,3);
    plot((5:10:T-5)/1e3,N)
    xlabel('Time (s)');
    ylabel('Count');
    gcaformat;
subplot(2,2,4);
     boxplotNB(1,C,'b',5)
     gcaformat;
     ylabel('Correlation');
%}

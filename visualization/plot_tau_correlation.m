
load('E:\Research_Projects\004_Propofol\data\experiments\scalp_EEG\raw\timeInformation.mat')
t0 = timeInfo.infusion_onset-timeInfo.object_drop;
data = load('E:\Research_Projects\004_Propofol\data\experiments\scalp_EEG\analyzed\Cz_multitaper_mean.mat');
pre = [];
for i = 1:length(t0)
    pre(:,i) = nanmedian(data.psd(:,data.time<t0(i),i),2);
    totalPower_Data(i) = sum(pre(:,i)*mean(diff(data.freq)));
end
freq = data.freq;

load('unitary_spectra_different_tau_decay.mat','tau','tau2');   
load('unitary_spectra_passive.mat');

figureNB(3.2,4.8);
[params,synFun] = synDetrend(f(f<100),mean(P(f<100,:),2)./mean(P(1,:)),0,'lorenz',[15e-3,1e-3,0,-1]);
subplot(2,1,1);
    plotwitherror(f,P,'Q','LineWidth',1,'color','k'); hold on;
    % plot(f,mean(P,2),'color','k','LineWidth',1); hold on;
    ylim([10^-16.5,10^-14.5]);
    yticks([1e-16,1e-15])
    set(gca,'yscale','log');
    set(gca,'xscale','log')
    xticks([1,10,100]);
    xticklabels([1,10,100]);
    xlabel('Frequency (Hz)');
    ylabel(['PSD (' char(956) 'V^2/Hz)'])
    xlim([1,100])
    plot(f,mean(P(1,:))*10.^synFun(f,[params(1:3),-Inf]),'--r');
    plot(f,mean(P(1,:))*10.^synFun(f,[params(1:2),-Inf,params(4)]),'--r');
    % plot(f,mean(P(1,:))*10.^synFun(f,params),'r');
    gcaformat;
    text(1.3,2e-15,'\tau_1','FontSize',6,'color','r')
    text(1.3,9e-17,'\tau_2','FontSize',6,'color','r')
subplot(2,1,2);
    errorbar(tau,mean(tau2,2),std(tau2,[],2),'.k','MarkerSize',10,'LineWidth',1)
    line([0,50],[0,50],'color','k')
    ylabel('EEG spectrum \tau_1 (ms)')
    xlim([7.5,37.5])
    ylim([5,45])
    ylabel('\tau_1 (ms)')
    xlabel('GABA_AR \tau_{decay} (ms)')
    gcaformat;

subplot(1,3,3)
    plotwitherror(freq,pre,'Q','LineWidth',1,'color','b'); hold on;
    plotwitherror(f,P*25e9,'Q','LineWidth',1,'color','k'); hold on;
    set(gca,'xscale','log')
    set(gca,'yscale','log')
    ylim([1e-7,1e2]); yticks([1e-6,1e-2,1e2])
    xlim([1,100]);
    xticks([1,10,100]);
    xticklabels([1,10,100]);
    xlabel('Frequency (Hz)')
    ylabel(['PSD (' char(956) 'V^2/Hz)'])
    gcaformat



for j = 1:600
    temp = resample(network2.results.dipoles(2:end,:,j),2e3,16e3);
    if(j==1)
        dp = zeros(size(temp,1),3,600);
    end
    dp(:,:,j) = temp;
end
Q = gpuArray(dp);
idcs = sa.cortex2K.in_from_cortex75K;
L1 = squeeze(sa.cortex75K.V_fem(49,:,:)); % czIDX = 49


h = waitbar(0);
V = gpuArray(zeros(length(idcs),600));
for j = 1:length(idcs)
    waitbar(j/length(idcs));
    L0 = L1(idcs(j),:);
    vz = sa.cortex75K.normals(idcs(j),:);
    [vx,vy] = getOrthBasis(vz);
    
    q = (Q(:,1,:).*vx+Q(:,2,:).*vy+Q(:,3,:).*vz)*1e-12; % nAum -> mAm

    W = squeeze(sum(L0.*q,2)*1e6); % V -> uV
    V(j,:) = var(W);
end
V = gather(V);
delete(h)




for j = 1:length(network2.neurons)
    nrn = network2.neurons(j).name;
    fid = fopen(['C:\Users\brake\Documents\temp\passive2\postsynaptic_network\connections\' nrn '.csv']);
    n = 0;
    tline = fgetl(fid);
    while ischar(tline)
      tline = fgetl(fid);
      n = n+1;
    end
    fclose(fid);
    N(j) = n;
end

clrs0 = clrsPT.lines(11);
% clrs = clrs0(G,:);
clrs0 = [1,0,0;0,0,1];
ei = cellfun(@(x)~isempty(strfind(x,'I_')),mTypes)+1;

figureNB(8.9,15);
ax = axes('Units','centimeters');
ax.Position = [6.64,4.66,1.62,3.28];
    scatter(N,mean(V),5,clrs0(ei,:),'filled'); hold on;
    scatter(mean(N),mean(mean(V)),15,[0,0,0],'filled','MarkerFaceColor','w','MarkerEdgeColor','k','LineWidth',0.5);
    set(gca,'yscale','log')
    xlabel('Synapse count')
    ylabel('Avg. EEG power (uV^2)')
    gcaformat
    xax = get(gca,'xaxis')
    xax.Exponent = 0
    xlim([0,25000])
    ylim([1e-15,1e-13]);
    set(gca,'yscale','linear');

    x0 = mean(N);
    y0 = mean(mean(V));
    xl = get(gca,'xlim');
    yl = get(gca,'ylim');
    ax.Units = 'normalized';
    xx =  @(x0) ax.Position(1)+ax.Position(3)*(x0-xl(1))/diff(xl);
    yy =  @(x0) ax.Position(2)+ax.Position(4)*(y0-yl(1))/diff(yl);
    A = annotation('arrow',[xx(x0+5e3+1e3),xx(x0+1e3)],[yy(y0),yy(y0)],'HeadWidth',5,'HeadLength',5);

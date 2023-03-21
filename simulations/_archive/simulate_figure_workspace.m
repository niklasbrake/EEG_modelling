load('E:\Research_Projects\004_Propofol\Modelling\head_models\sa_nyhead.mat')
X.vertices = sa.cortex75K.vc;
X.faces= sa.cortex75K.tri;

%%%%%%%%%%%% NYHM %%%%%%%%%%%%%%%%%

folder = 'E:\Research_Projects\004_Propofol\Modelling\neuron_simulations\data\simulations\cortical_column\baseline\20220518'
% folder = 'E:\Research_Projects\004_Propofol\Modelling\neuron_simulations\data\simulations\cortical_column\propofol';
F = dir(folder); F = F(3:end);
Q = zeros(32001,3,length(F));
n = length(F);
for i = 1:n
    data = load(fullfile(folder,F(i).name));
    Q(:,:,i) = data.Q';
end
Q = [Q(:,1,:),-Q(:,3,:),Q(:,2,:)];

% Correlated dipoles
idcs = sa.cortex2K.in_from_cortex75K;
N2 = length(idcs);
m = 1000;

p = zeros(16385,m);
h = waitbar(0); tic;
for i = 1:m
    waitbar(i/m,h,['Estimated time remaining: ' char(duration(0,0,toc/i*(m-i)))]);
    eeg = zeros(32001,1);
    idcs = randi(N2,n,1);
    for j = 1:n
        eeg = eeg + getEEG(Q(:,:,j),sa,idcs(j)); % nA.um -> mA.m and V -> uV
    end
    [p(:,i),f] = pmtm(eeg/n,2,[],1e3*16);
end

N = 25e10;
load('E:\Research_Projects\004_Propofol\Experiments\scalp_EEG\analyzed_data\pre_post_mean.mat')
freq = 0.5:0.5:512;
for i = 1:14
    pre(:,i) = exp(preP{i}(:,2));
end
pre = mean(pre,2);
load('E:\Research_Projects\004_Propofol\Modelling\neuron_simulations\data\summary_data\expected_PSD_baseline.mat')


SIG0 = mean(sum(p*mean(diff(f))));
p0 = mean(p,2)/SIG0;

SIG = @(rho) N*SIG0 + N*(N-1)*rho*SIG0;
p_lower = p0*SIG(0);
p_upper = p0*SIG(1);
p_est = p0*SIG(1e-6);

fig = figure('color','w');
    plot(f,p_est,'linewidth',1,'color','b'); hold on;
    plot(f,p_lower,'linewidth',1,'color','r'); hold on;
    plot(f,p_upper,'linewidth',1,'color','r'); hold on;
    % plot(dataPSD.data.freq,nanmean(nanmean(dataPSD.data.psd(:,1000:2000,:),2),3),'LineWidth',1,'color','k')
    plot(freq,pre,'LineWidth',1,'color','k')
    set(gca,'xscale','log');
    set(gca,'yscale','log');
    xlim([0.5,300]);
    gcaformat;
    xticks([1,10,100]);
    xticklabels([1,10,100]);
    xlabel('Frequency (Hz)');
    ylabel(['Power (' char(956) 'V^2/Hz)'])


p1 = interp1(f,p0,freq); p1 = p1(freq<300);

while(min(pre./p1)>=1)
    p1 = p1*1.1;
end


N = 25e10;
SIG = @(rho) N*SIG0 + N*(N-1)*rho*SIG0;
base = load('E:\Research_Projects\004_Propofol\Modelling\neuron_simulations\data\summary_data\expected_PSD_baseline.mat')
SIG0 = mean(sum(base.p*mean(diff(f))));
p0 = mean(base.p,2)/SIG0;

prop = load('E:\Research_Projects\004_Propofol\Modelling\neuron_simulations\data\summary_data\expected_PSD_propofol.mat')
SIG1 = mean(sum(prop.p*mean(diff(f))));
p1 = mean(prop.p,2)/SIG1;


fig = figure('color','w');
    plot(base.f,p0*SIG(8e-6),'linewidth',1,'color','k'); hold on;
    plot(prop.f,p1*SIG(8e-6),'linewidth',1,'color','r'); 
    set(gca,'xscale','log');
    set(gca,'yscale','log');
    xlim([0.5,300]);
    gcaformat;
    xticks([1,10,100]);
    xticklabels([1,10,100]);
    xlabel('Frequency (Hz)');
    ylabel(['Power (' char(956) 'V^2/Hz)'])


p0 = mean(base.p,2);
PSD = @(C) N*p0+N*p0.*2*(N-1)*C;


C_blank = ones(length(base.f),1);
C0 = C_blank*0;
C1 = C_blank*1;
C2 = C_blank*3e-6;
C3 = C_blank.*(exp(-(f-10).^2./2)+3e-6)./(1+3e-6);

fig = figure('color','w');
    plot(f,PSD(C0),'linewidth',1,'color','k'); hold on;
    plot(f,PSD(C1),'linewidth',1,'color','r'); 
    plot(f,PSD(C2),'linewidth',1,'color','b'); 
    plot(f,PSD(C3),'linewidth',1,'color','b'); 
    set(gca,'xscale','log');
    set(gca,'yscale','log');
    xlim([0.5,300]);
    gcaformat;
    xticks([1,10,100]);
    xticklabels([1,10,100]);
    xlabel('Frequency (Hz)');
    ylabel(['Power (' char(956) 'V^2/Hz)'])


fig = figure('color','w');
plot(f,(PSD(C3)-PSD(C0))./(PSD(C0)*2*(N-1)),'linewidth',1,'color','b'); 




fig = figure('color','w');
    plot(f,p0,'linewidth',1,'color','r'); hold on;
    set(gca,'xscale','log');
    set(gca,'yscale','log');
    xlim([0.5,300]);
    gcaformat;
    xticks([1,10,100]);
    xticklabels([1,10,100]);
    xlabel('Frequency (Hz)');
    ylabel(['Power (' char(956) 'V^2/Hz)'])




load('E:\Research_Projects\004_Propofol\Experiments\scalp_EEG\raw_data\time_series_all_channels.mat')
m = 1000;
eeg = zeros(32001,1);
idcs = randi(N2,size(Q,3),1);
for j = 1:n
    eeg = eeg + getEEG(Q(:,:,j),sa,idcs(j)); % nA.um -> mA.m and V -> uV
end
eeg0 = TimeDomainAligned(and(Time>-223,Time<-221),2,3);

fig=  figureNB(7,3);
axes('Position',[0,0.5,1,0.5]);
    plot(linspace(0,2,length(eeg0)),eeg0,'color','k')
    xlim([0,2])
    ylim([-5,5]*std(eeg0));
    axis off;
    line([0.5,0.5],[-33,-23]+10,'LineWidth',1,'color','k')
    text(0.475,-28+10,['1 ' char(956) 'V'],'FontSize',6,'HorizontalAlignment','right')
axes('Position',[0,0,1,0.5]);
    plot((1:length(eeg))/16*1e-3,eeg,'color','r')
    xlim([0,2])
    axis off;
    ylim([-7,7]*std(eeg));
    line([0.5,0.5],[-5,-4]*1e-6+1e-6,'LineWidth',1,'color','k')
    text(0.475,-4.5*1e-6+1e-6,'1 pV','FontSize',6,'HorizontalAlignment','right')



base = load('E:\Research_Projects\004_Propofol\Modelling\neuron_simulations\data\summary_data\expected_PSD_baseline.mat')
p0 = mean(base.p,2);
f = base.f;
load('E:\Research_Projects\004_Propofol\Experiments\_archive\FlavieCtrlsSpectraCzDropAligned20191210.mat')
y = nanmean(nanmean(SerialSpectraAmplCz(:,1:100,:),2),3);
mask = or(freq<55,freq>65);
figureNB;
    plot(freq,y.*mask,'LineWidth',1,'color','k')
    ylim([2e-1,5e1]*0.5)
    set(gca,'yscale','log');
    ylabel(['Power (' char(956) 'V^2/Hz)'])
    yyaxis right;
    plot(f,p0,'LineWidth',1,'color','r')
    set(gca,'yscale','log');
    ylim([2e-19,5e-17]*0.75)
    xlim([1,100]);
    set(gca,'xscale','log');
    xticks([1,10,100])

    ya = get(gca,'yaxis');
    ya(2).Visible = 'off';
    gcaformat
    xlabel('Frequency (Hz)');

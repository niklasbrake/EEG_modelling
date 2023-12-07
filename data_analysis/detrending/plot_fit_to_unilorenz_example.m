
dataFolder = 'E:\Research_Projects\004_Propofol\manuscript\Version3\Data';

load(fullfile(dataFolder,'data_time_information.mat'));
t0 = timeInfo.infusion_onset-timeInfo.object_drop;

load('E:\Research_Projects\004_Propofol\data\experiments\scalp_EEG\analyzed\Cz_multitaper_aligned.mat')
freq = aligned.freq;

[~,synFun] = fittingmodel('unilorenz');
[~,synFun2] = fittingmodel('exp2');

i = 1;

p0_eq6 = [20e-3,4e-3,-11.55,3.65];
p0_eq5 = [20e-3,-7.4,1.75];
sp = [10e-3,-7,2];
figureNB;
subplot(1,2,1);
    idcs = find(aligned.time<t0(i));
    y = nanmedian(aligned.psd(:,aligned.time<-1,i),2);
    % [px,synFun,full_model] = synDetrend(freq(freq<100),y(freq<100),3,'eq5',sp);
    plot(freq,y,'k','LineWidth',1);
    hold on;
    plot(freq,10.^synFun(freq,p0_eq5),'LineWidth',1,'color','b');
    plot(freq,10.^synFun2(freq,p0_eq6),'LineWidth',1,'color','r');
    set(gca,'xscale','log');
    set(gca,'yscale','log');
    xlim([0.5,100]);
    xticks([0.5,5,50]);
    xticklabels([0.5,5,50]);
    set(get(gca,'xaxis'),'MinorTick','off');
    set(get(gca,'yaxis'),'MinorTick','off');
    xlabel('Frequency (Hz)');
    ylim([1e-4,1e4])
    yticks(10.^[-2:2:4]);
    ylabel(['PSD (' char(956) 'V^2/Hz)']);
    gcaformat

p1_eq6 = [32e-3,4e-3,-14,3.8];
p1_eq5 = [32e-3,-1,2.2];
subplot(1,2,2);
    idcs = find(and(aligned.time>-10,aligned.time<0));
    y = nanmedian(aligned.psd(:,idcs,i),2);
    plot(freq,y,'k','LineWidth',1);
    hold on;
    plot(freq,10.^synFun(freq,p1_eq5),'LineWidth',1,'color','b');
    plot(freq,10.^synFun2(freq,p1_eq6),'LineWidth',1,'color','r');
    set(gca,'xscale','log');
    set(gca,'yscale','log');
    xlim([0.5,100]);
    xticks([0.5,5,50]);
    xticklabels([0.5,5,50]);
    set(get(gca,'xaxis'),'MinorTick','off');
    set(get(gca,'yaxis'),'MinorTick','off');
    xlabel('Frequency (Hz)');
    ylim([1e-4,1e4])
    yticks(10.^[-2:2:4]);
    ylabel(['PSD (' char(956) 'V^2/Hz)']);
    gcaformat





A = 10.^(p0_eq5(3))
lam = 10.^(p0_eq5(3))*exp(p0_eq5(2))

A = 10.^(p0_eq6(4))
lam = 10.^(p0_eq6(4))*exp(p0_eq6(3))

A = 10.^(p1_eq5(3))
lam = 10.^(p1_eq5(3))*exp(p1_eq5(2))

A = 10.^(p1_eq6(4))
lam = 10.^(p1_eq6(4))*exp(p1_eq6(3))
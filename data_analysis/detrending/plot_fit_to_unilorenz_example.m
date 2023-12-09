load(fullfile(dataFolder,'EEG_data','data_time_information.mat'));
t0 = timeInfo.infusion_onset-timeInfo.object_drop;

load(fullfile(dataFolder,'EEG_data','electrode2_Cz.mat'),'psd','freq','time');

% Compute Baseline spectrum
for i = 1:14
    preInfusion(:,i) = nanmedian(psd(:,and(time>=t0(i)-10,time<t0(i)),i),2);
    preLOC(:,i) = nanmedian(psd(:,and(time>=-10,time<0),i),2);
end

iNoise = find(and(freq>55,freq<65));
preInfusion(iNoise,:) = nan; preInfusion = fillgaps(preInfusion,5);
preLOC(iNoise,:) = nan; preLOC = fillgaps(preLOC,5);

[~,synFun] = fittingmodel('eq5');
[~,synFun2] = fittingmodel('eq6');

i = 1;

p0_eq6 = [20e-3,4e-3,-11.55,3.65];
p0_eq5 = [20e-3,-7.4,1.75];
sp = [10e-3,-7,2];
figureNB(7.7,4.4);
subplot(1,2,1);
    y = preInfusion(:,i);
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

p1_eq6 = [46e-3,4e-3,-13.7,3.66];
p1_eq5 = [46e-3,-11,2.2];
subplot(1,2,2);
    y = preLOC(:,i);
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



convert_parameters(p0_eq5,'eq5');
convert_parameters(p1_eq5,'eq5');
convert_parameters(p0_eq6,'eq6');
convert_parameters(p1_eq6,'eq6');

A = 10.^(p0_eq5(3))
lam = 10.^(p0_eq5(3))*exp(p0_eq5(2))

A = 10.^(p0_eq6(4))
lam = 10.^(p0_eq6(4))*exp(p0_eq6(3))

A = 10.^(p1_eq5(3))
lam = 10.^(p1_eq5(3))*exp(p1_eq5(2))

A = 10.^(p1_eq6(4))
lam = 10.^(p1_eq6(4))*exp(p1_eq6(3))
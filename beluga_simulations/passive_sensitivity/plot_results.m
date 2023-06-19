figureNB(18,6.3);


folder = 'E:\Research_Projects\004_Propofol\data\simulations\raw\passive_sensitivity\passive_sensitivity_analysis';

F = dir(folder);
F = F(3:end);

tauIRange = [5:5:30];
tauERange = [1:0.5:3.5];
erevRange = [-70:5:-45];
eFiringRange = round(10.^linspace(-1,1,6),2,'significant');
iFiringRange = round(10.^linspace(0,log10(20),6),2,'significant');
g_leak = round(10.^linspace(-4.5,-2.5,6),2,'significant');

clrs = clrsPT.sequential(10);
clrs = clrs(5:end,:);

P = {};
parValue = zeros(length(F),1);


for i = 1:5
    ax(i) = axes('Position',[0.3+0.09*i,0.67,0.06,0.3]);
end

for i = 1:5
    ax(i+5) = axes('Position',[0.3+0.09*i,0.17,0.06,0.3]);
end

f = (0.5:0.5:500)';
for i = 1:length(F)
    results = load(fullfile(folder,F(i).name));
    axes(ax(results.i))
    txt = strsplit(F(i).name,'_');
    colormap(clrs);
    switch txt{2}
        case 'tauI'
            set(gca,'CLim',[min(tauIRange),max(tauIRange)]);
            clr = interp1(tauIRange,clrs,results.pars.iSynParams.tau2);
            parValue(i) = results.pars.iSynParams.tau2;
            Ptemp = results.P;
            [params,synFun] = synDetrend(f,mean(Ptemp,2)./mean(Ptemp(1,:),2),0,'lorenz',[15e-3,1e-3,0,-1]);
            tauI(results.j) = params(1)*1e3
        case 'tauE'
            set(gca,'CLim',[min(tauERange),max(tauERange)]);
            clr = interp1(tauERange,clrs,results.pars.eSynParams.tau2);
            parValue(i) = results.pars.eSynParams.tau2;
            Ptemp = results.P;
            [params,synFun] = synDetrend(f,mean(Ptemp,2)./mean(Ptemp(1,:),2),0,'lorenz',[15e-3,1e-3,0,-1]);
            tauE(results.j) = params(2)*1e3
        case 'iFiring'
            set(gca,'CLim',[min(iFiringRange),max(iFiringRange)]);
            clr = interp1(iFiringRange,clrs,results.pars.iCellParams.firingRate);
            parValue(i) = results.pars.iCellParams.firingRate;
        case 'eFiring'
            set(gca,'CLim',[min(eFiringRange),max(eFiringRange)]);
            clr = interp1(eFiringRange,clrs,results.pars.eCellParams.firingRate);
            parValue(i) = results.pars.eCellParams.firingRate;
        case 'erev'
            set(gca,'CLim',[min(erevRange),max(erevRange)]);
            clr = interp1(erevRange,clrs,results.pars.biophys_pars.pas_mem_pars.erev_leak);
            parValue(i) = results.pars.biophys_pars.pas_mem_pars.erev_leak;
    end
    plotwitherror(f,results.P,'SE','color',clr);
    gcaformat;
    set(get(gca,'yaxis'),'visible','off');
    % colorbar;
    % title(txt{2});
    set(gca,'xscale','log');
    set(gca,'yscale','log');
    xlim([0.5,150]);
    ylim([1e-18,1e-15]);
    xticks([1,10,100]);
    xticklabels([1,10,100]);
    P{i} = sum(results.P.*0.5);
    I(i) = results.i;
    J(i) = results.j;
    paramChange{i} = txt{2};
    idcs = find(and(f>5,f<55));
    % [oofPre,oofFun] = getFOOOF(f(idcs),results.P(idcs,:),false);
    % alpha{i} = oofPre(:,2);
end

axL = axes();
axL.Position = ax(1).Position;
ylim([1e-18,1e-15]);
ylabel(['PSD (' char(956) 'V^2/Hz)'])
set(gca,'yscale','log');
set(get(gca,'xaxis'),'visible','off');
box off;
ylabel(['PSD (' char(956) 'V^2/Hz)']);
gcaformat;
axL.Position(1) = axL.Position(1)-0.01;
axL.Position(3) = 0.01;
axL.Position(4) = ax(1).Position(4);
axL.Position(2) = ax(1).Position(2);

for i = 1:5
    % i = iI(j);
    axes(ax(i+5));
    idcs = find(I==i);
    j = J(idcs);
    [x,jj] = sort(parValue(idcs));
    y = cat(1,P{idcs(jj)});
    % y = cat(2,alpha{idcs(jj)})';
    plotwitherror(x,y,'Q','color','k','LineWidth',1);
    xlabel(paramChange{idcs(1)});
    % ylim([0,1.5]);
    set(gca,'yscale','log');
    ylim([1e-15,1e-13]);
    set(get(gca,'yaxis'),'visible','off');
    gcaformat;
end

axL = axes();
axL.Position = ax(6).Position;
ylim([1e-15,1e-13]);
set(gca,'yscale','log');
ylabel(['Power (' char(956) 'V^2)']);
% ylim([0,1.5]);
% ylabel('Slope 5-55 Hz');
set(get(gca,'xaxis'),'visible','off');
box off;
gcaformat;
axL.Position(1) = axL.Position(1)-0.01;
axL.Position(3) = 0.01;
axL.Position(4) = ax(6).Position(4);
axL.Position(2) = ax(6).Position(2);


dataFolder = 'E:\Research_Projects\004_Propofol\manuscript\Version3\Data';
sims = load(fullfile(dataFolder,'simulation_passive_spectra.mat'));
load(fullfile(dataFolder,'simulation_passive_different_tau_decay.mat'),'GABAR_tau','spectra_fitted_tau');

[params,synFun] = synDetrend(sims.f(sims.f<100),mean(sims.P(sims.f<100,:),2)./mean(sims.P(1,:)),0,'lorenz',[15e-3,1e-3,0,-1]);

% figureNB;
% subplot(2,8,1);
axes('Position',[0.066,0.67,0.1,0.3]);
    plot(sims.f,sims.P,'color',[0.6,0.6,0.6,0.1])
    hold on;
    plot(sims.f,mean(sims.P,2),'k','LineWidth',1)
    ylim([10^-17,10^-14]);
    yticks([1e-16,1e-15])
    set(gca,'yscale','log');
    set(gca,'xscale','log')
    xticks([1,10,100]);
    xticklabels([1,10,100]);
    xlabel('Frequency (Hz)');
    ylabel(['PSD (' char(956) 'V^2/Hz)'])
    xlim([1,100])
    gcaformat;
axes('Position',[0.21,0.67,0.1,0.3]);
    plot(sims.f,mean(sims.P,2),'k','LineWidth',1);
    hold on;
    plot(sims.f,mean(sims.P(1,:))*10.^synFun(sims.f,[params(1:3),-Inf]),'--r');
    plot(sims.f,mean(sims.P(1,:))*10.^synFun(sims.f,[params(1:2),-Inf,params(4)]),'--r');
    % ylim([10^-16.5,10^-14.5]);
    ylim([10^-17,10^-14]);
    yticks([1e-16,1e-15])
    set(gca,'yscale','log');
    set(gca,'xscale','log')
    xticks([1,10,100]);
    xticklabels([1,10,100]);
    xlabel('Frequency (Hz)');
    % ylabel(['PSD (' char(956) 'V^2/Hz)'])
    yticklabels({});
    xlim([1,100])
    gcaformat;
    text(1.4,2e-15,'\tau_1','FontSize',7,'color','r')
    text(1.4,9e-17,'\tau_2','FontSize',7,'color','r')
axes('Position',[0.066,0.17,0.1,0.3]);
    plot(tauIRange,tauI,'k','LineWidth',1)
    xlim([5,20])
    ylabel('\tau_1 (ms)')
    xlabel('GABAR decay (ms)')
    gcaformat
axes('Position',[0.21,0.17,0.1,0.3]);
    plot(tauERange,tauE,'k','LineWidth',1)
    ylabel('\tau_2 (ms)')
    xlim([1,4])
    xlabel('AMPAR decay (ms)');
    gcaformat

load(fullfile(dataFolder,'simulation_passive_EEG_variance.mat'));
load(fullfile(dataFolder,'data_time_information.mat'))
t0 = timeInfo.infusion_onset-timeInfo.object_drop;
data = load(fullfile(dataFolder,'data_Cz_multitaper_meanRef.mat'));
pre = [];
for i = 1:length(t0)
    pre(:,i) = nanmedian(data.psd(:,data.time<t0(i),i),2);
end
freq = data.freq;
pre(and(freq>55,freq<65),:) = nan;
axes('Position',[0.886,0.17,0.1,0.8]);
    plotwitherror(freq,pre,'Q','LineWidth',1,'color','k'); hold on;
    plotwitherror(sims.f,sims.P*16e9,'Q','LineWidth',1,'color',clrsPT.diverging_CM(1,:)); hold on;
    set(gca,'xscale','log')
    set(gca,'yscale','log')
    ylim([1e-7,1e2]); yticks([1e-6,1e-2,1e2])
    xlim([1,100]);
    xticks([1,10,100]);
    xticklabels([1,10,100]);
    xlabel('Frequency (Hz)')
    ylabel(['PSD (' char(956) 'V^2/Hz)'])
    text(1.3,2.8e-4,sprintf('Poisson cortex\n(model)'),'FontSize',7,'color',clrsPT.diverging_CM(1,:));
    text(1.3,1,sprintf('Data'),'FontSize',7,'color','k');
    gcaformat




















dataFolder = 'E:\Research_Projects\004_Propofol\manuscript\Version3\Data';
mResults = load(fullfile(dataFolder,'simulation_avalanches_spectra.mat'));


clrs = clrsPT.sequential(10);
clrs = clrs(5:end,:);

figureNB(12.8,3.6);
subplot(1,3,1);
    for i = 1:size(P{2},2)
        plot(f,P{2}(:,i),'color',clrs(i,:),'LineWidth',1);
        hold on;
    end
    ylabel(['PSD (' char(956) 'V^2/Hz)']);
    xlabel('Frequency (Hz)')
    xlim([0.5,100]);
    ylim(10.^[-17,-14.5])
    % yticks(10.^[-16,-14]);
    xticks([1,10,100]);
    % xticklabels([1,10,100]);
    set(gca,'xscale','log')
    set(gca,'yscale','log')
    title(paramChange{2});
    gcaformat
subplot(1,3,2);
    for i = 1:size(P{4},2)
        plot(f,P{4}(:,i),'color',clrs(i,:),'LineWidth',1);
        hold on;
    end
    ylabel(['PSD (' char(956) 'V^2/Hz)']);
    xlabel('Frequency (Hz)')
    xlim([0.5,100]);
    ylim(10.^[-17,-14.5])
    % yticks(10.^[-16,-14]);
    xticks([1,10,100]);
    % xticklabels([1,10,100]);
    set(gca,'xscale','log')
    set(gca,'yscale','log')
    title(paramChange{4});
    gcaformat
subplot(1,3,3);
    for i = 1:size(P{5},2)
        plot(f,P{5}(:,i),'color',clrs(i,:),'LineWidth',1);
        hold on;
    end
    ylabel(['PSD (' char(956) 'V^2/Hz)']);
    xlabel('Frequency (Hz)')
    xlim([0.5,100]);
    ylim(10.^[-17,-14.5])
    % yticks(10.^[-16,-14]);
    xticks([1,10,100]);
    % xticklabels([1,10,100]);
    set(gca,'xscale','log')
    set(gca,'yscale','log')
    title(paramChange{5});
    gcaformat

figureNB(12.8,3.6);
subplot(1,3,1);
    for i = 1:size(P{1},2)
        plot(f,P{1}(:,i),'color',clrs(i,:),'LineWidth',1);
        hold on;
    end
    ylabel(['PSD (' char(956) 'V^2/Hz)']);
    xlabel('Frequency (Hz)')
    xlim([0.5,100]);
    ylim(10.^[-17,-14.5])
    % yticks(10.^[-16,-14]);
    xticks([1,10,100]);
    % xticklabels([1,10,100]);
    set(gca,'xscale','log')
    set(gca,'yscale','log')
    title(paramChange{1});
    gcaformat
subplot(1,3,2);
    % for i = 1:size(P{3},2)
    %     plot(f,P{3}(:,i),'color',clrs(i,:),'LineWidth',1);
    %     hold on;
    % end
    plot(data.f,mean(P_baseline,2),'color',clrs(1,:),'LineWidth',1);
    hold on;
    plot(data.f,mean(P_baseline_osc,2),'color',clrs(end,:),'LineWidth',1);
    ylabel(['PSD (' char(956) 'V^2/Hz)']);
    xlabel('Frequency (Hz)')
    xlim([0.5,100]);
    ylim(10.^[-17,-14.5])
    % yticks(10.^[-16,-14]);
    xticks([1,10,100]);
    % xticklabels([1,10,100]);
    set(gca,'xscale','log')
    set(gca,'yscale','log')
    title(paramChange{3});
    gcaformat
subplot(1,3,3);
    for i = 1:size(mResults.spectra,3)
        plot(mResults.f,mean(mResults.spectra(:,:,i),2),'color',clrs(i,:),'LineWidth',1);
        hold on;
    end
    ylabel(['PSD (' char(956) 'V^2/Hz)']);
    xlabel('Frequency (Hz)')
    xlim([0.5,100]);
    ylim(10.^[-17,-14.5])
    % yticks(10.^[-16,-14]);
    xticks([1,10,100]);
    % xticklabels([1,10,100]);
    set(gca,'xscale','log')
    set(gca,'yscale','log')
    title('branching');
    gcaformat






data = load(fullfile(dataFolder,'simulation_oscillations_spectra_0.1Hz.mat'));
P_tau_crit = data.P(:,:,2);
P_crit = data.P(:,:,4);
P_baseline = 1.2*data.P(:,:,6);
P_baseline_osc = data.P(:,:,7);
P_tau_osc = data.P(:,:,8);
P_tau = data.P(:,:,9);
P_crit_osc = data.P(:,:,10);

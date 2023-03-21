function plot_unitary_spectrum_and_tau
[sa,X] = network_simulation.getHeadModel;

load('E:\Research_Projects\004_Propofol\data\experiments\scalp_EEG\raw\timeInformation.mat')
t0 = timeInfo.infusion_onset-timeInfo.object_drop;

load('unitary_spectra_different_tau_decay.mat','tau','tau2');
load('E:\Research_Projects\004_Propofol\data\simulations\raw\passive\PSD.mat')
load('E:\Research_Projects\004_Propofol\data\simulations\raw\passive\totalVariance.mat')


data = load('E:\Research_Projects\004_Propofol\data\experiments\scalp_EEG\analyzed\Cz_multitaper_mean.mat');
pre = [];
for i = 1:length(t0)
    pre(:,i) = nanmedian(data.psd(:,data.time<t0(i),i),2);
end
freq = data.freq;
pre(and(freq>55,freq<65),:) = nan;

figureNB(12,6.2);
[params,synFun] = synDetrend(f(f<100),mean(P(f<100,:),2)./mean(P(1,:)),0,'lorenz',[15e-3,1e-3,0,-1]);
axes('Position',[0.17,0.61,0.16,0.3]);
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
    text(1.4,2e-15,'\tau_1','FontSize',7,'color','r')
    text(1.4,9e-17,'\tau_2','FontSize',7,'color','r')
axes('Position',[0.44,0.61,0.16,0.32]);
    errorbar(tau(1:end-1),mean(tau2(1:end-1,:),2),std(tau2(1:end-1,:),[],2),'.k','MarkerSize',7.5,'LineWidth',1)
    line([0,50],[0,50],'color','k')
    ylabel('EEG spectrum \tau_1 (ms)')
    xlim([7.5,32.5])
    ylim([5,45])
    ylabel('\tau_1 (ms)')
    xlabel('GABA_AR \tau_{decay} (ms)')
    gcaformat;
axes('Position',[0.17,0.21,0.43,0.21]);
    histogram(log10(mean(V)),'LineStyle','none','BinWidth',0.02,'FaceColor','k')
    xlim([-14.5,-13])
    ylabel('count')
    xlabel('Avg. power (uV^2)')
    xax = get(gca,'xaxis');
    xax.MinorTick = 'on';
    xax.MinorTickValues = [log10(linspace(10.^(-15),9*10.^(-15),11)),log10(linspace(10.^(-14),10.^(-13),11))];
    xticks([-14,-13])
    xticklabels({'10^{-14}','10^{-13}'})
    xlabel('Average single-neuron EEG power (uV^2)')
    gcaformat
axes('Position',[0.734,0.21,0.155,0.7]);
    plotwitherror(freq,pre,'Q','LineWidth',1,'color','b'); hold on;
    plotwitherror(f,P*25e9*3e6,'Q','LineWidth',1,'color','k'); hold on;
    set(gca,'xscale','log')
    set(gca,'yscale','log')
    ylim([1e-7,1e2]); yticks([1e-6,1e-2,1e2])
    xlim([1,100]);
    xticks([1,10,100]);
    xticklabels([1,10,100]);
    xlabel('Frequency (Hz)')
    ylabel(['PSD (' char(956) 'V^2/Hz)'])
    text(1.3,2.8e-4,sprintf('Poisson cortex\n(model)'),'FontSize',7,'color','k');
    text(1.3,1,sprintf('Data'),'FontSize',7,'color','b');
    gcaformat
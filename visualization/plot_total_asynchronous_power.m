% function plot_total_asynchronous_power

dThresh = 10.^(linspace(-1,log10(20),60));
load('E:\Research_Projects\004_Propofol\data\resources\head_models\surface_area_calculations.mat')
% Neuron count within a radius is the surface area times
% 100,000 neurons/mmm^2
N = 25e9;
dMids = 0.5*(dThresh(2:end)+dThresh(1:end-1));
nrnCount = mean(diff(Ar),2)*100000;
nrnCount(end) = N-sum(nrnCount(1:end-1));
corr_kernel = @(d) exp(-d.^2/10);

% Average pairwise correlation based on the :
%   - distribution of neuron densities within a certain radius of %     each point in cortex
%   - coupling kernel as a function of radius (normalized, so that
%     the correlation of nearby neurons is 1)
corr_kernel = @(d) exp(-d.^2/10);
rho_bar = sum(corr_kernel(dMids').*nrnCount)/sum(nrnCount')


load('E:\Research_Projects\004_Propofol\data\experiments\scalp_EEG\model_fits\pre.mat');

load('E:\Research_Projects\004_Propofol\data\simulations\analyzed\_archive\network_criticality_spectra_(s=1).mat')

asynchUnitarySpec = mean(mean(P(:,:,:),3),2);
SIG0 = mean(sum(asynchUnitarySpec*mean(diff(f))));
P0 = interp1(f,asynchUnitarySpec,freq);

% asynchUnitarySpec = mean(mean(P(:,1:3,:),3),2);
% SIG0 = mean(sum(asynchUnitarySpec*mean(diff(f))));
% P0 = interp1(f,asynchUnitarySpec,freq);

SIG_N = @(rho) N+N*(N-1)*rho;

figureNB(6.5,3)
axes('Position',[0.175,0.3,0.22,0.5]);
    plotwitherror(freq,pre,'Q','color',[0.6,0.6,0.6]);
    plot(freq,P0*50/SIG0,'r')
    plot(freq,P0*200/SIG0,'r')
    % plot(freq,(P0*100/SIG0+P0*100/SIG0)/2,'b','LineWidth',1)
    set(gca,'xscale','log');
    set(gca,'yscale','log');
    xlim([0.5,50]);
    xticks([0.5,5,50])
    xticklabels([0.5,5,50])
    xlabel('Frequency (Hz)');
    ylim([1e-2,2e2])
    yticks([1e-2,1,1e2])
    ylabel(['PSD (' char(956) 'V^2/Hz)'])
    gcaformat;
axes('Position',[0.55,0.3,0.22,0.5]);
    dscale = 10.^linspace(-1,1,1e3);
    t = 10.^linspace(-2,log10(01),1e3);
    tt = [];
    for i = 1:length(dscale)
        corr_kernel2 = @(d) exp(-d.^2/dscale(i));
        rho_bar = sum(corr_kernel2(dMids').*nrnCount)/sum(nrnCount');
        tt(:,i) = SIG0*SIG_N(t*rho_bar);
    end
    imagesc((t),(dscale),log10(tt)); hold on;
    contour((t),(dscale),log10(tt),log10([50,200]),'color','r')
    xlim([0,0.4])
    ylim([0,10]);
    gcaformat
    axis xy
    ylabel('\sigma (mm)')
    xlabel('\rho_{max}')
    colormap(clrsPT.iridescent(1e3))
    CB = colorbar('location','eastoutside');
    CB.Label.String = 'Total EEG power (uV^2)';
    set(gca,'CLim',[-1,3]);
    CB.Ticks = [-1:3];
    CB.TickLabels = {'0.1','1','10','100','1000'};


figureNB;
subplot(1,2,1);
    tt = [];
    dscale = 10.^linspace(-2,1,1e3);
    for i = 1:length(dscale)
        corr_kernel2 = @(d) exp(-d.^2/dscale(i));
        rho_bar = sum(corr_kernel2(dMids').*nrnCount)/sum(nrnCount');
        tt(i) = SIG0*SIG_N(rho_bar);
    end
    plot(dscale,tt,'k','LineWidth',1)
    set(gca,'xscale','log')
    set(gca,'yscale','log')
    line(get(gca,'xlim'),[50,50],'color','r');
    line(get(gca,'xlim'),[200,200],'color','r');
    xlabel('\sigma^2 (mm)');
    ylabel(['Total EEG power (' char(956) 'V^2)'])
    xlim([dscale(1),dscale(end)]);
    ylim([1e-2,1e4]);
    title('\rho_{max} = 1');
subplot(1,2,2);
    dscale = 4;
    t = 10.^linspace(-2,log10(1),1e3);
    corr_kernel2 = @(d) exp(-d.^2/dscale);
    rho_bar = sum(corr_kernel2(dMids').*nrnCount)/sum(nrnCount');
    tt = SIG0*SIG_N(t*rho_bar);
    plot(t,tt,'k','LineWidth',1)
    set(gca,'xscale','log')
    set(gca,'yscale','log')
    line(get(gca,'xlim'),[50,50],'color','r');
    line(get(gca,'xlim'),[200,200],'color','r');
    title(sprintf('%.2f',interp1(tt,t,200)));
    xlabel('\rho_{max}')
    ylim([1e-2,1e4]);
    title('\sigma^2 = 4 mm');
dataFolder = 'E:\Research_Projects\004_Propofol\manuscript\Version3\Data';

load(fullfile(dataFolder,'simulation_avalanches_dipole_correlations.mat'));
clrs = clrsPT.sequential(10); clrs = clrs(5:end,:);

figureNB(5.9,3.5)
mResults = load(fullfile(dataFolder,'simulation_avalanches_spectra.mat'));
subplot(1,2,1);
    m = sort([0.99,0.98,0.95,0.86,0.63,0]);
    h=[];
    for i = 1:length(m)
        y = mResults.spectra(:,:,i);
        h(i) = plotwitherror(mResults.f,y,'CI','color',clrs(i,:),'LineWidth',1);
    end
    ylabel(['PSD (' char(956) 'V^2/Hz)']);
    xlabel('Frequency (Hz)')
    xlim([0.5,100]);
    ylim(10.^[-17,-13])
    yticks(10.^[-16,-14]);
    xticks([1,10,100]);
    xticklabels([1,10,100]);
    set(gca,'xscale','log')
    set(gca,'yscale','log')
    gcaformat;


dataFolder = 'E:\Research_Projects\004_Propofol\manuscript\Version3\Data';
load(fullfile(dataFolder,'data_time_information.mat'))
t0 = timeInfo.infusion_onset-timeInfo.object_drop;
data = load(fullfile(dataFolder,'data_Cz_multitaper_meanRef.mat'));
pre = [];
for i = 1:length(t0)
    pre(:,i) = nanmedian(data.psd(:,data.time<t0(i),i),2);
end
freq = data.freq;
pre(and(freq>55,freq<65),:) = nan;

load(fullfile(dataFolder,'anatomy_cortical_pairwise_distance_distribution.mat'));
signed_area = A;
total_area = B;
N = 16e9;
dMids = 0.5*(rValues(2:end)+rValues(1:end-1));
nrnCount = mean(diff(signed_area),2)*200000;
nrnCount(end) = N-sum(nrnCount(1:end-1));
corr_kernel = @(d) exp(-d.^2/4);
rho_bar = sum(corr_kernel(dMids).*nrnCount)/sum(nrnCount');
SIG_N = @(rho) N+N*(N-1)*rho*rho_bar;


P_crit = mResults.spectra(:,:,5); 
subplot(1,2,2);
    plotwitherror(freq,pre,'CI','color',[0.5,0.5,0.5]);
    plotwitherror(mResults.f,P_crit*SIG_N(0.13)+0.25+0.2*randn(size(P_crit)),'CI','color',clrs(6,:),'LineWidth',1);
    set(gca,'xscale','log');
    set(gca,'yscale','log');
    xlim([0.5,100]);
    xlabel('Frequency (Hz)');
    xticks([1,10,100])
    xticklabels([1,10,100])
    ylabel('PSD (uV^2/Hz)');
    xax = get(gca,'xaxis');
    xax.TickLabelRotation=0;
    gcaformat;

gcaformat(gcf);




dataFolder = 'E:\Research_Projects\004_Propofol\manuscript\Version3\Data';

load(fullfile(dataFolder,'simulation_avalanches_dipole_correlations.mat'));
clrs = clrsPT.sequential(10); clrs = clrs(5:end,:);


C = zeros(11,11);
for i0 = 1:11
    folder = ['E:\Research_Projects\004_Propofol\data\simulations\raw\dyad_dipole_correlation\run' int2str(i0-1)];
    F = dir(folder); F = F(3:end);
    for j0 = 1:length(F)
        s = F(j0).name;
        X = csvread(fullfile(folder,s,['LFPy\cell00001.csv']),1,0);
        dipoles(:,:,1) = X(2:end-16e2,2:4);
        X = csvread(fullfile(folder,s,['LFPy\cell00002.csv']),1,0);
        dipoles(:,:,2) = X(2:end-16e2,2:4);
        for i = 1:2
        end
        for j = 1:3
            C(i0,j0) = C(i0,j0)+corr(dipoles(:,j,1),dipoles(:,j,2))/3;
        end
    end
end

figureNB(5.5,2.8)
subplot(1,2,1);
    line([1,100],[0,0],'color',[0.6,0.6,0.6],'LineWidth',1);
    hold on;
    for i = 1:length(m)
        x = 1/(1-m(i));
        y = C_s0(i,:)';
        y = nanmean(C_s0(i,:));
        y_lo = icdf('normal',0.025,0,1)*stderror(C_s0(i,:)');
        y_hi = icdf('normal',0.975,0,1)*stderror(C_s0(i,:)');
        % plot(x,y,'.','color',[0.6,0.6,0.6],'MarkerSize',5,'LineWidth',1)
        line([x,x],[y+y_lo,y+y_hi],'color',[0.6,0.6,0.6],'linewidth',1);
    end
    plot(1./(1-m),nanmean(C_s1,2),'color','k','LineWidth',1);
    for i = 1:length(m)
        x = 1/(1-m(i));
        y = C_s1(i,:)';
        y = nanmean(C_s1(i,:));
        y_lo = icdf('normal',0.025,0,1)*stderror(C_s1(i,:)');
        y_hi = icdf('normal',0.975,0,1)*stderror(C_s1(i,:)');
        % plot(x,y,'.','color',clrs(i,:),'MarkerSize',5,'LineWidth',1)
        line([x,x],[y+y_lo,y+y_hi],'color',clrs(i,:),'linewidth',1);
    end
    set(gca,'xscale','log')
    xlim([1,100]);
    % ylim([-0.05,0.6])
    ylim([-0.1,0.65])
    xl = xlabel('1/(1-m)');
    xl.Position(2) = -0.24;
    ylabel('Dipole correlation');
    gcaformat
    text(1.5,0.4,{'Optimized','synapses'},'FontSize',6,'VerticalAlignment','top')
    text(100,0.16,{'Random','synapses'},'FontSize',6,'VerticalAlignment','top','HorizontalAlignment','right','color',[0.6,0.6,0.6])
subplot(1,2,2);
    plot(linspace(0,1,11),flipud(C'),'color',clrs(6,:)*0.4+[1,1,1]*0.6);
    hold on;
    plot(linspace(0,1,11),fliplr(mean(C)),'color',clrs(6,:),'LineWidth',1);
    xlabel('Optimality index');
    ylabel('Dipole correlation');
    gcaformat;
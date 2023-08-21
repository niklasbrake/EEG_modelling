
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
corr_kernel = @(d) exp(-d.^2/5);
rho_bar = sum(corr_kernel(dMids).*nrnCount)/sum(nrnCount');
SIG_N = @(rho) N+N*(N-1)*rho*rho_bar;


dataFolder = 'E:\Research_Projects\004_Propofol\manuscript\Version3\Data';

load(fullfile(dataFolder,'simulation_avalanches_dipole_correlations.mat'));
clrs = clrsPT.sequential(10); clrs = clrs(5:end,:);


MM = round(M_S1+5e-3,3);


figureNB(5.9,3.5)
mResults = load(fullfile(dataFolder,'simulation_avalanches_spectra.mat'));
subplot(1,2,1);
    m = sort([0.99,0.98,0.95,0.86,0.63,0]);
    h=[];
    for i = 1:length(m)
        y = mResults.spectra(:,:,i);
        y = y*SIG_N(MM(i));
        h(i) = plotwitherror(mResults.f,y,'CI','color',clrs(i,:),'LineWidth',1);
    end
    ylabel(['PSD (' char(956) 'V^2/Hz)']);
    xlabel('Frequency (Hz)')
    xlim([0.5,100]);
    % ylim(10.^[-17,-13])
    % yticks(10.^[-16,-14]);
    xticks([1,10,100]);
    xticklabels([1,10,100]);
    set(gca,'xscale','log')
    set(gca,'yscale','log')
    gcaformat;





subplot(1,2,2);
    m = sort([0.99,0.98,0.95,0.86,0.63,0]);
    plotwitherror(freq,pre,'CI','color',[0.5,0.5,0.5]);
    h=[];
    for i = 4:length(m)
        y = mResults.spectra(:,:,i);
        y = y*SIG_N(MM(i)*0.65)+0.25+0.2*randn(size(y));
        h(i) = plotwitherror(mResults.f,y,'CI','color',clrs(i,:),'LineWidth',1);
    end
    ylabel(['PSD (' char(956) 'V^2/Hz)']);
    xlabel('Frequency (Hz)')
    xlim([0.5,100]);
    % ylim(10.^[-17,-13])
    % yticks(10.^[-16,-14]);
    xticks([1,10,100]);
    xticklabels([1,10,100]);
    set(gca,'xscale','log')
    set(gca,'yscale','log')
    gcaformat;

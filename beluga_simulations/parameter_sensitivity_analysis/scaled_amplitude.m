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


load('E:\Research_Projects\004_Propofol\data\simulations\raw\parameter_sensitivity_analysis2\fitted_spectra.mat')

Pm = splitapply(@(x) nanmedian(x,2),psd,pars(9,:));
P0 = splitapply(@(x) quantile(x,0.975,2),psd,pars(9,:));
P1 = splitapply(@(x) quantile(x,0.025,2),psd,pars(9,:));

P0 = sum(P0.*m,2)/sum(m);
P1 = sum(P1.*m,2)/sum(m);
Pm = sum(Pm.*m,2)/sum(m);

f = f(:);

clr = clrsPT.diverging_CM(1,:);
figureNB;
    plotwitherror(freq,pre,'M','LineWidth',1,'color','k'); hold on;
    F = fill([f;flipud(f)],[P0;flipud(P1)]*16e9,'k');
    F.FaceAlpha = 0.2;
    F.FaceColor = clr;
    F.EdgeColor = 'none';
    plot(f,Pm*16e9,'LineWidth',1,'color',clr)
    set(gca,'xscale','log')
    set(gca,'yscale','log')
    ylim([1e-7,1e2]); yticks([1e-6,1e-2,1e2])
    xlim([1,100]);
    xticks([1,10,100]);
    xticklabels([1,10,100]);
    ylim([1e-8,1e2]); 
    xlabel('Frequency (Hz)')
    ylabel(['PSD (' char(956) 'V^2/Hz)'])
    text(1.3,2.8e-4,sprintf('Poisson cortex\n(model)'),'FontSize',7,'color',clrsPT.diverging_CM(1,:));
    text(1.3,1,sprintf('Data'),'FontSize',7,'color','k');
    gcaformat
ax = gca
ax.Units = 'centimeters'
ax.Position(3:4) = [2,4.75];
[sa,X] = network_simulation_beluga.getHeadModel;
locations = randi(size(sa.cortex75K.vc,1),100,1); % Random location


for i = 1:14
    load(['E:\Research_Projects\005_Aperiodic_EEG\data\simulations\EI_crit\WT\run' int2str(i) '\LFPy\simulation_data.mat'])
    dp_WT(:,:,i) = dipoles(time>=100,:);
    load(['E:\Research_Projects\005_Aperiodic_EEG\data\simulations\EI_crit\PMS\run' int2str(i) '\LFPy\simulation_data.mat'])
    dp_PMS(:,:,i) = dipoles(time>=100,:);
end
[f,P_WT] = network_simulation_beluga.expectedEEGspectrum(dp_WT,sa,locations);
[f,P_PMS] = network_simulation_beluga.expectedEEGspectrum(dp_PMS,sa,locations);

idcs = find(and(f>5,f<55));
OOF_WT = getFOOOF(f(idcs)',P_WT(idcs,:),false);
OOF_PMS = getFOOOF(f(idcs)',P_PMS(idcs,:),false);

% idcs2 = randperm(100,12);
idcs2 = 1:size(P_WT,2);


clr1 = [210,133,65]/255;
clr2 = [70,115,157]/255;
figureNB;
subplot(1,2,1);
    plotwitherror(f,P_WT,'CI','color',clr1,'LineWidth',1);
    plotwitherror(f,P_PMS,'CI','color',clr2,'LineWidth',1);
    set(gca,'yscale','log')
    xlim([0,60])
    xlabel('Frequency (Hz)');
    ylabel(['PSD (' char(956) 'V^2/Hz)'])
    gcaformat;
subplot(1,2,2);
    bp = boxplotNB(1,OOF_WT(idcs2,2),clr1,10);
    bp.area.FaceAlpha = 0.7;
    bp.data.MarkerFaceAlpha = 1;
    bp=boxplotNB(2,OOF_PMS(idcs2,2),clr2,10);
    bp.area.FaceAlpha = 0.7;
    bp.data.MarkerFaceAlpha = 1;
    [~,p] = ttest2(OOF_WT(idcs2,2),OOF_PMS(idcs2,2))
    % ylim([-1.5,0])
    ylim([0.85,1.6])
    xticks([1,2])
    xlim([0.5,2.5]);
    xticklabels({'WT','PMS'});
    ylabel('Slope (10-50 Hz)')
    set(get(gca,'xaxis'),'visible','off')
    gcaformat
    text(1,0.8,'WT','HorizontalAlignment','center','FontSize',8,'color',clr1,'Rotation',0)
    text(2,0.8,'PMS','HorizontalAlignment','center','FontSize',8,'color',clr2,'Rotation',0)
    yticks([0.9:0.1:1.5])






figureNB(4.5,4.5);
subplot(1,2,1);
    bar(1,0.35,'FaceColor',clr1); hold on;
    bar(2,0.55,'FaceColor',clr2);
    ylabel('Avg. excitatory firing rate (Hz)')
    xticks([1,2]);
    xticklabels({'WT','Shank3B^{-/-}'})
    xl = get(gca,'xaxis');
    xl.TickLabelRotation = 45;
subplot(1,2,2);
    bar(1,1.25,'FaceColor',clr1); hold on;
    bar(2,0.95,'FaceColor',clr2);
    ylabel('Avg. inhibitory firing rate (Hz)')
    xticks([1,2]);
    xticklabels({'WT','Shank3B^{-/-}'})
    xl = get(gca,'xaxis');
    xl.TickLabelRotation = 45;
gcaformat(gcf)

figureNB(4,4.5);
    plotwitherror(f,P_WT*1e16,'CI','color',clr1,'LineWidth',1);
    plotwitherror(f,P_PMS*1e16,'CI','color',clr2,'LineWidth',1);
    set(gca,'yscale','log')
    xlim([0,60])
    xlabel('Frequency (Hz)');
    ylabel(['PSD (' char(956) 'V^2/Hz)'])
    gcaformat;

figureNB;
    bp = boxplotNB(1,OOF_WT(idcs2,2),clr1,10);
    bp.area.FaceAlpha = 0.7;
    bp.data.MarkerFaceAlpha = 1;
    bp=boxplotNB(2,OOF_PMS(idcs2,2),clr2,10);
    bp.area.FaceAlpha = 0.7;
    bp.data.MarkerFaceAlpha = 1;
    [~,p] = ttest2(OOF_WT(idcs2,2),OOF_PMS(idcs2,2))
    % ylim([-1.5,0])
    ylim([0.85,1.6])
    xticks([1,2])
    xlim([0.5,2.5]);
    xticklabels({'WT','PMS'});
    ylabel('Slope (10-50 Hz)')
    set(get(gca,'xaxis'),'visible','off')
    gcaformat
    text(1,0.8,'WT','HorizontalAlignment','center','FontSize',8,'color',clr1,'Rotation',0)
    text(2,0.8,'PMS','HorizontalAlignment','center','FontSize',8,'color',clr2,'Rotation',0)
    yticks([0.9:0.1:1.5])

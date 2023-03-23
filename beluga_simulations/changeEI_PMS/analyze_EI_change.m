[sa,X] = network_simulation_beluga.getHeadModel;
locations = randi(size(sa.cortex75K.vc,1),100,1); % Random location

load('E:\Research_Projects\005_Aperiodic_EEG\data\simulations\passive\WT\LFPy\simulation_data.mat')
[f,P_WT] = network_simulation_beluga.expectedEEGspectrum(dipoles(time>=100,:,:),sa,locations);

load('E:\Research_Projects\005_Aperiodic_EEG\data\simulations\passive\PMS\LFPy\simulation_data.mat')
[f,P_PMS] = network_simulation_beluga.expectedEEGspectrum(dipoles(time>=100,:,:),sa,locations);


figureNB;
    plotwitherror(f,P_WT,'CI','color','k','LineWidth',1);
    plotwitherror(f,P_PMS,'CI','color','r','LineWidth',1);
    set(gca,'yscale','log')
    xlim([0,60])
    xlabel('Frequency (Hz)');
    ylabel(['PSD (' char(956) 'V^2/Hz)'])
    gcaformat;


idcs = find(and(f>10,f<50));
OOF_WT = getFOOOF(f(idcs)',P_WT(idcs,:),false);
OOF_PMS = getFOOOF(f(idcs)',P_PMS(idcs,:),false);

idcs2 = randperm(100,12);

figureNB;
subplot(1,2,1);
    plotwitherror(f,P_WT,'CI','color','k','LineWidth',1);
    plotwitherror(f,P_PMS,'CI','color','r','LineWidth',1);
    set(gca,'yscale','log')
    xlim([0,60])
    xlabel('Frequency (Hz)');
    ylabel(['PSD (' char(956) 'V^2/Hz)'])
    gcaformat;
subplot(1,2,2);
    boxplotNB(1,-OOF_WT(idcs2,2),'k',10);
    boxplotNB(2,-OOF_PMS(idcs2,2),'r',10);
    [~,p] = ttest2(OOF_WT(idcs2,2),OOF_PMS(idcs2,2))
    ylim([-1.5,0])
    xticks([1,2])
    xlim([0.5,2.5]);
    xticklabels({'WT','PMS'});
    ylabel('Slope (10-50 Hz)')
    set(get(gca,'xaxis'),'visible','off')
    gcaformat
    text(1,-1.45,'WT','HorizontalAlignment','center','FontSize',8,'color','k','Rotation',45)
    text(2,-1.45,'PMS','HorizontalAlignment','center','FontSize',8,'color','r','Rotation',45)


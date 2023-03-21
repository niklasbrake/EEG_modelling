tauResults = load('E:\Research_Projects\004_Propofol\data\simulations\analyzed\analzye_simulations_tau.mat');
mResults = load('E:\Research_Projects\004_Propofol\data\simulations\analyzed\analzye_simulations_crit.mat');
load('E:\Research_Projects\004_Propofol\data\simulations\raw\test\analzye_simulations.mat')


delT = tauResults.P1(:,:,6)./tauResults.P1(:,:,1);
delM = mResults.P1(:,:,5)./mResults.P1(:,:,1);


fig = figureNB(18.3,6);
axes('Position',[0.07,0.15,0.23,0.75]);
    h(1) = plotwitherror(mResults.f,mResults.P1(:,:,1),'Q','LineWidth',1);
    hold on;
    h(2) = plotwitherror(mResults.f,mResults.P1(:,:,5),'Q','LineWidth',1);
    xlim([0.5,100]);
    ylim(10.^[-18,-13])
    set(gca,'xscale','log');
    set(gca,'yscale','log');
    title('Criticality changed (\tau = 10)')
    xlabel('Frequency (Hz)')
    ylabel(['PSD (' char(956) 'V^2/Hz)']);
    L = legend(h,{'m = 0','m = 0.98'}); L.ItemTokenSize = [10,10];
    L.Box = 'off';
axes('Position',[0.4,0.15,0.23,0.75]);
    h(1) = plotwitherror(tauResults.f,tauResults.P1(:,:,1),'Q','LineWidth',1);
    hold on;
    h(2) = plotwitherror(tauResults.f,tauResults.P1(:,:,6),'Q','LineWidth',1);
    xlim([0.5,100]);
    ylim(10.^[-18,-13])
    set(gca,'xscale','log');
    set(gca,'yscale','log');
    title('Synapses changed (m = 0)')
    xlabel('Frequency (Hz)')
    ylabel(['PSD (' char(956) 'V^2/Hz)']);
    L = legend(h,{'\tau = 10','\tau = 30'}); L.ItemTokenSize = [10,10];
    L.Box = 'off';
axes('Position',[0.73,0.15,0.23,0.75]);
    h(1) = plotwitherror(tauResults.f,tauResults.P1(:,:,1),'Q','LineWidth',1);
    h(2) = plotwitherror(f,P1(:,:,2),'Q','LineWidth',1);
    h(3) = plotwitherror(tauResults.f,tauResults.P1(:,:,1).*delM.*delT,'Q','LineWidth',1);
    xlim([0.5,100]);
    ylim(10.^[-18,-13])
    set(gca,'xscale','log');
    set(gca,'yscale','log');
    title('Both changed')
    xlabel('Frequency (Hz)')
    ylabel(['PSD (' char(956) 'V^2/Hz)']);
    L = legend(h,{'m = 0, \tau = 10','m = 0.98, \tau = 30','Prediction from multiplication'}); L.ItemTokenSize = [10,10];
    L.Box = 'off';
gcaformat(gcf)

% load('E:\Research_Projects\004_Propofol\data\simulations\raw\parameter_sensitivity_analysis\analyzed_results_all.mat')
% load('E:\Research_Projects\004_Propofol\data\simulations\raw\parameter_sensitivity_analysis\fitted_spectra.mat')

% tauI when gL is high
% idcs = find(pars(6,:)>1e-3);
idcs = 1:length(pars);
psd0 = psd(:,idcs);
idcs2 = interp1(linspace(0,1,5),1:5,P(4,idcs),'nearest','extrap');
p_gL_high = splitapply(@(x) mean(x,2),psd0,idcs2);

% tauI when gL is low
% idcs = find(pars(6,:)<1e-4);
idcs = 1:length(pars);
psd0 = psd(:,idcs);
idcs2 = interp1(linspace(0,1,5),1:5,P(5,idcs),'nearest','extrap');
p_gL_low = splitapply(@(x) mean(x,2),psd0,idcs2);


figureNB(5.3,3);
clrs = clrsPT.sequential(10); clrs = clrs(5:end,:);
ax(1) = axes('Position',[0.2,0.26,0.32,0.6]);
ax(2) = axes('Position',[0.62,0.26,0.32,0.6]);
for i = 1:5
    axes(ax(1));
        plot(f,p_gL_high(:,i),'color',clrs(i,:),'LineWidth',1);
        hold on;
    axes(ax(2));
        plot(f,p_gL_low(:,i),'color',clrs(i,:),'LineWidth',1);
        hold on;
end

axes(ax(1));
    set(gca,'xscale','log');
    set(gca,'yscale','log');
    xlabel('Frequency (Hz)')
    ylabel(['PSD (' char(956) 'V^2/Hz)'])
    xlim([1,100]);
    ylim([1e-17,2e-15])
    txt = text(2,6.6e-17,'\lambda_I min','FontSize',6,'color',clrs(1,:));
    txt = text(4,4.4e-16,'\lambda_I max','FontSize',6,'color',clrs(5,:));
    title('g_L high','FontSize',7,'Fontweight','normal')
axes(ax(2));
    set(gca,'xscale','log');
    set(gca,'yscale','log');
    xlabel('Frequency (Hz)')
    xlim([1,100]);
    ylim([1e-17,2e-15])
    yticks([]);
    title('g_L low','FontSize',7,'Fontweight','normal')
gcaformat(gcf)



figureNB(5.3,3);
subplot(1,2,1);
    idcs = find(pars(8,:)>1e-3);
    beta = oof_1_40(idcs,2);
    EI = pars(2,idcs)./pars(1,idcs);
    % EI = pars(6,idcs)./pars(5,idcs);
    % idcs2 = find(and(EI>0.1,EI<50));
    % EI = EI(idcs2); beta = beta(idcs2);
    h = plot(EI,beta,'.','color',[0,0,0],'MarkerSize',1);
    hold on
    FT = fitlm(log10(EI),beta);
    FT.Coefficients{2,4}
    corr(log10(EI)',beta)
    t = 10.^linspace(-5,5,1e3)';
    plot(t,FT.predict(log10(t)),'color','r','LineWidth',1);
    set(gca,'xscale','log');
    xlim([0.02,50])
    ylabel('Slope (1-40 Hz)');
    xticks([0.1,1,10])
    xticklabels({'10:1','1:1','1:10'})
    xlabel('E:I ratio')
    ylim([-0.5,1.5])
    txt1 = text(0.1,1,sprintf('\\rho = %.2f',corr(log10(EI)',beta)),'FontSize',7,'FontWeight','normal','Color','r');
    title('g_L high','FontSize',7,'Fontweight','normal')
    gcaformat;
subplot(1,2,2);
    idcs = find(pars(8,:)<1e-4);
    beta = oof_1_40(idcs,2);
    EI = pars(2,idcs)./pars(1,idcs);
    % EI = pars(6,idcs)./pars(5,idcs);
    % idcs2 = find(and(EI>0.1,EI<50));
    % EI = EI(idcs2); beta = beta(idcs2);
    h = plot(EI,beta,'.','color',[0,0,0],'MarkerSize',1);
    hold on
    FT = fitlm(log10(EI),beta);
    FT.Coefficients{2,4}
    
    t = 10.^linspace(-5,5,1e3)';
    plot(t,FT.predict(log10(t)),'color','r','LineWidth',1);
    set(gca,'xscale','log');
    xlim([0.02,50])
    xticks([0.1,1,10])
    xticklabels({'10:1','1:1','1:10'})
    xlabel('E:I ratio')
    ylim([-0.5,1.5])
    txt2 = text(0.1,1,sprintf('\\rho = %.2f',corr(log10(EI)',beta)),'FontSize',7,'FontWeight','normal','Color','r');
    title('g_L low','FontSize',7,'Fontweight','normal')
    gcaformat;

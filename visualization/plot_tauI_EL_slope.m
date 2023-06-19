load('E:\Research_Projects\004_Propofol\data\simulations\raw\parameter_sensitivity_analysis\fitted_spectra.mat')
clrs = clrsPT.sequential(10); clrs = clrs(5:end,:);

idcs2 = interp1(linspace(0,1,5),1:5,P(4,:),'nearest','extrap');
p_tauI = splitapply(@(x) mean(x,2),psd,idcs2);

idcs2 = interp1(linspace(0,1,5),1:5,P(5,:),'nearest','extrap');
p_EL = splitapply(@(x) mean(x,2),psd,idcs2);

figureNB(5.3,3);
subplot(1,2,1);
    for i = 1:5
        p = p_tauI(:,i);
        plot(f,10*log10(p/p(3)),'color',clrs(i,:),'LineWidth',1);
        hold on;
    end
    set(gca,'xscale','linear');
    xlabel('Frequency (Hz)')
    ylabel('Power (dB)');
    xlim([1,100]);
    xlim([1,40])
    xticks([1,20,40]);
    text(30,0,'\tau_I min','FontSize',7,'color',clrs(1,:));
    text(10,-6,'\tau_I max','FontSize',7,'color',clrs(5,:));
subplot(1,2,2);
    for i = 1:5
        p = p_EL(:,i);
        plot(f,10*log10(p/p(3)),'color',clrs(i,:),'LineWidth',1);
        hold on;
    end
    set(gca,'xscale','linear');
    xlabel('Frequency (Hz)')
    xlim([1,100]);
    xlim([1,40])
    xticks([1,20,40]);
    ylabel('Power (dB)');
    text(30,0,'E_L min','FontSize',7,'color',clrs(1,:));
    text(10,-6,'E_L max','FontSize',7,'color',clrs(5,:));
gcaformat(gcf)

% load('E:\Research_Projects\004_Propofol\data\simulations\raw\parameter_sensitivity_analysis\analyzed_results_all.mat')
% load('E:\Research_Projects\004_Propofol\data\simulations\raw\parameter_sensitivity_analysis\fitted_spectra.mat')
load('E:\Research_Projects\004_Propofol\data\simulations\raw\parameter_sensitivity_analysis2\fitted_spectra.mat')
load('E:\Research_Projects\004_Propofol\data\simulations\raw\parameter_sensitivity_analysis2\parameter_sensitivity_analysis2.mat')

% Define distributions of parameters
lambda_E = @(p) cdf('logn',p,log(0.5),1);
lambda_I = @(p) cdf('logn',p,log(2.5),1);
tauE = @(p) cdf('unif',p,1,3.5);
tauI = @(p) cdf('unif',p,5,20);
gE = @(p) cdf('unif',p,0.2e-3,2e-3);
gI = @(p) cdf('unif',p,0.2e-3,2e-3);
erev = @(p) cdf('unif',p,-75,-45);
gleak = @(p) cdf('unif',log10(p),-5,log10(0.005));
% Distribution of morphologies
neuronTypes = {'L23E_oi24rpy1';'L23I_oi38lbc1';'L23I_oi38lbc1';'L4E_53rpy1';'L4E_j7_L4stellate';'L4E_j7_L4stellate';'L4I_oi26rbc1';'L4I_oi26rbc1';'L5E_oi15rpy4';'L5E_j4a';'L5I_oi15rbc1';'L5I_oi15rbc1';'L6E_51_2a_CNG';'L6E_oi15rpy4';'L6I_oi15rbc1';'L6I_oi15rbc1'};
nrnAbundance = [26.8,3.2,4.3,9.5,9.5,9.5,5.6,1.5,4.9,1.3,0.6,0.8,14,4.6,1.9,1.9];
mTypeCount = accumarray(findgroups(neuronTypes),nrnAbundance);
mTypeP = cumsum(mTypeCount)/sum(mTypeCount);
mType = @(p) interp1(mTypeP,1:11,p,'next','extrap');
samplePars = @(P) [lambda_E(P(:,1)), lambda_I(P(:,2)), tauE(P(:,3)), ...
tauI(P(:,4)), gE(P(:,5)), gI(P(:,6)), erev(P(:,7)), gleak(P(:,8)), mType(P(:,9))];
P = samplePars(pars')';




EI = pars(2,:)./pars(1,:);
K = linspace(log10(min(EI)),log10(max(EI)),300);
b0 = 0.5*(K(2:end)+K(1:end-1));
h = histcounts(log10(EI),K);

K0 = quantile(EI,linspace(0,1,6));

clrs = clrsPT.sequential(10); clrs(1:5,:) = [];

figureNB(9,2.8);
axes('Position',[0.03, 0.27, 0.16, 0.78])
    for i = 1:5
        fill([b0,fliplr(b0)],max(h)*1.2*(i-1)+[h*0,fliplr(h)],'k','FaceAlpha',0.2,'EdgeColor','none');
        hold on;

        m0 = K0(i);
        m1 = K0(i+1);
        idcs = find(and(EI>=m0,EI<m1));
        h1 = histcounts(log10(EI(idcs)),K);
        fill([b0,fliplr(b0)],max(h)*1.2*(i-1)+[h1*0,fliplr(h1)],clrs(i,:),'EdgeColor','none');
        xlim(log10([0.05,50]))
        line(log10([0.05,50]),max(h)*1.2*(i-1).*[1,1],'color','k');
    end
    xticks([-1,0,1]);
    xticklabels({'10:1','1:1','1:10'});
    % xlabel('E:I ratio');
    xlabel('E:I ratio (\lambda_E:\lambda_I)')
    set(get(gca,'yaxis'),'visible','off')
    gcaformat;
axes('Position',[0.28, 0.27, 0.14, 0.68]);
    for i = 1:5
        m0 = K0(i);
        m1 = K0(i+1);
        idcs = find(and(EI>=m0,EI<m1));
        p = mean(psd(:,idcs),2);
        plot(f,10*log10(p./p(1)),'LineWidth',1,'color',clrs(i,:));
        xticks([1,20,40]);
        xlim([1,40]);
        xlabel('Frequency (Hz)');
        ylabel('Power (dB)');
        hold on;
    end
    gcaformat;
axes('Position',[0.55, 0.27, 0.18, 0.68])
    idcs = find(pars(8,:)<1e-4);
    beta = oof_1_40(idcs,2);
    EI = pars(2,idcs)./pars(1,idcs);
    h = plot(EI,beta,'.','color',[0,0,0],'MarkerSize',1);
    hold on
    FT = fitlm(log10(EI),beta);
    rho = corr(log10(EI)',beta);
    R = FT.Rsquared.Adjusted;
    p = FT.Coefficients{2,4};
    fprintf('\\rho = %.2f, R^2 = %.2f, p=%.2e\n',rho,R,p);
    t = 10.^linspace(-5,5,1e3)';
    plot(t,FT.predict(log10(t)),'color','r','LineWidth',1);
    set(gca,'xscale','log');
    xlim([0.02,50])
    xticks([0.1,1,10])
    xticklabels({'10:1','1:1','1:10'})
    xlabel('E:I ratio (\lambda_E:\lambda_I)')
    ylim([-0.5,1.5])
    txt2 = text(0.03,-0.45,sprintf('\\rho = %.2f',corr(log10(EI)',beta)),'FontSize',7,'FontWeight','normal','Color','r','VerticalAlignment','bottom');
    text(0.05,1.5,'g_L low','FontSize',7,'Fontweight','normal')
    yl = ylabel('Slope (1-40 Hz)');
    % yl.Position(1) = 0.003;
    gcaformat;
axes('Position',[0.8, 0.27, 0.18, 0.68])
    idcs = find(pars(8,:)>1e-3);
    beta = oof_1_40(idcs,2);
    EI = pars(2,idcs)./pars(1,idcs);
    h = plot(EI,beta,'.','color',[0,0,0],'MarkerSize',1);
    hold on
    FT = fitlm(log10(EI),beta);
    rho = corr(log10(EI)',beta);
    R = FT.Rsquared.Adjusted;
    p = FT.Coefficients{2,4};
    fprintf('\\rho = %.2f, R^2 = %.2f, p=%.2e\n',rho,R,p);
    t = 10.^linspace(-5,5,1e3)';
    plot(t,FT.predict(log10(t)),'color','r','LineWidth',1);
    set(gca,'xscale','log');
    xlim([0.02,50])
    xticks([0.1,1,10])
    xticklabels({'10:1','1:1','1:10'})
    xlabel('E:I ratio (\lambda_E:\lambda_I)')
    ylim([-0.5,1.5])
    txt1 = text(0.03,-0.45,sprintf('\\rho = %.2f',corr(log10(EI)',beta)),'FontSize',7,'FontWeight','normal','Color','r','VerticalAlignment','bottom');
    text(0.05,1.5,'g_L high','FontSize',7,'Fontweight','normal')
    gcaformat;



figureNB(7.3,10.5);
subplot(3,2,1);
    idcs = find(pars(8,:)>1e-3);
    % beta = oof_1_40(idcs,2);
    beta = oof_30_50(idcs,2);
    EI = pars(2,idcs)./pars(1,idcs);
    % EI = pars(6,idcs)./pars(5,idcs);
    % idcs2 = find(and(EI>0.1,EI<50));
    % EI = EI(idcs2); beta = beta(idcs2);
    h = plot(EI,beta,'.','color',[0,0,0],'MarkerSize',1);
    ylim([-1,4]);
    hold on
    FT = fitlm(log10(EI),beta);
    rho = corr(log10(EI)',beta);
    R = FT.Rsquared.Adjusted;
    p = FT.Coefficients{2,4};
    fprintf('\\rho = %.2f, R^2 = %.2f, p=%.2e\n',rho,R,p);
    t = 10.^linspace(-5,5,1e3)';
    plot(t,FT.predict(log10(t)),'color','r','LineWidth',1);
    set(gca,'xscale','log');
    xlim([0.02,50])
    yl = ylabel('Slope (30-50 Hz)');
    yl.Position(1) = 0.003;
    xticks([0.1,1,10])
    xticklabels({'10:1','1:1','1:10'})
    xlabel('E:I ratio (\lambda_E:\lambda_I)')
    txt1 = text(0.05,4,sprintf('\\rho = %.2f',corr(log10(EI)',beta)),'FontSize',7,'FontWeight','normal','Color','r','VerticalAlignment','top');
    title('g_L high','FontSize',7,'Fontweight','normal')
    gcaformat;
subplot(3,2,2);
    idcs = find(pars(8,:)<1e-4);
    % beta = oof_1_40(idcs,2);
    beta = oof_30_50(idcs,2);
    EI = pars(2,idcs)./pars(1,idcs);
    % EI = pars(6,idcs)./pars(5,idcs);
    % idcs2 = find(and(EI>0.1,EI<50));
    % EI = EI(idcs2); beta = beta(idcs2);
    h = plot(EI,beta,'.','color',[0,0,0],'MarkerSize',1);
    ylim([-1,4]);
    hold on
    FT = fitlm(log10(EI),beta);
    rho = corr(log10(EI)',beta);
    R = FT.Rsquared.Adjusted;
    p = FT.Coefficients{2,4};
    fprintf('\\rho = %.2f, R^2 = %.2f, p=%.2e\n',rho,R,p);
    t = 10.^linspace(-5,5,1e3)';
    plot(t,FT.predict(log10(t)),'color','r','LineWidth',1);
    set(gca,'xscale','log');
    xlim([0.02,50])
    xticks([0.1,1,10])
    xticklabels({'10:1','1:1','1:10'})
    xlabel('E:I ratio (\lambda_E:\lambda_I)')
    txt2 = text(0.05,4,sprintf('\\rho = %.2f',corr(log10(EI)',beta)),'FontSize',7,'FontWeight','normal','Color','r','VerticalAlignment','top');
    title('g_L low','FontSize',7,'Fontweight','normal')
    gcaformat;



subplot(3,2,3);
    idcs = find(pars(8,:)>1e-3);
    beta = oof_1_40(idcs,2);
    EI = pars(6,idcs)./pars(5,idcs);
    h = plot(EI,beta,'.','color',[0,0,0],'MarkerSize',1);
    hold on
    FT = fitlm(log10(EI),beta);
    rho = corr(log10(EI)',beta);
    R = FT.Rsquared.Adjusted;
    p = FT.Coefficients{2,4};
    fprintf('\\rho = %.2f, R^2 = %.2f, p=%.2e\n',rho,R,p);
    t = 10.^linspace(-5,5,1e3)';
    plot(t,FT.predict(log10(t)),'color','r','LineWidth',1);
    set(gca,'xscale','log');
    xlim([0.02,50])
    yl = ylabel('Slope (1-40 Hz)');
    yl.Position(1) = 0.003;
    xticks([0.1,1,10])
    xticklabels({'10:1','1:1','1:10'})
    xlabel('E:I ratio (g_E:g_I)')
    ylim([-0.5,1.5])
    txt1 = text(0.05,1.5,sprintf('\\rho = %.2f',corr(log10(EI)',beta)),'FontSize',7,'FontWeight','normal','Color','r','VerticalAlignment','top');
    title('g_L high','FontSize',7,'Fontweight','normal')
    gcaformat;
subplot(3,2,4);
    idcs = find(pars(8,:)<1e-4);
    beta = oof_1_40(idcs,2);
    EI = pars(6,idcs)./pars(5,idcs);
    h = plot(EI,beta,'.','color',[0,0,0],'MarkerSize',1);
    hold on
    FT = fitlm(log10(EI),beta);
    rho = corr(log10(EI)',beta);
    R = FT.Rsquared.Adjusted;
    p = FT.Coefficients{2,4};
    fprintf('\\rho = %.2f, R^2 = %.2f, p=%.2e\n',rho,R,p);
    t = 10.^linspace(-5,5,1e3)';
    plot(t,FT.predict(log10(t)),'color','r','LineWidth',1);
    set(gca,'xscale','log');
    xlim([0.02,50])
    xticks([0.1,1,10])
    xticklabels({'10:1','1:1','1:10'})
    xlabel('E:I ratio (g_E:g_I)')
    ylim([-0.5,1.5])
    txt2 = text(0.05,1.5,sprintf('\\rho = %.2f',corr(log10(EI)',beta)),'FontSize',7,'FontWeight','normal','Color','r','VerticalAlignment','top');
    title('g_L low','FontSize',7,'Fontweight','normal')
    gcaformat;

subplot(3,2,5);
    idcs = find(pars(8,:)>1e-3);
    beta = oof_30_50(idcs,2);
    EI = pars(6,idcs)./pars(5,idcs);
    h = plot(EI,beta,'.','color',[0,0,0],'MarkerSize',1);
    ylim([-1,4]);
    hold on
    FT = fitlm(log10(EI),beta);
    rho = corr(log10(EI)',beta);
    R = FT.Rsquared.Adjusted;
    p = FT.Coefficients{2,4};
    fprintf('\\rho = %.2f, R^2 = %.2f, p=%.2e\n',rho,R,p);
    t = 10.^linspace(-5,5,1e3)';
    plot(t,FT.predict(log10(t)),'color','r','LineWidth',1);
    set(gca,'xscale','log');
    xlim([0.02,50])
    yl = ylabel('Slope (30-50 Hz)');
    yl.Position(1) = 0.003;
    xticks([0.1,1,10])
    xticklabels({'10:1','1:1','1:10'})
    xlabel('E:I ratio (g_E:g_I)')
    txt1 = text(0.05,4,sprintf('\\rho = %.2f',corr(log10(EI)',beta)),'FontSize',7,'FontWeight','normal','Color','r','VerticalAlignment','top');
    title('g_L high','FontSize',7,'Fontweight','normal')
    gcaformat;
subplot(3,2,6);
    idcs = find(pars(8,:)<1e-4);
    beta = oof_30_50(idcs,2);
    EI = pars(6,idcs)./pars(5,idcs);;
    h = plot(EI,beta,'.','color',[0,0,0],'MarkerSize',1);
    ylim([-1,4]);
    hold on
    FT = fitlm(log10(EI),beta);
    rho = corr(log10(EI)',beta);
    R = FT.Rsquared.Adjusted;
    p = FT.Coefficients{2,4};
    fprintf('\\rho = %.2f, R^2 = %.2f, p=%.2e\n',rho,R,p);
    t = 10.^linspace(-5,5,1e3)';
    plot(t,FT.predict(log10(t)),'color','r','LineWidth',1);
    set(gca,'xscale','log');
    xlim([0.02,50])
    xticks([0.1,1,10])
    xticklabels({'10:1','1:1','1:10'})
    xlabel('E:I ratio (g_E:g_I)')
    txt2 = text(0.05,4,sprintf('\\rho = %.2f',corr(log10(EI)',beta)),'FontSize',7,'FontWeight','normal','Color','r','VerticalAlignment','top');
    title('g_L low','FontSize',7,'Fontweight','normal')
    gcaformat;




% Plot example fit
% i = randi(size(psd,1));
i = 4788;
figureNB(2.5,3);
    plot(f,psd(:,i),'color','k','LineWidth',1)
    hold on;
    plot(f,10.^oof_1_40(i,1)./f.^oof_1_40(i,2),'color','r','LineWidth',1)
    set(gca,'xscale','log')
    set(gca,'yscale','log')
    xlim([1,40])
    ylabel(['PSD (' char(956) 'V^2/Hz)']);;
    xlabel('Frequency (Hz)');
    xticks([1,40])
    gcaformat


EI = pars(2,:)./pars(1,:);
K = linspace(log10(min(EI)),log10(max(EI)),300);
b0 = 0.5*(K(2:end)+K(1:end-1));
h = histcounts(log10(EI),K);

K0 = quantile(EI,linspace(0,1,6));

clrs = clrsPT.sequential(10); clrs(1:5,:) = [];
figureNB;
for i = 1:5
    fill([b0,fliplr(b0)],max(h)*1.2*(i-1)+[h*0,fliplr(h)],'k','FaceAlpha',0.2,'EdgeColor','none');
    hold on;

    m0 = K0(i);
    m1 = K0(i+1);
    idcs = find(and(EI>=m0,EI<m1));
    h1 = histcounts(log10(EI(idcs)),K);
    fill([b0,fliplr(b0)],max(h)*1.2*(i-1)+[h1*0,fliplr(h1)],clrs(i,:),'EdgeColor','none');
    xlim(log10([0.05,50]))
    line(log10([0.05,50]),max(h)*1.2*(i-1).*[1,1],'color','k');
end
xticks([-1,0,1]);
xticklabels({'10:1','1:1','1:10'});
xlabel('E:I ratio');
set(get(gca,'yaxis'),'visible','off')
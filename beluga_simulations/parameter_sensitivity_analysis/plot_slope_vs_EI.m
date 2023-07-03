% load('E:\Research_Projects\004_Propofol\data\simulations\raw\parameter_sensitivity_analysis\analyzed_results_all.mat')
% load('E:\Research_Projects\004_Propofol\data\simulations\raw\parameter_sensitivity_analysis\fitted_spectra.mat')
load('E:\Research_Projects\004_Propofol\data\simulations\raw\parameter_sensitivity_analysis2\fitted_spectra.mat')

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
    FT = fitlm(log10(EI),beta)
    FT.Coefficients{2,4}
    corr(log10(EI)',beta)
    t = 10.^linspace(-5,5,1e3)';
    plot(t,FT.predict(log10(t)),'color','r','LineWidth',1);
    set(gca,'xscale','log');
    xlim([0.02,50])
    ylabel('Slope (1-40 Hz)');
    xticks([0.1,1,10])
    xticklabels({'10:1','1:1','1:10'})
    xlabel('E:I ratio (\lambda_E:\lambda_I)')
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
    FT = fitlm(log10(EI),beta)
    FT.Coefficients{2,4}
    
    t = 10.^linspace(-5,5,1e3)';
    plot(t,FT.predict(log10(t)),'color','r','LineWidth',1);
    set(gca,'xscale','log');
    xlim([0.02,50])
    xticks([0.1,1,10])
    xticklabels({'10:1','1:1','1:10'})
    xlabel('E:I ratio (\lambda_E:\lambda_I)')
    ylim([-0.5,1.5])
    txt2 = text(0.1,1,sprintf('\\rho = %.2f',corr(log10(EI)',beta)),'FontSize',7,'FontWeight','normal','Color','r');
    title('g_L low','FontSize',7,'Fontweight','normal')
    gcaformat;


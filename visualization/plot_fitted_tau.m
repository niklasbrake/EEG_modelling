load('E:\Research_Projects\004_Propofol\data\simulations\raw\parameter_sensitivity_analysis\fitted_spectra.mat')
f = f(2:end);
psd = psd(2:end,:);

k = 6;


idcs2 = interp1(linspace(0,1,k),1:k,P(3,:),'nearest','extrap');
p_tauE = splitapply(@(x) mean(x,2),psd,idcs2);
tau1 = zeros(50,k);
tauE = zeros(50,k);
for j = 1:k
    idcs = find(idcs2==j);
    m = length(idcs);
    for i = 1:50
        j0 = randi(m,m,1);
        tauE(i,j) = mean(pars(3,idcs(j0)));
        P_AMPA(:,i) = mean(psd(:,idcs(j0)),2);
    end
    [params,synFun] = synDetrend(f,P_AMPA./P_AMPA(1,:),0,'lorenz',[15e-3,1e-3,0,-1]);
    tau1(:,j) = params(:,2)*1e3;
    p_tauE(:,j) = mean(P_AMPA,2);
end

idcs2 = interp1(linspace(0,1,k),1:k,P(4,:),'nearest','extrap');
p_tauI = splitapply(@(x) mean(x,2),psd,idcs2);
tau2 = zeros(50,k);
tauI = zeros(50,k);
for j = 1:k
    idcs = find(idcs2==j);
    m = length(idcs);
    for i = 1:50
        j0 = randi(m,m,1);
        tauI(i,j) = mean(pars(4,idcs(j0)));
        P_GABA(:,i) = mean(psd(:,idcs(j0)),2);
    end
    [params,synFun] = synDetrend(f,P_GABA./P_GABA(1,:),0,'lorenz',[15e-3,1e-3,0,-1]);
    tau2(:,j) = params(:,1)*1e3;
    p_tauI(:,j) = mean(P_GABA,2);
end




clrs = clrsPT.sequential(k+5); clrs = clrs(5:end,:);



dataFolder = 'E:\Research_Projects\004_Propofol\manuscript\Version3\Data';
sims = load(fullfile(dataFolder,'simulation_passive_spectra.mat'));
sims.P = sims.P;
% sims.f = f;
idcs = find(sims.f<100);
p_fit = mean(sims.P(idcs,:)./mean(sims.P(1,:)),2);
[params,synFun] = synDetrend(sims.f(idcs),p_fit,0,'lorenz',[15e-3,2e-3,0,-1]);



figureNB(8.3,2.7);
B = 0.21;
axes('Position',[0.13,0.28,B,0.64]);
    plot(sims.f,mean(sims.P,2),'k','LineWidth',1);
    hold on;
    plot(sims.f,mean(sims.P(1,:))*10.^synFun(sims.f,[params(1:3),-Inf]),'--r');
    plot(sims.f,mean(sims.P(1,:))*10.^synFun(sims.f,[params(1:2),-Inf,params(4)]),'--r');
    ylim([10^-16.5,10^-14.5]);
    yticks([1e-16,1e-15])
    set(gca,'yscale','log');
    set(gca,'xscale','log')
    xticks([1,10,100]);
    xticklabels([1,10,100]);
    xlabel('Frequency (Hz)');
    ylabel(['PSD (' char(956) 'V^2/Hz)'])
    xlim([1,100])
    gcaformat;
    text(1.4,2e-15,'\tau_1','FontSize',7,'color','r')
    text(1.4,9e-17,'\tau_2','FontSize',7,'color','r')
    set(gca,'XTickLabelRotation',0)
A = 0.96-B;
axes('Position',[A,0.28,B,0.64]);
    plot(tauE,tau1,'.k','MarkerSize',5);
    FT = fitlm(tauE(:),tau1(:));
    hold on;
    t = linspace(0,50,1e3)';
    plot(t,FT.predict(t),'-r','LineWidth',1); 
    xlim([1,3.5]);
    xlabel('AMPAR \tau_E (ms)');
    ylabel('\tau_2 (ms)')
axes('Position',[(0.13+A)/2,0.28,B,0.64]);
   plot(tauI,tau2,'.k','MarkerSize',5);
    FT = fitlm(tauI(:),tau2(:));
    hold on;
    t = linspace(0,50,1e3)';
    plot(t,FT.predict(t),'-r','LineWidth',1);
    xlim([5,20]);
    xlabel('GABAR \tau_I (ms)');
    ylabel('\tau_1 (ms)')
gcaformat(gcf)
folder = 'E:\Research_Projects\004_Propofol\data\simulations\raw\passive_sensitivity\passive_sensitivity_analysis';

F = dir(folder);
F = F(3:end);

for i = 1:length(F)
    results = load(fullfile(folder,F(i).name));
    txt = strsplit(F(i).name,'_');
    % data.(txt{2}).
    switch txt{2}
        case 'tauI'
            parValue(i) = results.pars.iSynParams.tau2;
        case 'tauE'
            parValue(i) = results.pars.eSynParams.tau2;
            Ptemp = results.P;
            % [params,synFun] = synDetrend(f,mean(Ptemp,2)./mean(Ptemp(1,:),2),0,'lorenz',[15e-3,1e-3,0,-1]);
            % tauE(results.j) = params(2)*1e3
        case 'iFiring'
            parValue(i) = results.pars.iCellParams.firingRate;
        case 'eFiring'
            parValue(i) = results.pars.eCellParams.firingRate;
        case 'erev'
            parValue(i) = results.pars.biophys_pars.pas_mem_pars.erev_leak;
        case 'gleak'
            parValue(i) = results.pars.biophys_pars.pas_mem_pars.g_leak;
    end
    data.(txt{2}).P(:,:,results.j) = results.P;
    data.(txt{2}).par(results.j) = parValue(i);
    data.(txt{2}).pw(results.j,:) = 0.5*sum(results.P);
end

f = 0.5:0.5:500;
fn = fieldnames(data);
[a,b] = getSubplotSize(length(fn));


clrs = clrsPT.sequential(10); clrs = clrs(5:end,:);
figureNB;
for i = 1:length(fn)
    subplot(a,b,i)
    [~,I] = sort(data.(fn{i}).par);
    data.(fn{i}).par = data.(fn{i}).par(I);
    data.(fn{i}).P = data.(fn{i}).P(:,:,I);
    data.(fn{i}).pw = data.(fn{i}).pw(I,:);
    for j = 1:length(I)
        plotwitherror(f,data.(fn{i}).P(:,:,j),'CI','color',clrs(j,:));
    end
    title(fn{i});
    set(gca,'xscale','log');
    set(gca,'yscale','log');
    xlim([0.5,150]);
    ylim([1e-18,1e-15]);
    xticks([1,10,100]);
    xticklabels([1,10,100]);
    ylabel(['PSD (' char(956) 'V^2/Hz)'])
    set(gca,'yscale','log');
    set(get(gca,'xaxis'),'visible','off');
    box off;
end


for i = 1:length(fn)
    par = repmat(data.(fn{i}).par',[1,100]);
    x = mean(par,2);
    x = (x-min(x))/(max(x)-min(x));
    y = mean(data.(fn{i}).pw,2);

    FT = fitlm(x,y)
    beta(i) = abs(FT.Coefficients{2,1});
end

p_names = {'\lambda_E','E_{leak}','g_{leak}','\lambda_I','\tau_E','\tau_I'};

figureNB(5.7,3);
subplot(1,3,1);
    bar(beta)
    xticklabels(p_names);
subplot(1,3,2);
    plotwitherror(data.eFiring.par,data.eFiring.pw,'Q','LineWidth',1,'color','k');
    ylim([1e-15,1e-13]);
    set(gca,'yscale','log');
    set(gca,'xscale','log')
    xlabel('\lambda_E (Hz)');
    ylabel('Power (uV^2)')
subplot(1,3,3);
    plotwitherror(data.gleak.par,data.gleak.pw,'Q','LineWidth',1,'color','k');
    ylim([1e-15,1e-13]);
    set(gca,'yscale','log');
    set(gca,'xscale','log')
    xlabel('g_{leak} (S/cm^2)');
    yticklabels({});

figureNB(5.7,2.5);
axes('Position',[0.18,0.35,0.15,0.58]);
    plotwitherror(data.eFiring.par,data.eFiring.pw,'Q','LineWidth',1,'color','k');
    ylim([1e-15,1e-13]);
    set(gca,'yscale','log');
    set(gca,'xscale','log')
    % xlabel('\lambda_E (Hz)');
    xlabel('\lambda_E');
    ylabel('Power (uV^2)')
axes('Position',[0.38,0.35,0.15,0.58]);
    plotwitherror(data.iFiring.par,data.iFiring.pw,'Q','LineWidth',1,'color','k');
    ylim([1e-15,1e-13]);
    set(gca,'yscale','log');
    set(gca,'xscale','log')
    % xlabel('\lambda_I (Hz)');
    xlabel('\lambda_I');
    yticklabels({});
axes('Position',[0.58,0.35,0.15,0.58]);
    plotwitherror(data.erev.par,data.erev.pw,'Q','LineWidth',1,'color','k');
    ylim([1e-15,1e-13]);
    set(gca,'yscale','log');
    % xlabel('E_{leak} (mV)');
    xlabel('E_{leak}');
    yticklabels({});
axes('Position',[0.78,0.35,0.15,0.58]);
    plotwitherror(data.gleak.par,data.gleak.pw,'Q','LineWidth',1,'color','k');
    ylim([1e-15,1e-13]);
    set(gca,'yscale','log');
    set(gca,'xscale','log')
    % xlabel('g_{leak} (S/cm^2)');
    xlabel('g_{leak}');
    yticklabels({});
gcaformat(gcf)



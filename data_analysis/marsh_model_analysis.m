% function marsh_model_analysis(e0)

load('E:\Research_Projects\004_Propofol\manuscript\Version3\Data\data_time_information.mat');

[t,p,pum,tauFold,pBins,M] = marsh_model(1.6);

tauOrser = [19.4,29.6,44.5]/19.4;
pOrser = [0,0.5,2];
tauKita = [1,1.2,1.4,1.45,2];
pKita = [0.1,0.3,1,3,10];
pWhittington2 = [0,0.5,1,2,5,10];
tauWhittington2 = [17.5,17.5,19,22.5,50,75]/17.5;
pWhittington1 = [0,0.5,1,2,5,10];
tauWhittington1 = [17.5,17.5,19,48,58,75]/17.5;

% pWhittington1 = pOrser;
% tauWhittington1 = tauOrser;


pWhittington2 = pOrser;
tauWhittington2 = tauOrser;

[PP,I] = sort([pOrser,pKita,pWhittington1,pWhittington2]);
TT = [tauOrser,tauKita,tauWhittington1,tauWhittington2];
TT = TT(I)-1;

ftfun = fittype('a*x.^n./(k.^n + x.^n)');
FT = fit(PP(:),TT(:),ftfun,'StartPoint',[3,1,1],'Lower',[0,0,0],'Upper',[Inf,Inf,8])
confInt = predint(FT,PP);
%{

figureNB
subplot(1,4,1);
    % plotwitherror(tvec,tau_unscaled,'CI','color','k');
    plotwitherror(tRescaled,1e3*squeeze(params(1,:,:)),'CI','color','k');
    xlabel('LOC-aligned time (s)');
    ylabel(['\tau (ms)'])
    % xlim([-200,60])
    xlim([-1.5,0.5]);
    xticks([-1,0])
    xticklabels({'Infusion','LOC'})
subplot(1,4,2);
    for i = 1:14
        plot(t{i}./dt(i),p{i},'color',[0.6,0.6,0.6]);
        hold on;
    end
    xlabel('LOC-aligned time (s)');
    ylabel(['Effect-site concentration (' char(956) 'M)'])
    % xlim([-200,60])
    xlim([-1.5,0.5]);
    xticks([-1,0])
    xticklabels({'Infusion','LOC'})
subplot(1,4,3);
    scatter(pum(:),tauFold(:),2,[0.5,0.5,0.5],'filled','MarkerFaceAlpha',0.2);
    hold on;
    % plot(nanmedian(pVec2,2),nanmedian(tauSmooth,2),'k')
    plot(pBins,M,'color','k','LineWidth',1);
    set(gca,'xscale','log');
    xlabel(['Effect-site concentration (' char(956) 'M)']);
    ylabel('Estimated \tau (fold change)')
    ylabel('\tau (fold change)')
    xlim([0.01,10])
    line(get(gca,'xlim'),[1,1],'lineStyle','--','color','k')
    % plotPaperResults


ke0 = 10.^linspace(log10(0.26),log10(1.21),10);
clrs = clrsPT.sequential(length(ke0)+4);
clrs = clrs(4:end,:);
subplot(1,4,4);
for i = 1:length(ke0)
    e0 = ke0(i);
    [t,p,pum,tauFold,pBins,M] = marsh_model(e0);
    plot(pBins,M,'color',clrs(i,:)*0.7+0.3,'LineWidth',1);
    hold on;
end
set(gca,'xscale','log');
xlabel(['Effect-site concentration (' char(956) 'M)']);
ylabel('Estimated \tau (fold change)')
xlim([0.01,10])
ylim([0.5,3.5]);
line(get(gca,'xlim'),[1,1],'lineStyle','--','color','k')
plotPaperResults;


PP_test = linspace(0.01,100,1e3);
confInt = predint(FT,PP_test);
plot(PP_test,FT(PP_test)+1,'k','LineWidth',1,'LineStyle','--');
% F = fill([sort(PP_test),fliplr(sort(PP_test))],1+[confInt(:,1);flipud(confInt(:,2))],'b','LineStyle','none','FaceAlpha',0.1)

gcaformat(gcf);


%}


figureNB(8,3.4)
subplot(1,2,1);
    scatter(pum(:),tauFold(:),2,[0.5,0.5,0.5],'filled','MarkerFaceAlpha',0.2);
    hold on;
    % plot(nanmedian(pVec2,2),nanmedian(tauSmooth,2),'k')
    plot(pBins,M,'color','k','LineWidth',1);
    set(gca,'xscale','log');
    xl = xlabel(['Estimated effect-site propofol concentration (' char(956) 'M)']);
    xl.Position(1) = 12;
    ylabel('Estimated \tau (fold change)')
    ylabel('\tau_1 (fold change)')
    xlim([0.01,10])
    xticks([0.01,0.1,1,10]);
    xticklabels([0.01,0.1,1,10]);
    line(get(gca,'xlim'),[1,1],'lineStyle','--','color','k')
    ylim([0,4.5])
    % plotPaperResults
ke0 = 10.^linspace(log10(0.26),log10(1.21),10);
clrs = clrsPT.sequential(length(ke0)+4);
clrs = clrs(4:end,:);
subplot(1,2,2);
% for i = 1:length(ke0)
%     e0 = ke0(i);
%     [t,p,pum,tauFold,pBins,M] = marsh_model(e0);
%     plot(pBins,M,'color',clrs(i,:)*0.7+0.3,'LineWidth',1);
%     hold on;
% end
% plot(pBins,M,'color','k','LineWidth',2);
scatter(pum(:),tauFold(:),2,[0.5,0.5,0.5],'filled','MarkerFaceAlpha',0.2);
hold on;
set(gca,'xscale','log');
% xlabel(['Effect-site concentration (' char(956) 'M)']);
ylabel('\tau_1 (fold change)')
xlim([0.01,10])
xticks([0.01,0.1,1,10]);
xticklabels([0.01,0.1,1,10]);
ylim([0.5,3.5]);
line(get(gca,'xlim'),[1,1],'lineStyle','--','color','k')
plotPaperResults;
ylim([0,4.5])


PP_test = linspace(0.01,100,1e3);
confInt = predint(FT,PP_test);
plot(PP_test,FT(PP_test)+1,'k','LineWidth',1,'LineStyle','--');
% F = fill([sort(PP_test),fliplr(sort(PP_test))],1+[confInt(:,1);flipud(confInt(:,2))],'b','LineStyle','none','FaceAlpha',0.1)

gcaformat(gcf);
% end
function [t,p,pum,tauFold,pBins,M] = marsh_model(ke0);
    if(nargin==0)
        ke0 = 1.21;
        % ke0 = 0.26;% /min
    end
    V1 = 0.228*1e3; % ml/kg
    V2 = 0.463*1e3; % ml/kg
    V3 = 2.893*1e3; % ml/kg

    k10 = 0.119; %/min
    k12 = 0.112; %/min
    k13 = 0.0419; %/min
    k21 = 0.055; %/min
    k31 = 0.0033; %/min

    A = [-k12-k13-k10, V2/V1*k21, V3/V1*k31, 0; ...
            V1/V2*k12,-k21,0, 0; ...
            V1/V3*k13,0,-k31, 0; ...
            ke0,0,0,-ke0];
    B = [1/V1;0;0;0];

    load('E:\Research_Projects\004_Propofol\manuscript\Version3\Data\data_time_information.mat');
    load('E:\Research_Projects\004_Propofol\manuscript\Version3\Data\data_fitted_params_rescaled_time.mat')
    timeInfo.infusion_onset(1) = timeInfo.object_drop(1)-200;
    tau = squeeze(params(1,:,:))*1e3;

    % load('E:\Research_Projects\004_Propofol\data\experiments\scalp_EEG\model_fits\rescaled_20230701.mat')
    % tau = squeeze(pars(:,1,:))*1e3;
    % tRescaled =  linspace(-1.5,0.5,200);

    fits = load('E:\Research_Projects\004_Propofol\data\experiments\scalp_EEG\model_fits\20230701\param_save.mat');
    tau = squeeze(fits.pars(:,1,:))*1e3;

    pVec2 = zeros(length(tRescaled),14);
    for j = 1:14
        dt = 1e-3;
        tmax = timeInfo.end_time(j)/60;
        T = 0:dt:tmax;
        N = length(T);

        x = zeros(4,N);
        u = zeros(N,1);
        u(and(T>timeInfo.infusion_onset(j)/60,T<timeInfo.object_drop(j)/60)) = 1e3; % 1mg/kg/min


        for i = 1:N-1
            x(:,i+1) = x(:,i) + dt*(A*x(:,i) + B*u(i));
        end
        t{j} = (60*T-timeInfo.object_drop(j));
        p{j} = x(4,:)';
        dtau = abs(timeInfo.infusion_onset(j)-timeInfo.object_drop(j))/60;
        tVec = (T-timeInfo.object_drop(j)/60)/dtau;
        pVec2(:,j) = interp1(tVec(:),p{j},tRescaled(:),'linear','extrap');
    end

    tauBL = nanmean(nanmean(tau(tRescaled<-1,:)));
    tauSmooth = smoothdata(tau,'gaussian',0.01);
    tauSmooth(tauSmooth==0)=nan;

    pBins = 10.^linspace(-3,max(pVec2(:)),100);
    dh = min(diff(pBins));
    pEdges = [pBins(1)-dh,pBins+dh/2];

    M = zeros(size(pBins));
    SD = zeros(size(pBins));
    for i = 1:length(pEdges)-1
        idcs = find(and(pVec2(:)>pEdges(i),pVec2(:)<=pEdges(i+1)));
        if(length(idcs)<5)
            M(i) = nan;
            SD(i) = nan;
        else
            M(i) = nanmean(tauSmooth(idcs));
            SD(i) = stderror(tauSmooth(idcs));
        end
    end

    % Convert from ug/ml to uM
    pBins = pBins/((178.27*1e6)*0.022)*1e6;
    pum = pVec2/((178.27*1e6)*0.022)*1e6;
    for i = 1:length(p)
       p{i} = p{i}/((178.27*1e6)*0.022)*1e6;
    end

    tauFold = tauSmooth/tauBL;
    M = M/tauBL;
end
function plotPaperResults()
    clrs = clrsPT.lines(4);

    tauOrser = [19.4,29.6,44.5]/19.4;
    pOrser = [0,0.5,2];
    plot(pOrser,tauOrser,'*','MarkerSize',5,'LineWidth',0.75,'color',clrs(1,:));

    tauKita = [1,1.2,1.4,1.45,2];
    pKita = [0.1,0.3,1,3,10];
    plot(pKita,tauKita,'*','MarkerSize',5,'LineWidth',0.75,'color',clrs(2,:));


    pWhittington2 = [0,0.5,1,2,5,10];
    tauWhittington2 = [17.5,17.5,19,22.5,50,75]/17.5;
    plot(pWhittington2,tauWhittington2,'*','MarkerSize',5,'LineWidth',0.75,'color',clrs(3,:));

    pWhittington1 = [0,0.5,1,2,5,10];
    tauWhittington1 = [17.5,17.5,19,48,58,75]/17.5;
    plot(pWhittington1,tauWhittington1,'*','MarkerSize',5,'LineWidth',0.75,'color',clrs(4,:));
end
folder = 'E:\Research_Projects\004_Propofol\data\simulations\raw\critical_dipole_correlation\S=1';

load(fullfile(folder,'m=0.00.mat'));
S1{1} = nanmean(C,3);
load(fullfile(folder,'m=0.63.mat'));
S1{2} = nanmean(C,3);
load(fullfile(folder,'m=0.86.mat'));
S1{3} = nanmean(C,3);
load(fullfile(folder,'m=0.91.mat'));
S1{4} = nanmean(C,3);
load(fullfile(folder,'m=0.95.mat'));
S1{5} = nanmean(C,3);
load(fullfile(folder,'m=0.98.mat'));
S1{6} = nanmean(C,3);
load(fullfile(folder,'m=0.99.mat'));
S1{7} = nanmean(C,3);

folder = 'E:\Research_Projects\004_Propofol\data\simulations\raw\critical_dipole_correlation\S=0';

load(fullfile(folder,'m=0.00.mat'));
S0{1} = nanmean(C,3);
load(fullfile(folder,'m=0.63.mat'));
S0{2} = nanmean(C,3);
load(fullfile(folder,'m=0.86.mat'));
S0{3} = nanmean(C,3);
% load(fullfile(folder,'m=0.91.mat'));
S0{4} = nan(size(S0{3}));
load(fullfile(folder,'m=0.95.mat'));
S0{5} = nanmean(C,3);
load(fullfile(folder,'m=0.98.mat'));
S0{6} = nanmean(C,3);
load(fullfile(folder,'m=0.99.mat'));
S0{7} = nanmean(C,3);


mValues = [0,0.63,0.86,0.91,0.95,0.98,0.99];

clrs = clrsPT.sequential(11); clrs = clrs(5:end,:);

A = cat(3,S1{:});
M_S1 = squeeze(nanmean(nanmean(A,1),2));

figureNB(7,7)
subplot(2,2,1);
plot(1./(1-mValues),M_S1,'k','LineWidth',1);
hold on;
for i = 1:length(mValues)
    m1 = 1./(1-mValues(i));
    
    M = nanmean(S1{i});
    s1 = quantile(M,0.05);
    s2 = quantile(M,0.95);
    plot(m1,nanmean(M),'.','color',clrs(i,:),'MarkerSize',10);
    hold on;
    line(m1*[1,1],[s1,s2],'LineWidth',1,'color',clrs(i,:));

    M = nanmean(S0{i});
    s1 = quantile(M,0.05);
    s2 = quantile(M,0.95);
    plot(m1,nanmean(M),'.','color',[0.6,0.6,0.6],'MarkerSize',10);
    hold on;
    line(m1*[1,1],[s1,s2],'LineWidth',1,'color',[0.6,0.6,0.6]);

end
set(gca,'xscale','log')
xlim([0.5,200])
xticks([1,10,100])
ylim([-0.02,0.23])

line(get(gca,'xlim'),[0,0],'LineWidth',1,'color',[0.6,0.6,0.6])

ylabel('Dipole correlation')
xlabel('1/(1-m)');
gcaformat;


folder = 'E:\Research_Projects\004_Propofol\data\simulations\raw\critical_dipole_correlation\m=0.98';

S = [0:0.2:1];

load(fullfile('E:\Research_Projects\004_Propofol\data\simulations\raw\critical_dipole_correlation\S=0','m=0.98.mat'));
% load(fullfile(folder,'dipole_correlation_S=0.0.mat'));
S_sub{1} = nanmean(C,3);
load(fullfile(folder,'dipole_correlation_S=0.2.mat'));
S_sub{2} = nanmean(C,3);
load(fullfile(folder,'dipole_correlation_S=0.4.mat'));
S_sub{3} = nanmean(C,3);
load(fullfile(folder,'dipole_correlation_S=0.6.mat'));
S_sub{4} = nanmean(C,3);
load(fullfile(folder,'dipole_correlation_S=0.8.mat'));
S_sub{5} = nanmean(C,3);
% load(fullfile(folder,'dipole_correlation_S=1.0.mat'));
load(fullfile('E:\Research_Projects\004_Propofol\data\simulations\raw\critical_dipole_correlation\S=1','m=0.98.mat'));
S_sub{6} = nanmean(C,3);


A = cat(3,S_sub{:});
M_Sub = squeeze(nanmean(nanmean(A,1),2));

subplot(2,2,2);

% plot(S,M_Sub,'color',clrs(end-1,:),'LineWidth',1);
hold on;
for i = 1:length(S)
    M = nanmean(S_sub{i});
    s1 = quantile(M,0.05);
    s2 = quantile(M,0.95);
    plot(S(i),nanmean(M),'.','color',clrs(end-1,:),'MarkerSize',10);
    hold on;
    line(S(i)*[1,1],[s1,s2],'LineWidth',1,'color',clrs(end-1,:));
end
xlim([0,1]);
ylim([-0.02,0.23])
ylabel('Dipole correlation')
xlabel('Optimality Index');
gcaformat;

S2 = permute(S,[1,3,2]);
B = ones(100,100).*S2;
[X,Y] = prepareCurveData(B(:),A(:));
FT = polyfit(X,Y,3);
t = linspace(0,1,1e3);
% plot(t,polyval(FT,t),'color',clrs(end-1,:),'LineWidth',1);
rts = roots([FT(1:end-1),FT(end)-0.15*nanmean(M)]);
% rts(imag(rts)==0);
plot(t,interp1(S,M_Sub,t,'cubic'),'color',clrs(end-1,:),'LineWidth',1);

% fun = fittype('A/(1+exp(-(x-C)/B))');
% FT = fit(X,Y,fun,'StartPoint',[1,2,0.5]);
% plot(t,FT(t),'color',clrs(end-1,:),'LineWidth',1);
interp1(M_Sub,S,0.3*nanmean(M),'linear');


folder = 'E:\Research_Projects\004_Propofol\data\simulations\raw\critical_dipole_correlation\PSD';
F = dir(folder); F = F(3:end);
for i = 1:length(F)
    load(fullfile(folder,F(i).name));
    psd(:,:,i) = P/2;
end

dataFolder = 'E:\Research_Projects\004_Propofol\manuscript\Nature Communications\_final_submission\_data'
% Get baseline EEG spectrum from propofol cohort
data = load(fullfile(dataFolder,'electrode2_Cz.mat'));
data.baseline = squeeze(nanmedian(data.psd(:,data.tRescaled<-1,:),2));

load(fullfile(dataFolder,'anatomy_cortical_pairwise_distance_distribution.mat'));
signed_area = A;
total_area = B;
N = 16e9;
dMids = 0.5*(rValues(2:end)+rValues(1:end-1));
nrnCount = mean(diff(signed_area),2)*200000;
nrnCount(end) = N-sum(nrnCount(1:end-1));
corr_kernel = @(d) exp(-d.^2/5);
rho_bar = sum(corr_kernel(dMids).*nrnCount)/sum(nrnCount');
SIG_N = @(rho) N+N*(N-1)*rho*rho_bar;



clrs = clrsPT.sequential(10); clrs = clrs(5:end,:);

MM = round(M_S1+5e-3,3);
MM(1) = 0;

subplot(2,2,3);
    m = sort([0.99,0.98,0.95,0.86,0.63,0]);
    h=[];
    for i = 1:length(m)
        y = psd(:,:,i);
        % y = y*SIG_N(MM(i));
        h(i) = plotwitherror(f,y,'CI','color',clrs(i,:),'LineWidth',1);
    end
    ylabel(['PSD (' char(956) 'V^2/Hz)']);
    xlabel('Frequency (Hz)')
    xlim([0.5,100]);
    ylim(10.^[-16.5,-13.5])
    yticks(10.^[-16,-15,-14])
    xticks([1,10,100]);
    xticklabels([1,10,100]);
    set(gca,'xscale','log')
    set(gca,'yscale','log')
    gcaformat;


rho_max = 0.02;
scaling_factor = rho_max/M_S1(end-1);
interp1(M_Sub,S,rho_max,'linear')

subplot(2,2,4);
    plotwitherror(data.freq,data.baseline,'M','color',[0.5,0.5,0.5]);
    h=[];
    for i = 1:length(m)
        y = psd(:,:,i);
        y = y*SIG_N(scaling_factor*MM(i));%0.25+0.2*randn(size(y));
        h(i) = plotwitherror(f,y,'CI','color',clrs(i,:),'LineWidth',1);
    end
    ylabel(['PSD (' char(956) 'V^2/Hz)']);
    xlabel('Frequency (Hz)')
    xlim([0.5,100]);
    ylim(10.^[-2,2])
    xticks([1,10,100]);
    xticklabels([1,10,100]);
    set(gca,'xscale','log')
    set(gca,'yscale','log')
    gcaformat;



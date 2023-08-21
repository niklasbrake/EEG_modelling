folder = 'E:\Research_Projects\004_Propofol\data\simulations\raw\critical_dipole_correlation\S=1';

load(fullfile(folder,'m=0.00'));
S1{1} = nanmean(C,3);
load(fullfile(folder,'m=0.63'));
S1{2} = nanmean(C,3);
load(fullfile(folder,'m=0.86'));
S1{3} = nanmean(C,3);
load(fullfile(folder,'m=0.91'));
S1{4} = nanmean(C,3);
load(fullfile(folder,'m=0.95'));
S1{5} = nanmean(C,3);
load(fullfile(folder,'m=0.98'));
S1{6} = nanmean(C,3);
load(fullfile(folder,'m=0.99'));
S1{7} = nanmean(C,3);

folder = 'E:\Research_Projects\004_Propofol\data\simulations\raw\critical_dipole_correlation\S=0';

load(fullfile(folder,'m=0.00'));
S0{1} = nanmean(C,3);
load(fullfile(folder,'m=0.63'));
S0{2} = nanmean(C,3);
load(fullfile(folder,'m=0.86'));
S0{3} = nanmean(C,3);
% load(fullfile(folder,'m=0.91'));
S0{4} = nan(size(S0{3}));
load(fullfile(folder,'m=0.95'));
S0{5} = nanmean(C,3);
load(fullfile(folder,'m=0.98'));
S0{6} = nanmean(C,3);
load(fullfile(folder,'m=0.99'));
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

line(get(gca,'xlim'),[0,0],'LineWidth',1,'color',[0.6,0.6,0.6])

ylabel('Dipole correlation')
xlabel('1/(1-m)');
gcaformat;
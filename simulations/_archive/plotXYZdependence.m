% load('E:\Research_Projects\Propofol\Modelling\neuron_simulations\simulation_20220310_13h54.mat')
% load('E:\Research_Projects\Propofol\Modelling\neuron_simulations\simulation_20220310_16h22.mat')
load('E:\Research_Projects\Propofol\Modelling\neuron_simulations\simulation_20220311_13h00.mat')
N = size(synPos,1);

dt = mean(diff(t));

fig = figureNB;
fig.Position = [0,0,11.4,12]
ax = axes('Position',[0.07,0.09,0.4,0.86]); hold on;
for i = 1:length(morphology.segments)
    x = morphology.coordinates(morphology.segments{i},3:5);
    plot(x(:,1),-x(:,2),'k','parent',ax); hold on;
end
set(gca,'DataAspectRatio',[1,1,1])
axis off;

CM = lines(N);
QN = sqrt(Qx.^2+Qy.^2+Qz.^2);
Xphi = atan(synPos(:,1)./synPos(:,2));
Qphi = zeros(size(Xphi));
Qphi_all = zeros(size(Xphi,1),length(t));
for i = 1:N
    [~,iM] = max(QN(i,:));
    Qphi(i) = atan(Qx(i,iM)./Qy(i,iM));
    Qphi_all(i,:) = atan(Qx(i,:)./Qy(i,:));
end

for i = 1:N
    plot(synPos(i,1),-synPos(i,2),'.','MarkerSize',5,'color',CM(i,:),'parent',ax)
    [~,iM] = max(QN(i,:));
    subplot(3,2,2)
        plot(synPos(i,1),Qx(i,iM),'.','color',CM(i,:)); hold on;
    subplot(3,2,6)
        plot(-synPos(i,2),Qy(i,iM),'.','color',CM(i,:)); hold on;
    subplot(3,2,4)
        plot(synPos(i,3),Qz(i,iM),'.','color',CM(i,:)); hold on;
end

subplot(3,2,2)
    xlabel('x position (um)')
    ylabel('nA um ms'); gcaformat;
    xlim([-300,500]); ylim([-30,30]);
subplot(3,2,6)
    xlabel('z position (um)')
    ylabel('nA um ms'); gcaformat;
    xlim([-300,500]); ylim([-30,30]);
subplot(3,2,4)
    xlabel('y position (um)')
    ylabel('nA um ms'); gcaformat;
    xlim([-300,500]); ylim([-30,30]);




fig = figureNB;
fig.Position = [0,0,4.6,2.6];
y = QN(1,:);
plot(y,'k','Linewidth',1); hold on;
[am,im] = max(y);
plot(im,am,'.','MarkerSize',20,'color',[1,0.5,0.5])
xlim([0,450])
axis off


fig = figureNB;
fig.Position = [0,0,4.6,2.6];
y = QN(1,:);
[am,im] = max(y);
fill([1:length(y),length(y):-1:1],[y,zeros(1,length(y))],[1,0.5,0.5],'LineStyle','none'); hold on;
plot(y,'k','Linewidth',1); hold on;
xlim([0,450])
axis off


% ft = fittype(@(a,c,d,x) a*(x>5).*(exp(-(x-5)*c)-exp(-(x-5)*d)));
% ft2 = fittype(@(a,c,x) a*(x>5).*(x-5).*(exp(-(x-5)*c)));
% for i = 1:size(QN,1)
%     % FT = fit(t(:),QN(i,:)',ft,'StartPoint',[80,0.6,0.7],'Lower',[0,0,0])
%     FT = fit(t(:),QN(i,:)',ft2,'StartPoint',[4,0.6],'Lower',[0,0])
%     C(i,:) = coeffvalues(FT);
%     cla
%     plot(t(:),QN(i,:)); hold on;
%     plot(FT);
%     drawnow;
% end

% fig = figureNB;
% qiv = quiver(0*Qx(:,1),0*Qx(:,1),Qx(:,1),Qy(:,1));
% xlim([-2.2,2.2]);
% ylim([-2.2,2.2]);
% set(gca,'DataAspectRatio',[1,1,1]);
% for i  = 1:length(t)
%     qiv.UData = Qx(:,i);
%     qiv.VData = Qy(:,i);
%     drawnow;
% end


% bins = linspace(0,pi/2,20);
% ctrs = 0.5*(bins(1:end-1)+bins(2:end));
% idcs = discretize(abs(Qphi),bins);
% x = length(ctrs);
% y = zeros(length(ctrs),length(t));
% for i = 1:length(ctrs)
%     i0 = find(idcs==i);
%     y(i,:) = nanmean(QN(i0,:));
% end
% imagesc(t,ctrs,y)


figure
subplot(1,2,1);
    plot(t,QN(1,:),'LineWidth',1);
    xlabel('Time (ms)')
    gcaformat
subplot(1,2,2);
    [p,f] = pspectrum([QN(1,:),zeros(1,1e5)],1e3/dt,'FrequencyLimits',[0.5,500],'FrequencyResolution',1);
    plot(f,p,'LineWidth',1)
    set(gca,'xscale','log');
    set(gca,'yscale','log');
    xticks([1,10,100]);
    xticklabels({1,10,100});
    xlabel('Frequency (Hz)');
    ylabel('log power');
    yticks([]);
    xlim([0.5,512])
    gcaformat


figure;
for i = 1:size(QN,1)
    [p(:,i),f] = pspectrum([QN(i,:),zeros(1,1e5)],1e3/dt,'FrequencyLimits',[0.5,500],'FrequencyResolution',1);
end
    plot(f,p,'LineWidth',0.5)
    set(gca,'xscale','log');
    set(gca,'yscale','log');
    xticks([1,10,100]);
    xticklabels({1,10,100});
    xlabel('Frequency (Hz)');
    ylabel('log power');
    yticks([]);
    xlim([0.5,512])
    gcaformat


bins = linspace(0,pi/2,20);
ctrs = 0.5*(bins(1:end-1)+bins(2:end));
idcs = discretize(abs(Qphi),bins);
x = length(ctrs);
y = zeros(length(ctrs),length(f));
for i = 1:length(ctrs)
    i0 = find(idcs==i);
    y(i,:) = nanmean(p(:,i0),2);
end
imagesc(f,ctrs,log(y))

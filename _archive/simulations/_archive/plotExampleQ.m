load('E:\Research_Projects\Propofol\Modelling\neuron_simulations\simulation_20220310_13h54.mat')
N = size(synPos,1);

dt = mean(diff(t));

fig = figureNB;
fig.Position = [0,0,11.4,12]
ax = axes('Position',[0.07,0.09,0.4,0.86]); hold on;
for i = 1:length(morphology.segments)
    x = morphology.coordinates(morphology.segments{i},3:5);
    plot(x(:,1),x(:,2),'k','parent',ax); hold on;
end
set(gca,'DataAspectRatio',[1,1,1])
axis off;

CM = lines(N);
QN = sqrt(Qx.^2+Qy.^2+Qz.^2);
for i = 1:N
        plot(synPos(i,1),synPos(i,2),'.','MarkerSize',5,'color',CM(i,:),'parent',ax)
    [~,iM] = max(QN(i,:));
    subplot(3,2,2)
        plot(synPos(i,1),dt*sum(Qx(i,:)),'.','color',CM(i,:)); hold on;
    subplot(3,2,6)
        plot(synPos(i,2),dt*sum(Qy(i,:)),'.','color',CM(i,:)); hold on;
    subplot(3,2,4)
        plot(synPos(i,3),dt*sum(Qz(i,:)),'.','color',CM(i,:)); hold on;
end

subplot(3,2,2)
    xlabel('x position (um)')
    ylabel('nA um ms'); gcaformat;
    xlim([-200,800]); ylim([-30,30]);
subplot(3,2,6)
    xlabel('z position (um)')
    ylabel('nA um ms'); gcaformat;
    xlim([-200,800]); ylim([-30,30]);
subplot(3,2,4)
    xlabel('y position (um)')
    ylabel('nA um ms'); gcaformat;
    xlim([-200,800]); ylim([-30,30]);




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
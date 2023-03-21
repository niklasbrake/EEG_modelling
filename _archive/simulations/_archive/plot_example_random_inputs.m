data = load('E:\Research_Projects\004_Propofol\Modelling\neuron_simulations\data\simulations\86955_simulation01.mat');
W = getEEG(data);

dt = mean(diff(data.t));

V = data.simulatenous.V;
CMb = [linspace(0,1,500)',linspace(0,1,500)',ones(500,1)];
CMr = [ones(500,1),linspace(1,0,500)',linspace(1,0,500)'];
% CMr = [ones(500,1),ones(500,1),linspace(1,0,500)'];
CM = [CMb;CMr];


cAxis = [linspace(min([-75,min(V(:))]),-65,500),linspace(-65.1,min([max(V(:)),-50]),500)];
% cAxis = linspace(-70,-50,1e3);
% CM = parula(1e3);
clr = CM(interp1(cAxis,1:1e3,V(:,1),'nearest'),:);

fig = figureNB;
fig.Position = [0,0,25,12]; fig.Color = 'k';
ax1 = axes('Position',[0.067,0.05,0.35,0.94]);
for i = 1:length(data.morphology.Vx)
    x = [data.morphology.Vx(i,:)',data.morphology.Vy(i,:)'];
    h(i) = plot(x(:,1),-x(:,2),'color','w','LineWidth',1); hold on;
end

% idx1 = interp1(cAxis,1:1e3,-75,'nearest','extrap');
cRange = [cAxis(1),cAxis(end)];

% colormap(CM(idx1:end,:));
colormap(CM);
C = colorbar('Location','northoutside');
C.Position = [0.04,0.68,0.125,0.04];
C.Ticks = ([cRange(1):5:cRange(2)]-cRange(1))/range(cRange);
C.TickLabels = [cRange(1):5:cRange(2)];
C.Color = 'w';
C.Label.String = 'Membrane Voltage (mV)';

set(gca,'DataAspectRatio',[1,1,1])
set(gca,'FontSize',12);
axis off;

ax2 = axes('Position',[0.5,0.25,0.42,0.44]);
h2 = plot(data.t,W,'color','y','LineWidth',1);
yl = get(gca,'ylim');
h2 = plot(data.t,W*nan,'color','y','LineWidth',1);
set(gca,'color','none');
ylim(yl);
% yticks([-1:0.5:1]);
set(get(gca,'xaxis'),'color','w'); xlabel('Time (ms)');
set(get(gca,'yaxis'),'color','w'); ylabel('EEG (pV)');
set(gca,'tickdir','out')
box off;
set(gca,'FontSize',18);
xlim([data.t(1),data.t(end)]);


v = VideoWriter(fullfile('E:\Research_Projects\004_Propofol\Modelling\neuron_simulations\figures\videos','Vm_random_synaptic_input_propofol.avi'));
v.FrameRate = 30;
open(v);


str = text(0,600,['Time: ' num2str(data.t(1),2)],'color','w');
for j = 1:16:size(V,2)-15
    clr = CM(interp1(cAxis,1:1e3,V(:,j),'nearest','extrap'),:);
    for i = 1:length(data.morphology.Vx)
        h(i).Color = clr(i,:);
    end
    h2.YData(j:j+15) = W(j:j+15);
    drawnow;
    writeVideo(v,getframe(fig));
end

close(v);

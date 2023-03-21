data = load('C:\Users\brake\Documents\temp\passive2\data.mat');
network2 = data.network.importResults;


idcs = randi(2e3,20,1);
[f,P] = network2.expectedEEGspectrum([],sa,idcs);


data = load('E:\Research_Projects\004_Propofol\data\experiments\scalp_EEG\analyzed\Cz_multitaper_mean.mat');
load('E:\Research_Projects\004_Propofol\data\experiments\scalp_EEG\raw\timeInformation.mat')
t0 = timeInfo.infusion_onset-timeInfo.object_drop;
pre = [];
for i = 1:length(t0)
    pre(:,i) = nanmedian(data.psd(:,data.time<t0(i),i),2);
    totalPower_Data(i) = sum(pre(:,i)*mean(diff(data.freq)));
end
freq = data.freq;


folder = 'E:\Research_Projects\004_Propofol\data\resources\cortical_column_Hagen\segmentations';
F = dir(folder); F = F(3:end);


layers = [0,-300,-300,-750,-1000,-1335];
x0(1,:) = [1150,layers(2)-20];
x0(2,:) = [-100,layers(2)-100];
x0(3,:) = [760,layers(4)];
x0(4,:) = [-40,layers(4)];
x0(5,:) = [1075,layers(4)];
x0(6,:) = [220,layers(5)];
x0(7,:) = [-300,layers(5)];
x0(8,:) = [1250,layers(5)];
x0(9,:) = [950,layers(6)-50];
x0(10,:) = [550,layers(6)+50];
x0(11,:) = [-130,layers(6)];

x_min = Inf;
x_max = -Inf;
y_min = Inf;
y_max = -Inf;

figureNB(8.9,4);
axes('Position',[0.02,0.08,0.38,0.83]);
LY = [0; -78; -500; -767; -975; -1250]*1.2;
for i = 1:length(LY)
    line([-1e3,2e3],[1,1]*LY(i),'color',[0.75,0.75,0.75]);
end

for k = 1:size(x0,1)
    morphData = load(fullfile(folder,F(k).name));
    if(isempty(strfind(F(k).name,'I_')))
        clr = clrsPT.qualitative_CM.red;
    else
        clr = clrsPT.qualitative_CM.blue;
    end
    for i = 1:length(morphData.area)
        x = morphData.x(i,:);
        y = morphData.y(i,:);
        z = morphData.z(i,:);
        L = vecnorm(diff([x;y;z]'));
        r = morphData.area(i)/2*pi/L;
        % line(x+x0(k,1),z+x0(k,2),'color',clr,'LineWidth',2./(1+10*exp(-(r-1)/2)));
        x_min = min(x_min,min(x+x0(k,1)));
        x_max = max(x_max,max(x+x0(k,1)));
        y_min = min(y_min,min(z+x0(k,2)));
        y_max = max(y_max,max(z+x0(k,2)));
    end
    set(gca,'DataAspectRatio',[1,1,1])
end
xlim([x_min,x_max]);
ylim([y_min,y_max]);
axis off;


[params,synFun] = synDetrend(f(f<100),mean(P(f<100,:),2)./mean(P(1,:)),0,'lorenz',[15e-3,1e-3,0,-1]);

axes('Position',[0.5,0.79,0.44,0.15]);
    plot(t,sum(eeg,2),'color','k');
    ylim([-4*std(sum(eeg,2)),4*std(sum(eeg,2))])
    xlim([-100,1000]);
    line([-50,-50],[-5,5]*1e-6,'color','k','LineWidth',1)
    text(-70,0,'10 nV','VerticalAlignment','middle','HorizontalAlignment','right','FontSize',7)
    line([0,100],[-6,-6]*1e-6,'color','k','LineWidth',1)
    text(50,-7e-6,'100 ms','VerticalAlignment','top','HorizontalAlignment','center','FontSize',7)   
    axis off;
axes('Position',[0.56,0.26,0.14,0.38]);
    plotwitherror(f,P,'Q','LineWidth',1,'color','k'); hold on;
    % plot(f,mean(P,2),'color','k','LineWidth',1); hold on;
    ylim([10^-16.5,10^-14.5]);
    yticks([1e-16,1e-15])
    set(gca,'yscale','log');
    set(gca,'xscale','log')
    xticks([1,10,100]);
    xticklabels([1,10,100]);
    xlabel('Frequency (Hz)');
    ylabel(['PSD (' char(956) 'V^2/Hz)'])
    xlim([1,100])
    plot(f,mean(P(1,:))*10.^synFun(f,[params(1:3),-Inf]),'--r');
    plot(f,mean(P(1,:))*10.^synFun(f,[params(1:2),-Inf,params(4)]),'--r');
    % plot(f,mean(P(1,:))*10.^synFun(f,params),'r');
    gcaformat;
    title('Unitary spectrum','fontweight','normal','fontsize',7)
axes('Position',[0.83,0.26,0.14,0.38]);
    plotwitherror(freq,pre,'Q','LineWidth',1,'color','b'); hold on;
    plotwitherror(f,P*25e9,'Q','LineWidth',1,'color','k'); hold on;
    set(gca,'xscale','log')
    set(gca,'yscale','log')
    ylim([1e-7,1e2]); yticks([1e-6,1e-2,1e2])
    xlim([1,100]);
    xticks([1,10,100]);
    xticklabels([1,10,100]);
    xlabel('Frequency (Hz)')
    ylabel(['PSD (' char(956) 'V^2/Hz)'])
    gcaformat

gcaformat(gcf);
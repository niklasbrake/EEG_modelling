
t = 1:11
load('C:\Users\brake\Documents\temp\correlations_0.98_subopt.mat')
m = sort([0.99,0.98,0.95,0.86,0.63,0]);
clrs = clrsPT.sequential(length(m)+4);
clrs = clrs(5:end,:);
clrs0 = clrs(5,:);
clrs = [linspace(0.6,clrs0(1),length(t));linspace(0.6,clrs0(2),length(t));linspace(0.6,clrs0(3),length(t))]';

subplot(1,3,3);
    t = pi*(1-linspace(0,1,11));
    plot(t,mean(C),'color','k','LineWidth',1);
    hold on;
    for i = 1:length(t)
        x = t(i);
        y = mean(C(:,i));
        y_lo = icdf('normal',0.05,0,1)*stderror(C(:,i));
        y_hi = icdf('normal',0.95,0,1)*stderror(C(:,i));
        plot(x,y,'.','color',clrs(i,:),'MarkerSize',10,'LineWidth',1)
        line([x,x],[y+y_lo,y+y_hi],'color',clrs(i,:),'linewidth',1);
    end
    ylim([-0.05,0.4])
    % xlim([0,pi])
    xlim([-0.2,pi+0.2])
    xticks([0,pi/2,pi])
    xticklabels({'0',[char(960) '/2'],char(960)})
    ylabel('Dipole correlation')
    xlabel('Embedding perturbation')
    gcaformat

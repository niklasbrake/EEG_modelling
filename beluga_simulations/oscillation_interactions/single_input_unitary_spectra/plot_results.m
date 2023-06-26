
load('E:\Research_Projects\004_Propofol\data\simulations\raw\peak_trend_sensitivity\trend_peak_interaction2.mat')
f = 0.1:0.1:500;

rhythms = load('E:\Research_Projects\004_Propofol\data\simulations\raw\peak_trend_sensitivity\rhythm_char.mat');
rhythms = load('E:\Research_Projects\004_Propofol\data\simulations\raw\peak_trend_sensitivity\rhythm_spectra.mat');
rhythms.psd = rhythms.psd([1,5,3,4,2]);


P3 = zeros(length(f),5,5);
P3(:,:,1) = mean(P2(:,:,1:5),2);
P3(:,:,2) = mean(P2(:,:,21:25),2);
P3(:,:,3) = mean(P2(:,:,11:15),2);
P3(:,:,4) = mean(P2(:,:,16:20),2);
P3(:,:,5) = mean(P2(:,:,6:10),2);



figureNB(18,11);
for i = 1:5
    for j = 2:5
        subplot(5,5,5*(i-1)+j)
        if(j==1)
            plot(f,P3(:,1,1),'color',[0.6,0.6,0.6],'LineWidth',1);
            hold on;
            plot(f,P3(:,j,i),'color','k','LineWidth',1);
        else
            plot(f,P3(:,1,i),'color','k','LineWidth',1);
            hold on;
            plot(f,P3(:,j,i),'color','r','LineWidth',1);
        end
        set(gca,'xscale','log')
        set(gca,'yscale','log')
        xlim([0.1,100])
        xticks([1,100]);
        ylim([1e-17,1e-15]);
        % if(i==2)
        %     ylim([1e-17,1e-13]);
        % else
        %     ylim([1e-17,1e-15]);
        % end
        if(i<5)
            xticklabels({});
        else
            xticklabels([1,100]);
        end
        % yticks([]);
        gcaformat;
        set(gca,'color','none');
        drawnow;
        xticks([]);
        yticks([]);
        % box on;
        if(i==5);
            xticks([0.1,1,10,100]);
            xticklabels([0.1,1,10,100]);
            xlabel('Frequency (Hz)');
            xax = get(gca,'xaxis');
            xax.TickLabelRotation = 0;
        end
    end
end
% set(gcf,'color','none')


for i = 1:5
    subplot(5,5,5*(i-1)+1)
    cla
    % plot(rhythms.f,mean(rhythms.psd{1},2),'color',[0.6,0.6,0.6],'LineWidth',0.5);
    hold on;
    plot(rhythms.f,mean(rhythms.psd{i},2),'color','k','LineWidth',1);
    gcaformat;
    set(gca,'color','none');
    set(gca,'xscale','log');
    set(gca,'yscale','log');
    xlim([0.1,100]);
    ylim([1e-6,1e0]);

    xticks([]);
    yticks([]);
    % box on;

    if(i==5);
        xticks([0.1,1,10,100]);
        xticklabels([0.1,1,10,100]);
        xlabel('Frequency (Hz)');
        xax = get(gca,'xaxis');
        xax.TickLabelRotation = 0;
    end
end


% figureNB(15,12);
% for i = 1:5
%     subplot(5,5,5*(i-1)+1)
%     plot(rhythms.t,rhythms.eeg{i},'color','k','LineWidth',0.5);
%     xlim([1.8e4,2.3e4])
%     ylim([0.5,1.5]);
%     axis off;
% end
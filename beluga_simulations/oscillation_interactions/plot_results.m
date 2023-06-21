[sa,X] = network_simulation_beluga.getHeadModel;

locations = randi(size(sa.cortex75K.vc,1),100,1); % Random location

%{
data = load('E:\Research_Projects\004_Propofol\data\simulations\raw\aperiodic_sensitivity\analyzed_results.mat');
location = 64030;
eeg = network_simulation_beluga.getEEG(dipoles,sa,location);
eeg1 = resample(eeg,1e3,16e3);
[psd,f2] = pmtm(detrend(eeg1,'constant'),2,[],1e3);
idcs(1,:) = (data.parameters(1,:) == 10);
idcs(2,:) = (data.parameters(2,:) == -58.5);
idcs(3,:) = (data.parameters(3,:) == 0.5);
idcs(4,:) = (data.parameters(4,:) == 0.98);
for i = 1:4
    K = unique(data.parameters(i,:));
    for k = 1:2
        temp = setdiff(1:4,i);
        idcs0 = find(prod(idcs(temp,:)));
        idcs2 = find(data.parameters(i,:)==K(k));
        P_crit(:,i) = mean(psd(:,intersect(idcs0,idcs2)),2);
    end
end
P_crit = P_crit(:,[4,1,2,3]);
%}


% load('E:\Research_Projects\004_Propofol\data\simulations\analyzed\oscillation_interactions.mat');
load('E:\Research_Projects\004_Propofol\data\simulations\raw\peak_trend_sensitivity\peak_trend_sensitivity.mat')
f = 0.1:0.1:500;

rhythms = load('E:\Research_Projects\004_Propofol\data\simulations\raw\peak_trend_sensitivity\rhythm_char.mat');



P3 = zeros(length(f),4,5);
P3(:,:,1) = mean(P2(:,:,1:4),2);
% P3(:,:,2) = interp1(f2,P_crit,f,'linear','extrap');
P3(:,:,2) = mean(P2(:,:,17:20),2);
P3(:,:,3) = mean(P2(:,:,9:12),2);
P3(:,:,4) = mean(P2(:,:,13:16),2);
P3(:,:,5) = mean(P2(:,:,5:8),2);


figureNB(15,12);
for i = 1:5
    for j = 1:4
        subplot(5,5,5*(i-1)+j+1)
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
            xticks([1,10,100]);
            xticklabels([1,10,100]);
            xlabel('Frequency (Hz)');
        end
    end
end
% set(gcf,'color','none')

for i = 1:5
    subplot(5,5,5*(i-1)+1)
    plot(rhythms.f,mean(rhythms.psd{1},2),'color',[0.6,0.6,0.6],'LineWidth',0.5);
    hold on;
    plot(rhythms.f,mean(rhythms.psd{i},2),'color','k','LineWidth',0.5);
    gcaformat;
    set(gca,'color','none');
    set(gca,'xscale','log');
    set(gca,'yscale','log');
    xlim([0.5,100]);
    ylim([1e-6,1e0]);

    xticks([]);
    yticks([]);
    % box on;

    if(i==5);
        xticks([1,10,100]);
        xticklabels([1,10,100]);
        xlabel('Frequency (Hz)');
    end
end


figureNB(15,12);
for i = 1:5
    subplot(5,5,5*(i-1)+1)
    plot(rhythms.t,rhythms.eeg{i},'color','k','LineWidth',0.5);
    xlim([1.8e4,2.3e4])
    ylim([0.5,1.5]);
    axis off;
end
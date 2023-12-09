
load(fullfile(dataFolder,'simulations','trend_peak_interactions','parameter_changes.mat'));

f = 0.1:0.1:500;

rhythms = load(fullfile(dataFolder,'simulations','trend_peak_interactions','rhythm_spectra.mat'));
rhythms.psd = rhythms.psd([1,5,3,4,2]);

P3 = zeros(length(f),5,5);
P3(:,:,1) = mean(P2(:,:,1:5),2);
P3(:,:,2) = mean(P2(:,:,21:25),2);
P3(:,:,3) = mean(P2(:,:,11:15),2);
P3(:,:,4) = mean(P2(:,:,16:20),2);
P3(:,:,5) = mean(P2(:,:,6:10),2);


[full_model,synFun] = fittingmodel('eq1');
p(:,1) = [0.01,0.002,-0.35,-0.2];
p(:,2) = [0.03,0.002,-0.16,-0.41];
p(:,3) = [0.01,0.002,-0.07,-0.6];
p(:,4) = [0.01,0.0025,-1,-0.02];
p(:,5) = [0.011,0.0025,-0.25,-0.28];



figureNB(18,11);
for i = 1:5
    for j = 1:5
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

figureNB(18*0.9,11);
fs = 100; t = 0:1/fs:1000; y = 1+0.1*sin(2*pi*t*2);
[psd,f] = pmtm(detrend(y,'constant'),2,[],fs);

X{5} = csvread('../../rhythm1.csv');
X{3} = csvread('../../rhythm2.csv');
X{4} = csvread('../../rhythm3.csv');
X{2} = csvread('../../rhythm4.csv');
X{1} = [X{2}(:,1),1+0.1*randn(size(X{2}(:,1)))];

for i = 1:5
    subplot(5,2,2*(i-1)+2)
    if(i~=2)
        plot(rhythms.f/5,mean(rhythms.psd{i},2),'color','k','LineWidth',1);
    else
        plot(rhythms.f,mean(rhythms.psd{i},2),'color','k','LineWidth',1);
    end
    if(i==5)
        plot(f,psd,'color','k','LineWidth',1);
    end
    gcaformat;
    set(gca,'color','none');
    set(gca,'xscale','log');
    set(gca,'yscale','log');
    xlim([0.1,100]);
    ylim([1e-8,1e0]);

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

    subplot(5,2,2*(i-1)+1)
        plot(X{i}(:,1),X{i}(:,2),'color','k');
        axis off;
end



function figureS5(dataFolder)

if(nargin<1)
    error('Path to data required as input argument. Data can be downloaded from link in README file.');
end 

load(fullfile(dataFolder,'simulations','trend_peak_interactions','parameter_changes.mat'));
f = 0.1:0.1:500;

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

clrs = clrsPT.lines(5);

figureNB(18,11);
for i = 1:5
    for j = 1:5
        subplot(5,6,6*(i-1)+j)
        if(j==1)
            % plot(f,P3(:,1,1),'color',[0.6,0.6,0.6],'LineWidth',1);
            plot(f,P3(1,j,1)*10.^synFun(f,p(:,1)),'color',[0.6,0.6,0.6]);
            hold on;
            plot(f,P3(:,j,i),'color',clrs(j,:),'LineWidth',1);
        else
            plot(f,P3(1,1,1)*10.^synFun(f,p(:,1)),'color',[0.6,0.6,0.6],'LineStyle','--');
            hold on;
            plot(f,P3(1,j,1)*10.^synFun(f,p(:,j)),'color',[0.6,0.6,0.6]);
            plot(f,P3(:,j,i),'color',clrs(j,:),'LineWidth',1);
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
        ylabel('log power')
        gcaformat;
        set(gca,'color','none');
        drawnow;
        xticks([]);
        yticks([]);

        if(i==5);
            xticks([0.1,1,10,100]);
            xticklabels([0.1,1,10,100]);
            xlabel('Frequency (Hz)');
            xax = get(gca,'xaxis');
            xax.TickLabelRotation = 0;
        end

        subplot(5,6,6*(i-1)+6)
        plot(f,10*log10(P3(:,j,i)./P3(1,j,1)./10.^synFun(f(:),p(:,j))),'color',clrs(j,:));
        hold on;
        set(gca,'xscale','log')
        xlim([0.1,100])
        xticks([1,100]);
        gcaformat;
        set(gca,'color','none');
        drawnow;
        xticks([]);
        ylabel('Power (dB)');
        ylim([-1,6]);

        if(i==5);
            xticks([0.1,1,10,100]);
            xticklabels([0.1,1,10,100]);
            xlabel('Frequency (Hz)');
            xax = get(gca,'xaxis');
            xax.TickLabelRotation = 0;
        end
    end
end
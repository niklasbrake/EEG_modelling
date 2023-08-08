for i = 1:14
    [time,freq,psd1(:,:,i),psd2(:,:,i)] = toy_model();
end


figureNB(14,8);
subplot(2,1,1);
    plot(t2-70,Y,'k');
    xlabel('LOC-aligned time (s)')
    xlim([-70,30]);
    yticks([]);
    ylabel('Toy model (neural dynamics)')
% subplot(2,1,2);
    % plot(T-70,eeg+20*randn(size(eeg)))
    % xlabel('LOC-aligned time (s)')
    % xlim([-70,30]);
    % yticks([]);
    % ylabel('Toy model (EEG)')

figureNB(14,8);
subplot(2,1,1);
    plot(T-70,eeg,'k');
    xlabel('LOC-aligned time (s)')
    xlim([-70,30]);
    yticks([]);
    ylabel('Toy model (EEG)')
subplot(2,2,3)
    P0 = log10(nanmedian(psd2./nanmean(psd2(:,time<20,:),2),3));
    imagesc(time-70,freq,P0)
    % colorbar
    ylim([0.5,50])
    axis xy;
    xlabel('LOC-aligned time (s)');
    ylabel('Frequency (Hz)')
subplot(2,2,4);
    plot(time-70,10*nanmean(P0(freq<=4,:)),'k');
    xlabel('LOC-aligned time (s)')
    ylabel('Delta power (dB)');
    xlim([-70,30]);
colormap(flipud(clrsPT.iridescent(1e3)))
gcaformat(gcf);

figureNB(14,8);
subplot(2,1,1);
    plot(t2-70,Y,'k');
    xlabel('LOC-aligned time (s)')
    xlim([-70,30]);
    yticks([]);
    ylabel('Toy model (neural dynamics)')
subplot(2,2,3)
    P1 = log10(nanmedian(psd1./nanmean(psd1(:,time<20,:),2),3));
    imagesc(time-70,freq,P1);
    % colorbar
    ylim([0.5,50])
    axis xy;
    xlabel('LOC-aligned time (s)');
    ylabel('Frequency (Hz)')
subplot(2,2,4);
    plot(time-70,10*nanmean(P1(freq<=4,:)),'k');
    xlabel('LOC-aligned time (s)')
    ylabel('Delta power (dB)');
    xlim([-70,30]);
colormap(flipud(clrsPT.iridescent(1e3)))
gcaformat(gcf);





figureNB(14,9);
axes('Units','centimeters','Position',[1.82, 4.67, 10.84, 2.73])
    plot(T-70,eeg,'k');
    xlabel('LOC-aligned time (s)')
    xlim([-70,30]);
    yticks([]);
    ylabel('Toy model (EEG)')
axes('Units','centimeters','Position',[1.82, 0.88, 4.68, 2.73])
    P0 = log10(nanmedian(psd2./nanmean(psd2(:,time<20,:),2),3));
    imagesc(time-70,freq,P0)
    % colorbar
    ylim([0.5,50])
    axis xy;
    xlabel('LOC-aligned time (s)');
    ylabel('Frequency (Hz)')
axes('Units','centimeters','Position',[7.98, 0.88, 4.68, 2.73]);
    plot(time-70,10*nanmean(P0(freq<=4,:)),'k');
    xlabel('LOC-aligned time (s)')
    ylabel('Delta power (dB)');
    xlim([-70,30]);
axes('Units','centimeters','Position',[1.82, 7.9, 10.84, 1])
    A = @(t) 1 + 1./(1+exp(-(t-60)/6));
    plot(t2-70,A(t2)-1,'k','LineWidth',1);
    ylim([-0.1,1.1]);
    xlabel('LOC-aligned time (s)')
    ylabel('\tau (ms)');
    yticks([0,0.5,1]);
    yticklabels([15,27.5,40]);
    yyaxis right;
    ylabel('A');
    yticks([0,0.5,1]);
    xlim([-70,30]);
    yticklabels([1,1.5,2]);
    yax = get(gca,'yaxis');
    yax(2).Color = 'k';
colormap(flipud(clrsPT.iridescent(1e3)))
gcaformat(gcf);


figureNB(18.3,10);
subplot(2,1,1);
    A = @(t) 1 + 1./(1+exp(-(t-60)/6));
    plot(t2-70,A(t2)-1,'k','LineWidth',1);
    xlabel('LOC-aligned time (s)')
    ylabel('\tau (ms)');
    yticks([0,0.5,1]);
    yticklabels([15,27.5,40]);
    yyaxis right;
    ylabel('A');
    yticks([0,0.5,1]);
    xlim([-70,30]);
    yticklabels([1,1.5,2]);
    yax = get(gca,'yaxis');
    yax(2).Color = 'k';
subplot(2,1,2);
    plot(T-70,eeg,'k','LineWidth',1)
    xlabel('LOC-aligned time (s)')
    xlim([-70,30]);
    yticks([]);
    ylabel('Toy model (EEG)')
    gcaformat(gcf);

function [time,freq,psd1,psd2] = toy_model()
  dt = 1/1024;

    tcov = 0:dt:0.4;
    m = length(tcov);
    ecov = @(tau) exp(-tcov(:)/tau);

    T = 0:dt:100;


    t2 = 0:dt:100.4; t2 = t2(:);

    tau = @(t) (20+30./(1+exp(-(t-60)/10)))*1e-3;
    A = @(t) 1 + 1./(1+exp(-(t-60)/10));
    B = @(t) 1 + 1.25./(1+exp(-2*(t-70)));

    Y = sin(2*pi*t2).*B(t2) + 0.5*randn(size(t2));
    eeg = zeros(length(T),1);

    for i = 1:length(T);
        eeg(i) = A(T(i))*sum(Y(i:i+m-1).*ecov(tau(T(i))));
    end

    [b0,a0] = butter(1,0.1*dt,'high');

    figureNB;
    subplot(2,2,1)
        [freq,time,psd1] = eegfft(t2,Y,4,3);
        imagesc(time,freq,log10(psd1./nanmedian(psd1(:,time<10),2)));
        colorbar
        set(gca,'CLim',[0,1])
        ylim([0.5,50])
        axis xy
    subplot(2,2,2);
        [freq,time,psd2] = eegfft(t2,eeg,4,3);
        imagesc(time,freq,log10(psd2./nanmedian(psd2(:,time<10),2)));
        colorbar
        set(gca,'CLim',[0,1])
        ylim([0.5,50])
        axis xy
    subplot(2,1,2);
        plot(T,eeg,'k');
        drawnow;
end
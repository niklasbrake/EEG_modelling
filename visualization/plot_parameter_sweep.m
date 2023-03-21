load('E:\Research_Projects\004_Propofol\data\simulations\raw\parameter_search\power_spectra_all.mat')

idcs0 = find((dgam==1).*(dlam==1).*(dtau==1))
P_baseline = P(:,:,idcs0);
tau_baseline = tau(idcs0);


fig = figureNB(20,6);
for j = 1:length(tau_baseline)
    t = tau_baseline(j);
    idcs = setdiff(find(tau==tau_baseline(j)),idcs0);
    [~,idcs2] = sortrows([dtau(idcs)',dgam(idcs)',dlam(idcs)']);
    idcs = idcs(idcs2);
    for i = 1:length(idcs)
        ax = axes('units','centimeters');
        ax.Position = [0.1+i,0.1+j,0.75,0.75];
        plot(f,nanmean(P_baseline(:,:,j),2),'k','LineWidth',0.5); hold on;
        plot(f,nanmean(P(:,:,idcs(i)),2),'r','LineWidth',0.5); hold on;
        set(gca,'xscale','log');
        set(gca,'yscale','log');
        xlim([0.5,100]);
        xticks([1,10,100]);
        xticklabels([1,10,100]);
        xlabel('Frequency (Hz)')
        ylabel(['PSD (' char(956) 'V^2/Hz)'])
        ylim([1e-18,1e-13]);
        if(i>1)
            ylabel('');
            yticklabels({});
        end
        if(j>1)
            xlabel('');
            xticklabels({});
        end
        % axis off;
        gcaformat;
        drawnow;
        set(gca,'FontSize',2)
        set(gca,'LineWidth',0.25)
        text(20,1e-13,sprintf('\\tau = %d\n\\Delta\\tau =%d\n\\Delta\\gamma =%d\n\\Delta\\lambda =%.1f',tau(idcs(i)),dtau(idcs(i)),dgam(idcs(i)),dlam(idcs(i))),'FontSize',2,'HorizontalAlignment','left','VerticalAlignment','top');
    end
end

filename = 'E:\Research_Projects\004_Propofol\data\simulations\raw\parameter_search\parameter_search.svg';
print(filename,fig,'-dsvg');
network = network_simulation_beluga('C:/Users/brake/Documents/temp');
M = 2e3;
dh = 4;
tmax = 300e3;


network.branchNo = 0;
[ids,ts,ei] = network.simulatespikes_critplane(M,tmax);
N0 = computeAvalancheSize(ids,ts,dh);

network.branchNo = 0.98;
[ids,ts,ei] = network.simulatespikes_critplane(M,tmax);
N1 = computeAvalancheSize(ids,ts,dh);

% bins = 10.^linspace(0,log10(2*M),20);
bins = linspace(0,500,500);
c0 = histcounts(N0,bins,'Normalization','pdf');
c1 = histcounts(N1,bins,'Normalization','pdf');

figureNB;
subplot(1,3,1);
    plot(log(bins(1:end-1)),log(c0),'.-r');
    hold on;
    plot(log(bins(1:end-1)),log(c1),'.-k')
    xlabel('log avalanche size')
    ylabel('log P')
    L = legend('Poisson (m=0)','m=0.98')
    L.ItemTokenSize = [10,10];
    L.Box = 'off';
    % ylim([1e-4,1e-1])
    gcaformat
    xlim(log([bins(1),bins(end)]));


function N = computeAvalancheSize(ids,ts,dh)
    [h,bins] = histcounts(ts,'BinWidth',dh);
    t = bins(2:end)-dh/2;

    isAv = 0;
    counter = 1;
    t0 = [];
    for i = 1:length(h)
        if(h(i)>0)
            if(isAv==0)
                t0(counter,1) = t(max(i-1,1));
            end
            isAv = 1;
        else
            if(isAv == 1)
                t0(counter,2) = t(i);
                counter = counter+1;
            end
            isAv = 0;
        end
    end

    N = zeros(size(t0,1),1);
    for j = 1:size(t0,1)
        N(j) = sum((ts>=t0(j,1)).*(ts<=t0(j,2)));
        % N(j) = length(unique(ids(find((ts>=t0(j,1)).*(ts<t0(j,2))))));
    end
end
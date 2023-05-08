N = 250;
r = 0.1*ones(N,1);
V = 0.1*(1-0.1);
gam = icdf('normal',r,0,1);
fun = @(x,S) S-mvncdf(gam(1:2),0,[1,x;x,1])+r(1)*r(2);

rhoValues = [0,0.01,0.02,0.05,0.1,0.2];
bins = [-0.5:80.5];

figureNB;
for i = 1:length(rhoValues)
    S = rhoValues(i)*(0.1*(1-0.1));
    l0 = fzero(@(x) fun(x,S),[-0.999,0.999]);

    L = ones(N,N);
    L = l0*L+(1-l0)*eye(N);
    x = mvnrnd(gam(:),L,5e4)>0;
    ss = sum(x,2);
    h(:,i) = histcounts(ss,'BinEdges',bins,'Normalization','pdf');

    [ts,ids] = find(x);
    switch i
    case 1
        subplot(2,2,1);
        raster(ids,ts,gcf);
        xlim([1,200]);
        ylabel('Neuron index');
        yticks([50:50:250]);
        set(get(gca,'yaxis'),'visible','on');
        axis ij;
        xticks([100,200]);
    case 4
        subplot(2,2,2);
        raster(ids,ts,gcf);
        xlim([1,200]);
        ylabel('Neuron index');
        yticks([50:50:250]);
        set(get(gca,'yaxis'),'visible','on');
        axis ij;
        xticks([100,200]);
    case 5
        subplot(2,2,3);
        raster(ids,ts,gcf);
        xlim([1,200]);
        ylabel('Neuron index');
        yticks([50:50:250]);
        set(get(gca,'yaxis'),'visible','on');
        axis ij;
        xticks([100,200]);
    end
end
subplot(2,2,4)
    plot([0:80],h,'LineWidth',1)
    ylabel('Frequency');
    xlabel('Syncrhonous spikes');


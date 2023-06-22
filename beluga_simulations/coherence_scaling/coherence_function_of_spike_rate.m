x = csvread('C:\Users\brake\Documents\GitHub\EEG_modelling\beluga_simulations\oscillation_interactions\rhythm4.csv');
plot(x(:,1),x(:,2));


Y = cumsum(x(:,2))/sum(x(:,2));



nValues = 10.^[0:8];
figureNB;
clrs = clrsPT.sequential(length(nValues));
for i=  1:length(nValues)
    N = nValues(i);
    ts = interp1(Y,x(:,1),rand(N,1));
    H1 = histcounts(ts,[0:10:1e4]);
    ts = interp1(Y,x(:,1),rand(N,1));
    H2 = histcounts(ts,[0:10:1e4]);

    eeg = detrend([H1(:),H2(:)],'constant');
    P0 = sum(pmtm(eeg,2,[],1e3),2);
    [P1,f] = pmtm(sum(eeg,2),2,[],1e3);

    C(:,i) = (P1-P0)./P0;
    plot(f,C(:,i),'color',clrs(i,:));
    hold on;
    set(gca,'xscale','log')
    drawnow;
end

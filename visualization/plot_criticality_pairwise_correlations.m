folder = 'E:\Research_Projects\004_Propofol\data\simulations\raw\correlations';
F = dir(folder); F = F(3:end);
for i=  1:length(F)
    data{i} = csvread(fullfile(folder,F(i).name));
end
m = cellfun(@(x)str2num(x(3:end-4)),{F(:).name});
[m,I] =sort(m);
fr = cellfun(@(x)length(x),data)/(30000^2)*100;
fr = fr(I);

subplot(1,3,2);
    plot(1./(1-m(:)),fr(:),'.-k','MarkerSize',15,'LineWidth',1)
    ylabel('Correlated synapse pairs (STTC>0.25)')
    xlabel('1/(1-m)')
    set(gca,'xscale','log')
    ylim([0,1])
    yticks([0,0.5,1])
    yticklabels({'0%','0.1%','1%'})
    gcaformat
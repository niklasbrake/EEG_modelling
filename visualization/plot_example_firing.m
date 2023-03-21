folder = 'C:\Users\brake\Documents\temp\action_potentials2\LFPy';

F = dir(folder);
F = F(3:end);

for i = 1:length(F)
    data = csvread(fullfile(folder,F(i).name),1,0);
    V(:,i) = data(:,5);
    t = data(:,1);
end


fid = fopen('C:\Users\brake\Documents\temp\action_potentials2\postsynaptic_network\mTypes.txt')
morphs = textscan(fid,'%s');
morphs = morphs{1};
[~,morphs] = cellfun(@(x)fileparts(x),morphs,'UniformOutput',false);
fclose(fid);

ei = [];
N = [];
Y = cell(2,1);
for i = 1:size(V,2)
    [y,x] = findpeaks(V(:,i),'MinPeakHeight',-10);
    N(i) = length(x);
    ei(i) = isempty(strfind(morphs{i},'E_'));
    if(ei(i))
        Y{1} = [Y{1},N(i)/5];
    else
        Y{2} = [Y{2},N(i)/5];
    end
end


[~,idcs] = unique(morphs);
idcs(3) = idcs(3)+1;

fig = figureNB(18.3,10);

xh = 0.165;
yh = 0.2;
dx = (0.95-xh-0.05)/3;
dy = (0.95-yh-0.1)/2;

for i = 1:12
    y = 2-floor((i-1)/4);
    x = mod(i-1,4);
    axes('Position',[0.08+x*dx,0.1+dy*y,xh,yh]);
    if(i<12)
        plot(t,V(:,idcs(i)),'k');
        xlim([1e3,6e3]);
        xticks([1e3:1e3:6e3]);
        xticklabels([0:5]);
        xlabel('Time (s)');
        ylabel('Voltage (mV)');
        title(morphs{idcs(i)},'interpret','latex');
    end
end
boxplotNB(1,N(ei==0)/5,'r',5)
hold on;
boxplotNB(2,N(ei==1)/5,'b',5)
set(gca,'yscale','log')
xlim([0.5,2.5])
ylabel('Frequency (Hz)')
xticks([1,2]);
xticklabels({'E','I'});
ylim([0.1,100])
yticks([0.1,1,10,100])
gcaformat(fig);




figureNB
plot(network.results.dipoles(:,:,1)+linspace(-100,100,3))

masterPath = 'E:\Research_Projects\004_Propofol\Modelling\neuron_simulations\data\simulations\dipole_orientation';
network = network_simulation(masterPath);

M = 500;
network = network.initialize_postsynaptic_network(M);

network.branchingIdx = 0;
EI = rand(1,M)<0.85;
network = network.setsynapsecount(M);
network.save_presynaptic_network((1:M)',150*ones(M,1),EI);


load('E:\Research_Projects\004_Propofol\Modelling\neuron_simulations\data\cortical_column_Hagen\segment_areas.mat');

SD = zeros(M,3);
for j = 1:M
    mData = nrnSegs.(network.neurons(j).mType);
    sa = mData.area;
    seg = randperm(length(sa),1);
    pos = [mean(mData.x,2),mean(mData.y,2),mean(mData.z,2)];
    segDir = pos./vecnorm(pos,2,2);
    SD(j,:) = segDir(seg,:);
    filename = fullfile(network.postNetwork,'connections',[network.neurons(j).name '.csv']);
    fid = fopen(filename,'w');
    fprintf(fid,'%s,%s\n',int2str(seg),int2str(j));
    fclose(fid);
end
save(fullfile(masterPath,'synapse_compartments.mat'),'SD');

network.tmax = 240;
network = network.simulate;
network = network.importResults;

d = network.results.dipoles;

dMag = squeeze(max(vecnorm(d,2,2)));
dMag./max(dMag)*3

PD = squeeze(nanmedian(d./vecnorm(d,2,2)))';
PD = [PD(:,1),PD(:,3),-PD(:,2)];
PD(find(EI==0),:) = -PD(find(EI==0),:);


clrs = [230, 25, 75;
    60, 180, 75;
    255, 225, 25;
    0, 130, 200;
    245, 130, 48;
    70, 240, 240;
    240, 50, 230;
    250, 190, 212;
    0, 128, 128;
    220, 190, 255;
    170, 110, 40;
    255, 250, 200;
    128, 0, 0;
    170, 255, 195;
    0, 0, 128;
    128, 128, 128;
    255, 255, 255;
    0, 0, 0]/255;


[theta0,phi0] = cart2sph(SD(:,1),SD(:,2),SD(:,3));
[theta1,phi1] = cart2sph(PD(:,1),PD(:,2),PD(:,3));
SD_pol = [theta0,phi0];
PD_pol = [theta1,phi1];

% gIdcs = findgroups({network.neurons.mType});
gIdcs = cellfun(@(x)isempty(strfind(x,'E')),{network.neurons.mType});
clrs = [1,0,0;0,0,1];
% clrs = clrsPT.lines(length(unique(gIdcs)));
clrs = clrs(gIdcs+1,:);

fig = figureNB(5,2.75);
ax(1) = axes('Position',[0.16,0.3,0.3,0.55]);
ax(2) = axes('Position',[0.6,0.3,0.3,0.55]);
for i = 1:2
    axes(ax(i));
    scatter(SD_pol(:,i),PD_pol(:,i),2*dMag,clrs,'filled'); hold on;
    FT = fitlm(SD_pol(:,i),PD_pol(:,i),'Weights',dMag);
    % plot([-1,1],FT.predict([-1;1]),'color','k');    R2(i) = FT.Rsquared.Adjusted;
end
axes(ax(1));
    gcaformat;
    title('Azimuth','FontSize',7);
    ylabel('Dipole orientation')
    xlim([-pi,pi]); xticks([-pi,0,pi]); xticklabels({'-\pi','0','\pi'});
    ylim([-pi,pi]); yticks([-pi,0,pi]); yticklabels({'-\pi','0','\pi'});
    line([-pi,pi],[-pi,pi]);
    % text(pi,-pi,sprintf('R^2=%.2f',R2(1)),'FontSize',6,'HorizontalAlignment','right','VerticalAlignment','bottom')
axes(ax(2));
    gcaformat;
    title('Elevation','FontSize',7);
    xlabel('Synapse orientation')
    xlim([-pi,pi]/2); xticks([-pi/2,0,pi/2]); xticklabels({'-\pi/2','0','\pi/2'})
    ylim([-pi,pi]/2); yticks([-pi/2,0,pi/2]); yticklabels({'-\pi/2','0','\pi/2'})
    line([-pi,pi]/2,[-pi,pi]/2);
    % text(pi/2,-pi/2,sprintf('R^2=%.2f',R2(2)),'FontSize',6,'HorizontalAlignment','right','VerticalAlignment','bottom')
gcaformat(fig);